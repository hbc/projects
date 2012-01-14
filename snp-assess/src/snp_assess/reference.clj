;; Create reference set of variations for control lanes from known
;; mixed populations

(ns snp-assess.reference
  (:import [org.biojava3.core.sequence.io FastaReaderHelper]
           [org.broadinstitute.sting.utils.variantcontext
            VariantContextBuilder Allele])
  (:use [clojure.java.io]
        [clojure.string :only [join]]
        [clojure.algo.generic.functor :only [fmap]])
  (:require [clj-yaml.core :as yaml]))

;; External file interaction

(defn get-fasta-seq-map [in-file]
  "Parse FASTA input file to in-memory map of ids to sequences."
  (let [seq-map (FastaReaderHelper/readFastaDNASequence (file in-file) true)]
    (reduce #(assoc %1 %2 (->> %2 (.get seq-map) (.toString))) {} (.keySet seq-map))))

;; Generate VCF output file components using GATK API

(defn convert-to-vc [contig pos base-freqs]
  "Convert base frequency information into a VCF VariantContext."
  (letfn [(to-alleles [bases]
            (for [[i [base _]] (map-indexed vector bases)]
              (Allele/create base (= i 0))))
          (to-freqs [bases]
            (join "," (rest (map second bases))))]
    (let [ordered-bases (sort-by second > (vec base-freqs))]
      (-> (VariantContextBuilder. contig contig pos (+ 1 pos) (to-alleles ordered-bases))
          (.attributes (-> {}
                           (#(if (> (count ordered-bases) 1)
                               (assoc % "AF" (to-freqs ordered-bases))
                               %))))
          (.make)))))

;; Organize input FASTA and YAML frequency description into frequency
;; of bases at each position.

(defn gen-ref-at-pos [ref-bases percents]
  "Generate reference bases and frequencies at a position."
  (letfn [(add-base-percent [coll key]
            (let [base (get ref-bases key)
                  freq (get percents key)]
              (assoc coll base (+ (get coll base 0)
                                  freq))))]
    (fmap #(/ % 100.0)
          (reduce add-base-percent {} (keys percents)))))

(defn per-pos-seqs [seqs]
  "Lazy list of maps with bases at each position in a set of sequences."
  (for [i (range (-> seqs vals first count))]
    (reduce #(assoc %1 %2 (str (nth (get seqs %2) i)))
            {} (keys seqs))))

(defn gen-ref
  "Generate VCF formatted reference sequence for defined population."
  [seqs percents]
  {:pre [(= 100.0 (apply + (vals percents)))
         (= 1 (count (set (map count (vals seqs)))))]}
  (map-indexed vector
               (map #(gen-ref-at-pos % percents) (per-pos-seqs seqs))))

(defn -main [ref-fasta ref-config]
  (let [seqs (get-fasta-seq-map ref-fasta)
        config (-> ref-config slurp yaml/parse-string)
        percents (into {} (for [[k v] (:reference config)] [(name k) v]))]
    (gen-ref seqs percents)))
