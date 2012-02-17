;; Create reference set of variations for control lanes from known
;; mixed populations

(ns snp-assess.reference
  (:import [org.biojava3.core.sequence.io FastaReaderHelper]
           [org.broadinstitute.sting.utils.variantcontext Allele
            VariantContextBuilder]
           [org.broadinstitute.sting.utils.codecs.vcf StandardVCFWriter
            VCFHeader VCFInfoHeaderLine VCFHeaderLineCount VCFHeaderLineType]
           [net.sf.picard.reference ReferenceSequenceFileFactory]
           [net.sf.picard.sam CreateSequenceDictionary])
  (:use [clojure.java.io]
        [clojure.string :only [join]]
        [clojure.algo.generic.functor :only [fmap]]
        [ordered.map :only [ordered-map]]
        [bcbio.variation.variantcontext :only [parse-vcf]])
  (:require [clj-yaml.core :as yaml]
            [fs.core :as fs]))

;; External file interaction

(defn get-fasta-seq-map [in-file]
  "Parse FASTA input file to in-memory map of ids to sequences."
  (let [seq-map (FastaReaderHelper/readFastaDNASequence (file in-file) true)]
    (reduce #(assoc %1 %2 (->> %2 (.get seq-map) (.toString))) {} (.keySet seq-map))))

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

;; Generate VCF output file components using GATK API

(defn convert-to-vc
  "Convert base frequency information into a VCF VariantContext."
  [contig pos base-freqs & {:keys [depth aa-finder]}]
  (letfn [(to-alleles [bases]
            (for [[i [base _]] (map-indexed vector bases)]
              (Allele/create base (= i 0))))
          (to-freqs [bases]
            (join "," (rest (map second bases))))
          (aa-changes [bases]
            (join "," (map #(aa-finder pos %) (keys bases))))]
    (let [ordered-bases (sort-by second > (vec base-freqs))]
      (-> (VariantContextBuilder. contig contig (+ 1 pos) (+ 1 pos) (to-alleles ordered-bases))
          (.attributes (-> {}
                           (#(if (> (count ordered-bases) 1)
                               (assoc % "AF" (to-freqs ordered-bases))
                               %))
                           (#(if-not (nil? depth)
                               (assoc % "DP" depth)
                               %))
                           (#(if-not (nil? aa-finder)
                               (assoc % "AA_CHANGE" (aa-changes ordered-bases))
                               %))))
          (.make)))))

(defn- to-seq-dict [in-fasta]
  (let [idx-file (str in-fasta ".dict")]
    (if-not (fs/exists? idx-file)
      (.instanceMain (CreateSequenceDictionary.)
                     (into-array [(str "R=" in-fasta)
                                  (str "O=" idx-file)])))
    (.getSequenceDictionary
     (ReferenceSequenceFileFactory/getReferenceSequenceFile (file in-fasta)))))

(defn- get-header-with-attrs
  "Retrieve VCFHeader with allele frequency (AF) and Depth (DP) info line"
  []
  (VCFHeader. #{(VCFInfoHeaderLine. "AF" VCFHeaderLineCount/A
                                    VCFHeaderLineType/Float "Allele Frequency")
                (VCFInfoHeaderLine. "DP" 1 
                                    VCFHeaderLineType/Integer "Total Depth")
                (VCFInfoHeaderLine. "AA_CHANGE" VCFHeaderLineCount/A
                                    VCFHeaderLineType/String
                                    "Amino acid change caused by variants")}))

(defn write-vcf-ref [seqs percents ref-fasta out-file]
  "Write output reference file in VCF format"
  (let [ref-name (ffirst (sort-by second > (vec percents)))]
    (with-open [writer (StandardVCFWriter. (file out-file) (to-seq-dict ref-fasta))]
      (.writeHeader writer (get-header-with-attrs))
      (doseq [vc (map (fn [[i x]] (convert-to-vc ref-name i x))
                      (gen-ref seqs percents))]
        (.add writer vc)))))

(defn write-vcf-calls
  "Write output calls in VCF format"
  [vcs out-file ref-file]
  (with-open [writer (StandardVCFWriter. (file out-file) (to-seq-dict ref-file))]
    (.writeHeader writer (get-header-with-attrs))
    (doseq [vc vcs]
      (.add writer vc))))

;; Retrieve reference data information from VCF input

(defn- frequency-attributes [vc]
  "Retrieve frequency values from VCF attributes."
  (let [raw (-> vc :attributes (get "AF"))
        raw-list (map #(Float/parseFloat %) (into [] (cond
                                                      (string? raw) [raw]
                                                      (nil? raw) []
                                                      :else raw)))]
    (cons (- 1.0 (apply + raw-list)) raw-list)))

(defn- vc-to-freqs [vc]
  "Convert a variant context into list of bases and associated frequencies."
  (let [alleles (map #(.getBaseString %) (into [] (.getAlleles (:vc vc))))
        freqs (frequency-attributes vc)]
    (map vector alleles
         (map (partial * 100.0) freqs))))

(defn read-vcf-ref
  "Read reference information into ordered map of positions + base and frequency.
   Optionally filter by a maximum frequency to include."
  ([in-file]
     (reduce #(assoc %1 [(:chr %2) (:start %2) (:base %2)] (:freq %2)) (ordered-map)
              (flatten
               (for [vc (parse-vcf in-file)]
                 (for [[allele freq] (vc-to-freqs vc)]
                   {:chr (:chr vc) :start (- (:start vc) 1) :base allele :freq freq})))))
  ([in-file max-freq]
     (into (ordered-map) (filter #(-> % val (<= max-freq))
                                 (read-vcf-ref in-file)))))

(defn -main [ref-fasta ref-config out-file]
  (let [seqs (get-fasta-seq-map ref-fasta)
        config (-> ref-config slurp yaml/parse-string)
        percents (into {} (for [[k v] (:reference config)] [(name k) v]))]
    (write-vcf-ref seqs percents ref-fasta out-file)))
