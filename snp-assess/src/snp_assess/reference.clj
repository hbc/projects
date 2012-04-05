(ns snp-assess.reference
  "Create reference variation set for control from known mixed populations."
  (:import [org.biojava3.core.sequence.io FastaReaderHelper FastaReader
            GenericFastaHeaderParser DNASequenceCreator]
           [org.biojava3.core.sequence.compound AmbiguityDNACompoundSet]
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
        [bcbio.variation.variantcontext :only [parse-vcf get-vcf-source]]
        [snp-assess.off-target :only [parse-pos-line]])
  (:require [clj-yaml.core :as yaml]
            [fs.core :as fs]))

;; External file interaction

(defn get-fasta-seq-map [in-file]
  "Parse FASTA input file to in-memory map of ids to sequences."
  (let [seq-map (.process (FastaReader. (file in-file) (GenericFastaHeaderParser.)
                                        (DNASequenceCreator.
                                         (AmbiguityDNACompoundSet/getDNACompoundSet))))]
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
            (map-indexed (fn [i [base _]]
                           (Allele/create (str base) (= i 0)))
                         bases))
          (to-freqs [bases]
            (join "," (rest (map second bases))))
          (aa-changes [bases]
            (join "," (map #(aa-finder pos %) (rest (keys bases)))))]
    (let [ordered-bases (sort-by second > (vec base-freqs))]
      (-> (VariantContextBuilder. contig contig (+ 1 pos) (+ 1 pos) (to-alleles ordered-bases))
          (.attributes (-> {}
                           (#(if (> (count ordered-bases) 1)
                               (assoc % "AF" (to-freqs ordered-bases))
                               %))
                           (#(if-not (nil? depth)
                               (assoc % "DP" depth)
                               %))
                           (#(if (and (not (nil? aa-finder))
                                      (> (count ordered-bases) 1))
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
                                    "Amino acid change caused by a variant")}))

(defn write-vcf-calls
  "Write output calls in VCF format"
  [vcs out-file ref-file]
  (with-open [writer (StandardVCFWriter. (file out-file) (to-seq-dict ref-file))]
    (.writeHeader writer (get-header-with-attrs))
    (doseq [vc vcs]
      (.add writer vc))))

(defn write-vcf-ref-percents
  "Write control reference file in VCF format using mixture of inputs."
  [seqs config ref-fasta out-file]
  (let [percents (into {} (for [[k v] (:reference config)] [(name k) v]))
        ref-name (ffirst (sort-by second > (vec percents)))]
    (write-vcf-calls (map (fn [[i x]] (convert-to-vc ref-name i x))
                          (gen-ref seqs percents))
                     out-file ref-fasta)))

(defn write-vcf-ref-variable
  "Write control reference in VCF using defined constant and variable regions."
  [seqs config ref-fasta out-file]
  {:pre [(= 1 (count (vals seqs)))]}
  (letfn [(parse-freq-file [freqs in-file]
            (with-open [rdr (reader in-file)]
              (reduce (fn [coll [_ pos base freq]]
                        (assoc coll pos [base (/ freq 100.0)]))
                      freqs (->> rdr
                                 line-seq
                                 (map parse-pos-line)))))
          (seq-by-index [seqs defined]
            (->> (map-indexed (fn [i x] [i (str x)]) (first (vals seqs)))
                 (drop-while #(< (first %) (apply min (keys defined))))
                 (take-while #(< (first %) (apply max (keys defined))))))
          (get-base-freqs [ref-base [defined-base defined-freq]]
            (let [ambig-map {"R" ["A" "G"] "Y" ["C" "T"] "M" ["C" "A"]
                             "W" ["A" "T"] "S" ["G" "C"] "K" ["T" "G"]}
                  final-ref-base (if-not (contains? ambig-map ref-base) ref-base
                                         (first (remove #(= defined-base %)
                                                        (get ambig-map ref-base))))]
              (cond
               (= defined-freq 1.0) {defined-base defined-freq}
               (nil? defined-freq) {final-ref-base 1.0}
               :else {defined-base defined-freq
                      final-ref-base (- 1.0 defined-freq)})))]
    (let [defined (-> {}
                      (parse-freq-file (:constant config))
                      (parse-freq-file (:variable config)))]
      (write-vcf-calls (map (fn [[i ref-base]] (convert-to-vc (first (keys seqs)) i
                                                              (get-base-freqs ref-base
                                                                              (get defined i))))
                            (seq-by-index seqs defined))
                       out-file ref-fasta))))

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
  ([in-file ref-file]
     (with-open [vcf-source (get-vcf-source in-file ref-file)]
       (reduce #(assoc %1 [(:chr %2) (:start %2) (:base %2)] (:freq %2)) (ordered-map)
               (flatten
                (for [vc (parse-vcf vcf-source)]
                  (for [[allele freq] (vc-to-freqs vc)]
                    {:chr (:chr vc) :start (- (:start vc) 1) :base allele :freq freq}))))))
  ([in-file ref-file max-freq]
     (into (ordered-map) (filter #(-> % val (<= max-freq))
                                 (read-vcf-ref in-file ref-file)))))

(defn -main [ref-fasta ref-config out-file]
  (let [seqs (get-fasta-seq-map ref-fasta)
        config (-> ref-config slurp yaml/parse-string)]
    (case (:method config)
      "variable" (write-vcf-ref-variable seqs config ref-fasta out-file)
      (write-vcf-ref-percents seqs config ref-fasta out-file))))
