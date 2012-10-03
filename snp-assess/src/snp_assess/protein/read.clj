(ns snp-assess.protein.read
  "Read based assessment of protein changes. Handles multiple variant changes
   at a single amino acid as well as phasing between connected amino acids."
  (:import [net.sf.samtools SAMFileReader SAMFileReader$ValidationStringency]
           [org.biojava3.core.sequence DNASequence]
           [org.broadinstitute.sting.utils.variantcontext VariantContextBuilder]
           [org.broadinstitute.sting.utils.codecs.vcf
            VCFHeader VCFInfoHeaderLine VCFHeaderLineCount VCFHeaderLineType])
  (:use [clojure.java.io]
        [bcbio.run.broad :only [index-bam]]
        [bcbio.variation.variantcontext :only [parse-vcf get-vcf-iterator write-vcf-w-template]]
        [snp-assess.protein.calc :only [calc-aa-change]]
        [snp-assess.reference :only [get-variants-by-pos]])
  (:require [clojure.string :as string]
            [clj-yaml.core :as yaml]
            [bcbio.run.itx :as itx]))

(defn- mapped-reads-by-pos
  "Lazy list of base and position for a mapped read."
  [rec map-seq]
  (flatten
   (for [block (.getAlignmentBlocks rec)]
     (let [read-start (dec (.getReadStart block))
           ref-start (dec (.getReferenceStart block))]
       (for [i (range (.getLength block))]
         {:contig (.getReferenceName rec)
          :ref-i (+ ref-start i)
          :read-i (+ read-start i)
          :read-seq (str (nth map-seq (+ read-start i)))})))))

(defn- at-read-end?
  "Identify bases located within kmer distance of front and end of reads"
  [map-seq kmer-size x]
  (let [half-kmer (/ (dec kmer-size) 2)]
    (or (< (:read-i x) half-kmer)
        (>= (:read-i x) (- (count map-seq) half-kmer)))))

(defn- reformat-variant-to-change
  "Convert variant position information to format for aa calculation.
   Returns nil for non-analyzed positions and pre-filtered calls."
  [variant-map pos]
  (when-let [cur-vars (get variant-map ((juxt :contig :ref-i) pos))]
    (when (contains? (set (cons (:ref cur-vars) (:alt-alleles cur-vars))) (:read-seq pos))
      {:position (:ref-i pos)
       :majority (:ref cur-vars)
       :new (:read-seq pos)})))

(defn- call-aa
  [prot-map changes]
  {:position (apply max (map :position changes))
   :aa (apply calc-aa-change (cons prot-map changes))})

(defn- aa-on-read
  "Retrieve amino acid changes predicted by current read."
  [rec prot-map count-map variant-map kmer-size]
  (let [map-seq (.getReadString rec) 
        orig-seq (if (.getReadNegativeStrandFlag rec)
                   (-> (DNASequence. map-seq)
                       .getReverseComplement
                       .getSequenceAsString)
                   map-seq)
        cnt (get count-map orig-seq 1)]
    (->> (mapped-reads-by-pos rec map-seq)
         (remove (partial at-read-end? map-seq kmer-size))
         (map (partial reformat-variant-to-change variant-map))
         (remove nil?)
         (group-by #(:aa-pos (get prot-map (:position %))))
         vals
         (filter #(= 3 (count %)))
         (map (partial call-aa prot-map))
         (map #(assoc % :count cnt)))))

(defn- get-read-counts
  "Map of original read sequence to counts from input YAML dump file."
  [count-file]
  (reduce (fn [coll x]
            (assoc coll (:seq x) (:count x)))
          {} (-> count-file slurp yaml/parse-string)))

(defn calc-aa-from-reads
  "Prepare amino acid changes by position based on raw read data."
  [bam-file call-file ref-file prot-map kmer-size & {:keys [count-file]}]
  
  (SAMFileReader/setDefaultValidationStringency SAMFileReader$ValidationStringency/LENIENT)
  (let [count-map (get-read-counts count-file)
        variant-map (get-variants-by-pos call-file ref-file)]
    (with-open [bam-rdr (SAMFileReader. (file bam-file) (file (index-bam bam-file)))
                bam-iter (.iterator bam-rdr)]
      (reduce (fn [coll x]
                (let [counts (get coll (:position x) {})]
                  (assoc coll (:position x)
                         (assoc counts (:aa x)
                                (+ (:count x) (get counts (:aa x) 0))))))
              {} (flatten (map #(aa-on-read % prot-map count-map variant-map kmer-size)
                               (iterator-seq bam-iter)))))))

(defn- add-aa-to-vc
  "Add amino acid change information to the current variant context."
  [change-map vc]
  (letfn [(get-change-freqs [i]
            (when-let [changes (get change-map i)]
              (let [total (apply + (vals changes))]
                {:freqs (when (> (count changes) 1)
                          (->> (map (fn [[x n]] {:aa x :freq (/ n total)}) changes)
                               (sort-by :freq >)
                               rest))
                 :total total})))
          (maybe-add-aa-change [attrs freqs]
            (if freqs
              (-> attrs
                  (assoc "AA_CHANGE" (string/join "," (map :aa freqs)))
                  (assoc "AA_AF" (string/join "," (map #(format "%.2f" (* 100.0 (:freq %)))
                                                       freqs))))
              attrs))]
    (if-let [changes (get-change-freqs (dec (:start vc)))]
      (-> (VariantContextBuilder. (:vc vc))
          (.attributes (-> (:attributes vc)
                           (maybe-add-aa-change (:freqs changes))
                           (assoc "AA_DP" (:total changes))))
          .make)
      (:vc vc))))

(defn- add-aa-info-to-header
  [_ header]
  (doto header
    (.addMetaDataLine (VCFInfoHeaderLine. "AA_CHANGE" VCFHeaderLineCount/UNBOUNDED
                                          VCFHeaderLineType/String
                                          "Amino acid change caused by minority variants"))
    (.addMetaDataLine (VCFInfoHeaderLine. "AA_AF" VCFHeaderLineCount/UNBOUNDED
                                          VCFHeaderLineType/Float
                                          "Amino acid change frequencies"))
    (.addMetaDataLine (VCFInfoHeaderLine. "AA_DP" 1 VCFHeaderLineType/Integer
                                          "Number of reads contributing to amino acid calls"))))

(defn annotate-calls-w-aa
  "Add grouped amino acid calls to a VCF file of call information."
  [bam-file call-file ref-file prot-map kmer-size & {:keys [count-file]}]
  (let [out-file (itx/add-file-part call-file "protein")]
    (when (itx/needs-run? out-file)
      (let [change-map (calc-aa-from-reads bam-file call-file ref-file prot-map kmer-size
                                           :count-file count-file)]
        (with-open [vcf-iter (get-vcf-iterator call-file ref-file)]
          (write-vcf-w-template call-file {:out out-file}
                                (map (partial add-aa-to-vc change-map) (parse-vcf vcf-iter))
                                ref-file
                                :header-update-fn add-aa-info-to-header))))
    out-file))
