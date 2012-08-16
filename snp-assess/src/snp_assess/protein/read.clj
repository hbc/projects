(ns snp-assess.protein.read
  "Read based assessment of protein changes. Handles multiple variant changes
   at a single amino acid as well as phasing between connected amino acids."
  (:import [net.sf.samtools SAMFileReader SAMFileReader$ValidationStringency]
           [org.biojava3.core.sequence DNASequence])
  (:use [clojure.java.io]
        [bcbio.run.broad :only [index-bam]]
        [snp-assess.reference :only [get-variants-by-pos]])
  (:require [clj-yaml.core :as yaml]))

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
          :read-seq (str (nth map-seq (+ read-start i)))})))))

(defn- is-minority-variant?
  [variant-map pos]
  (when-let [cur-vars (get variant-map ((juxt :contig :ref-i) pos))]
    (contains? (set (:alt-alleles cur-vars)) (:read-seq pos))))

(defn- aa-on-read
  "Retrieve amino acid changes predicted by current read."
  [rec prot-map count-map variant-map]
  (let [map-seq (.getReadString rec) 
        orig-seq (if (.getReadNegativeStrandFlag rec)
                   (-> (DNASequence. map-seq)
                       .getReverseComplement
                       .getSequenceAsString)
                   map-seq)]
    (->> (mapped-reads-by-pos rec map-seq)
         (filter (partial is-minority-variant? variant-map))
         (group-by #(:aa-pos (get prot-map (:ref-i %)))))
    {:count (get count-map orig-seq 1)}))

(defn- get-read-counts
  "Map of original read sequence to counts from input YAML dump file."
  [count-file]
  (reduce (fn [coll x]
            (assoc coll (:seq x) (:count x)))
          {} (-> count-file slurp yaml/parse-string)))

(defn calc-aa-from-reads
  "Prepare amino acid changes by position based on raw read data."
  [bam-file call-file ref-file prot-map & {:keys [count-file]}]
  (index-bam bam-file)
  (SAMFileReader/setDefaultValidationStringency SAMFileReader$ValidationStringency/LENIENT)
  (let [count-map (get-read-counts count-file)
        variant-map (get-variants-by-pos call-file ref-file)]
    (with-open [bam-rdr (SAMFileReader. (file bam-file))
                bam-iter (.iterator bam-rdr)]
      (doseq [x (take 20 (iterator-seq bam-iter))]
        (aa-on-read x prot-map count-map variant-map)))))