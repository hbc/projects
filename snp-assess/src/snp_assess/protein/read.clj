(ns snp-assess.protein.read
  "Read based assessment of protein changes. Handles multiple variant changes
   at a single amino acid as well as phasing between connected amino acids."
  (:import [net.sf.samtools SAMFileReader SAMFileReader$ValidationStringency])
  (:use [clojure.java.io]
        ;[bcbio.run.broad :only [index-bam]]
        ))

(defn calc-aas-from-reads
  "Prepare amino acid changes by position based on raw read data."
  [bam-file ref-file prot-map & {:keys [count-file]}]
  ;(index-bam bam-file)
  (SAMFileReader/setDefaultValidationStringency SAMFileReader$ValidationStringency/LENIENT)
  (with-open [bam-rdr (SAMFileReader. (file bam-file))
              bam-iter (.iterator bam-rdr)]
    (doseq [x (take 5 (iterator-seq bam-iter))]
      (println x))))