(ns snp-assess.annotate
  "Annotate minority variant calls with useful metrics for identifying potential false positives."
  (:require [bcbio.run.broad :as broad]
            [bcbio.run.itx :as itx]))

(defn add-mv-annotations
  "Add minority variant annotations to variant calls"
  [in-vcf align-bam ref-file]
  (let [file-info {:out-vcf (itx/add-file-part in-vcf "annotate")}
        args ["-R" ref-file
              "-I" align-bam
              "--variant" in-vcf
              "-o" :out-vcf
              "-A" "MeanNeighboringBaseQuality"
              "-A" "FisherStrand"]]
    (broad/index-bam align-bam)
    (broad/run-gatk "VariantAnnotator" args file-info {:out [:out-vcf]})
    (:out-vcf file-info)))