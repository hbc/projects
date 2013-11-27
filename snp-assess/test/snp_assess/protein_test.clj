(ns snp-assess.protein-test
  (:use [midje.sweet]
        [snp-assess.protein.calc]
        [snp-assess.protein.read]
        [snp-assess.annotate])
  (:require [me.raynes.fs :as fs]
            [bcbio.run.itx :as itx]))

(background
 (around :facts
         (let [data-dir (str (fs/file "test" "data"))
               ref-file (str (fs/file data-dir "protein" "hxb2-ref.fa"))
               mut-file (str (fs/file data-dir "protein" "known_integrase_mutations.csv"))
               prot-map (prep-protein-map {:files ref-file :known mut-file
                                           :frame-offset 0 :aa-offset 49})
               bam-file (str (fs/file data-dir "align" "S6-test.bam"))
               count-file (str (fs/file data-dir "align" "S6-test-counts.yaml"))
               call-file (str (fs/file data-dir "align" "S6-calls.vcf"))
               call-out-file (itx/add-file-part call-file "protein")
               ann-out-file (itx/add-file-part call-out-file "annotate")]
           (doseq [x [call-out-file ann-out-file]]
             (itx/remove-path x))
           ?form)))

(facts "Convert reference sequence into map of codons and known amino acid changes."
  (get prot-map 0) => {:offset 0 :codon [\G \C \C] :aa-pos 49 :known {}}
  (get prot-map 2) => {:offset 2, :codon [\G \C \C] :aa-pos 49 :known {}}
  (-> prot-map (get 319) :known (get "S")) => ["EVG30" "RAL15"]
  (-> prot-map (get 313) :known) => {"Y" ["EVG15"]})

(facts "Calculate amino acid changes caused by variations."
  (calc-aa-change {0 {:offset 0 :codon [\G \C \C] :aa-pos 1 :known {}}}
                  {:position 0 :new "C"}) => "A1P"
  (calc-aa-change {2 {:offset 2 :codon [\G \C \C] :aa-pos 1 :known {}}}
                  {:position 2 :new "C"}) => "A1A"
  (calc-aa-change {0 {:offset 0 :codon [\G \G \G] :aa-pos 1 :known {}}
                   2 {:offset 2 :codon [\G \G \G] :aa-pos 1 :known {}}}
                  {:position 0 :new "T"} {:position 2 :new "A"}) => "G1*"
  (calc-aa-change {0 {:offset 0 :codon [\G \C \C]
                      :aa-pos 1 :known {"P" ["EVG15"]}}}
                  {:position 0 :new "C"}) => "A1P_EVG15"
  (calc-aa-change {0 {:offset 0 :codon [\G \C \C]
                      :aa-pos 1 :known {"P" ["EVG30" "RAL15"]}}}
                  {:position 0 :new "C"}) => "A1P_EVG30_RAL15")


(facts "Generate amino acid changes based on input reads."
  (let [kmer-size 13
        aa-changes (calc-aa-from-reads bam-file call-file ref-file
                                       prot-map kmer-size :count-file count-file)]
    (get aa-changes 32) => {"G59R" 13, "G59G" 97}
    (annotate-calls-w-aa bam-file call-file ref-file prot-map
                         kmer-size :count-file count-file) => call-out-file
    (add-mv-annotations call-out-file bam-file ref-file) => ann-out-file))
