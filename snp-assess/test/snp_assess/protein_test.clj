(ns snp-assess.reference-test
  (:use [midje.sweet]
        [snp-assess.protein.calc]
        [snp-assess.protein.read])
  (:require [fs.core :as fs]))

(background
 (around :facts
         (let [data-dir (str (fs/file "test" "data"))
               ref-file (str (fs/file data-dir "protein" "hxb2-ref.fa"))
               mut-file (str (fs/file data-dir "protein" "known_integrase_mutations.csv"))
               prot-map (prep-protein-map {:files ref-file :known mut-file
                                           :frame-offset 0 :aa-offset 49})
               bam-file (str (fs/file data-dir "align" "S6-test.bam"))
               count-file (str (fs/file data-dir "align" "S6-test-counts.yaml"))]
           ?form)))

(fact "Convert reference sequence into map of codons and known amino acid changes."
  (get prot-map 0) => {:offset 0 :codon [\G \C \C] :aa-pos 49 :known {}}
  (get prot-map 2) => {:offset 2, :codon [\G \C \C] :aa-pos 49 :known {}}
  ;;prot-map => nil
  (-> prot-map (get 319) :known (get "S")) => ["EVG30" "RAL15"]
  (-> prot-map (get 313) :known) => {"Y" ["EVG15"]})

(fact "Calculate amino acid changes caused by variations."
  (calc-aa-change {0 {:offset 0 :codon [\G \C \C] :aa-pos 1 :known {}}} 0 "C") => "A1P"
  (calc-aa-change {2 {:offset 2 :codon [\G \C \C] :aa-pos 1 :known {}}} 2 "C") => "A1A"
  (calc-aa-change {2 {:offset 2 :codon [\T \G \G] :aa-pos 1 :known {}}} 2 "A") => "W1*"
  (calc-aa-change {0 {:offset 0 :codon [\G \C \C]
                      :aa-pos 1 :known {"P" ["EVG15"]}}} 0 "C") => "A1P_EVG15"
  (calc-aa-change {0 {:offset 0 :codon [\G \C \C]
                      :aa-pos 1 :known {"P" ["EVG30" "RAL15"]}}} 0 "C") => "A1P_EVG30_RAL15")

(fact "Generate amino acid changes based on input reads."
  (calc-aas-from-reads bam-file ref-file prot-map :count-file count-file) => nil)
