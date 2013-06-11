(ns snp-assess.reference-test
  (:use [midje.sweet]
        [snp-assess.reference])
  (:require [fs.core :as fs]))

(tabular
 (fact "Correctly generate reference bases at a position."
   (gen-ref-at-pos ?refs ?freqs) => ?expected)
 ?refs ?freqs ?expected
 {:a "G" :b "C"} {:a 90.0 :b 10.0} {"C" 0.1 "G" 0.9}
 {:a "G" :b "G"} {:a 90.0 :b 10.0} {"G" 1.0}
 {:a "G" :b "G" :c "G"} {:a 80 :b 5.0 :c 15} {"G" 1.0}
 {:a "G" :b "G" :c "A"} {:a 80 :b 5.0 :c 15} {"G" 0.85 "A" 0.15})

(fact "Generate map of per-sequence base positions"
  (let [seq-map (per-pos-seqs {:a "GATC" :b "AAAA" :c "TTTT"})]
    (count seq-map) => 4
    (first seq-map) => {:a "G" :b "A" :c "T"}
    (second seq-map) => {:a "A" :b "A" :c "T"}))

(fact "Generate VCF output from reference information."
  (let [alt-vc (convert-to-vc "test" 0 {"G" 0.9 "A" 0.1})
        ref-vc (convert-to-vc "test" 0 {"G" 1.0})]
    (.getChr alt-vc) => "test"
    (.getStart ref-vc) => 1
    (-> alt-vc (.getAttributes) (.get "AF")) => "0.10"
    (.getAttributes ref-vc) => {}
    (-> alt-vc (.getReference) (.getBaseString)) => "G"
    (.isBiallelic alt-vc) => true
    (.isBiallelic ref-vc) => false))

(fact "Parse known frequency information from reference VCF"
  (let [vcf-file (str (fs/file "test" "data" "count_data" "known_freqs.vcf"))
        ref-file (str (fs/file "test" "data" "reference" "hxb2.fa"))
        freqs (read-vcf-ref vcf-file ref-file)
        filter-freqs (read-vcf-ref vcf-file ref-file 5.1)]
    (first freqs) => (contains [["HXB2" 225 "G"] (roughly 61.1)])
    (count freqs) => 6
    (last freqs) => [["HXB2" 227 "T"] 100.0]
    (count filter-freqs) => 2
    (first filter-freqs) => (contains [["HXB2" 226 "A"] (roughly 5.0)])))
