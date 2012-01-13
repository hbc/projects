(ns snp-assess.reference-test
  (:use [midje.sweet]
        [snp-assess.reference]))

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
