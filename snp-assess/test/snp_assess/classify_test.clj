(ns snp-assess.classify-test
  (:use [midje.sweet]
        [snp-assess.classify])
  (:require [fs]))

(facts "Read variant positions from file"
  (let [vrn-file (fs/join (fs/cwd) "test" "data" "coverage_pos" "pos.tsv")]
    (read-vrn-pos vrn-file 5.0) => #{["HXB2_IUPAC_93-5" 951 "A"] ["HXB2_IUPAC_93-5" 953 "T"]}
    (read-vrn-pos vrn-file 3.0) => #{["HXB2_IUPAC_93-5" 951 "A"]}))
