(ns snp-assess.classify-test
  (:use [midje.sweet]
        [snp-assess.classify]
        [snp-assess.config])
  (:require [fs]))

(let [data-dir (fs/join (fs/cwd) "test" "data")
      vrn-file (fs/join data-dir "coverage_pos" "pos.tsv")
      data-file (fs/join data-dir "raw" "raw_variations.tsv")]
  (facts "Read variant positions from file"
    (read-vrn-pos vrn-file 5.0) => #{["HXB2_IUPAC_93-5" 951 "A"] ["HXB2_IUPAC_93-5" 953 "T"]}
    (read-vrn-pos vrn-file 3.0) => #{["HXB2_IUPAC_93-5" 951 "A"]})

  (facts "Classification data from file"
    (prep-classifier-data data-file vrn-file default-config)))
