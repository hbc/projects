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
    (let [data (prep-classifier-data data-file vrn-file default-config)]
      (count data) => 3 ; 3 positions of interest in file
      (map last (first data)) => [:neg :neg :neg :neg :pos] ; classifications
      (count (second data)) => 2 ; minority variants at second position
      (ffirst data) => (contains [0.228 :neg])
      )))

(facts "Remap raw data for classification"
  (finalize-raw-data ["notused" 20 1.0E-4 100] :test default-config) =>
    (contains [(roughly 0.516) (roughly 9.0E-4) 0.4 :test])
  (finalize-raw-data ["notused" 30 1.0E-2 50] :test2 default-config) =>
    (contains [(roughly 0.8387) (roughly 0.0999) 0.2 :test2]))
