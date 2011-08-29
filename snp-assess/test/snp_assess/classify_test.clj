(ns snp-assess.classify-test
  (:use [midje.sweet]
        [snp-assess.classify]
        [snp-assess.config])
  (:require [fs]))

(let [data-dir (fs/join (fs/cwd) "test" "data")
      vrn-file (fs/join data-dir "coverage_pos" "pos.tsv")
      data-file (fs/join data-dir "raw" "raw_variations.tsv")
      config (assoc default-config :classification {:max-pct 5.0 :naive-min-score 1.0})]
  (facts "Read variant positions from file"
    (read-vrn-pos vrn-file 5.0) => #{["HXB2_IUPAC_93-5" 951 "A"] ["HXB2_IUPAC_93-5" 953 "T"]}
    (read-vrn-pos vrn-file 3.0) => #{["HXB2_IUPAC_93-5" 951 "A"]})

  (facts "Classification data from file"
    (let [data (prep-classifier-data data-file vrn-file config)]
      (count data) => 8 ; Total minority variants
      (count (first data)) => 4 ; items at each position
      (map last data) => [0 0 0 0 1 0 0 0] ; classifications
      (first data) => (contains [0.228 0])))

  (facts "Build a full classifier from file data"
    (let [c (train-classifier data-file vrn-file config)]
      (map float (.coefficients c)) => [0.0 0.0 0.0 0.0 0.125]
      (.toString c) => (contains "Linear Regression Model")
      ((classifier-checker c config) 20 1.0E-3 200) => false)))

(facts "Remap raw data for classification"
  (finalize-raw-data ["notused" 20 1.0E-4 100] :test default-config) =>
    (contains [(roughly 0.516) (roughly 0.0909) 0.4 :test])
  (finalize-raw-data ["notused" 30 1.0E-2 50] :test2 default-config) =>
    (contains [(roughly 0.8387) (roughly 1.0) 0.2 :test2]))
