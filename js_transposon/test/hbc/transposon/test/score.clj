;; Tests for normalizing and scoring combined position counts

(ns hbc.transposon.test.score
  (:use [hbc.transposon.score]
        [midje.sweet])
  (:require [clojure.string :as string]
            [incanter.core :as icore])
  (:import [java.io StringReader]))

(let [ds (icore/dataset [:chr :pos :one :two :three :seq]
                        [["chr1" 10  0  50 10 "GG"]
                         ["chr1" 20 10 200 10 "CC"]])
      config {:experiments [{:name "one" :expnorm "auto"}
                            {:name "two" :expnorm "auto"}
                            {:name "three" :expnorm "auto"}]
              :algorithm {:rownorm "auto"}}
      config2 {:experiments [{:name "one"}
                            {:name "two"}
                            {:name "three"}]
               :algorithm {:rownorm ""}}]
  (fact "Normalize a CSV file by reads and experiment counts."
    (let [ncounts (normalize-counts ds config :base 100)
          ncounts-raw (normalize-counts ds config2)
          n2counts (normalize-pos-ratios ncounts config)
          n2counts-raw (normalize-pos-ratios ncounts-raw config2)]
      (icore/sel ncounts :cols :one) => [0.0 100.0]
      (icore/sel ncounts :cols :two) => [20.0 80.0]
      (icore/col-names ncounts) => [:chr :pos :one :two :three :seq]
      (icore/sel n2counts :cols :one) => [0.0 1.0]
      (icore/sel n2counts :cols :two) => [0.4 0.8]
      (icore/col-names n2counts) => [:chr :pos :one :two :three :seq]
      (icore/sel ncounts-raw :cols :one) => [0 10]
      (icore/sel n2counts-raw :cols :one) => [0 10]
      (fact "Calculate dataset statistics"
        (let [stats (summarize-count-statistics n2counts)]
          (count stats) => 3
          (first stats) => (contains [:one 0.5 0.5]))))))

(let [ds (icore/dataset [:chr :pos :one :two :three :seq]
                        [["chr1" 10 1.0 0.0 1.0 "GG"]
                         ["chr1" 20 1.0 0.5 1.0 "GG"]])
      config {:experiments [{:name "one" :controls ["three"]}
                            {:name "two"}
                            {:name "three"}]}]
  (fact "Filter dataset to remove single experiment only rows."
    (let [fds (filter-by-multiple ds config)]
      (icore/nrow fds) => 1
      (icore/sel fds :rows 0) => ["chr1" 20 1.0 0.5 1.0 "GG"])))

(let [ds (icore/dataset [:chr :pos :one :two :seq]
                        [["chr1" 10 0.5 1.0 "GG"]
                         ["chr1" 20 1.0 0.5 "GG"]])
      config {:experiments [{:name "one" :controls ["two"]}
                            {:name "two" :controls []}]}
      config2 {:experiments [{:name "one"}
                             {:name "two" :controls []}]
               :controls ["two"]}]
  (fact "Filter dataset to remove control dominated rows."
    (let [fds (filter-by-controls ds config)
          fds2 (filter-by-controls ds config2)]
      (icore/nrow fds) => 1
      (icore/nrow fds2) => 1
      (icore/sel fds :rows 0) => ["chr1" 20 1.0 0.5 "GG"]
      (icore/sel fds2 :rows 0) => ["chr1" 20 1.0 0.5 "GG"])))

(let [ds (icore/dataset [:chr :pos :one :two :three :seq]
                        [["chr1" 10  0  50  5 "GG"]
                         ["chr1" 20 10 200 20 "CC"]])
      config {:experiments [{:sample "one"}
                            {:sample "two"}
                            {:sample "three"}]}]
  (fact "Filter dataset to use other processed samples as background controls"
    (let [fds (filter-by-control-samples ds config)]
      (icore/sel fds :rows 0) => ["chr1" 10 0.0 50 0.0 "GG"]
      (icore/sel fds :rows 1) => ["chr1" 20 0.0 200 20 "CC"])))