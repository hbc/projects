;; Tests for normalizing and scoring combined position counts

(ns hbc.transposon.test.score
  (:use [hbc.transposon.score]
        [midje.sweet])
  (:require [clojure.string :as string]
            [incanter.core :as icore])
  (:import [java.io StringReader]))

(let [ds (icore/dataset [:chr :pos :one :two :three :seq]
                        [["chr1" 10  0  50 10 "GG"]
                         ["chr1" 20 10 200 10 "CC"]])]
  (fact "Normalize a CSV file by reads and experiment counts."
    (let [ncounts (normalize-counts ds :base 100)
          n2counts (normalize-pos-ratios ncounts)]
      (icore/sel ncounts :cols :one) => [0.0 100.0]
      (icore/sel ncounts :cols :two) => [20.0 80.0]
      (icore/col-names ncounts) => [:chr :pos :one :two :three :seq]
      (icore/sel n2counts :cols :one) => [0.0 1.0]
      (icore/sel n2counts :cols :two) => [0.4 0.8]
      (icore/col-names n2counts) => [:chr :pos :one :two :three :seq]
      (fact "Calculate dataset statistics"
        (let [stats (summarize-count-statistics n2counts)]
          (count stats) => 3
          (first stats) => (contains [:one 0.5 0.5]))))))
