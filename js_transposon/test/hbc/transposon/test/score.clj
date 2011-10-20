;; Tests for normalizing and scoring combined position counts

(ns hbc.transposon.test.score
  (:use [hbc.transposon.score]
        [midje.sweet])
  (:require [clojure.string :as string]
            [incanter.core :as icore])
  (:import [java.io StringReader]))

(fact "Normalize a CSV file to reads per million counts."
  (let [ncounts (normalize-counts "test.csv" :base 100)]
    (icore/sel ncounts :cols :one) => [0.0 100.0]
    (icore/sel ncounts :cols :two) => [20.0 80.0]
    (icore/col-names ncounts) => [:chr :pos :one :two :three :seq])
  (against-background
    (icore/get-input-reader "test.csv") => (StringReader. (string/join "\n"
      ["chr,pos,one,two,three,seq"
       "chr1,10,0,50,10,GG"
       "chr1,20,10,200,10,CC"]))))
