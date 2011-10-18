(ns hbc.transposon.test.core
  (:use [hbc.transposon.core]
        [midje.sweet])
  (:require [clojure.java.io :as io]
            [clojure.string :as string])
  (:import [java.io StringReader]))

(fact "Read locations from custom tab-delimited file."
  (let [data (read-custom-positions "")]
    (first data) => {:strand "+" :space "chr1" :start 98328731 :count 32
                     :seq "TAGGTCTAGAGAAGTAGCTGAGTGTTCTATATCCAGATCATCAG"}
    (count data) => 2)
  (against-background
    (io/reader anything) => (StringReader. (string/join "\n"
    ["HWI-ST782:118:C07NCACXX:7:2308:6744:99171 1:N:0:TAGCTT\t+\tchr1\t98328731\tTAGGTCTAGAGAAGTAGCTGAGTGTTCTATATCCAGATCATCAG\tEAH4.77?ECC>@C?;;@((66.5;B@CCA>>3;,55:?5>@A5\t0\t32\t98328731\tTAGGTCTAGAGAAGTAGCTGAGTGTTCTATATCCAGATCATCAG"
     "HWI-ST782:118:C07NCACXX:7:1208:5008:100804 1:N:0:TAGCTT\t+\tchr1\t95519114\tAGGGTACCATTACCGTCCCTGGACATGCACGTCAGCAACACCTC\tFHGGHHIJJJJJJJHHHFFFEDEEEDDDDDDDDDDDDDDDDDDD\t0\t71\t95519114\tAGGGTACCATTACCGTCCCTGGACATGCACGTCAGCAACACCTC"]))))
