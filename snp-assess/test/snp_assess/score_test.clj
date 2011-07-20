(ns snp-assess.score-test
  (:use [midje.sweet]
        [snp-assess.score]
        [snp-assess.config :only [default-config]]))

(facts "Minority variant detection and analysis."
  (minority-variants
   {"A" 10 "C" 1  "G" 19 "T" 20} default-config) => [["G" 19/50] ["A" 1/5]]
  (minority-variants
   {"A" 10 "C" 5  "G" 15 "T" 20} default-config) => [["G" 3/10] ["A" 1/5] ["C" 1/10]])

;.;. FAIL at (NO_SOURCE_FILE:1)
;.;.     Expected: [["C" 1/3]]
;.;.       Actual: #<workflow$buffer$fn__368 cascalog.workflow$buffer$fn__368@5043cc83>
(facts "Retrieve minority frequencies from list of bases"
  (minor-target-freq [["G"] ["G"] ["C"]]) => [["C" 1/3]])
