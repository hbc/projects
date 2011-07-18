(ns snp-assess.score-test
  (:use [midje.sweet]
        [snp-assess.score]
        [snp-assess.config :only [default-config]]))

(facts "minority-variants"
  (minority-variants
   {"A" 12 "C" 10 "G" 15 "T" 200} default-config) => [["A" 0.5] ["G" 1.0]])
