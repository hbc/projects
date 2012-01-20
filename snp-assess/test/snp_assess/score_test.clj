(ns snp-assess.score-test
  (:use [midje.sweet]
        [snp-assess.score]
        [snp-assess.config :only [default-config]]))

(facts "Minority variant detection and analysis."
  (minority-variants
   {"A" 10 "C" 1  "G" 19 "T" 300} default-config) => [["G" 19/330] ["A" 1/33]]
  (minority-variants
   {"A" 10 "C" 5  "G" 15 "T" 20} default-config) => [["G" 3/10] ["A" 1/5] ["C" 1/10]])

(facts "Retrieve minority frequencies from list of bases"
  (minor-target-freq [["G"] ["G"] ["C"]] default-config) => [["C" 1/3]]
  (minor-target-freq (concat (repeat 300 ["A"]) (repeat 5 ["G"]) (repeat 1 ["C"]))
                     default-config) => [["G" 5/306]]
  (minor-target-freq (concat (repeat 50 ["A"]) (repeat 5 ["G"]) (repeat 5 ["C"]))
                     default-config) => [["G" 1/12] ["C" 1/12]])

(facts "Test for presence of a variant"
  (has-variant? "G" nil (concat (repeat 500 ["A"]) (repeat 5 ["G"]) (repeat 1 ["C"]))
                default-config) => true
  (has-variant? "G" 5/480 (concat (repeat 500 ["A"]) (repeat 5 ["G"]) (repeat 1 ["C"]))
                default-config) => true
  (has-variant? "G" 5/530 (concat (repeat 500 ["A"]) (repeat 5 ["G"]) (repeat 1 ["C"]))
                default-config) => true
  (has-variant? "G" 1/1000 (concat (repeat 500 ["A"]) (repeat 5 ["G"]) (repeat 1 ["C"]))
                default-config) => false
  (has-variant? "G" nil [["G"] ["G"] ["C"]] default-config) => false)

(facts "Sample for minimum coverage by random removal."
  (random-min-coverage "G" nil [["G"] ["G"] ["C"]] default-config) => nil
  (random-min-coverage "G" 5/500 (concat (repeat 500 ["A"]) (repeat 5 ["G"]) (repeat 1 ["C"]))
                       default-config) => (roughly 400 300))
