;; Adjustable filtering and analysis paramters
(ns snp-assess.config)

(def default-config
  {:kmer-range [1e-5 0.10]
   :qual-range [4.0 35.0]
   :map-score-range [0.0 250.0]
   :random-coverage-step 25
   :random-coverage-sample 10
   :min-score 1.1
   :min-freq 0.0035})
