;; Adjustable filtering and analysis paramters
(ns snp-assess.config)

(def default-config
  {:kmer-range [1e-5 0.10]
   :qual-range [4.0 35.0]
   :map-score-range [0.0 250.0]})
