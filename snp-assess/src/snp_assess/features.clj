(ns snp-assess.features
  "Create features from a set of metrics for input to classification algorithms.
  Handles normalization for linear regression analysis and feature expansion for
  random forest separation.")

;; ## Linear regression

(defn min-max-norm
  "Perform min-max normalization, truncated larger values at max and min."
  [score [minv maxv]]
  (let [trunc-score-max (if (< score maxv) score maxv)
        trunc-score (if (> trunc-score-max minv) trunc-score-max minv)]
    (/ (- trunc-score minv) (- maxv minv))))

(defn normalize-params
  "Prep quality score, k-mer and map score with min/max normalization."
  [qual kmer-pct map-score config]
  [(min-max-norm qual (:qual-range config))
   (min-max-norm kmer-pct (:kmer-range config))
   (min-max-norm map-score (:map-score-range config))])

;; ## Random forest
