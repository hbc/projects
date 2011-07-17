;; Score and assess reads

(ns snp-assess.score
  (:use [cascalog.api]))

(defn min-max-norm [score [minv maxv]]
  "Perform min-max normalization, truncated larger values at max and min."
  (let [trunc-score-max (if (< score maxv) score maxv)
        trunc-score (if (> trunc-score-max minv) trunc-score-max minv)]
    (/ (- trunc-score minv) (- maxv minv))))

(defn combine-score-w-config [config]
  (defmapop combine-score [kmer-pct qual map-score]
    (+ (min-max-norm kmer-pct (:kmer-range config))
       (min-max-norm qual (:qual-range config))
       (min-max-norm map-score (:map-score-range config)))))
