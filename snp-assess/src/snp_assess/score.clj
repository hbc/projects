;; Score and assess reads

(ns snp-assess.score
  (:use [cascalog.api]))

(defn min-max-norm [score [minv maxv]]
  "Perform min-max normalization, truncated larger values at max and min."
  (let [trunc-score-max (if (< score maxv) score maxv)
        trunc-score (if (> trunc-score-max minv) trunc-score-max minv)]
    (/ (- trunc-score minv) (- maxv minv))))

(defn score-calc [kmer-pct qual map-score config]
  "Calculate read score given kmer, quality and mapping scores."
  (+ (min-max-norm kmer-pct (:kmer-range config))
     (min-max-norm qual (:qual-range config))
     (min-max-norm map-score (:map-score-range config))))

(defn read-passes? [kmer-pct qual map-score config]
  (>= (score-calc kmer-pct qual map-score config)
      (:min-score config)))

(defn minority-variants [base-counts config]
  "List of minority variant bases and frequencies.")

(defn score-calc-cascalog [config]
  "Prepare cascalog ready function for calculating scores."
  (defmapop score-calc [kmer-pct qual map-score]
    (score-calc kmer-pct qual map-score config)))

(defn read-filter-cascalog [config]
  "Prepare cascalog ready function to filter reads with configuration."
  (deffilterop read-filter [kmer-pct qual map-score]
    (read-passes? kmer-pct qual map-score config)))

(defbufferop minor-target-freq [read-scores]
  "Return frequency of minority variants at the position.
   read-score is: base, score.")
