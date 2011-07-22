;; Score and assess reads

(ns snp-assess.score
  (:use [cascalog.api]
        [clojure.contrib.math]))

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
  "List of minority variant bases and frequencies.
  The highest frequency base, the majority variant, is removed,
  along with any that are below the minimum configured frequency."
  (let [total (apply + (vals base-counts))]
    (->> (interleave (keys base-counts)
                     (map #(/ % total) (vals base-counts)))
         (partition 2)
         (sort-by second >)
         rest
         (filter #(>= (second %) (:min-freq config))))))

(defn minor-target-freq [read-bases config]
  "Return frequency of minority variants at the position.
  read-bases are a list of called bases, returns list of base/frequency"
  (->> read-bases
       (map first)
       frequencies
       (#(minority-variants % config))))

(defn has-variant? [base read-bases config]
  "Determine if the read bases correspond to the expected variant."
  (let [minor-freqs (minor-target-freq read-bases config)]
    (and (>= (count minor-freqs) 1)
         (= base (str (ffirst minor-freqs))))))

(defn random-min-coverage [base read-bases config]
  "Randomly remove reads to determine minimum coverage to detect a variant base."
  (letfn [(remove-random [bases step]
            (drop step (shuffle bases)))]
    (loop [cur-bases read-bases, cur-count nil]
      (if (has-variant? base cur-bases config)
        (recur (remove-random cur-bases (:random-coverage-step config))
               (count cur-bases))
        cur-count))))

;; Create data bins for visualization; adapted from John Lawrence Aspden:
;; https://github.com/johnlawrenceaspden/hobby-code/blob/master/rule_of_succession.clj

(defn ceil-places [number decimals]
  "Get the ceiling at the specified decimal place
  http://stackoverflow.com/questions/5072492/
  how-do-i-trim-the-decimal-of-a-number-using-clojure-or-jython"
  (let [factor (expt 10 decimals)]
    (bigdec (/ (ceil (* factor number)) factor))))

(defn histogram-bins [vals bins & [max-val]]
  (letfn [(intervals [n max] (partition 2 1 (range 0.0 max (/ max n))))
          (between [xs a b] (/ (count (filter #(and (>= % a) (< % b)) xs))
                               (count xs)))]
    (let [sort-vals (sort (flatten vals))
          max (if (nil? max-val) (ceil-places (last sort-vals) 1) max-val)]
      [(map first (intervals bins max))
       (map (fn [[a b]] (between sort-vals a b)) (intervals bins max))])))

;; Cascalog ready generating functions -- need configuration dictionary

(defn score-calc-cascalog [config]
  "Prepare cascalog ready function for calculating scores."
  (defmapop score-calc-c [kmer-pct qual map-score]
    (score-calc kmer-pct qual map-score config)))

(defn read-filter-cascalog [config]
  "Prepare cascalog ready function to filter reads with configuration."
  (deffilterop read-filter [kmer-pct qual map-score]
    (read-passes? kmer-pct qual map-score config)))

(defn minor-target-cascalog [config]
  "Cascalog function for getting minor target frequencies given positional calls."
  (defbufferop minor-target [read-bases]
    (map #(-> % second float)
         (minor-target-freq read-bases config))))

(defn min-coverage-cascalog [config]
  "Cacalog function that calculates minimum coverage to detect a variant"
  (defbufferop min-coverage [info]
    (let [exp-base (ffirst info)
          read-bases (map second info)]
       (repeatedly (:random-coverage-sample config)
                   #(random-min-coverage exp-base read-bases config)))))
