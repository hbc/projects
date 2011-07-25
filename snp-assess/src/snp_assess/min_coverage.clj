;; Determine distribution of minimum coverage required to detect
;; variations at different frequencies. This helps estimate reads
;; required to effectively detect minority variations.

(ns snp-assess.min-coverage
  (:use [clojure.string :only [split]]
        [cascalog.api]
        [incanter.core :only [save]]
        [incanter.charts :only [xy-plot add-lines box-plot add-box-plot]]
        [snp-assess.config :only [default-config]]
        [snp-assess.score :only [min-coverage-cascalog read-filter-cascalog
                                 histogram-bins]]
        [snp-assess.core :only [snpdata-from-hfs]]
        [snp-assess.off-target :only [pos-from-hfs]])
  (:require [cascalog [ops :as ops]])
  (:gen-class))

;; Coverage for minority variant detection

(defn min-coverage [snpdata var-positions min-coverage-fn filter-fn]
  "Query for the minimum coverage necessary to detect minority variant."
  (??<- [?chr ?pos ?var-base ?var-freq ?coverage]
        (snpdata ?chr ?pos ?base ?qual ?kmer-pct ?map-score)
        (var-positions ?chr ?pos ?var-base ?var-freq)
        (filter-fn ?kmer-pct ?qual ?map-score)
        (min-coverage-fn ?var-base ?var-freq ?base :> ?coverage)))

(defn min-coverage-plot [data-dir pos-dir]
  "Plot minimum coverage for detecting minority variants."
  (let [freq-cov (min-coverage (snpdata-from-hfs data-dir)
                               (pos-from-hfs pos-dir)
                               (min-coverage-cascalog default-config)
                               (read-filter-cascalog default-config))
        cov-by-freq (reduce (fn [cov-map [freq cov]]
                              (assoc cov-map freq (cons cov (get cov-map freq))))
                            {} (map #(drop 3 %) freq-cov))
        want-freqs (set (map str [40.0 33.5 6.5 5.0 1.5 1.0 0.5]))
        plot-info (->> cov-by-freq
                       (filter #(contains? want-freqs (str (first %))))
                       (sort-by first))
        [freq x] (first plot-info)
        plot (box-plot x :series-label freq
                       :legend true :title "Minimum coverage for variant detection"
                       :x-label "" :y-label "Passing reads")]
    (doseq [[freq x] (rest plot-info)]
      (add-box-plot plot x :series-label freq))
    (save plot "min-coverage-frequencies.png")
    (println plot-info)))

;; Overall coverage for an experiment

(defn coverage-dist [snpdata]
  "Distribution of coverage at each position."
  (??<- [?chr ?pos ?count]
        (snpdata ?chr ?pos ?base _ _ _)
        (ops/count ?count)))

(defn filtered-coverage-dist [snpdata filter-fn]
  "Distribution of filtered coverage at each position."
  (??<- [?chr ?pos ?count]
        (snpdata ?chr ?pos ?base ?qual ?kmer-pct ?map-score)
        (filter-fn ?kmer-pct ?qual ?map-score)
        (ops/count ?count)))

(defn coverage-dist-plot [data-dir]
  (letfn [(counts-only [xs] (map last xs))]
    (let [cov (counts-only
               (coverage-dist (snpdata-from-hfs data-dir)))
          cov-filter (counts-only
                      (filtered-coverage-dist (snpdata-from-hfs data-dir)
                                              (read-filter-cascalog default-config)))
          total-reads (/ (apply + cov) 1e6)
          num-bins 100.0
          [cov-x cov-y] (histogram-bins cov num-bins)
          [f-cov-x f-cov-y] (histogram-bins cov-filter num-bins)]
      (doto (xy-plot cov-x cov-y :series-label "raw" :legend true
                     :title (format "Coverage: %.1f million reads" total-reads)
                     :x-label "Coverage" :y-label "")
        (add-lines f-cov-x f-cov-y :series-label "filter")
        (save "coverage-distribution.png")))))

(defn -main [data-dir pos-dir]
  (coverage-dist-plot data-dir)
  (min-coverage-plot data-dir pos-dir))
