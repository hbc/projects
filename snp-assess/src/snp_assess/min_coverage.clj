;; Determine distribution of minimum coverage required to detect
;; variations at different frequencies. This helps estimate reads
;; required to effectively detect minority variations.

(ns snp-assess.min-coverage
  (:use [clojure.string :only [split]]
        [cascalog.api]
        [incanter.core :only [save]]
        [incanter.charts :only [xy-plot add-lines]]
        [snp-assess.config :only [default-config]]
        [snp-assess.score :only [min-coverage-cascalog read-filter-cascalog
                                 histogram-bins]]
        [snp-assess.core :only [snpdata-from-hfs]]
        [snp-assess.off-target :only [pos-from-hfs]])
  (:require [cascalog [ops :as ops]])
  (:gen-class))

(defn min-coverage [snpdata var-positions min-coverage-fn filter-fn]
  (??<- [?chr ?pos ?var-base ?var-freq ?coverage]
        (snpdata ?chr ?pos ?base ?qual ?kmer-pct ?map-score)
        (var-positions ?chr ?pos ?var-base ?var-freq)
        ;(filter-fn ?kmer-pct ?qual ?map-score)
        (min-coverage-fn ?var-base ?base :> ?coverage)))

(defn min-coverage-plots [data-dir pos-dir]
  (let [freq-cov (min-coverage (snpdata-from-hfs data-dir)
                               (pos-from-hfs pos-dir)
                               (min-coverage-cascalog default-config)
                               (read-filter-cascalog default-config))
        cov-by-freq (reduce (fn [cov-map [freq cov]]
                              (assoc cov-map freq (cons cov (get cov-map freq))))
                            {} (map #(drop 3 %) freq-cov))]
    (let [num-bins 100.0
          hist-info (sort-by first
                             (map (fn [[ x cs]] (cons x (histogram-bins cs num-bins)))
                                  cov-by-freq))
          [freq x y] (first hist-info)
          plot (xy-plot x y :series-label freq
                        :legend true :title "Minimum detection coverage"
                        :x-label "Passing reads" :y-label "")]
        (doseq [[freq x y] (rest hist-info)]
          (add-lines plot x y :series-label freq))
        (doto plot
          (save "min-coverage-frequencies.png"))
      (println hist-info))))

(defn -main [data-dir pos-dir]
  (min-coverage-plots data-dir pos-dir))
