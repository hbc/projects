;; Determine distribution of minimum coverage required to detect
;; variations at different frequencies. This helps estimate reads
;; required to effectively detect minority variations.

(ns snp-assess.min-coverage
  (:use [clojure.string :only [split]]
        [cascalog.api]
        [incanter.core :only [save dataset]]
        [incanter.pdf :only [save-pdf]]
        [incanter.io]
        [incanter.charts :only [xy-plot add-lines box-plot add-box-plot
                                scatter-plot add-points set-x-label
                                set-y-label set-title]]
        [snp-assess.score :only [min-coverage-cascalog read-filter-cascalog
                                 histogram-bins]]
        [snp-assess.core :only [snpdata-from-hfs load-config]]
        [snp-assess.off-target :only [pos-from-hfs]]
        [snp-assess.classify :only [raw-data-frequencies]])
  (:require [cascalog [ops :as ops]]
            [fs.core :as fs])
  (:gen-class))

;; Coverage for minority variant detection

(defn min-coverage [snpdata var-positions min-coverage-fn filter-fn]
  "Query for the minimum coverage necessary to detect minority variant."
  (??<- [?chr ?pos ?var-base ?var-freq ?coverage]
        (snpdata ?chr ?pos ?base ?qual ?kmer-pct ?map-score)
        (var-positions ?chr ?pos ?var-base ?var-freq)
        (filter-fn ?kmer-pct ?qual ?map-score)
        (min-coverage-fn ?var-base ?var-freq ?base :> ?coverage)))

(defn min-coverage-plot [data-dir pos-dir image-dir config]
  "Plot minimum coverage for detecting minority variants."
  (let [freq-cov (min-coverage (snpdata-from-hfs data-dir)
                               (pos-from-hfs pos-dir)
                               (min-coverage-cascalog config)
                               (read-filter-cascalog config))
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
    (save plot (str (fs/file image-dir "min-coverage-frequencies.png")))
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

(defn coverage-dist-plot [data-dir image-dir config]
  (letfn [(counts-only [xs] (map last xs))]
    (let [cov (coverage-dist (snpdata-from-hfs data-dir))
          cov-filter (filtered-coverage-dist (snpdata-from-hfs data-dir)
                                             (read-filter-cascalog config))
          total-reads (/ (apply + (counts-only cov)) 1e6)
          num-bins 100.0
          [cov-x cov-y] (histogram-bins (counts-only cov) num-bins)
          [f-cov-x f-cov-y] (histogram-bins (counts-only cov-filter) num-bins)]
      (doto (xy-plot cov-x cov-y :series-label "raw" :legend true
                     :title (format "Coverage: %.1f million reads" total-reads)
                     :x-label "Coverage" :y-label "")
        (add-lines f-cov-x f-cov-y :series-label "filter")
        (save (str (fs/file image-dir "coverage-distribution.png"))))
      (doto (scatter-plot (map second cov) (counts-only cov) :series-label "raw" :legend true
                          :title "Coverage by position"
                          :x-label "Position" :y-label "Coverage")
        (add-points (map second cov-filter) (counts-only cov-filter) :series-label "filter")
        (save (str(fs/file image-dir "coverage-by-pos.png")))))))

(defn coverage-by-pos-plot [title data-file image-dir config]
  "Plot of read coverage by position."
  (let [raw-freqs (raw-data-frequencies data-file (assoc config :min-freq 0.0))
        cov-x (map #(-> % :position second) raw-freqs)
        cov-y (map #(/ (:total %) 1e6) raw-freqs)]
    (println "Ready to plot")
    (doto (xy-plot cov-x cov-y)
      (set-x-label "Position")
      (set-y-label "Coverage (million reads)")
      (set-title (str (if (get config :unique-only false)
                                 "Unique reads: "
                                 "Total coverage: ")
                               title))
      (save-pdf (str (fs/file image-dir (format "%s-coverage.pdf" (fs/base-name data-file true))))))
    {:x cov-x :y cov-y}))

(defn write-raw-coverage [cov out-dir data-files]
  (let [out-file (str (fs/file out-dir (format "%s-coverage.csv"
                                               (fs/base-name (first data-files) true))))
        ds (dataset (cons "pos" (map #(fs/base-name % true) data-files))
                    (map (fn [i] (cons (-> cov first :x (nth i))
                                       (map #(-> % :y (nth i)) cov)))
                         (range (-> cov first :x count))))]
    (save ds out-file)))

(defn -main [work-dir config-file & data-files]
  (let [image-dir (str (fs/file work-dir "images"))
        config (load-config config-file)]
    (if-not (fs/exists? image-dir)
      (fs/mkdirs image-dir))
    (let [data-files (partition 2 data-files)
          cov (map (fn [[title fname]]
                     (coverage-by-pos-plot title fname image-dir config))
                   data-files)]
      (write-raw-coverage cov image-dir (map second data-files)))
    ;(coverage-dist-plot data-dir image-dir config)
    ;(min-coverage-plot data-dir pos-dir image-dir config)
    ))
