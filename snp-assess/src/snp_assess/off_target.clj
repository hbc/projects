;; Examine the range of off-target variations in deep sequencing
;; experiments. The goal is to estimate the lowest level of possible
;; detection by examining the distribution of background variants.

(ns snp-assess.off-target
  (:use [clojure.string :only [split]]
        [cascalog.api]
        [incanter.core :only [save]]
        [incanter.charts :only [xy-plot add-lines]]
        [snp-assess.config :only [default-config]]
        [snp-assess.score :only [minor-target-cascalog read-filter-cascalog
                                 histogram-bins]]
        [snp-assess.core :only [snpdata-from-hfs]])
  (:require [cascalog [ops :as ops]])
  (:gen-class))

(defn off-target-freqs [snpdata no-var-positions minority-freq-fn]
  (??<- [?freq]
        (snpdata ?chr ?pos ?base _ _ _)
        (no-var-positions ?chr ?pos _ _)
        (minority-freq-fn ?base :> ?freq)))

(defn off-target-freqs-filter [snpdata no-var-positions minority-freq-fn filter-fn]
  (??<- [?freq]
        (snpdata ?chr ?pos ?base ?qual ?kmer-pct ?map-score)
        (no-var-positions ?chr ?pos _ _)
        (filter-fn ?kmer-pct ?qual ?map-score)
        (minority-freq-fn ?base :> ?freq)))

(defmapop parse-pos-line [line]
  (let [[space pos base freq] (take 4 (split line #"\t"))]
    [space (Integer/parseInt pos) base (Float/parseFloat freq)]))

(defn pos-from-hfs [dir]
  "Retrieve space, pos, expected base and frequency from HDFS files"
  (let [source (hfs-textline dir)]
    (<- [?chr ?pos ?base ?freq]
        (source ?line)
        (parse-pos-line ?line :> ?chr ?pos ?base ?freq))))

(defn off-target-plots [data-dir pos-dir]
  (let [freq (off-target-freqs (snpdata-from-hfs data-dir)
                               (pos-from-hfs pos-dir)
                               (minor-target-cascalog default-config))
        filter-freq (off-target-freqs-filter (snpdata-from-hfs data-dir)
                                             (pos-from-hfs pos-dir)
                                             (minor-target-cascalog default-config)
                                             (read-filter-cascalog default-config))]
    (let [num-bins 100.0
          [freq-x freq-hist] (histogram-bins freq num-bins)
          [freq-filter-x freq-filter-hist] (histogram-bins filter-freq num-bins)]
      (doto (xy-plot freq-x freq-hist :series-label "raw"
                     :legend true :title "Off-target"
                     :title "Off-target" :x-label "Frequency" :y-label "")
        (add-lines freq-filter-x freq-filter-hist :series-label "filtered")
        (save "off-target-frequencies.png")))
    (println (histogram-bins freq 100.0))
    (println (histogram-bins filter-freq 100.0))))

(defn -main [data-dir pos-dir]
  (off-target-plots data-dir pos-dir))
