;; Examine the range of off-target variations in deep sequencing
;; experiments. The goal is to estimate the lowest level of possible
;; detection by examining the distribution of background variants.

(ns snp-assess.off-target
  (:use [cascalog.api]
        [clojure.string :only [split]]
        [snp-assess.config :only [default-config]]
        [snp-assess.score :only [minor-target-cascalog read-filter-cascalog]]
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
    (println (flatten freq))
    (println (flatten filter-freq))))

(defn -main [data-dir pos-dir]
  (off-target-plots data-dir pos-dir))
