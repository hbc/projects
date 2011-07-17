;; Assess variation in experiments with very deep sequencing
;; coverage, to help in setting parameters for variation calling.

(ns snp-assess.core
  (:use [clojure.string :only [split]]
        [cascalog.api]
        [snp-assess.config :only [default-config]]
        [snp-assess.score :only [combine-score-w-config]])
  (:require [cascalog [ops :as ops]])
  (:gen-class))

;; Summary statistics for positions of interest

(defn calc-snpdata-stats [snpdata targets score-combine-fn]
  (??<- [?chr ?pos ?base ?count ?avg-score ?avg-kmer-pct ?avg-qual ?avg-map ?type]
        (snpdata ?chr ?pos ?base ?qual ?kmer-pct ?map-score)
        (targets ?chr ?pos ?base ?type)
        (ops/count ?count)
        (ops/avg ?kmer-pct :> ?avg-kmer-pct)
        (ops/avg ?qual :> ?avg-qual)
        (ops/avg ?map-score :> ?avg-map)
        (score-combine-fn ?kmer-pct ?qual ?map-score :> ?score)
        (ops/avg ?score :> ?avg-score)))

;; Parsing variation data files

(defmapop parse-snpdata-line [line]
  (let [[space pos base qual kmer-pct map-score] (split line #"\t")]
    [space (Integer/parseInt pos) base (Integer/parseInt qual)
     (Float/parseFloat kmer-pct) (Integer/parseInt map-score)]))

(defn snpdata-from-hfs [dir]
  "Parse a directory of variation data files from HDFS"
  (let [source (hfs-textline dir)]
    (<- [?chr ?pos ?base ?qual ?kmer-pct ?map-score]
        (source ?line)
        (parse-snpdata-line ?line :> ?chr ?pos ?base ?qual ?kmer-pct ?map-score))))

;; Parsing target files of positions to query

(defmapop parse-pos-line [line]
  (let [[space pos base type] (take 4 (split line #"\t"))]
    [space (Integer/parseInt pos) base type]))

(defn target-pos-from-hfs [dir]
  "Parse a directory of chromosome position targets to be assessed."
  (let [source (hfs-textline dir)]
    (<- [?chr ?pos ?base ?type]
        (source ?line)
        (parse-pos-line ?line :> ?chr ?pos ?base ?type))))

(defn target-snpdata-stats [data-dir pos-dir]
  "Output individual base call statistics at specified positions."
  (let [out (calc-snpdata-stats (snpdata-from-hfs data-dir)
                                (target-pos-from-hfs pos-dir)
                                (combine-score-w-config default-config))
        sort-out (sort-by #(vec (map % [0 1 2])) out)]
    (doseq [row sort-out]
      (println (apply format "| %s | %s | %s | %5s | %.1f | %.1e | %.1f | %.1f | %10s |" row)))))

(defn -main [data-dir pos-dir]
  (target-snpdata-stats data-dir pos-dir))
