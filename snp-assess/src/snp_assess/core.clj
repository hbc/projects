;; Assess variation in experiments with very deep sequencing
;; coverage, to help in setting parameters for variation calling.

(ns snp-assess.core
  (:use [clojure.string :only [split]]
        [cascalog.api])
  (:require [cascalog [ops :as ops]])
  (:gen-class))

;; Summary statistics for positions of interest

(defn calc-snpdata-stats [snpdata targets]
  (??<- [?chr ?pos ?base ?count ?avg-kmer-pct ?avg-qual ?avg-map ?type]
        (snpdata ?chr ?pos ?base ?qual ?kmer-pct ?map-score)
        (targets ?chr ?pos ?base ?type)
        (ops/count ?count)
        (ops/avg ?kmer-pct :> ?avg-kmer-pct)
        (ops/avg ?qual :> ?avg-qual)
        (ops/avg ?map-score :> ?avg-map)))

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
                                (target-pos-from-hfs pos-dir))
        sort-out (sort-by #(vec (map % [0 1 2])) out)]
    (doseq [row sort-out]
      (println (apply format "| %s | %s | %s | %5s | %.1e | %.1f | %.1f | %10s |" row)))))

(defn -main [data-dir pos-dir]
  (target-snpdata-stats data-dir pos-dir))
