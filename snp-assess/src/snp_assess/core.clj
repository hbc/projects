;; Assess variation in experiments with very deep sequencing
;; coverage, to help in setting parameters for variation calling.

(ns snp-assess.core
  (:use [clojure.java.io]
        [clojure.string :only [split]]
        [cascalog.api]
        [snp-assess.score :only [score-calc-cascalog]])
  (:require [cascalog [ops :as ops]]
            [clj-yaml.core :as yaml]))

;; Summary statistics for positions of interest

(defn calc-snpdata-stats [snpdata targets score-calc-fn]
  (??<- [?chr ?pos ?base ?count ?avg-score ?avg-kmer-pct ?avg-qual ?avg-map ?type]
        (snpdata ?chr ?pos ?base ?qual ?kmer-pct ?map-score)
        (targets ?chr ?pos ?base ?type)
        (ops/count ?count)
        (ops/avg ?kmer-pct :> ?avg-kmer-pct)
        (ops/avg ?qual :> ?avg-qual)
        (ops/avg ?map-score :> ?avg-map)
        (score-calc-fn ?kmer-pct ?qual ?map-score :> ?score)
        (ops/avg ?score :> ?avg-score)))

;; Parsing variation data files

(defn parse-snpdata-line-tuple [line]
  (let [parts (split line #"\t")
        [space pos base num qual kmer-pct map-score]
        (condp = (count parts)
          6 (concat (take 3 parts) ["1"] (drop 3 parts))
          7 parts)]
    [space (Integer/parseInt pos) base (Integer/parseInt num) (Integer/parseInt qual)
     (Float/parseFloat kmer-pct) (Integer/parseInt map-score)]))

(defn parse-snpdata-line [line]
  (zipmap [:space :pos :base :num :qual :kmer-pct :map-score]
          (parse-snpdata-line-tuple line)))

(defn snpdata-from-hfs [dir]
  "Parse a directory of variation data files from HDFS"
  (let [source (hfs-textline dir)]
    (<- [?chr ?pos ?base ?qual ?kmer-pct ?map-score]
        (source ?line)
        (parse-snpdata-line-tuple ?line :> ?chr ?pos ?base ?num ?qual ?kmer-pct ?map-score))))

;; Parsing target files of positions to query

(defn parse-pos-line [line]
  (let [[space pos base type] (take 4 (split line #"\t"))]
    [space (Integer/parseInt pos) base type]))

(defn target-pos-from-hfs [dir]
  "Parse a directory of chromosome position targets to be assessed."
  (let [source (hfs-textline dir)]
    (<- [?chr ?pos ?base ?type]
        (source ?line)
        (parse-pos-line ?line :> ?chr ?pos ?base ?type))))

(defn target-snpdata-stats [data-dir pos-dir config]
  "Output individual base call statistics at specified positions."
  (let [out (calc-snpdata-stats (snpdata-from-hfs data-dir)
                                (target-pos-from-hfs pos-dir)
                                (score-calc-cascalog config))
        sort-out (sort-by #(vec (map % [0 1 2])) out)]
    (doseq [row sort-out]
      (println (apply format "| %s | %s | %s | %5s | %.1f | %.1e | %.1f | %.1f | %10s |" row)))))

(defn load-config [in-file]
  (-> in-file slurp yaml/parse-string :params))

(defn -main [data-dir pos-dir config-file]
  (target-snpdata-stats data-dir pos-dir (load-config config-file)))
