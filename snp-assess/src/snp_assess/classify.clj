;; Create classifier for evaluating reads associated with
;; low percentage variations. Build a linear classifier that
;; distinguishes real reads from false positives based on quality
;; score, map score and k-mer frequency.

(ns snp-assess.classify
  (:use [clojure.java.io]
        [clj-ml.utils :only [serialize-to-file deserialize-from-file]]
        [snp-assess.core :only [parse-snpdata-line]]
        [snp-assess.off-target :only [parse-pos-line]])
  (:require [fs]))

(defn read-vrn-pos [pos-file max-pct]
  "Prepare set of chromosome, position base for expected variations."
  (with-open [rdr (reader pos-file)]
    (->> rdr
         line-seq
         (map parse-pos-line)
         (filter #(<= (last %) max-pct))
         (map #(take 3 %))
         set)))

(defn build-classifier [data-file pos-file config]
  "Build classifier for reads based on variant/non-variant positions"
  (let [positives (read-vrn-pos pos-file (:max-pct config))]))

(defn prepare-classifier [data-file pos-file work-dir config]
  "High level work to get classifier, included serialization to a file."
  (let [out-dir (fs/join work-dir "classifier")
        classifier-file (fs/join out-dir "")]
    (if-not (fs/exists? out-dir)
      (fs/mkdirs out-dir))
    (if (fs/exists? classifier-file)
      (deserialize-from-file classifier-file)
      (let [classifier (build-classifier data-file pos-file config)]
        (serialize-to-file classifier classifier-file)
        classifier))))

(defn -main [data-file pos-file work-dir]
  (let [config {:max-pct 5.0}]
    (prepare-classifier data-file pos-file work-dir config)))
