;; Create classifier for evaluating reads associated with
;; low percentage variations. Build a linear classifier that
;; distinguishes real reads from false positives based on quality
;; score, map score and k-mer frequency.

(ns snp-assess.classify
  (:use [clj-ml.utils :only [serialize-to-file deserialize-from-file]])
  (:require [fs]))

(defn build-classifier [data-dir pos-dir]
  "Build classifier for reads based on variant/non-variant positions")

(defn prepare-classifier [data-dir pos-dir work-dir]
  "High level work to get classifier, included serialization to a file."
  (let [out-dir (fs/join work-dir "classifier")
        classifier-file (fs/join out-dir "")]
    (if-not (fs/exists? out-dir)
      (fs/mkdirs out-dir))
    (if (fs/exists? classifier-file)
      (deserialize-from-file classifier-file)
      (let [classifier (build-classifier data-dir pos-dir)]
        (serialize-to-file classifier classifier-file)
        classifier))))

(defn -main [data-dir pos-dir work-dir]
  (prepare-classifier data-dir pos-dir work-dir))
