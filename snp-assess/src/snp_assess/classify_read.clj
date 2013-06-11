(ns snp-assess.classify-read
  "Classify file of reads, providing output compatible with Topcoder grading"
  (:use [clojure.java.io]
        [clojure.string :only [split]]
        [clj-ml.data :only [make-instance]]
        [clj-ml.classifiers :only [classifier-classify]]
        [snp-assess.core :only [load-config]]
        [snp-assess.features :only [metrics-to-features]]
        [snp-assess.classify :only [prepare-classifier get-dataset]])
  (:require [fs.core :as fs]))

(defn classifier-score [classifier config]
  "Calculate read score using pre-built classifier."
  (let [ds (get-dataset [])]
    (fn [qual kmer-pct map-score]
      (let [[nq nk nm] (metrics-to-features qual kmer-pct map-score
                                            (assoc-in config [:classication :classifier]
                                                      [:regression :linear]))]
        (classifier-classify classifier (make-instance ds {:qual nq :kmer nk
                                                           :map-score nm :c -1}))))))

(defn parse-read-line [line]
  "Extract quality, kmer and map-score from input file line."
  (let [[qual kmer map-score] (take 3 (split line #" "))]
    [(Integer/parseInt qual) (Float/parseFloat kmer) (Integer/parseInt map-score)]))

(defn scale-score [x]
  "Scale output score to match expected integer range (0-1e9) for TopCoder"
  (let [max-scale-score 1e9
        max-score 2.0]
    (int (* (/ (min (max 0.0 x) max-score)
               max-score)
            max-scale-score))))

(defn scores-for-read-file [in-file scorer]
  "Provide a score for each set of read characteristics in the input file."
  (let [out-file (str (fs/file (fs/parent in-file)
                               (str (fs/name in-file) "-scores.txt")))]
    (with-open [rdr (reader in-file)
                wrtr (writer out-file)]
      (->> rdr
           line-seq
           rest
           (map parse-read-line)
           (map #(apply scorer %))
           (map scale-score)
           (map #(.write wrtr (str % "\n")))
           doall))
    out-file))

(defn -main [read-file work-dir config-file]
  (let [config (-> (load-config config-file)
                   (assoc :verbose true))
        c (prepare-classifier nil nil work-dir config)
        scorer (classifier-score c config)]
    (scores-for-read-file read-file scorer)))
