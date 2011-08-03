;; Create classifier for evaluating reads associated with
;; low percentage variations. Build a linear classifier that
;; distinguishes real reads from false positives based on quality
;; score, map score and k-mer frequency.

(ns snp-assess.classify
  (:use [clojure.java.io]
        [clj-ml.utils :only [serialize-to-file deserialize-from-file]]
        [snp-assess.core :only [parse-snpdata-line]]
        [snp-assess.score :only [minority-variants naive-read-passes?]]
        [snp-assess.off-target :only [parse-pos-line]]
        [snp-assess.config :only [default-config]])
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

(defn get-vrn-class [vrn-data positives]
  "Determine variation data class based on position and base.")

(defn- minority-vrns-from-raw [pos-data config]
  "Retrieve minority variants at a position given raw data."
  (letfn [(count-bases [[id xs]] [(last id) (count xs)])
          (good-read? [read] (apply naive-read-passes? (conj (vec (take-last 3 read)) config)))
          (filter-reads [[id xs]] [id (filter good-read? xs)])
          (below-thresh? [[_ freq]] (<= (* freq 100.0)
                                        (-> config :classification :max-pct)))]
    (->> pos-data
         (map filter-reads)
         (map count-bases)
         flatten
         (apply hash-map)
         (#(minority-variants % config))
         (filter below-thresh?)
         (map first)
         set)))

(defn data-from-pos [pos-data positives config]
  "Retrieve classification data at a particular read position."
  (let [want-bases (minority-vrns-from-raw pos-data config)
        want-keys (filter #(contains? want-bases (last %)) (keys pos-data))]
    (doseq [cur-id want-keys]
      (let [class (if (contains? positives cur-id) :pos :neg)]
        (println cur-id (first (get pos-data cur-id)) class)))
    ))

(defn prep-classifier-data [data-file pos-file config]
  "Retrieve classification data based on variant/non-variant positions"
  (let [positives (read-vrn-pos pos-file (-> config :classification :max-pct))]
    (with-open [rdr (reader data-file)]
      (->> rdr
           line-seq
           (map parse-snpdata-line)
           (partition-by #(take 2 %))
           (map (fn [xs] (group-by #(take 3 %) xs)))
           (map #(data-from-pos % positives config))
           vec))))

(defn build-classifier [data-file pos-file config]
  (let [class-data (prep-classifier-data data-file pos-file config)]))

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
  (prepare-classifier data-file pos-file work-dir default-config))
