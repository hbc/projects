;; Create classifier for evaluating reads associated with
;; low percentage variations. Build a linear classifier that
;; distinguishes real reads from false positives based on quality
;; score, map score and k-mer frequency.

(ns snp-assess.classify
  (:use [clojure.java.io]
        [clj-ml.utils :only [serialize-to-file deserialize-from-file]]
        [clj-ml.data :only [make-dataset dataset-set-class make-instance]]
        [clj-ml.classifiers :only [make-classifier classifier-train
                                   classifier-evaluate classifier-classify]]
        [snp-assess.core :only [parse-snpdata-line]]
        [snp-assess.score :only [minority-variants naive-read-passes?
                                 normalize-params roughly-freq?]]
        [snp-assess.off-target :only [parse-pos-line]]
        [snp-assess.config :only [default-config]])
  (:require [fs]))

;; Prepare clasifier data: list of normalized parameters (quality,
;; kmer and mapping scores) and naive classifier based on simple
;; linear combination

(defn read-vrn-pos [pos-file max-pct]
  "Prepare set of chromosome, position base for expected variations."
  (letfn [(pos-map [xs]
            (reduce #(assoc %1 (take 3 %2) (last %2)) {} xs))]
    (with-open [rdr (reader pos-file)]
      (->> rdr
           line-seq
           (map parse-pos-line)
           (filter #(<= (last %) max-pct))
           pos-map))))

(defn vrn-data-plus-config [read-data config]
  "Retrieve parameter data for classification with config appended"
  (conj (vec (take-last 3 read-data)) config))

(defn- minority-vrns-from-raw [pos-data config]
  "Retrieve minority variants at a position given raw data."
  (letfn [(count-bases [[id xs]] [(last id) (count xs)])
          (good-read? [read] (apply naive-read-passes? (vrn-data-plus-config read config)))
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

(defn finalize-raw-data [raw-data class config]
  "Prepare normalized raw data for input into classifiers."
  (let [normal-data (apply normalize-params (vrn-data-plus-config raw-data config))]
    (conj normal-data class)))

(defn data-from-pos [pos-data positives config]
  "Retrieve classification data at a particular read position."
  (let [want-bases (minority-vrns-from-raw pos-data config)
        want-keys (filter #(contains? want-bases (last %)) (keys pos-data))
        num-vals (count (vrn-data-plus-config (first (vals pos-data)) config))]
    (for [cur-id want-keys]
      (let [class (if (contains? positives cur-id) 1 0)]
        (map #(finalize-raw-data % class config) (get pos-data cur-id))))))

(defn raw-reads-by-pos [rdr config]
  "Lazy generator of read statistics at each position in data file."
  (->> rdr
       line-seq
       (map parse-snpdata-line)
       (partition-by #(take 2 %))))

(defn- random-sample-negatives [data]
  "Randomly sample negative examples to make total positives. Since there are
   a larger number of negative examples, this prevents learning to just call all
   negatives in the classifier."
  (let [positives (filter #(== 1 (last %)) data)
        negatives (take (count positives) (shuffle (filter #(== 0 (last %)) data)))]
    (concat positives negatives)))

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
           flatten
           (partition 4)
           random-sample-negatives
           vec))))

;; Do the work of classification, with a prepared set of data inputs

(defn- get-dataset [data]
  "Weka dataset ready for classification or training."
  (let [header [:qual :kmer :map-score :c]]
    (make-dataset "ds" header data {:class :c})))

(defn train-classifier [data-file pos-file config]
  "Manage retrieving data and training the classifier."
  (let [class-data (prep-classifier-data data-file pos-file config)
        ds (get-dataset class-data)
        c (-> (make-classifier :regression :linear) (classifier-train ds))]
    c))

(defn prepare-classifier [data-file pos-file work-dir config]
  "High level work to get classifier, included serialization to a file."
  (let [out-dir (fs/join work-dir "classifier")
        classifier-file (fs/join out-dir "build.bin")]
    (if-not (fs/exists? out-dir)
      (fs/mkdirs out-dir))
    (if (fs/exists? classifier-file)
      (deserialize-from-file classifier-file)
      (let [classifier  (train-classifier data-file pos-file config)]
        (serialize-to-file classifier classifier-file)
        classifier))))

;; Call variants given a trained read classifier:
;; - Walk through each position in the input file.
;;    - Filter reads based on the pre-built linear classifier
;;    - Retrieve percentage of called at a position + classification
;;    - Check known variants and determine if correct,

(defn call-vrns-at-pos [reads score-calc config]
  "With read parameters at a position, get map of bases and percent present."
  (letfn [(sum-by-base [xs]
            (reduce (fn [m [b s]] (assoc m b (+ s (get m b 0.0)))) {} xs))
          (percents [xs]
            (let [total (apply + (vals xs))]
              (reduce (fn [m [k v]] (assoc m k (* 100.0 (/ v total)))) {} xs)))
          (remove-low [min-freq orig]
            (let [min-freq-percent (* 100.0 min-freq)
                  want (->> orig
                            (filter #(> (second %) min-freq-percent))
                            (map first))]
              (select-keys orig want)))]
      (let [position (take 2 (first reads))
            base-calls (->> reads
                            (map #(conj % (apply score-calc (drop 3 %))))
                            (map #(list (nth % 2) (last %)))
                            sum-by-base
                            percents
                            (remove-low (:min-freq config)))]
        [position base-calls])))

(defn classifier-checker [classifier config]
  "Calculate probability of read inclusion using pre-built classifier."
  (let [ds (get-dataset [])]
    (fn [qual kmer-pct map-score]
      (let [[nq nk nm] (normalize-params qual kmer-pct map-score config)]
        (classifier-classify classifier (make-instance ds {:qual nq :kmer nk
                                                           :map-score nm :c -1}))))))

(defn assign-position-type [pos base-counts known-vrns config]
  "Given read calls, determine if true/false and positive/negative, type."
  (letfn [(annotate-base [pos base freq known-vrns]
            {:pos pos :base base :freq freq
             :exp-freq (get known-vrns [(first pos) (second pos) base])})
          (minority-base [bases]
            {:pre (== 1 (count (filter #(not (nil? (:exp-freq %))) bases)))}
            (first (filter #(not (nil? (:exp-freq %))) bases)))
          (minority-matches? [bases config]
            (let [base (minority-base bases)]
              (roughly-freq? (:freq base) (:exp-freq base) config)))
          (finalize [bases e-base-count call]
            (let [freq (if (== 0 e-base-count) 100.0
                           (:exp-freq (minority-base bases)))]
              {:class call :freq freq :calls base-counts}))]
    (let [all-bases ["A" "C" "G" "T"]
          ready-bases (map (fn [b] (annotate-base pos b (get base-counts b) known-vrns))
                           all-bases)
          e-bases (map :base (filter #(not (nil? (:exp-freq %))) ready-bases))
          c-bases (map :base (filter #(not (nil? (:freq %))) ready-bases))]
      (finalize ready-bases (count e-bases)
                (case (count e-bases)
                      0 (cond
                         (== (count c-bases) 1) :true-positive
                         (> (count c-bases) 1) :false-positive
                         (== (count c-bases) 0) :false-negative)
                      1 (cond
                         (== (count c-bases) 1) :false-negative
                         (> (count c-bases) 2) :false-negative
                         (minority-matches? ready-bases config) :true-positive
                         :else :false-negative))))))

(defn assess-classifier [data-file pos-file classifier config]
  "Determine rates of true/false positive/negatives with trained classifier"
  (let [classifier-score (classifier-checker classifier config)
        known-vrns (read-vrn-pos pos-file 101.0)]
    (with-open [rdr (reader data-file)]
      (->> (raw-reads-by-pos rdr config)
           (map (fn [xs] (call-vrns-at-pos xs classifier-score config)))
           (map (fn [[pos bases]] (assign-position-type pos bases known-vrns config)))
           (map println)
           vec))))

(defn -main [data-file pos-file work-dir]
  (let [config default-config
        c (prepare-classifier data-file pos-file work-dir config)]
    (println c)
    (println
     (assess-classifier data-file pos-file c config))))
