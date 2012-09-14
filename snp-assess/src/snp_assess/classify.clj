(ns snp-assess.classify
  "Create classifier for evaluating reads associated with
   low percentage variations. Build a linear classifier that
   distinguishes real reads from false positives based on quality
   score, map score and k-mer frequency."
  (:use [clojure.java.io]
        [clj-ml.utils :only [serialize-to-file deserialize-from-file]]
        [clj-ml.data :only [make-dataset dataset-set-class make-instance]]
        [clj-ml.classifiers :only [make-classifier classifier-train
                                   classifier-evaluate classifier-classify]]
        [snp-assess.core :only [parse-snpdata-line load-config load-run-config
                                parse-pos-line]]
        [snp-assess.score :only [minority-variants naive-read-passes?
                                 roughly-freq?]]
        [snp-assess.features :only [metrics-to-features]]
        [snp-assess.classify-eval :only [write-assessment print-vrn-summary
                                         roc-summarize-assessment]]
        [snp-assess.reference :only [read-vcf-ref]])
  (:require [clojure.string :as string]
            [fs.core :as fs]))

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

(defn read-ref [pos-file ref-file max-pct]
  "Read reference file of variations; handling flat and VCF formats."
  (let [first-line (with-open [rdr (reader pos-file)]
                     (first (line-seq rdr)))]
    (cond
     (= (.indexOf first-line "##fileformat=VCF") 0) (read-vcf-ref pos-file ref-file max-pct)
     :else (read-vrn-pos pos-file max-pct))))

(defn get-vrn-data [read-data]
  "Retrieve parameter data for classification."
  (vec ((juxt :qual :kmer-pct :map-score) read-data)))

(defn- minority-vrns-from-raw [pos-data config]
  "Retrieve minority variants at a position given raw data."
  (letfn [(count-bases [[id xs]] [(last id) (count xs)])
          (good-read? [read] (apply naive-read-passes? (conj (get-vrn-data read)
                                                             config)))
          (filter-reads [[id xs]] [id (filter good-read? xs)])]
    (->> pos-data
         (map filter-reads)
         (map count-bases)
         flatten
         (apply hash-map)
         (#(minority-variants % config))
         (#(apply hash-map (flatten %))))))

(defn finalize-raw-data [raw-data class config]
  "Prepare normalized raw data for input into classifiers."
  (let [normal-data (apply metrics-to-features (conj (get-vrn-data raw-data) config))]
    (conj (vec normal-data) class)))

(defn is-low-minority-vrn? [cur-id class base-freqs config]
  (let [freq (* 100.0 (get base-freqs (last cur-id)))
        min-freq (* 100.0 (:min-freq config))]
    (case class
      (0 :b) (and (<= freq (-> config :classification :max-neg-pct))
                  (>= min-freq))
      (1 :a) (and (<= freq (-> config :classification :max-pos-pct))
                  (>= min-freq)))))

(defn data-from-pos [pos-data positives config]
  "Retrieve classification data at a particular read position."
  (let [want-bases (minority-vrns-from-raw pos-data config)
        want-keys (filter #(contains? want-bases (last %)) (keys pos-data))]
    (reduce (fn [coll cur-id]
              (let [group (case [(contains? positives cur-id)
                                 (get-in config [:classification :group])]
                            [true :category] :a
                            [false :category] :b
                            [true :numerical] 1
                            [false :numerical] 0)]
                (if-not (is-low-minority-vrn? cur-id group want-bases config) coll
                        (concat coll
                                (map #(finalize-raw-data % group config) (get pos-data cur-id))))))
            [] want-keys)))

(defn raw-reads-by-pos [rdr config & {:keys [is-training?]}]
  "Lazy generator of read statistics at each position in data file."
  (let [assess-range (if is-training?
                       (get-in config [:classification :train-bases]
                               (get-in config [:classification :assess-bases]))
                       (get-in config [:classification :assess-bases]))
        raw-out (->> rdr
                     line-seq
                     (map parse-snpdata-line)
                     (partition-by (juxt :space :pos)))]
    (if (nil? assess-range) raw-out
        (->> raw-out
             (drop-while #(< (-> % first :pos) (first assess-range)))
             (take-while #(< (-> % first :pos) (second assess-range)))))))

(defn- random-sample-negatives [config data]
  "Randomly sample negative examples to match total positives. Since there are
   a larger number of negative examples, this prevents learning to just call all
   negatives in the classifier."
  (let [positives (filter #(contains? #{1 :a} (last %)) data)
        negatives (take (count positives) (shuffle (filter #(contains? #{0 :b} (last %)) data)))]
    (concat positives negatives)))

(defn downsample-at-pos
  "Downsample reads to a given percentage of total input reads at a position."
  [pct reads]
  (if (nil? pct)
    reads
    (take (Math/round (* pct (count reads))) (shuffle reads))))

(defn prep-classifier-data [data-file pos-file ref-file config & {:keys [downsample]}]
  "Retrieve classification data based on variant/non-variant positions"
  (let [positives (read-ref pos-file ref-file
                            (-> config :classification :max-pos-pct))]
    (with-open [rdr (reader data-file)]
      (->> (raw-reads-by-pos rdr config :is-training? true)
           (map (partial downsample-at-pos downsample))
           (map (fn [xs] (group-by (juxt :space :pos :base) xs)))
           (map #(data-from-pos % positives config))
           (apply concat)
           (#(if-not (get-in config [:classification :random-sample] true) %
                     (random-sample-negatives config %)))
           vec))))

;; Do the work of classification, with a prepared set of data inputs

(defn- unique-keywords [data]
  "Retrieve unique numerical keywords to match items in data."
  (vec (map keyword (map str (range (count data))))))

(defn get-dataset [data config]
  "Weka dataset ready for classification or training."
  (let [category? (= :category (get-in config [:classification :group]))
        header (conj (unique-keywords (drop-last (first data))) (if category? {:c [:a :b]} :c))]
    (make-dataset "ds" header data {:class :c})))

(defn train-classifier [data-file pos-file ref-file config & {:keys [downsample]}]
  "Manage retrieving data and training the classifier."
  (let [class-data (prep-classifier-data data-file pos-file ref-file config :downsample downsample)
        ds (get-dataset class-data config)
        c (-> (apply make-classifier (conj (get-in config [:classification :classifier])
                                           (get-in config [:classification :options] {})))
              (classifier-train ds))]
    c))

(defn prepare-classifier
  "High level work to get classifier, including serialization to a file."
  [data-file pos-file ref-file work-dir config & {:keys [downsample always-prep?]}]
  (let [out-dir (str (fs/file work-dir "classifier"))
        classifier-file (str (fs/file out-dir
                                      (format "%s-classifier%s.bin"
                                              (-> config :classification :classifier last name)
                                              (if downsample (format "-%.1f" downsample) ""))))]
    (if-not (fs/exists? out-dir)
      (fs/mkdirs out-dir))
    (if (and (not always-prep?) (fs/exists? classifier-file))
      (deserialize-from-file classifier-file)
      (let [classifier (train-classifier data-file pos-file ref-file config :downsample downsample)]
        (serialize-to-file classifier classifier-file)
        classifier))))

;; Call variants given a trained read classifier:
;; - Walk through each position in the input file.
;;    - Filter reads based on the pre-built linear classifier
;;    - Retrieve percentage of called at a position + classification
;;    - Check known variants and determine if correct,

(defn call-vrns-at-pos [reads passes? config]
  "With read parameters at a position, get map of bases and percent present."
  (letfn [(percents [xs]
            (let [total (apply + (vals xs))]
              (list
               (if (> total 0)
                 (reduce (fn [m [k v]] (assoc m k (* 100.0 (/ v total)))) {} xs)
                 {})
               total)))
          (percents-by-base [reads]
            (->> reads
                 (map #(repeat (if (get config :unique-only false) 1 (get % :num 1))
                               (:base %)))
                 flatten
                 frequencies
                 percents))
          (organize-filter-calls [min-freq [orig total]]
            (let [min-freq-percent (* 100.0 min-freq)
                  want (->> orig
                            (filter #(>= (second %) min-freq-percent))
                            (map first))
                  filter-calls (select-keys orig want)
                  back-filter-freq (->> orig
                                        (remove (fn [[k v]] (contains? (set want) k)))
                                        (sort-by second <)
                                        first
                                        second)]
              {:calls (select-keys orig want)
               :total total
               :back-filter-freq back-filter-freq}))]
    (let [position ((juxt :space :pos) (first reads))
          clean-reads (remove #(= "N" (:base %)) reads)
          raw-percents (first (percents-by-base clean-reads))
          vrns-at-pos (->> clean-reads
                           (filter (fn [xs]
                                     (or (> (get raw-percents (:base xs))
                                            (get-in config [:classification :max-neg-pct]))
                                         (apply passes?
                                                ((juxt :qual :kmer-pct :map-score) xs)))))
                           percents-by-base
                           (organize-filter-calls (:min-freq config)))]
      (assoc vrns-at-pos :position position))))

(defn classifier-checker [classifier config]
  "Calculate probability of read inclusion using pre-built classifier."
  (fn [qual kmer-pct map-score]
    (let [data (metrics-to-features qual kmer-pct map-score config)
          category? (= :category (get-in config [:classification :group]))
          default (if category? :b -1)
          ds-instance (make-instance (get-dataset [(conj (vec data) default)] config)
                                                    (-> (zipmap (unique-keywords data) data)
                                                        (assoc :c default)))
          score (classifier-classify classifier ds-instance)]
      (if category?
        (== score 0.0)
        (> score (-> config :classification :pass-thresh))))))

(defn assign-position-type [pos base-counts total known-vrns config]
  "Given read calls, determine if true/false and positive/negative, type."
  (letfn [(annotate-base [pos base freq known-vrns]
            {:pos pos :base base :freq freq
             :exp-freq (get known-vrns [(first pos) (second pos) base])})
          (is-match? [config base]
            (or (and (nil? (:freq base)) (nil? (:exp-freq base)))
                (roughly-freq? (:freq base) (:exp-freq base) config)))
          (called-match? [bases config]
            (every? (partial is-match? config) bases))
          (low-freq-match? [bases config]
            (let [sum-freqs (apply + (remove nil? (map :exp-freq bases)))]
              (if (< sum-freqs 99.0)
                (is-match? config (first (remove #(nil? (:exp-freq %)) bases)))
                false)))
          (finalize [bases call]
            (let [freq (when-not (nil? call)
                         (apply min (remove nil? (map :exp-freq bases))))]
              {:class call :freq freq :calls base-counts :total-reads total :pos pos}))]
    (let [all-bases ["A" "C" "G" "T"]
          ready-bases (map (fn [b] (annotate-base pos b (get base-counts b) known-vrns))
                           all-bases)
          e-bases (map :base (filter #(not (nil? (:exp-freq %))) ready-bases))
          c-bases (map :base (filter #(not (nil? (:freq %))) ready-bases))]
      (finalize ready-bases
                (case (count e-bases)
                  0 nil
                  1 (cond
                     (called-match? ready-bases config) :true-negative
                     (low-freq-match? ready-bases config) :true-positive
                     (> (count c-bases) 1) :false-positive
                     (== (count c-bases) 0) :false-negative)
                  (cond
                   (called-match? ready-bases config) :true-positive
                   :else :false-negative))))))

(defn raw-data-frequencies [data-file config]
  "Retrieve base frequencies based on raw input data, without filtering."
  (with-open [rdr (reader data-file)]
    (->> (raw-reads-by-pos rdr config)
         (map (fn [xs] (call-vrns-at-pos xs (fn [_ _ _] true) config)))
         doall)))

(defn assess-classifier
  "Determine rates of true/false positive/negatives with trained classifier"
  [data-file pos-file ref-file classifier config & {:keys [downsample]}]
  (let [classifier-passes? (classifier-checker classifier config)
        known-vrns (read-ref pos-file ref-file 101.0)]
    (with-open [rdr (reader data-file)]
      (->> (raw-reads-by-pos rdr config)
           (map (partial downsample-at-pos downsample))
           (map (fn [xs] (call-vrns-at-pos xs classifier-passes? config)))
           (map (fn [x] (assign-position-type (:position x) (:calls x)
                                              (:total x) known-vrns config)))
           (remove #(nil? (:class %)))
           (map (fn [xs] (when (:verbose config) (println xs)) xs))
           vec))))

(defn add-classification-info
  "Add detailed classification information to configuration."
  [config]
  (reduce (fn [coll [k v]] (assoc-in coll [:classification k] v))
          config
          (case (get-in config [:classification :algorithm])
            "random-forest" {:classifier [:decision-tree :random-forest]
                             :random-sample false
                             :options {:num-features-to-consider 14
                                       :num-trees-in-forest 120
                                       :random-seed 1}
                             :group :category}
            {:classifier [:regression :linear]
             :random-sample true
             :group :numerical})))

(defn pipeline-prep-classifier
  [data-file pos-file ref-file work-dir config & {:keys [evaluate? always-prep?
                                                         downsample]}]
  (let [c (prepare-classifier data-file pos-file ref-file work-dir config
                              :downsample downsample :always-prep? always-prep?)]
    (println c)
    (when evaluate?
      (let [a (assess-classifier data-file pos-file ref-file c config
                                 :downsample downsample)]
        (write-assessment a data-file work-dir)
        (print-vrn-summary a)))
    c))

(defn prep-config-files
  [run-config-file param-config-file work-dir]
  {:run-config (load-run-config run-config-file work-dir)
   :config (-> (load-config param-config-file)
                   (assoc :verbose true)
                   add-classification-info)})

(defn prep-classifier-info
  [run-config config work-dir]
  (let [classify-exp (first (filter :classify (:experiments run-config)))
        data-file (:files classify-exp)
        pos-file (get-in run-config [:ref :control])
        ref-file (get-in run-config [:ref :files])]
    {:c (pipeline-prep-classifier data-file pos-file ref-file work-dir config
                                  :evaluate? (:evaluate classify-exp))
     :ref-file ref-file
     :pos-file pos-file
     :data-file data-file}))

(defn -main [run-config-file param-config-file work-dir]
  (let [{:keys [config run-config]} (prep-config-files run-config-file param-config-file work-dir)
        {:keys [c ref-file pos-file data-file]} (prep-classifier-info run-config config work-dir)
        roc-classes [:true-positive :false-positive :true-negative :false-negative]
        out-file (file work-dir "classifier" "downsample.csv")]
    (with-open [wtr (writer out-file)]
      (.write wtr (str (string/join "," (concat ["downsample" "rep"] (map name roc-classes)))
                       "\n"))
      (doseq [pct (get-in run-config [:downsample :percents])]
        (doseq [i (range (get-in run-config [:downsample :replicates]))]
          (let [a (assess-classifier data-file pos-file ref-file c config
                                     :downsample pct)
                a-sum (roc-summarize-assessment a)]
            (.write wtr (str (string/join "," (concat [pct i] (map #(get a-sum % 0) roc-classes)))
                             "\n"))
            (.flush wtr)))))))
