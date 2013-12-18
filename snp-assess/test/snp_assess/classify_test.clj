(ns snp-assess.classify-test
  (:use [midje.sweet]
        [snp-assess.core :only [load-config]]
        [snp-assess.classify]
        [snp-assess.features :only [metrics-to-features]]
        [snp-assess.classify-eval :only [summarize-assessment roc-summarize-assessment]])
  (:require [me.raynes.fs :as fs]))

;; Top level definitions useful for all tests
(background
 (around :facts
         (let [data-dir (str (fs/file "test" "data"))
               vrn-file (str (fs/file data-dir "reference" "reference.vcf"))
               ref-file (str (fs/file data-dir "reference" "hxb2.fa"))
               data-file (str (fs/file data-dir "raw" "raw_variations.tsv"))
               config-file (str (fs/file "config" "classification.yaml"))
               config (-> (load-config config-file)
                          (assoc-in [:classification :naive-min-score] 1.2)
                          (assoc-in [:classification :max-neg-pct] 5.0)
                          (assoc-in [:classification :pass-thresh] 0.0)
                          (assoc-in [:classification :classifier] [:regression :linear])
                          (assoc-in [:classification :group] :numerical))]

           ?form)))

(facts "Read variant positions from file"
  (read-ref vrn-file ref-file 5.1) => (just {["HXB2" 951 "T"] (roughly 5.0),
                                             ["HXB2" 953 "A"] (roughly 1.0)})
  (read-ref vrn-file ref-file 3.0) => (just {["HXB2" 953 "A"]
                                             (roughly 1.0)}))

(facts "Classification data from file"
  (let [data (prep-classifier-data data-file vrn-file ref-file config)]
    (count data) => 6                   ; Total minority variants
    (count (first data)) => 4           ; items at each position
    (map last data) => [1 1 1 0 0 0]    ; classifications
    (first data) => (contains [0.228 1])))

(facts "Provide random downsampling of classification data."
  (let [prep-fn (partial prep-classifier-data data-file vrn-file ref-file config)]
    (count (first (prep-fn))) => 4
    (count (first (prep-fn :downsample 0.5))) => #(contains? #{0 4} %)))

(facts "Build a full classifier from file data"
  (let [c (train-classifier data-file vrn-file ref-file config)]
    (map float (.coefficients c)) => [0.0 0.0 0.0 0.0 0.5]
    (.toString c) => (contains "Linear Regression Model")
    ((classifier-checker c config) 20 1.0E-3 200) => true))

;.;. There is an inevitable reward for good deeds. -- Ming Fu Wu
(facts "Random forest classification"
  (let [rf-config (-> config
                      (assoc-in [:classification :algorithm] "random-forest")
                      add-classification-info)
        c (train-classifier data-file vrn-file ref-file rf-config)
        assess-data (assess-classifier data-file vrn-file ref-file c rf-config) => nil]
    (.toString c) => (contains "Random forest")
    (roc-summarize-assessment assess-data) => {:true-positive 1 :false-positive 1 :false-negative 1}
    (map :class assess-data) => [:true-positive :false-positive :false-negative]))

(facts "Assess a classifier based on variant calling ability"
  (let [c (train-classifier data-file vrn-file ref-file config)
        get-class (fn [p r e] (:class (assign-position-type p r 0 e config)))
        assess-data (assess-classifier data-file vrn-file ref-file c config)
        expect-constant {["chr" 10 "G"] 100.0}
        expect-low {["chr" 10 "G"] 95.0 ["chr" 10 "A"] 5.0}]
    (get-class ["chr" 10] {"G" 99.0} expect-constant) => :true-negative
    (get-class ["chr" 10] {"G" 95.0 "A" 5.0} expect-constant) => :false-positive
    (get-class ["chr" 10] {"G" 95.0 "A" 5.0} expect-low) => :true-positive
    (get-class ["chr" 10] {"G" 83.0 "A" 17.0} expect-low) => :false-negative
    (get-class ["chr" 10] {"G" 90.0 "A" 5.0 "C" 3.0} expect-low) => :false-negative
    (get-class ["chr" 10] {"G" 99.0} expect-low) => :false-negative
    (map :class assess-data) => [:false-negative :false-positive :false-negative]
    (first (summarize-assessment assess-data)) => (just [(roughly 5.0)
                                                         {:false-negative 1} []])))

(facts "Remap raw data for classification"
  (let [config (-> (load-config config-file)
                   (assoc-in [:classification :classifier] [:regression :linear]))]
    (finalize-raw-data {:qual 20 :kmer-pct 1.0E-4 :map-score 100} :test config) =>
      (contains [(roughly 0.45) (roughly 9.0E-4) 0.4 :test])
    (finalize-raw-data {:qual 30 :kmer-pct 1.0E-2 :map-score 50} :test2 config) =>
      (contains [(roughly 0.7) (roughly 0.0999) 0.2 :test2])))

(facts "Convert input metrics into classification features."
  (let [metrics [25 0.002 175]]
    (letfn [(test-get-features [x]
              (apply metrics-to-features (conj metrics
                                               (assoc-in config [:classification :classifier] x))))]
      (test-get-features [:regression :linear]) => (just [0.575 (roughly 0.0199) 0.7])
      (count (test-get-features [:decision-tree :fast-random-forest])) => 14)))

(facts "Read raw data collapsed by counts per unique read."
  (let [count-file (str (fs/file data-dir "count_data" "raw_variations_count.tsv"))
        config (-> (load-config config-file)
                   (assoc-in [:classification :assess-bases] nil)
                   (assoc :min-freq 0.0))]
    (let [raw-freqs (raw-data-frequencies count-file config)]
      (count raw-freqs) => 3
      (-> raw-freqs first :position) => ["HXB2" 5]   ; read position
      (-> raw-freqs first :total) => 71169 ; number of total reads
      (-> raw-freqs first :calls (get "G")) => (roughly 99.9747)
      (-> raw-freqs first :calls (get "T")) => (roughly 0.009835))))

