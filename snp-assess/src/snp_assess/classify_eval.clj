;; Evaluate classifications of minority variants.
;; These work at the high level of variants called, following read
;; selection via linear regression classification. Human readable
;; summaries allow comparison of classification strategies.

(ns snp-assess.classify-eval
  (:use [clojure.java.io]
        [clojure.string :only [split join]]
        [incanter.stats :only [median]])
  (:require [fs.core :as fs]
            [clj-yaml.core :as yaml]
            [doric.core :as doric]))

;; Provide summary details about read calls based on frequency

(defn summarize-assessment [data]
  "Summarize assessment data to provide reflection of true/false positives"
  (letfn [(summarize-calls [info]
            (reduce (fn [m x]
                      (assoc m (:class x) (inc (get m (:class x) 0))))
                    {} info))
          (collect-freqs [info]
            (->> info
                 (filter #(contains? #{:true-positive :true-negative} (:class %)))
                 (map #(-> % :calls vals sort first))))]
    (->> data
         (group-by :freq)
         (map (fn [[freq xs]] [freq (summarize-calls xs) (collect-freqs xs)])))))

(defn print-vrn-summary [data]
  "Print high level summary of correct and incorrect expected bases by frequency"
  (letfn [(sum-counts [counts want] (apply + (vals (select-keys counts want))))
          (all-counts [xs want] (apply + (map #(sum-counts % want) xs)))
          (str-w-percent [x xs]
            (if (= 0 (apply + xs))
              "0.0"
              (format "%s (%.1f%%)" x (* 100.0 (/ x (apply + xs))))))]
    (println (doric/table [{:name :percent :format #(format "%.1f" %)}
                           :correct
                           :wrong
                           {:name :median :format #(format "%.1f" %)}]
                          (for [[f counts efs] (sort-by first > (summarize-assessment data))]
                            {:percent f
                             :correct (sum-counts counts [:true-positive :true-negative])
                             :wrong (sum-counts counts [:false-positive :false-negative])
                             :median (median efs)})))
    (println (doric/table [:category :correct :wrong]
                          (for [[name s e] [["100%" 99.9 100.0] ["5+%" 5.0 99.9] ["0-5%" 0.0 5.0]]]
                            (let [c (->> (summarize-assessment data)
                                         (filter (fn [[x _ _]] (and (> x s) (<= x e))))
                                         (map second))
                                  y (all-counts c [:true-positive :true-negative])
                                  n (all-counts c [:false-positive :false-negative])]
                              {:category name
                               :correct (str-w-percent y [y n])
                               :wrong (str-w-percent n [y n])}))))))

;; Write calls to output file for querying and adjusting parameters:
;; incorrect false-positives and low frequency true-positives

(defn in-depth-calls [data]
  "Retrieve calls of interest for in-depth analysis")

(defn- classifier-check-fname [data-file work-dir ext]
  "Retrieve location of a classifier related file; ext is the unique file extension."
  (str (fs/file work-dir "classifier" (-> data-file
                                          fs/base-name
                                          (split #"\.")
                                          first
                                          (#(format "%s%s" % ext))))))

(defn write-indepth-calls [data data-file work-dir]
  "Write details on correct/wrong regions of interest for by-hand evaluation."
  (let [in-depth-fname (classifier-check-fname data-file work-dir "-indepth.tsv")]
    (if-not (fs/file? in-depth-fname)
      (with-open [wrtr (writer in-depth-fname)]
        (doseq [cur-indepth (in-depth-calls data)]
          (.write wrtr (format "%s\n" (join "\t" cur-indepth))))))))

;; Manage reading and writing details about each position to YAML files

(defn write-assessment [data data-file work-dir]
  "Write assessment information to YAML output file"
  (let [dump-file (classifier-check-fname data-file work-dir ".yaml")]
    (spit dump-file (yaml/generate-string data))
    dump-file))

(defn read-assessment
  "Read assessment details using standard dump file name from data file."
  ([dump-file] (->> dump-file
                    slurp
                    yaml/parse-string
                    (map (fn [x] (assoc x :class (keyword (:class x)))))))
  ([data-file work-dir] (read-assessment (data-file work-dir))))

(defn -main [data-file work-dir]
  (let [a (read-assessment data-file work-dir)]
    (print-vrn-summary a)))
