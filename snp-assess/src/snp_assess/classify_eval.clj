;; Evaluate classifications of minority variants.
;; These work at the high level of variants called, following read
;; selection via linear regression classification. Human readable
;; summaries allow comparison of classification strategies.

(ns snp-assess.classify-eval
  (:use [clojure.string :only [split]]
        [incanter.stats :only [median]])
  (:require [fs]
            [clj-yaml.core :as yaml]))

(defn summarize-assessment [data]
  "Summarize assessment data to provide reflection of true/false positives"
  (letfn [(summarize-calls [info]
            (reduce (fn [m x]
                      (assoc m (:class x) (inc (get m (:class x) 0))))
                    {} info))
          (collect-freqs [info]
            (->> info
                 (filter #(= (:class %) :true-positive))
                 (map #(-> % :calls vals sort first))))]
      (->> data
        (group-by :freq)
        (map (fn [[freq xs]] [freq (summarize-calls xs) (collect-freqs xs)])))))

(defn print-vrn-summary [data]
  "Print high level summary of correct and incorrect expected bases by frequency"
  (letfn [(sum-counts [counts want] (apply + (vals (select-keys counts want))))
          (all-counts [xs want] (apply + (map #(sum-counts % want) xs)))]
    (println "| percent | correct | wrong | median |")
    (doseq [[freq counts exp-freqs] (sort-by first > (summarize-assessment data))]
      (println (format "| %s | %s | %s | %.1f |" freq
                       (sum-counts counts [:true-positive])
                       (sum-counts counts [:false-positive :false-negative])
                       (median exp-freqs))))
    (println "|     | correct | wrong |")
    (doseq [[name s e] [["100%" 100.0 100.0] ["5+%" 5.1 99.9] ["0-5%" 0 5.0]]]
      (let [c (->> (summarize-assessment data)
                   (filter (fn [[x _ _]] (and (>= x s) (<= x e))))
                   (map second))
            y (all-counts c [:true-positive])
            n (all-counts c [:false-positive :false-negative])] 
        (println (format "| %s | %s (%.1f%%) | %s (%.1f%%) |" name
                         y (* 100.0 (/ y (+ y n)))
                         n (* 100.0 (/ n (+ y n)))))))))

(defn- dump-fname [data-file work-dir]
  "Retrieve location of a dump YAML File from the data file and work directory."
  (fs/join work-dir "classifier" (-> data-file
                                     fs/basename
                                     (split #"\.")
                                     first
                                     (#(format "%s.%s" % "yaml")))))

(defn write-assessment [data data-file work-dir]
  "Write assessment information to YAML output file"
  (let [dump-file (dump-fname data-file work-dir)]
    (spit dump-file (yaml/generate-string data))
    dump-file))

(defn read-assessment
  "Read assessment details using standard dump file name from data file."
  ([dump-file] (->> dump-file
                    slurp
                    yaml/parse-string
                    (map (fn [x] (assoc x :class (keyword (:class x)))))))
  ([data-file work-dir] (read-assessment (dump-fname data-file work-dir))))

(defn -main [data-file work-dir]
  (let [a (read-assessment data-file work-dir)]
    (print-vrn-summary a)))
