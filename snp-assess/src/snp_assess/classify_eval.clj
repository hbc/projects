;; Evaluate classifications of minority variants.
;; These work at the high level of variants called, following read
;; selection via linear regression classification. Human readable
;; summaries allow comparison of classification strategies.

(ns snp-assess.classifier-eval
  (:use [clojure.string :only [split]])
  (:reqires [fs]
            [clj-yaml.core :as yaml]))

(defn summarize-assessment [data]
  "Summarize assessment data to provide reflection of true/false positives"
  (letfn [(summarize-calls [info]
            (reduce (fn [m x]
                      (assoc m (:class x) (inc (get m (:class x) 0))))
                    {} info))
          (ratios [xs]
            (let [total (apply + (vals xs))]
              (if (> total 0)
                (reduce (fn [m [k v]] (assoc m k (* 100.0 (/ v total)))) {} xs)
                {})))]
      (->> data
        (group-by :freq)
        (map (fn [[freq xs]] [freq (summarize-calls xs)]))
        (map (fn [[freq xs]] [freq xs (ratios xs)])))))

(defn print-vrn-summary [data]
  "Print high level summary of correct and wrong by frequency"
  (letfn [(sum-counts [counts want] (+ (vals (select-keys counts want))))]
    (doseq [[freq counts _] (sort-by first > (summarize-assessment data))]
      (println (format "| %s | %s | %s | " freq
                       (sum-counts counts [:true-positive])
                       (sum-counts counts [:false-positive :false-negative]))))))

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
  [data-file work-dir] (read-assessment (dump-fname data-file work-dir))
  "Read assessment details using standard dump file name from data file."
  [dump-file]
  "Read assessment details from YAML dump file."
  (-> dump-file
      slurp
      yaml/parse-string))

(defn -main [data-file work-dir]
  (let [a (read-assessment data-file work-dir)]
    (print-vrn-summary a)))
