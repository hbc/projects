(ns hbc.transposon.core
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]))

(defn- close-to-last [config]
  "Predicate checking if a value is within range of the last, used to partition."
  (let [range (-> config :algorithm :distance)
        last (atom 0)]
    (fn [x]
      (let [check (+ @last range)]
        (reset! last (:start x))
        (>= check (:start x))))))

(defn- combine-by-position-space [data config]
  "Combine set of groups pre-configured to be on the same chromosome/contig space."
  (let [close-to-last? (close-to-last config)]
    (->> data
         (sort-by :start)
         (partition-by close-to-last?)
         )))

(defn combine-by-position [data config]
  "Combine read by position, grouping nearby locations together."
  (apply concat
   (map #(combine-by-position-space (second %) config)
        (group-by :space data))))

;; Read from custom tab-delimited input files

(defn read-custom-positions [file-name]
  "Read custom position file into location details on each line."
  (with-open [rdr (io/reader file-name)]
    (doall (map
            (fn [[_ strand space start seq _ _ n _ _]]
              {:strand strand :space space
               :pos (if (= strand "+") (Integer/parseInt start)
                        (+ (Integer/parseInt start) (count seq)))
               :seq seq :count (Integer/parseInt n)})
            (csv/read-csv rdr :separator \tab)))))
