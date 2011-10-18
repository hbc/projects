(ns hbc.transposon.core
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]))

;; Read from custom tab-delimited input files

(defn read-custom-positions [file-name]
  "Read custom position file into location details on each line."
  (with-open [rdr (io/reader file-name)]
    (doall (map
            (fn [[_ strand space start seq _ _ count _ _]]
              {:strand strand :space space :start (Integer/parseInt start)
               :seq seq :count (Integer/parseInt count)})
            (csv/read-csv rdr :separator \tab)))))
