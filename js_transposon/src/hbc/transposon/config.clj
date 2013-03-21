(ns hbc.transposon.config
  "Load experimental information for YAML input files plus optional Excel spreadsheet."
  (:require [clojure.java.io :as io]
            [clj-yaml.core :as yaml]
            [me.raynes.fs :as fs]
            [dk.ative.docjure.spreadsheet :as excel]))

(defn- read-exps
  "Read experimental information from an input excel file."
  [excel-file base-dir]
  (letfn [(find-file [x]
            (first (filter fs/exists?
                           [(str (io/file base-dir x))
                            (str (io/file base-dir (str x ".final.txt")))
                            (str (fs/file (fs/parent base-dir) x))
                            (str (fs/file (fs/parent base-dir) (str x ".final.txt")))
                            x (str x ".final.txt")])))
          (find-lineage [x]
            (cond
             (.startsWith x "Gr") "Gr"
             (.startsWith x "B") "B"
             (.startsWith x "T") "T"
             :else x))
          (process-row [items]
            (let [[fname sample lineage timepoint & args] (map excel/read-cell items)
                  in-file (find-file fname)]
              (when (nil? in-file)
                (throw (Exception. (format "Could not find input file for %s in %s"
                                           fname base-dir))))
              {:sample sample
               :name in-file
               :lineage (find-lineage (str lineage))
               :timepoint (str timepoint)
               :expnorm "auto"}))]
    (->> (excel/load-workbook excel-file)
         excel/sheet-seq
         first
         excel/row-seq
         (map excel/cell-seq)
         (map process-row))))

(defn do-load
  "Load configuration, potentially adding experimental details from input Excel."
  [base-dir yaml-file excel-file]
  (let [config (if yaml-file
                 (-> yaml-file slurp yaml/parse-string)
                 {:algorithm {:rownorm "" :distance 25}
                  :dir {:orig base-dir :out (if (.endsWith base-dir "merge")
                                              base-dir
                                              (str (io/file base-dir "merge")))}})
        ready-excel-file (if (and excel-file (not (fs/absolute? excel-file)))
                           (first (filter fs/exists?
                                          (map #(str (io/file % excel-file))
                                               [base-dir (fs/parent base-dir)])))
                           excel-file)]
    (if excel-file
      (if (nil? ready-excel-file)
        (throw (Exception. (str "Could not find input excel file: " excel-file)))
        (assoc config :experiments (read-exps ready-excel-file base-dir)))
      config)))