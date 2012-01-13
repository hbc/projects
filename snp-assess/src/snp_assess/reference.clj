;; Create reference set of variations for control lanes from known
;; mixed populations

(ns snp-assess.reference
  (:import [org.biojava3.core.sequence.io FastaReaderHelper])
  (:use [clojure.java.io])
  (:require [clj-yaml.core :as yaml]))

(defn get-fasta-seq-map [in-file]
  "Parse FASTA input file to in-memory map of ids to sequences."
  (let [seq-map (FastaReaderHelper/readFastaDNASequence (file in-file) true)]
    (reduce #(assoc %1 %2 (->> %2 (.get seq-map) (.toString))) {} (.keySet seq-map))))

(defn- generate-reference-at-pos [ref-bases percents]
  "Generate reference bases and frequencies at a position.")

(defn generate-reference
  "Generate VCF formatted reference sequence for defined population."
  [seqs percents]
  {:pre [(= 100.0 (apply + (vals percents)))
         (= 1 (count (set (map count (vals seqs)))))]})

(defn -main [ref-fasta ref-config]
  (let [seqs (get-fasta-seq-map ref-fasta)
        config (-> ref-config slurp yaml/parse-string :reference)]
    (generate-reference seqs config)))
