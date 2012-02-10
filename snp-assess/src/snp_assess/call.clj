(ns snp-assess.call
  "Perform variant calling on patient samples with pre-built read classifier
   Writes output as VCF file with information on allele frequencies, read
   depth and variant effects."
  (:use [clojure.java.io]
        [snp-assess.config :only [default-config]]
        [snp-assess.classify :only [prepare-classifier classifier-checker
                                    raw-reads-by-pos call-vrns-at-pos]]
        [snp-assess.reference :only [convert-to-vc write-vcf-calls]])
  (:require [fs.core :as fs]))

(defn vrns-by-pos
  "Lazy sequence of variation calls at each position."
  [rdr classifier config]
  (letfn [(vrns-as-vc [[[contig pos] freqs count]]
            (convert-to-vc contig pos freqs :depth count))]
    (let [classifier-passes? (classifier-checker classifier config)]
      (->> (raw-reads-by-pos rdr config)
           (map (fn [xs] (call-vrns-at-pos xs classifier-passes? config)))
           (map vrns-as-vc)))))

(defn write-calls-as-vcf
  "Call bases using supplied read classifier, writing as VCF output."
  [read-file ref-file classifier config]
  (let [out-file (str (fs/file (fs/parent read-file)
                               (str (fs/name read-file) "-calls.vcf")))]
    (with-open [rdr (reader read-file)]
      (write-vcf-calls (vrns-by-pos rdr classifier config) out-file ref-file))))

(defn -main [read-file ref-file work-dir]
  (let [config (-> default-config
                   (assoc :verbose true))
        c (prepare-classifier nil nil work-dir config)]
    (write-calls-as-vcf read-file ref-file c config)))
