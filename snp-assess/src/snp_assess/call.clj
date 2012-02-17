(ns snp-assess.call
  "Perform variant calling on patient samples with pre-built read classifier
   Writes output as VCF file with information on allele frequencies, read
   depth and variant effects."
  (:use [clojure.java.io]
        [snp-assess.config :only [default-config]]
        [snp-assess.classify :only [prepare-classifier classifier-checker
                                    raw-reads-by-pos call-vrns-at-pos]]
        [snp-assess.reference :only [convert-to-vc write-vcf-calls]]
        [snp-assess.protein :only [calc-aa-change prep-protein-map]])
  (:require [fs.core :as fs]
            [clj-yaml.core :as yaml]))

(defn vrns-by-pos
  "Lazy sequence of variation calls at each position."
  [rdr classifier aa-finder config]
  (letfn [(vrns-as-vc [[[contig pos] freqs count]]
            (convert-to-vc contig pos freqs :depth count :aa-finder aa-finder))]
    (let [classifier-passes? (classifier-checker classifier config)]
      (->> (raw-reads-by-pos rdr config)
           (map (fn [xs] (call-vrns-at-pos xs classifier-passes? config)))
           (map vrns-as-vc)))))

(defn write-calls-as-vcf
  "Call bases using supplied read classifier, writing as VCF output."
  [read-file ref-file classifier aa-finder config]
  (let [out-file (str (fs/file (fs/parent read-file)
                               (str (fs/name read-file) "-calls.vcf")))]
    (with-open [rdr (reader read-file)]
      (write-vcf-calls (vrns-by-pos rdr classifier aa-finder config)
                       out-file ref-file))))

(defn- read-run-config
  "Read input configuration file, mapping relative paths to absolute."
  [run-config-file work-dir]
  (letfn [(add-full-path [coll keys]
            (loop [ks keys
                   x coll]
              (if (empty? ks)
                x
                (recur (rest ks)
                       (assoc x (first ks)
                              (str (fs/file work-dir (get x (first ks)))))))))]
    (-> (-> run-config-file slurp yaml/parse-string)
        (#(assoc % :ref (add-full-path (:ref %) [:files :known])))
        (#(assoc % :experiments (map (fn [x] (add-full-path x [:files]))
                                     (:experiments %)))))))

(defn -main [run-config-file work-dir]
  (let [run-config (read-run-config run-config-file work-dir)
        _ (println run-config)
        config (-> default-config
                   (assoc :verbose true))
        c (prepare-classifier nil nil work-dir config)
        aa-finder (partial calc-aa-change (prep-protein-map (-> config :ref)))]
    (doseq [x (:experiments config)]
      (write-calls-as-vcf (:files x) (-> config :ref :files)
                          c aa-finder config))))
