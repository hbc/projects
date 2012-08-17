(ns snp-assess.call
  "Perform variant calling on patient samples with pre-built read classifier
   Writes output as VCF file with information on allele frequencies, read
   depth and variant effects."
  (:use [clojure.java.io]
        [snp-assess.core :only [load-config]]
        [snp-assess.classify :only [pipeline-prep-classifier classifier-checker
                                    raw-reads-by-pos call-vrns-at-pos
                                    add-classification-info]]
        [snp-assess.reference :only [convert-to-vc write-vcf-calls]]
        [snp-assess.protein.calc :only [prep-protein-map]]
        [snp-assess.protein.read :only [annotate-calls-w-aa]])
  (:require [fs.core :as fs]
            [clj-yaml.core :as yaml]))

(defn vrns-by-pos
  "Lazy sequence of variation calls at each position."
  [rdr classifier config]
  (letfn [(vrns-as-vc [x]
            (let [[contig pos] (:position x)]
              (convert-to-vc contig pos (:calls x) :depth (:total x)
                             :back-filter-freq (:back-filter-freq x))))]
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
      (write-vcf-calls (vrns-by-pos rdr classifier config)
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
        (#(assoc % :ref (add-full-path (:ref %) [:files :known :control])))
        (#(assoc % :experiments (map (fn [x] (add-full-path x [:files :align :count]))
                                     (:experiments %)))))))

(defn -main [run-config-file param-config-file work-dir]
  (let [run-config (read-run-config run-config-file work-dir)
        config (-> (load-config param-config-file)
                   (assoc :verbose true)
                   add-classification-info)
        classify-exp (first (filter :classify (:experiments run-config)))
        c (pipeline-prep-classifier (:files classify-exp)
                                    (get-in run-config [:ref :control])
                                    (get-in run-config [:ref :files])
                                    work-dir
                                    config :evaluate? (:evaluate classify-exp))
        prot-map (prep-protein-map (:ref run-config))]
    (doseq [x (:experiments run-config)]
      (let [ref-file (-> run-config :ref :files)
            mv-call-file (write-calls-as-vcf (:files x) ref-file c config)]
        (annotate-calls-w-aa (:align x) mv-call-file ref-file prot-map
                             :count-file (:count x))))))
