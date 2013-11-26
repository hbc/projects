(ns snp-assess.call
  "Perform variant calling on patient samples with pre-built read classifier
   Writes output as VCF file with information on allele frequencies, read
   depth and variant effects."
  (:use [clojure.java.io]
        [snp-assess.annotate :only [add-mv-annotations]]
        [snp-assess.classify :only [classifier-checker
                                    raw-reads-by-pos call-vrns-at-pos
                                    add-classification-info
                                    prep-config-files prep-classifier-info]]
        [snp-assess.reference :only [convert-to-vc write-vcf-calls]]
        [snp-assess.protein.calc :only [prep-protein-map]]
        [snp-assess.protein.read :only [annotate-calls-w-aa]])
  (:require [me.raynes.fs :as fs]
            [clj-yaml.core :as yaml]
            [bcbio.run.itx :as itx])
  (:gen-class))

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
    (when (itx/needs-run? out-file)
      (with-open [rdr (reader read-file)]
        (write-vcf-calls (vrns-by-pos rdr classifier config)
                         out-file ref-file)))
    out-file))

(defn id-variants [run-config-file param-config-file work-dir]
  (let [{:keys [config run-config]} (prep-config-files run-config-file param-config-file work-dir)
        {:keys [c ref-file pos-file data-file]} (prep-classifier-info run-config config work-dir)
        prot-map (prep-protein-map (:ref run-config))]
    (doseq [x (:experiments run-config)]
      (let [ref-file (-> run-config :ref :files)
            mv-call-file (write-calls-as-vcf (:files x) ref-file c config)
            aa-call-file (annotate-calls-w-aa (:align x) mv-call-file ref-file prot-map
                                              (:kmer-size config) :count-file (:count x))]
        (add-mv-annotations aa-call-file (:align x) ref-file)))))

(def ^{:doc "Mapping of command line arguments to main functions" :private true}
  altmain-map
  {:snp-call id-variants})

(defn -main [& args]
  (if-let [alt-fn (get altmain-map (keyword (first args)))]
    (apply alt-fn (rest args))
    (apply id-variants args)))
