(ns snp-assess.detection
  "Determine limit of detection for variant calling based on background false positives."
  (:use [clojure.java.io]
        [snp-assess.classify :only [prep-config-files prep-classifier-info
                                    classifier-checker read-ref
                                    raw-reads-by-pos downsample-at-pos]]))

(defn background-freqs-at-pos
  "Return background frequencies at this constant position, with downsampling."
  [reads passes? config run-config])

(defn find-background-freqs
  "Determine frequency of background mutations to set cutoffs that minimize false positives"
  [data-file pos-file ref-file classifier config run-config]
  (let [passes? (classifier-checker classifier config)
        known-vrns (read-ref pos-file ref-file 101.0)]
    (with-open [rdr (reader data-file)]
      (->> (raw-reads-by-pos rdr config)
           (map #(background-freqs-at-pos % passes? config run-config))))))

(defn -main [run-config-file param-config-file work-dir]
  (let [{:keys [config run-config]} (prep-config-files run-config-file param-config-file work-dir)
        {:keys [c ref-file pos-file data-file]} (prep-classifier-info run-config config work-dir)]
    (find-background-freqs data-file pos-file ref-file c config run-config)))