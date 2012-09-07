(ns snp-assess.detection
  "Determine limit of detection for variant calling based on background false positives."
  (:use [clojure.java.io]
        [incanter.stats :only [quantile]]
        [snp-assess.classify :only [prep-config-files prep-classifier-info
                                    classifier-checker read-ref
                                    raw-reads-by-pos call-vrns-at-pos downsample-at-pos]]))

(defn background-freqs-at-pos
  "Return background frequencies at this constant position, with downsampling."
  [reads passes? config run-config]
  (println (first reads))
  (mapcat (fn [downsample]
            (->> (call-vrns-at-pos (downsample-at-pos downsample reads) passes?
                                   (assoc config :min-freq 0.0))
                 :calls
                 (map second)
                 (sort >)
                 rest))
          (mapcat (partial repeat (get-in run-config [:downsample :replicates]))
                  (get-in run-config [:downsample :lod-percents]))))

(defn find-background-freqs
  "Determine frequency of background mutations to set cutoffs that minimize false positives"
  [data-file pos-file ref-file classifier config run-config]
  (let [passes? (classifier-checker classifier config)
        constant-is (set (->> (read-ref pos-file ref-file 101.0)
                              (filter #(= 100.0 (second %)))
                              (map first)
                              (map second)))]
    (with-open [rdr (reader data-file)]
      (let [probs [0.8 0.85 0.9 0.925 0.95 0.975 0.99 0.999 0.9999 0.99999 1.0]
            backgrounds (->> (raw-reads-by-pos rdr config)
                             (filter #(contains? constant-is (:pos (first %))))
                             (mapcat #(background-freqs-at-pos % passes? config run-config)))]
        (zipmap probs (quantile backgrounds :probs probs))))))

(defn -main [run-config-file param-config-file work-dir]
  (let [{:keys [config run-config]} (prep-config-files run-config-file param-config-file work-dir)
        {:keys [c ref-file pos-file data-file]} (prep-classifier-info run-config config work-dir)]
    (println (find-background-freqs data-file pos-file ref-file c config run-config))))
