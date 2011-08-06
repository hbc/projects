(defproject snp-assess "0.0.1-SNAPSHOT"
  :description "Deep sequence variation assessment with Hadoop."
  :dependencies [[org.clojure/clojure "1.2.1"]
                 [cascalog "1.8.0"]
                 [incanter "1.2.3" :exclusions [swank-clojure]]
                 [com.leadtune/clj-ml "0.1.2"]
                 [fs "0.8.1"]]
  :dev-dependencies [[org.apache.hadoop/hadoop-core "0.20.2-dev"]
                     [midje "1.1.1"]]
  :run-aliases {:snp-data snp-assess.core
                :off-target snp-assess.off-target
                :min-coverage snp-assess.min-coverage
                :classify snp-assess.classify}
  :aot [snp-assess.core])
