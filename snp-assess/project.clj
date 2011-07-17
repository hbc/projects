(defproject snp-assess "0.0.1-SNAPSHOT"
  :description "Deep sequence variation assessment with Hadoop."
  :dependencies [[org.clojure/clojure "1.2.1"]
                 [cascalog "1.8.0-SNAPSHOT"]
                 [incanter "1.2.3" :exclusions [swank-clojure]]]
  :dev-dependencies [[org.apache.hadoop/hadoop-core "0.20.2-dev"]]
  :run-aliases {:snp-data snp-assess.core}
  :aot [snp-assess.core])
