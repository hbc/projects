(defproject snp-assess "0.0.1-SNAPSHOT"
  :description "Deep sequence variation assessment for populations."
  :dependencies [[org.clojure/clojure "1.5.1"]
                 [org.clojure/algo.generic "0.1.0" :exclusions [org.clojure/clojure]]
                 [cascalog "1.8.5" :exclusions [org.clojure/clojure log4j]]
                 [incanter/incanter-core "1.5.1" :exclusions [org.clojure/clojure junit]]
                 [incanter/incanter-charts "1.5.1" :exclusions [org.clojure/clojure junit]]
                 [incanter/incanter-pdf "1.5.1" :exclusions [org.clojure/clojure junit]]
                 [bcbio.variation "0.1.1"]]
  :profiles {:dev
             {:dependencies
              [[org.apache.hadoop/hadoop-core "0.20.2-dev" :exclusions [commons-logging org.slf4j/slf4j-api
                                                                        org.slf4j/slf4j-log4j12 log4j
                                                                        commons-codec]]
               [midje "1.4.0" :exclusions [org.clojure/clojure ordered joda-time
                                           org.clojure/math.combinatorics]]]}}
  :plugins [[lein-midje "2.0.0-SNAPSHOT"]]
  :repositories {"biojava" "http://www.biojava.org/download/maven/"}
  :aliases {"snp-data" ["run" "-m" "snp-assess.core"]
            "snp-off-target" ["run" "-m" "snp-assess.off-target"]
            "snp-min-coverage" ["run" "-m" "snp-assess.min-coverage"]
            "snp-classify" ["run" "-m" "snp-assess.classify"]
            "snp-classify-eval" ["run" "-m" "snp-assess.classify-eval"]
            "snp-classify-read" ["run" "-m" "snp-assess.classify-read"]
            "snp-call" ["run" "-m" "snp-assess.call"]
            "snp-lod" ["run" "-m" "snp-assess.detection"]
            "snp-reference" ["run" "-m" "snp-assess.reference"]}
  ;;:jvm-opts ["-Xmx2g"]
  )
