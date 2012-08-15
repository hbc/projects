(defproject snp-assess "0.0.1-SNAPSHOT"
  :description "Deep sequence variation assessment for populations."
  :dependencies [[org.clojure/clojure "1.4.0"]
                 [org.clojure/algo.generic "0.1.0" :exclusions [org.clojure/clojure]]
                 [cascalog "1.8.5" :exclusions [org.clojure/clojure]]
                 [incanter/incanter-core "1.3.0" :exclusions [org.clojure/clojure]]
                 [incanter/incanter-charts "1.3.0" :exclusions [org.clojure/clojure]]
                 [incanter/incanter-pdf "1.3.0" :exclusions [org.clojure/clojure]]
                 [nz.ac.waikato.cms.weka/weka-stable "3.6.6"]
                 [org.clojars.chapmanb/fast-random-forest "0.98"]
                 [com.leadtune/clj-ml "0.2.2" :exclusions [lt/weka hr.irb/fastRandomForest
                                                           org.clojure/clojure]]
                 [fs "1.1.2" :exclusions [org.clojure/clojure]]
                 [clj-yaml "0.3.1"]
                 [ordered "1.0.0" :exclusions [org.clojure/clojure]]
                 [doric "0.7.0-SNAPSHOT"]
                 [bcbio.variation "0.0.1-SNAPSHOT"]]
  :profiles {:dev
             {:dependencies
              [[org.apache.hadoop/hadoop-core "0.20.2-dev" :exclusions [commons-logging org.slf4j/slf4j-api
                                                                        org.slf4j/slf4j-log4j12 log4j]]
               [midje "1.4.0" :exclusions [org.clojure/clojure]]]}}
  :plugins [[lein-midje "2.0.0-SNAPSHOT"]]
  :repositories {"biojava" "http://www.biojava.org/download/maven/"}
  :aliases {"snp-data" ["run" "-m" "snp-assess.core"]
            "snp-off-target" ["run" "-m" "snp-assess.off-target"]
            "snp-min-coverage" ["run" "-m" "snp-assess.min-coverage"]
            "snp-classify" ["run" "-m" "snp-assess.classify"]
            "snp-classify-eval" ["run" "-m" "snp-assess.classify-eval"]
            "snp-classify-read" ["run" "-m" "snp-assess.classify-read"]
            "snp-call" ["run" "-m" "snp-assess.call"]
            "snp-reference" ["run" "-m" "snp-assess.reference"]}
  :jvm-opts ["-Xmx2g"]
  :aot [snp-assess.core])
