(defproject hbc.transposon "0.0.4"
  :description "Track and analyze transposon insertion points over time."
  :dependencies [[org.clojure/clojure "1.5.1"]
                 [org.clojure/data.csv "0.1.2"]
                 [org.clojure/tools.cli "0.2.2"]
                 [clj-yaml "0.4.0"]
                 [dk.ative/docjure "1.6.0"]
                 [me.raynes/fs "1.4.0"]
                 [incanter/incanter-core "1.4.1"]
                 [incanter/incanter-charts "1.4.1"]
                 [incanter/incanter-io "1.4.1"]
                 [lonocloud/synthread "1.0.3" :exclusions [org.clojure/clojure]]]
  :min-lein-version "2.0.0"
  :profiles {:dev {:dependencies [[midje "1.5-RC1" :exclusions [org.clojure/clojure]]]}}
  :plugins [[lein-midje "3.0-RC1"]]
  :javac-options ["-nowarn" "-target" "1.6" "-source" "1.6"]
  :main hbc.transposon.core
  :aliases {"merge" ["run" "-m" "hbc.transposon.merge"]
            "score" ["run" "-m" "hbc.transposon.score"]})