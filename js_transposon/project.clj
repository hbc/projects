(defproject hbc.transposon "0.0.1-SNAPSHOT"
  :description "Track and analyze transposon insertion points over time."
  :dependencies [[org.clojure/clojure "1.3.0"]
                 [org.clojure/data.csv "0.1.0"]
                 [clj-yaml "0.3.1"]
                 [fs "0.11.1"]
                 [incanter "1.3.0-SNAPSHOT"]]
  :dev-dependencies [[midje "1.3.0" :exclusions [org.clojure/clojure]]]
  :run-aliases {:merge hbc.transposon.merge
                :score hbc.transposon.score})
