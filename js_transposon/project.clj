(defproject hbc.transposon "0.0.1-SNAPSHOT"
  :description "Track and analyze transposon insertion points over time."
  :dependencies [[org.clojure/clojure "1.3.0"]
                 [org.clojure/data.csv "0.1.0"]
                 [clj-yaml "0.3.1"]
                 [fs "0.10.1"]]
  :dev-dependencies [[midje "1.3-alpha4"]]
  :run-aliases {:merge hbc.transposon.core})
