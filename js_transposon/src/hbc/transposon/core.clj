(ns hbc.transposon.core
  "Top level access to functionality from pre-build jars."
  (:require [clojure.string :as string]
            [hbc.transposon.merge]
            [hbc.transposon.score])
  (:gen-class))

(def ^{:doc "Mapping of commandline argument to main functions"
       :private true}
  altmain-map
  {:merge hbc.transposon.merge/-main
   :score hbc.transposon.score/-main})

(defn -main [& args]
  (if-let [alt-fn (get altmain-map (keyword (first args)))]
    (apply alt-fn (rest args))
    (println "Unexpected sub-command. Need one of:" (string/join ", "(map name (keys altmain-map))))))