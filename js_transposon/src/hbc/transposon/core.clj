;; Organize sequence identified transposon insertion sites across
;; multiple experiments.

(ns hbc.transposon.core
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [clojure.string :as string]
            [clj-yaml.core :as yaml]))

;; Collapse read locations into related groups by position

(defn output-combined [out-fname exps combo]
  "Write combined data to output CSV file."
  (letfn [(output-a-pos [pos-data]
            (let [all (-> pos-data vals first)]
              (concat
               [(:space all)
                (string/join ";" (:pos all))]
               (map #(get-in pos-data [% :count] 0) exps)
               [(string/join ";" (:seq all))])))]
    (with-open [out-h (io/writer out-fname)]
      (csv/write-csv out-h (concat
                            [(concat ["chr" "pos"] exps ["seq"])]
                            (map output-a-pos combo)))
      out-h)))

(defn combine-locations [data]
  "Combine data locations, grouped by experiment at each insertion position"
  (letfn [(collapse-an-exp [[k v] all-info]
            [k (assoc all-info :count (apply + (map :count v)))])
          (collapse-a-pos [pos-data]
            (let [all-info {:space (-> pos-data first :space)
                             :pos (set (map :pos pos-data))
                             :seq (set (map :seq pos-data))}]
              (apply hash-map (flatten
                               (map #(collapse-an-exp % all-info) (group-by :exp pos-data))))))]
    (map collapse-a-pos data)))

(defn- close-to-last [config]
  "Predicate checking if a value is within range of the last, used to partition."
  (let [range (-> config :algorithm :distance)
        last (atom 0)]
    (fn [x]
      (let [check (+ @last range)]
        (reset! last (:start x))
        (>= check (:start x))))))

(defn- combine-by-position-space [data config]
  "Combine set of groups pre-configured to be on the same chromosome/contig space."
  (let [close-to-last? (close-to-last config)]
    (->> data
         (sort-by :start)
         (partition-by close-to-last?))))

(defn combine-by-position [data config]
  "Combine read by position, grouping nearby locations together."
  (apply concat
   (map #(combine-by-position-space (second %) config)
        (group-by :space data))))

;; Read from custom tab-delimited input files

(defn read-custom-positions [file-name]
  "Read custom position file into location details on each line."
  (with-open [rdr (io/reader file-name)]
    (doall (map
            (fn [[_ strand space start seq _ _ n _ _]]
              {:strand strand :space space
               :pos (if (= strand "+") (Integer/parseInt start)
                        (+ (Integer/parseInt start) (count seq)))
               :seq seq :count (Integer/parseInt n)})
            (csv/read-csv rdr :separator \tab)))))

;; Read configuration details from input YAML file

(defn -main [config-file]
  (let [config (-> config-file slurp yaml/parse-string)]
    (println config)))
