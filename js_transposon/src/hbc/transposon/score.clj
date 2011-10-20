;; Assess reliability of merged insertions sites by normalizing and
;; assigning quality scores.

(ns hbc.transposon.score
  (:use [incanter.io :only [read-dataset]])
  (:require [incanter.core :as icore]))

(defn normalize-counts
  "Normalize counts to standard metric based on totals in each experiment.
   By default normalizes to total count in a column scaled to 1 million reads."
  [in-file & {:keys [base ignore]
              :or {base 1e6
                   ignore #{:chr :pos :seq}}}]
  (letfn [(normalize [x total]
            (* (/ x total) base))
          (maybe-normalize [ds col]
            (if-not (contains? ignore col)
              (icore/transform-col ds col normalize
                                   (float (apply + (icore/$ col ds))))
              ds))]
    (loop [ds (read-dataset in-file :header true)
           cols (icore/col-names ds)]
      (if-let [col (first cols)]
        (recur (maybe-normalize ds col) (rest cols))
        ds))))

(defn normalize-pos-ratios [ds]
  "Normalize counts as the percentage of the max at a position.
   This normalizes by row, in contrast to normalize-counts
   which handles columns.")
