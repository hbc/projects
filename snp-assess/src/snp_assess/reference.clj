;; Create reference set of variations for control lanes from known
;; mixed populations

(ns snp-assess.reference)

(defn -main [ref-fasta ref-config]
  (let [a (read-assessment data-file work-dir)]
    (print-vrn-summary a)))
