(ns snp-assess.protein.calc
  "Assess influence of variations on protein sequence."
  (:import [org.biojava3.core.sequence.transcription TranscriptionEngine$Builder]
           [org.biojava3.core.sequence DNASequence])
  (:use [clojure.java.io]
        [clojure.string :only (split join)]
        [snp-assess.reference :only (get-fasta-seq-map)]))

;; ## Calculate amino acid changes caused by variants

(defn- translate-codon
  "Translate codon to amino acid sequence using BioJava."
  [codon]
  (let [tx-engine (-> (TranscriptionEngine$Builder.)
                      (.trimStop false)
                      .build)]
    (str (.translate tx-engine (DNASequence. (join codon))))))

(defn- update-codon-w-changes
  [kw codon-info changes prot-map]
  (reduce (fn [coll x]
            (if-let [new (get x kw)]
              (assoc coll :codon
                     (assoc (:codon coll) (:offset (get prot-map (:position x)))
                            new))
              coll))
          codon-info changes))

(defn calc-aa-change
  "Calculate amino acid change caused by a variation as a string representation."
  [prot-map & changes]
  {:pre [(every? #(contains? prot-map (:position %)) changes)]}
  (let [codon-info (update-codon-w-changes :majority
                                           (get prot-map (:position (first changes)))
                                           changes prot-map)
        orig-aa (translate-codon (:codon codon-info))
        new-aa (translate-codon (:codon (update-codon-w-changes :new codon-info changes
                                                                prot-map)))]
    (str orig-aa
         (:aa-pos codon-info)
         new-aa
         (if-let [known-names (get (:known codon-info) new-aa nil)]
           (str "_" (join "_" known-names))))))

;; ## Parse sequence information into map of positions to codons

(defn- ref-to-codons
  "Convert reference sequence into codons."
  [seq-map frame-offset]
  {:pre [(= 1 (count seq-map))]}
  (let [seq (first (vals seq-map))]
    (->> (cons (take frame-offset seq)
              (partition 3 3 "" (drop frame-offset seq)))
         (remove empty?)
         (map vec))))

(defn- add-codon-to-map
  "Add position to codon information for the current codon."
  [codon cur-bp-pos cur-aa-pos known-muts coll]
  (letfn [(add-codon-i-to-map [coll i]
            (assoc coll (+ cur-bp-pos i) {:offset i
                                          :codon codon
                                          :aa-pos cur-aa-pos
                                          :known (get known-muts cur-aa-pos {})}))]
    (reduce add-codon-i-to-map coll (range (count codon)))))

(defn read-known-mutations
  "Read positions of known mutations to associate with variations."
  [known-file]
  (letfn [(add-mut [coll [str-pos aa name]]
            (let [pos (Integer/parseInt str-pos)
                  cur (get coll pos {})]
              (assoc coll pos (assoc cur aa (cons name (get cur aa []))))))]
    (with-open [rdr (reader known-file)]
      (reduce add-mut {}
              (->> rdr
                   line-seq
                   (map #(split % #",")))))))

(defn prep-protein-map
  "Create a map of positions in a sequence to protein changes.
   ref-info is a map with the FASTA input file, offsets to use
   and known mutations."
  [ref-info]
  (let [codons (ref-to-codons (get-fasta-seq-map (:files ref-info))
                              (:frame-offset ref-info))
        known-muts (read-known-mutations (:known ref-info))]
    (loop [cur-codons codons
           cur-bp-pos 0
           cur-aa-pos (:aa-offset ref-info)
           coll {}]
      (if (empty? cur-codons)
        coll
        (recur (rest cur-codons) (+ (count (first cur-codons)) cur-bp-pos) (+ 1 cur-aa-pos)
               (add-codon-to-map (first cur-codons) cur-bp-pos cur-aa-pos
                                 known-muts coll))))))
