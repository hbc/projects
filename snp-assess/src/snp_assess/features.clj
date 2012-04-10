(ns snp-assess.features
  "Create features from a set of metrics for input to classification algorithms.
  Handles normalization for linear regression analysis and feature expansion for
  random forest separation."
  (:use [incanter.core :only [log pow exp]]
        [incanter.stats :only [euclidean-distance manhattan-distance]]))

;; ## Top-level functionality

(defmulti metrics-to-features
  "Convert input metrics into remapped features for classification."
  (fn [qual kmer-pct map-score config] (get-in config [:classification :classifier])))

;; ## Linear regression

(defn min-max-norm
  "Perform min-max normalization, truncated larger values at max and min."
  [score [minv maxv]]
  (let [trunc-score-max (if (< score maxv) score maxv)
        trunc-score (if (> trunc-score-max minv) trunc-score-max minv)]
    (/ (- trunc-score minv) (- maxv minv))))

(defmethod metrics-to-features [:regression :linear]
  [qual kmer-pct map-score config]
  "Prep quality score, k-mer and map score with min/max normalization."
  [(min-max-norm qual (:qual-range config))
   (min-max-norm kmer-pct (:kmer-range config))
   (min-max-norm map-score (:map-score-range config))])

;; ## Random forest

(defn- rf-distance-features
  "Provide distance based features of inputs in a unit cube.
  Normalizes the 3 input features to a point in a 1 by 1 by 1 cube,
  and then calculates distance metrics of this point compared to the
  start, middle and end of the cube."
  [qual k_hat map-score]
  (let [p [(/ (- qual 3) 32) (/ (+ k_hat 15) 15) (/ map-score 249)]
        p_start [0 0 0]
        p_middle [0.5 0.5 0.5]
        p_end [1 1 1]]
    (reduce (fn [all p_cur] (concat all
                                    [(exp (- (euclidean-distance p p_cur)))
                                     (- (log (manhattan-distance p p_cur)))]))
            [] [p_start p_middle p_end])))

(defn- get-rf-features
  [qual kmer-pct map-score _]
  "Provide feature expansion for random forest classification.
  This is a direct port of the work done by Chun-Sung Ferng on the TopCoder
  competition. Produces 14 feature metrics from input kmers. The first 8 metrics
  are log or square transformations on inputs multiplied together. The remaining
  6 metrics are similarity comparisons between the 3 metrics and start/median/end
  reference points in 3-space."
  (let [k_bar (log kmer-pct)
        k_hat (+ k_bar 20)
        m_hat (/ map-score 8)]
    (concat
     [(* (pow qual 2) k_hat (pow m_hat 2))
      (* (pow qual 2) (pow k_hat 2) m_hat)
      (* k_hat m_hat)
      (* qual k_hat m_hat)
      (* qual (pow k_hat 2) m_hat)
      (* qual (pow k_hat 2) (pow m_hat 2))
      (/ (* qual m_hat) (pow k_bar 2))
      (/ (* qual m_hat) (- k_bar))]
     (rf-distance-features qual k_hat map-score))))

(defmethod metrics-to-features [:decision-tree :fast-random-forest]
  [& args]
  (apply get-rf-features args))
  
(defmethod metrics-to-features [:decision-tree :random-forest]
  [& args]
  (apply get-rf-features args))
