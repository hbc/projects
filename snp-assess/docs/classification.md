## Classification

Low-frequency variations at each position in a HIV population contain both
true variants and false changes as a result of sequencing error. To
differentiate these cases, we trained machine learning classifiers using
known variations in our mixed population with three evaluation metrics:

- quality score: The Phred score of sequencing quality at a base, assigned by
  the sequencer.
- mapping score: The alignment score of a read, assigned by the Novoalign aligner.
- k-mer frequency: Frequency of the 13bp region surrounding a
  position. Very low-frequency k-mers are more likely to be sequencing
  artifacts.

We trained a linear classifier to distinguish true and false reads based on
these metrics. Following a TopCoder crowdsourcing competition to improve
classification, we expanded this approach to include a random forest
classifier. Below are plots of detected low frequency variants before and after
filtering. Classification reduced incorrect variants resulting from sequencing
variations to a low background frequency below our limit of detection.
