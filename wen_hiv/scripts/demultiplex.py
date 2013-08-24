"""
Script to perform demultiplexing of the reads

We are expecting it to go like this:

R1: possible_contaminant_sequence_barcode_HIV/MLV_genomic_sequence
R2: possible_contaminant_sequence_barcode_adapter_genomic_sequence

The plan:

- Take first 8 bases from each end
- use Bio.pairwise2 module to align each end to each of the barcodes
Use a large gap penalty, we want no gaps. 1 for correct, 1 for mismatch.
- take the best matching barcode as the correct barcode if it is at
least n-1 of match where n is the barcode length
- if none match then read is unclassified (unclassified)
- if both match the same barcode then the read is that barcode (full_evidence)
- if only one matches a barcode then the read is that barcode (half_evidence)
- if each end matches a different barcode, read is unclassified (ambiguous)

pairwise2.align.localxx("barcode", "sequence")
for a in pairwise2.align.localms("ACCGT", "ACG", 1, -1, 20):
give like -20 for a gap, no gaps
-1 for a mismatch
1 for a correct

"""
