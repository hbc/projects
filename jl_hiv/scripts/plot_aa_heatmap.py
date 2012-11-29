#!/usr/bin/env python
"""Provide a heatmap of amino acid changes across a set of patient samples.

Usage:
  plot_aa_heatmap.py <sample_dir>
"""
import os
import sys
import glob

import vcf
import pandas
import pandas.rpy.common as com
import rpy2.robjects as rpy

def main(sample_dir):
    max_freq = 1.5
    out_file = os.path.join(sample_dir, "patient_sample_heatmap.pdf")
    freqs = {}
    for fname in glob.glob(os.path.join(sample_dir, "*.vcf")):
        sample_name = os.path.splitext(os.path.basename(fname))[0]
        with open(fname) as in_handle:
            freqs[sample_name] = read_aa_vals(in_handle, max_freq)
    df = pandas.DataFrame(freqs)
    rpy.r.assign("out_file", out_file)
    rpy.r.assign("df", com.convert_to_r_dataframe(df))
    rpy.r.assign("col_names", df.columns.tolist())
    rpy.r('''
    library(pheatmap)
    library(RColorBrewer)
    colnames(df) <- col_names
    pdf(out_file, width=2.5)
    pheatmap(as.matrix(df),
             main="Amino acid change\n    population frequencies",
             cluster_rows=FALSE, cluster_cols=FALSE, border_color=NA,
             color=brewer.pal(9,"Greys")[4:9])
    dev.off()
    ''')
    #color=colorRampPalette(brewer.pal(9,"Blues")[5:9])(5))

def read_aa_vals(in_handle, max_freq):
    vals = {}
    for rec in vcf.VCFReader(in_handle):
        if rec.INFO.has_key("AA_AF"):
            af = rec.INFO["AA_AF"]
            change = rec.INFO["AA_CHANGE"]
            if isinstance(af, list):
                af = af[0]
                change = change.split(",")[0]
            change = change.split("_")[0]
            vals[int(change[1:-1])] = min(max_freq, float(af))
    return vals

if __name__ == "__main__":
    main(*sys.argv[1:])
    
