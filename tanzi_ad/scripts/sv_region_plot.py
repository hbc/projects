"""Prepare coverage plots of regions of interest from structural variant calls.
"""
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import pybedtools
import seaborn as sns
import yaml

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do

def main(toplot_file, final_dir):
    with open(toplot_file) as in_handle:
        toplot = yaml.safe_load(in_handle)
    for cur in toplot:
        out_file = os.path.join(os.path.dirname(toplot_file), "%s_%s_%s-coverage.pdf" % (
            cur["chrom"], cur["start"], cur["end"]))
        if not utils.file_exists(out_file):
            with file_transaction({}, out_file) as tx_out_file:
                with PdfPages(tx_out_file) as pdf_out:
                    sns.despine()
                    for batch, batch_items in sorted(cur["batches"].items()):
                        df = None
                        n_un = 0
                        n_a = 0
                        for descr, status in sorted(batch_items.items(),
                                                    key=lambda xs: (xs[1], xs[0])):
                            if status == "affected":
                                n_a += 1
                            else:
                                n_un += 1
                            name = "%s:%s" % (descr, status)
                            cur_df = calc_coverage_df(descr, cur["chrom"], cur["start"], cur["end"],
                                                      final_dir, out_file, name)
                            if df is None:
                                df = cur_df
                            else:
                                df = pd.merge(df, cur_df, on="position")
                        a_pal = sns.color_palette("coolwarm", 2 * n_a + 1)[-n_a:]
                        un_pal = sns.color_palette("coolwarm", 2 * n_un + 1)[:n_un]
                        sns.set_palette(a_pal + un_pal)
                        plot = df.plot(x="position", kind="line")
                        title = "%s [%s:%s-%s (%.1fkb)]" % (batch,
                                                            cur["chrom"], cur["start"], cur["end"],
                                                            (cur["end"] - cur["start"]) / 1000.0)
                        plot.set_title(title)
                        plot.set_ylim(0)
                        plot.get_xaxis().set_major_formatter(
                            matplotlib.ticker.FuncFormatter(lambda x, p: '{:,}'.format(int(x))))
                        plot.set_xlabel("")
                        plot.set_ylabel("coverage")
                        plot.get_xaxis().set_major_locator(matplotlib.ticker.MaxNLocator(6))
                        pdf_out.savefig(plot.get_figure())
                        plt.close()

def calc_coverage_df(descr, chrom, start, end, final_dir, out_base, name):
    pad = 100
    cov_file = "%s-%s.bed" % (os.path.splitext(out_base)[0], descr)
    coord_file = os.path.join(os.path.dirname(cov_file), "coords-%s_%s_%s.bed" % (chrom, start, end))
    with open(coord_file, "w") as out_handle:
        out_handle.write("%s\t%s\t%s\n" % (chrom, start - pad, end + pad))
    if not utils.file_exists(cov_file):
        bam_file = os.path.join(final_dir, descr, "%s-ready.bam" % descr)
        coords = "%s:%s-%s" % (chrom, start - pad, end + pad)
        cmd = "samtools view -b {bam_file} {coords} | bedtools coverage -abam - -b {coord_file} -d > {cov_file}"
        do.run(cmd.format(**locals()), "calculate coverage")
    print cov_file
    assert os.path.exists(cov_file)
    with open(cov_file) as in_handle:
        df = pd.read_table(in_handle, header=None, names=["chrom", "start", "end", "offset", name])
        df["position"] = df["start"] + df["offset"]
    return df[["position", name]]

if __name__ == "__main__":
    main(*sys.argv[1:])
