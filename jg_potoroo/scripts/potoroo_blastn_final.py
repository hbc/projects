from rkinf.toolbox import blastn, fasta
from bcbio import utils
from bcbio.utils import file_exists
from rkinf.cluster import start_cluster, stop_cluster
from rkinf.log import setup_logging, logger
from functools import partial
from rkinf.utils import build_results_dir, append_stem, flatten
from Bio import SeqIO
import yaml
import sys
import csv
import os
import pandas as pd
import itertools
import rpy2.robjects as robjects
import rpy2.robjects.vectors as vectors
from IPython.parallel import require

TO_KEEP = ["qseqid", "sseqid", "evalue", "length", "pident", "sstart",
           "ssend"]

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    setup_logging(config)
    start_cluster(config)

    # after the cluster is up, import a view to it
    from rkinf.cluster import view
    in_file = config.get("query")
    org_names = [x["name"] for x in config["refs"]]

    curr_files = in_file

    for stage in config["run"]:
        if stage == "blastn":
            logger.info("Running %s on %s." % (stage, curr_files))
            blastn_config = config["stage"][stage]
            refs = config["refs"]
            args = zip(*itertools.product([curr_files], refs,
                                          [blastn_config], [config]))
            blastn_results = view.map(blastn.run, *args)
            curr_files = blastn_results

        if stage == "annotate":
            logger.info("Running %s on %s." % (stage, curr_files))
            # annotate the data frames
            args = zip(*itertools.product(curr_files, ["sseqid"],
                                          org_names))
            annotated = view.map(_annotate_df, *args)
            curr_files = annotated

        if stage == "combine":
            out_fname = os.path.join(os.path.dirname(curr_files[0]),
                                                     append_stem(in_file,
                                                                 "combined"))
            logger.info("Combining %s into %s." % (curr_files, out_fname))
            org_names = [x["name"] for x in config["refs"]]
            #       combined = _make_combined_csv(curr_files, org_names, out_fname)

    stop_cluster()


def _make_combined_csv(fasta_file, input_files, org_names, out_file=None):
    """
    takes a list of output files from blastn and a fasta file and attaches
    the columns jagesh wants all into one big data frame and writes it out
    """
    suffixes = map(lambda x: "_" + x, org_names)

    # columns to keep
    TO_KEEP = list(flatten(["qseqid", map(lambda x: "sseqid" + x, suffixes),
                            map(lambda x: "evalue" + x, suffixes),
                            map(lambda x: "length" + x, suffixes),
                            map(lambda x: "pident" + x, suffixes),
                            map(lambda x: "sstart" + x, suffixes),
                            map(lambda x: "send" + x, suffixes)]))

    # read inputs as tables and merge into one big table
    inputs = zip(suffixes, map(pd.read_table, input_files))
    dfs = [inp[1].rename(columns=lambda name: name + inp[0]) for inp in inputs]
    d = {}
    for x in suffixes:
        d["qseqid" + x] = "qseqid"
    renamed = [x.rename(columns=d) for x in dfs]

    merged = reduce(lambda x, y: pd.merge(x, y, on="qseqid"),
                    renamed[1:], renamed[0])

    df_subset = merged[TO_KEEP]

    # add the sequence from the fasta file
    seqs = pd.DataFrame(fasta.apply_seqio(
        fasta_file, lambda x: {'qseqid': x.id, 'seq': str(x.seq)},
        "fasta"))

    merged = pd.merge(df_subset, seqs, on="qseqid")
    merged.to_csv(out_file, index=False, sep="\t")
    return out_file


@require(os)
def _annotate_df(in_file, join_column, organism, out_file=None):
    from rkinf.log import logger
    from rkinf.utils import append_stem
    from rpy2 import robjects
    ORG_TO_ENSEMBL = {"opossum": {"gene_ensembl": "mdomestica_gene_ensembl",
                                  "gene_symbol": "hgnc_symbol"},
                      "mouse": {"gene_ensembl": "mmusculus_gene_ensembl",
                                "gene_symbol": "mgi_symbol"},
                      "human": {"gene_ensembl": "hsapiens_gene_ensembl",
                                "gene_symbol": "hgnc_symbol"},
                      "taz": {"gene_ensembl": "sharrisii_gene_ensembl",
                              "gene_symbol": "hgnc_symbol"}}

    if organism not in ORG_TO_ENSEMBL:
        logger.error("organism not supported")
        exit(1)

    logger.info("Annotating %s." % (organism))
    if not out_file:
        out_file = append_stem(in_file, "annotated")
    if os.path.exists(out_file):
        return out_file
    # use biomaRt to annotate the data file
    r = robjects.r
    r.assign('join_column', join_column)
    r.assign('in_file', in_file)
    r.assign('out_file', out_file)
    r.assign('ensembl_gene', ORG_TO_ENSEMBL[organism]["gene_ensembl"])
    r.assign('gene_symbol', ORG_TO_ENSEMBL[organism]["gene_symbol"])
    r('''
    library(biomaRt)
    ensembl = useMart("ensembl", dataset = ensembl_gene)
    d = read.table(in_file, header=TRUE)
    a = getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id",
                gene_symbol, "description"),
                filters=c("ensembl_transcript_id"), values=d[,join_column],
                mart=ensembl)
    m = merge(d, a, by.x=join_column, by.y="ensembl_transcript_id")
    write.table(m, out_file, quote=FALSE, row.names=FALSE, sep="\t")
    ''')

    return out_file


def _combine_and_write(dataframes, out_file):
    from rkinf.log import logger
    import pandas as pd
    from bcbio.utils import file_exists
    logger.info("Writing combined file to %s." % (out_file))
    if file_exists(out_file):
        return out_file
    merged = pd.concat(dataframes)
    df_subset = merged[TO_KEEP]
    df_subset.to_csv(out_file, index=False, sep="\t")
    return out_file


def _merge_dataframes(dataframes, on):
    import pandas as pd
    """ returns a single dataframe of a list of dataframes merged on
    a column """
    return reduce(lambda x, y: pd.merge(x, y, on=on),
                  dataframes[1:], dataframes[0])


if __name__ == "__main__":
    main(*sys.argv[1:])
