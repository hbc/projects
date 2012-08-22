from rkinf.toolbox import blastn, fasta
from bcbio import utils
from rkinf.cluster import start_cluster, stop_cluster
from rkinf.log import logger, setup_logging
from functools import partial
from rkinf.utils import build_results_dir, append_stem, flatten
from Bio import SeqIO
import yaml
import sys
import csv
import os
import pandas as pd


@utils.memoize_outfile("-combined.tsv")
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


def _make_length_dict(input_files):
    """ make a dictionary of the start and end nucleotides on query
     sequences that make up the longest match to any of the matched
     databases """
    d = {}

    for f in input_files:
        with open(f) as in_handle:
            reader = csv.reader(in_handle, delimiter="\t")
            header = reader.next()
            for line in reader:
                linedict = dict(zip(header, line))
                if linedict["qseqid"] not in d:
                    d[linedict["qseqid"]] = linedict
                else:
                    curr = d[linedict["qseqid"]]
                    if(abs(int(curr["qend"]) - int(curr["qstart"])) <
                       abs(int(linedict["qend"]) - int(linedict["qstart"]))):
                        d[linedict["qseqid"]] = linedict
    return d


def _trim(fasta_file, input_files):
    # make dictionary of all sequence ids in them
    # if a sequence file isnt in the dictionary add it and its begin/end/length
    # if it is already and the length in the dictionary is shorter than the
    d = _make_length_dict(input_files)

    def trim_function(x, d):
        new_seq = x
        entry = d[x.id]
        start = int(entry["qstart"])
        end = int(entry["qend"])
        new_seq.seq = x.seq[start:end]
        return new_seq

    out_file = append_stem(fasta_file, "trimmed")
    out_handle = open(out_file, "w")

    def output_writer(x):
        return SeqIO.write(x, out_handle, "fasta")

    map(output_writer,
        fasta.apply_seqio(fasta_file, partial(trim_function, d=d), "fasta"))

    return out_file


def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    setup_logging(config)
    start_cluster(config)

    # after the cluster is up, import a view to it
    from rkinf.cluster import view
    in_file = config.get("query")

    # de-parallelize for now
    blast_results = []
    for stage in config["run"]:
        if config["stage"][stage]["program"] == "blastn":
            blastn_config = config["stage"][stage]
            blast_results = [blastn.run(in_file, ref, blastn_config, config) for
                             ref in config["refs"]]

    for identity in config["min_identity"]:
        filtered_results = []
        for blast_result in blast_results:
            filtered_results.append(blastn.filter_results_by_length(
                blast_result, identity))

        fasta_hits = set()
        for filtered_result in filtered_results:
            fasta_hits.update(blastn.get_id_of_hits(filtered_result))

        def in_set_predicate(x):
            return x.id in fasta_hits

        outfile = os.path.join(build_results_dir(blastn_config, config),
                               append_stem(os.path.basename(in_file),
                                           str(identity) + "_filt"))

        fasta_filtered = fasta.filter_fasta(in_file,
                                            in_set_predicate,
                                            outfile)

        trimmed = _trim(fasta_filtered, filtered_results)
        org_names = [x["name"] for x in config["refs"]]
        logger.info(trimmed)
        logger.info(filtered_results)
        logger.info(org_names)
        combined = _make_combined_csv(trimmed, filtered_results, org_names)

    stop_cluster()

if __name__ == "__main__":
    main(*sys.argv[1:])
