"""Combine and rank multiple differential expression experiments.

Provides logic to organize a set of top differential expression results
into a combined ranked list based on:

- Number of times differential expression identifies a gene
- Ranking within experiments
"""
import os
import csv
import collections

from bcbio.utils import file_exists

def identify_top_ranked(exp_info, config):
    out_file = "{0}diffexp-top.tsv".format(os.path.commonprefix([x[1] for x in exp_info]))
    print out_file
    if not file_exists(out_file):
        gene_info = _collect_gene_info(exp_info)
        top_genes = _organize_gene_list(gene_info["counts"], gene_info["exps"], config)
        with open(out_file, "w") as out_handle:
            _write_top_genes(out_handle, top_genes, gene_info["exps"],
                             gene_info["targets"], gene_info["attributes"])
    return out_file

def _collect_gene_info(exp_info):
    """Consolidate differential expression results by genes identified.
    """
    exps = collections.defaultdict(list)
    attributes = {}
    counts = collections.defaultdict(list)
    targets = collections.defaultdict(list)
    for exp_name, fname in exp_info:
        for i, (gene, target, attrs) in enumerate(_read_annotated_file(fname)):
            attributes[gene] = attrs
            exps[gene].append(exp_name)
            counts[gene].append(i)
            targets[gene].append(target)
    return {"exps": exps, "attributes": attributes,
            "counts": counts, "targets": targets}

def _read_annotated_file(fname):
    """Iterator over annotated differential expression targets.
    """
    with open(fname) as in_handle:
        reader = csv.reader(in_handle, dialect="excel-tab")
        reader.next() # header
        for parts in reader:
            target = parts[0]
            gene = parts[5]
            attrs = parts[6:]
            yield gene, target, attrs

def _organize_gene_list(counts, exps, config):
    """Order genes by priority, emphasizing multiple experiments and rankings.
    """
    to_rank = []
    for gene, positions in counts.iteritems():
        # Rank by:
        # 1. Number of times a gene seen
        # 2. Number of experiments a gene is present in
        # 3. Best position of the gene in rankings
        to_rank.append((len(positions), len(set(exps[gene])), -min(positions), gene))
    to_rank.sort(reverse=True)
    for show in to_rank[:10]:
        print show
    return [x[-1] for x in to_rank[:int(config["algorithm"]["final_targets"])]]

def _write_top_genes(out_handle, genes, exps, targets, attributes):
    """Output tab delimited file of top ranked genes.
    """
    writer = csv.writer(out_handle, dialect="excel-tab")
    writer.writerow(["gene", "biotype", "experiments", "targets",
                     "genesymbol", "description"])
    for gene in genes:
        biotype, name, description = attributes[gene]
        cur_targets = ";".join(sorted(list(set(targets[gene]))))
        cur_exps = ";".join(sorted(list(set(exps[gene]))))
        writer.writerow([gene, biotype, cur_exps, cur_targets,
                         name, description])
