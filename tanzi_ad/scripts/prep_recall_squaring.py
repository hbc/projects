"""Prepare full sample set for squaring off in defined regions of interest.
"""
import glob
import gzip
import os

import yaml

def main():
    vcf_dir = "/n/hsphS10/hsphfs1/chb/projects/tanzi_ad/data/recall_variants/freebayes"
    cram_glob = "/n/hsphS10/hsphfs2/tanzi_recalled/*alz-*/final/{name}/{name}-ready.cram"
    vbed_file = "/n/hsphS10/hsphfs1/chb/projects/tanzi_ad/data/ad_loci/square/alz_square_regions.bed"
    out_file = "alz-square.yaml"
    samples = [add_cram(s, cram_glob) for s in samples_to_vcf(vcf_dir)]
    print len(samples)
    out = {"upload": {"dir": "../final"},
           "details": [sample_to_config(s, vbed_file) for s in samples]}
    with open(out_file, "w") as out_handle:
        yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)

def sample_to_config(sample, vbed_file):
    return {"description": sample["name"],
            "analysis": "variant2",
            "files": sample["cram"],
            "vrn_file": sample["vcf"],
            "genome_build": "GRCh37",
            "metadata": {"batch": "talz"},
            "algorithm": {"aligner": False, "mark_duplicates": False,
                          "recalibrate": False, "realign": False, "variantcaller": False,
                          "jointcaller": "bcbio-variation-recall", "effects": "vep",
                          "variant_regions": vbed_file},
            "resources": {"vep": {"options": ["--pick"]}}}

def add_cram(sample, cram_glob):
    cram_files = glob.glob(cram_glob.format(name= sample["name"]))
    if len(cram_files) > 1:
        redo_cram_files = [x for x in cram_files if x.find("redo") > 0]
        if redo_cram_files:
            cram_files = redo_cram_files
        else:
            cram_files = cram_files[:1]
    assert len(cram_files) == 1, (sample, cram_files)
    cram_file = cram_files[0]
    assert os.path.exists(cram_file + ".crai"), (sample, cram_file)
    sample["cram"] = cram_file
    return sample

def samples_to_vcf(vcf_dir):
    out = []
    for fname in glob.glob(os.path.join(vcf_dir, "*.vcf.gz")):
        with gzip.open(fname) as in_handle:
            for line in in_handle:
                if line.startswith("#CHROM"):
                    parts = line.split("\t")
                    format_i = parts.index("FORMAT")
                    for s in parts[format_i + 1:]:
                        out.append({"vcf": fname, "name": s.strip()})
                    break
    return out

if __name__ == "__main__":
    main()
