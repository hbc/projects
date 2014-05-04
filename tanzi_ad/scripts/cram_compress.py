import glob

from cluster_helper.cluster import cluster_view

def cram_compress(in_bam):
    import os
    import subprocess
    from bcbio import utils
    from bcbio.distributed.transaction import file_transaction
    print in_bam
    ref_file = "/n/hsphS10/hsphfs1/chb/biodata/genomes/Hsapiens/GRCh37/seq/GRCh37.fa"
    out_file = "%s.cram" % os.path.splitext(in_bam)[0]
    jvm_opts = "-Xms3g -Xmx6g"
    if not utils.file_exists(out_file):
        print "cramming", out_file
        with file_transaction(out_file) as tx_out_file:
            cmd = ("cramtools {jvm_opts} cram "
                   "--input-bam-file {in_bam} "
                   "--capture-all-tags "
                   "--ignore-tags 'BD:BI' "
                   "--reference-fasta-file {ref_file} "
                   "--lossy-quality-score-spec '*8' "
                   "--output-cram-file {tx_out_file}")
            subprocess.check_call(cmd.format(**locals()), shell=True)
    if not utils.file_exists(out_file + ".crai"):
        print "indexing", out_file + ".crai"
        with file_transaction(out_file + ".crai") as tx_out_file:
            tx_in_file = os.path.splitext(tx_out_file)[0]
            utils.symlink_plus(out_file, tx_in_file)
            cmd = ("cramtools {jvm_opts} index "
                   "--input-file {tx_in_file}")
            subprocess.check_call(cmd.format(**locals()), shell=True)
    if os.path.exists(in_bam) and utils.file_exists(out_file):
       os.remove(in_bam)
    return out_file

def main():
    files = (sorted(glob.glob("/n/hsphS10/hsphfs2/tanzi_recalled/*alz-*/final/*/*-ready.bam")) +
             sorted(glob.glob("/n/hsphS10/hsphfs2/tanzi_recalled/*alz-*/final/*/*-ready.cram")))
    with cluster_view(scheduler="lsf", queue="hsph_bioinfo", num_jobs=20,
                      cores_per_job=6, extra_params={"mem": "7"}) as view:
        view.map(cram_compress, files)

if __name__ == "__main__":
    main()
