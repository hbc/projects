PROJDIR = "/n/scratch00/hsph/projects/js_trinity/"
TMPDIR = PROJDIR + "/tmp"
PICARD="~/opt/lib/java"
KHMER="~/opt/share/khmer/scripts"

quality_check = segment {
  fastqc.using(optargs:"-t 8")
}

oases = segment {
}

digital_normalization = segment {
  khmer_digital_normalization.using(khmer:$KHMER, optargs:"-N 4 -x 4e9")
}

kmer_histogram = segment {
}

filter_vector = segment {
  bwa_index.using(genome: "$PROJDIR/meta/bwa/genome", vector: "$PROJDIR/meta/vector_sequence_orig.fa") +
    bwa_aln.using(genome: "$PROJDIR/meta/bwa/genome", optargs: "-t 8") +
    bwa_samse.using(genome: "$PROJDIR/meta/bwa/genome") +
    samsort.using(optargs: "SORT_ORDER=queryname TMP_DIR=$TMPDIR") +
    samfilter.using(filter: "excludeAligned", optargs: "WRITE_READS_FILES=false TMP_DIR=$TMPDIR") +
    sam2fastq.using(optargs :"RC=false TMP_DIR=$TMPDIR") +
    tagdust.using(contam: "$PROJDIR/meta/vectorends.fa") +
    tagdust.using(contam: "$PROJDIR/meta/contam_second_round.fa")
}

filter_other = segment {
  tagdust.using(contam: "`cat $PROJDIR/meta/contaminants.fa`") +
    cutadapt.using(contam: "$PROJDIR/meta/for_cutadapt.fa", optargs:"-m 35 -q 20") +
    trim_polya(optargs:"-m 35 -q 20") +
    sickle(pair:"se", qtype:"illumina", optargs:"-q 20 -l 35")
}

analyze_kmer = segment {
  
}

Bpipe.run {
}
