PROJDIR="/n/scratch00/hsph/projects/js_trinity"
//BAIT=PROJDIR + "/meta/vector_sequence.fa"
//SCRAM=PROJDIR + "../meta/vector_sequence_scrambled.fa"
PICARD="~/opt/lib/java"
VECTOR="$PROJDIR/meta/vector_sequence_orig.fa"
//CONTAM = PROJDIR + "/meta/contaminants.fa"
TMPDIR = PROJDIR + "/tmp"
CONTAM = PROJDIR + "/meta/contam_second_round.fa"

@Transform("count")
count_fastq_reads = {
  exec "grep '@' $input | wc -l > $output"
  forward input
} 

@Filter("cut")
filter_contams_with_cutadapt = {
  doc title: "Filter ends of reads with cutadapt.",
desc: """Perform another filtering step using cutadapt, this will trim the
ends of the reads that may not have been matched during BWA mapping and
Tagdust mapping. For the vector sequences it will search both ends of
the read. For Illumina adaptor sequences and Nextera Transposon adaptor
sequences that could exist on the 3' end due to read-through, only the 3'
end is looked at.

Right now the contaminant file is a little ad-hoc, fiddled with the
sequences until it seemed to be filtering out most of the reads
that had vector contamination at one or both ends. This is an area of
improvement.

Removes reads that, after trimming, are < 35 bases. Also trims off
ends of reads where quality < 20.
""",
author: "rory.kirchner@gmail.com"

  exec "cutadapt `cat $PROJDIR/meta/for_cutadapt.fa` -m 35 -q 20 $input > $output" 
}

@Filter("nopolya")
trim_polya = {
  doc title: "Trims off poly-A from 3' ends of reads"
  exec "cutadapt -a AAAAAAAAAAAAAAAAAAAAAAAAAA -m 35 -q 20 $input > $output"
}

@Filter("sickled")
sickle_se = {
  doc title: "Trim ends of reads using sickle."
  exec "sickle se -f $input -o $output -t illumina -q 20 -l 35"
}

@Filter("dusted")
filter_contams_with_tagdust = {
  doc title: "Filter out possible contaminants with Tagdust.",
desc: """ Filters out possible contaminants listed in $CONTAM,
which includes Nextera transposon adaptors, Illumina paired end
adaptor sequences and parts of the vector. """,
author: "rory.kirchner@gmail.com"
  
  exec "tagdust -s -f 0.01 -o ${output} $CONTAM ${input}"
}

make_bwa_library_from_vector = {
  exec "bwa index -p $PROJDIR/meta/bwa/vector -a is $VECTOR"
}


@Transform("sam")
map_vector_with_bwa = {
  doc title: "Map the reads to the contaminating vector.",
desc: """ Maps reads to the vector that the sequences we are interested in
was cloned into, reads mapping to the vector will be discarded later on.""",
author: "rory.kirchner@gmail.com"  
  
  exec "bwa aln -t 4 $PROJDIR/meta/bwa/vector $input > ${input}.bwa"
  exec """
bwa samse $PROJDIR/meta/bwa/vector ${input}.bwa $input >
$output
"""
}

@Filter("sort")
sort_samfile = {
  exec """
java -jar $PICARD/SortSam.jar INPUT=$input OUTPUT=$output SORT_ORDER=queryname TMP_DIR=$TMPDIR
"""
}

@Filter("novector")
remove_aligned = {
  exec """
java -jar $PICARD/FilterSamReads.jar INPUT=$input OUTPUT=$output FILTER=excludeAligned
WRITE_READS_FILES=false TMP_DIR=TMPDIR
"""
}

@Transform("fastq")
sam_to_fastq = {
  exec """
java -jar $PICARD/SamToFastq.jar INPUT=$input FASTQ=$output RC=false TMP_DIR=TMPDIR
"""
}

sanity_check_filtering = {
  doc title: "Dirty check to see how effective the filtering steps have been.",
desc: """Takes a bunch of subsequences of the original vector and counts their
occurances in the reads, as a sanity check for how effective the filtering
has been.""",
author: "rory.kirchner@gmail.com"
  
  exec "echo '0-1k of vector, forward'"
  exec "grep -i AAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCCAATACGCAAACC $input | wc -l"
  exec "echo '1-2k of vector, forward'"
  exec "grep -i ACAATAACCCTGATAAATGCTTCAATAATATTGAAAAACGCGCGAATTGCAA $input | wc -l"
  exec "echo '2-3k of vector, forward'"
  exec "grep -i GCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAA $input | wc -l"
  exec "echo '3-4k of vector, forward'"
  exec "grep -i TGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGA $input | wc -l"
  forward input
}

Bpipe.run {
  /*"530_%_s" * [count_fastq_reads + map_vector_with_bwa + sort_samfile + remove_aligned +
    sam_to_fastq + filter_contams_with_tagdust + filter_contams_with_cutadapt +
    count_fastq_reads]*/
  //"530_%_s" * [filter_contams_with_tagdust]
  //"530_%_s" * [sickle_se]
  "530_%_s" * [trim_polya]
}
