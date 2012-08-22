#!/bin/bash
PROJID="js_trinity"
PROJDIR="/n/scratch00/hsph/projects/js_trinity"
BAIT="$PROJDIR/meta/vector_sequence.fa"
SCRAM="$PROJDIR/meta/vector_sequence_scrambled.fa"
VECTOR="$PROJDIR/meta/vector_sequence_orig.fa"

BWA_DISABLE=1
TAGDUST_DISABLE=1
CUTADAPT_DISABLE=1

if [ -z $BWA_DISABLE ]
then
    # make a genome file for bwa to use out of the vector
    bwa index -p $PROJDIR/meta/bwa/vector -a is $VECTOR
    echo "Mapping reads to vector to identify contaminants."
    for fname in "$@"
    do
	prefix=`basename $fname .fastq`
	bsub -P $PROJID -q hsph -o %J.out -e %J.err -J bwa "bwa aln -t 4 $PROJDIR/meta/bwa/vector $fname > ${prefix}_mapped_to_vector.bwa; bwa samse $PROJDIR/meta/bwa/vector ${prefix}_mapped_to_vector.bwa $fname > ${prefix}_mapped_to_vector.sam"
    done
fi

if [ -z $REMOVE_UNMAPPED_READS_FROM_SAM ]
then
    
    
    
    

if [ -z $TAGDUST_DISABLE ]
then
    echo "Removing likely contaminant kmers with tagdust."
    for fname in "$@"
    do
	prefix=`basename $fname .fastq`
	bsub -K -P $PROJID -o %J.out -e %J.err -q hsph -J tagdust tagdust -s -f 0.01 -o $PROJDIR/data/${prefix}.dusted.fastq $BAIT $fname &
	bsub -K -P $PROJID -o %J.out -e %J.err -q hsph -J tagdust tagdust -s -f 0.01 -o $PROJDIR/data/${prefix}.dusted.scramble.fastq $SCRAM $fname &
    done
    echo "Waiting for jobs to finish."
    wait
    echo "Jobs finished"
fi

if [ -z $CUTADAPT_DISABLE ]
then
    echo "Trimming contaminant ends of reads with cutadapt."
    CONTAM=`python cutadapt_generator.py ../meta/contaminants.fa`
    for fname in "$@"
    do
	echo "cutadapt -m 30 $CONTAM $fname"
    done
fi


echo "Processing complete."
