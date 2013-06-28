//VARIABLES
RESULTSDIR="/n/scratch00/hsph/projects/kc_rrbs/results"
QUALITY=30
REFERENCEGENOME="/n/scratch00/hsph/biodata/genomes/Mmusculus/mm9/iGenomes/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa"
TMPDIR="/n/scratch00/hsph/tmp"


// ANALYSES


// Align
@Transform("sam")
rrbsmap = {
exec 	"""
		bsmap -a $input -d $REFERENCEGENOME -D C-CGG -o $output
		"""	
}

// transform sam to bam 
@Transform("bam")
sam_to_bam = {
exec 	"""
		java -Xmx2g -Djava.io.tmpdir=${TMPDIR} -jar /n/HSPH/local/share/java/picard/SamFormatConverter.jar INPUT=$input OUTPUT=$output		
		"""
}

// sort bam 
@Filter("coordsorted")
sort = {
exec 	"""
		java -Xmx2g -Djava.io.tmpdir=${TMPDIR} -jar /n/HSPH/local/share/java/picard/SortSam.jar INPUT=$input OUTPUT=$output SORT_ORDER=coordinate
		"""
}


//
@Transform("bam.bai")
index_bam = {
exec	"""
		samtools index $input
		"""
		forward input
}


Bpipe.run {"S%.fastq" * [rrbsmap + sam_to_bam + sort + index_bam]}
 
