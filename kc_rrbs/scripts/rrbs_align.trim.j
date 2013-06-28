//VARIABLES
RESULTSDIR="/n/scratch00/hsph/projects/kc_rrbs/results"
DATADIR="/n/scratch00/hsph/projects/kc_rrbs/data/fastq/trimmed/"
QUALITY=30
REFERENCEGENOMEDIR="/n/scratch00/hsph/biodata/genomes/Mmusculus/mm9/bismark/UCSC/"
TMPDIR="/n/scratch00/hsph/tmp"
QUANTMETHSCRIPT="/n/scratch00/hsph/projects/kc_rrbs/scripts/quant_meth_methylkit.r"

// ANALYSES

// Trim and FastQC
@Transform("trimmed.fq")
trim_galore = {
exec 	"""
	trim_galore --rrbs --fastqc --fastqc_args "--outdir ${RESULTSDIR}/fastQC" --quality ${QUALITY} $input
	"""
}


// Align
@Transform("fq_bismark.sam")
bismarkalign = {
exec 	"""
	bismark -n 1 -l 50 $REFERENCEGENOMEDIR $input
	"""	
}


// sort sam 
@Filter("coordsorted")
sortsam = {
exec 	"""
		java -Xmx2g -Djava.io.tmpdir=${TMPDIR} -jar /n/HSPH/local/share/java/picard/SortSam.jar INPUT=$input OUTPUT=$output SORT_ORDER=coordinate
		"""
}


//quantitate methylation with methylkit, sam files will be parsed and CpG C/T conversions counted for each individual sample
@Transform("_CpG.txt")
quantmeth = {
exec	"""
		$QUANTMETHSCRIPT $input $DATADIR
		"""
}


Bpipe.run {"S%.fastq" * [trim_galore + bismarkalign + sortsam + quantmeth]}
