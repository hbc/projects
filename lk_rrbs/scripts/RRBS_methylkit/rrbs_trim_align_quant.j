////////////////////////////////////////////////////////////////////////////////////////////////
//VARIABLES

//DIRECTORIES
BASEDIR="/n/home08/jhutchin/projects/lk_rrbs/data/temp" //the directory with the fastq files you would like to process
TMPDIR="/n/scratch00/hsph/tmp"
SCRIPTDIR="/n/home08/jhutchin/projects/lk_rrbs/scripts/RRBS_methylkit" //directory where you have place the scripts
PICARDDIR="/n/HSPH/local/share/java/picard" //directory where the Picard tools are located

//TRIM VARIABLES
QUALITY=30 //trim bases with phred quality scores lower than this
ADAPTER="AGATCGGAAGAGC" //adapter to trim, if unknown, use the first 13bp of the Illumina adapter 'AGATCGGAAGAGC' and check the FASTQC overrepresented sequences for adapters to verify
MINTRIMMEDLENGTH=25

//BISMARK ALIGNER VARIABLES
BUILD="hg19" //genome build
DIRECTIONVAR="non_directional" //options are directional or non_directional
REFERENCEGENOMEDIR="/n/hsphS10/hsphfs1/chb/biodata/genomes/Hsapiens/hg19/bismark/UCSC" //bismark prepared genome

//METHYLKIT CpG QUANTITATION VARIABLES
MINIMUMCOVERAGE=10 //minimum read coverage to call a methylation status for a base
MINIMUMQUALITY=20 //minimum phred quality score to call a methylation status for a base

if ( DIRECTIONVAR=="non_directional") {
    TRIM_GALORE_DIRECTIONVAR="--non_directional"
} else {
    TRIM_GALORE_DIRECTIONVAR=""
}



////////////////////////////////////////////////////////////////////////////////////////////////
// ANALYSES

//setupdirectories
setupdirs = {
exec	"""mkdir -p ${BASEDIR}/fastqc/${input}/pretrim/"""
exec	"""mkdir -p ${BASEDIR}/fastqc/${input}/posttrim/"""
forward input
}

//run fastqc on untrimmed
fastqc = {
exec	"""
	fastqc --o ${BASEDIR}/fastqc/${input}/pretrim/ $input
	"""
forward input
}


// Trim and FastQC
@Transform("trimmed.fq")
trim_galore = {
exec 	"""
	trim_galore --rrbs ${TRIM_GALORE_DIRECTIONVAR} --fastqc --fastqc_args "--outdir ${BASEDIR}/fastqc/${input}/posttrim" --adapter ${ADAPTER} --length ${MINTRIMMEDLENGTH} --quality ${QUALITY} $input
	"""
}

// Align
@Transform("fq_bismark.sam")
bismarkalign = {
exec 	"""
	bismark -n 1 -l 50 --$DIRECTIONVAR ${REFERENCEGENOMEDIR}/ $input
	"""	
}

// sort sam 
@Filter("coordsorted")
sortsam = {
exec 	"""
		java -Xmx2g -Djava.io.tmpdir=${TMPDIR} -jar ${PICARDDIR}/SortSam.jar INPUT=$input OUTPUT=$output SORT_ORDER=coordinate
		"""
}


//quantitate methylation with methylkit, sam files will be parsed and CpG C/T conversions counted for each individual sample
quantmeth = {
exec	"""
		${SCRIPTDIR}/knitr_quant_meth_methylkit.r $input $BASEDIR $BUILD $SCRIPTDIR $MINIMUMCOVERAGE $MINIMUMQUALITY
		"""
}

compile_results = {
exec	"""
		bash ${SCRIPTDIR}/compile_results.sh $input.md
		"""
}


Bpipe.run {"%.fastq" * [setupdirs + fastqc + trim_galore + bismarkalign + sortsam + quantmeth + compile_results]}
