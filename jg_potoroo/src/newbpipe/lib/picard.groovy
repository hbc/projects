PICARD="~opt/lib/java"

@Filter("sort")
samsort = {
  exec "java -jar $PICARD/SortSam.jar INPUT=$input OUTPUT=$output $optargs"
}

@Filter("filt")
samfilter = {
  exec """java -jar $PICARD/FilterSamReads.jar INPUT=$input OUTPUT=$output
FILTER=$filter $optargs"""
}

@Transform("fastq")
sam2fastq = {
  exec """java -jar $PICARD/SamToFastq.jar INPUT=$input FASTQ=$output $optargs"""
}
