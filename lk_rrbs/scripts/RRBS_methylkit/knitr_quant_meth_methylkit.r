#!/n/sw/R-2.15.0/bin/Rscript

library(knitr)
args=commandArgs()
print(args)
samfilename=args[6] # coordinate sorted sam file from bismark alignment, passed at bpipe command line, specified in bpipe variables
sampleID=sub(".trimmed.fq_bismark.coordsorted.sam", "", basename(samfilename))
dataDir=args[7] # directory with the sam file, passed at bpipe command line, specified in bpipe variables
methquantDir=paste(dataDir, "methylation_quantitation_results", sep="/")
build=args[8] # genomic build, passed at bpipe command line, specified in bpipe variables
scriptDir=args[9] # directory containing this script and the R script to knit into markdown document, passed at bpipe command line, specified in bpipe variables
minimumcoverage=args[10]
minimumquality=args[11]

script_to_knit=file.path(scriptDir, "quant_meth_methylkit.rmd")
markdown_output=file.path(dataDir, "quant_meth_methylkit.md")
knit(input=script_to_knit, output=paste(sampleID, "quant_meth_methlykit.md", sep="_"), envir=environment())


q()
