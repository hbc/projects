args<-NULL;
args<-commandArgs(trailingOnly=TRUE);
test.data<-read.csv(file=args[1],header=F);
fisher.test(test.data);
