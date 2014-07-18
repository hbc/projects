core<-read.csv(file="core_genome.txt",header=F)
rownames<-core[,1]
core<-core[2:length(core[1,])]
core<-t(core)
core.means<-NULL;core.lb<-NULL;core.ub<-NULL;for (f in 1:length(core[1,])) {core.means<-c(core.means,mean(core[,f]));core.lb<-c(core.lb,quantile(core[,f],probs=0.025,names=FALSE));core.ub<-c(core.ub,quantile(core[,f],probs=0.975,names=FALSE))}
library(Hmisc)
pdf(file="core_genome.pdf")
errbar(1:length(core.means),core.means,yplus=core.ub,yminus=core.lb,xlab="Number of sampled genomes",ylab="Number of core CDSs")
y<-core.means
x<-1:length(core[1,])
core.model<-nls(y~K*exp(-x/T)+O,start=c(K=2000,T=2,O=1700),control = list(maxiter = 10000))
core.modelsum<-summary(core.model)
core.modelsum
try(confint(core.model))
points(x,core.modelsum$parameters[1,1]*exp(-x/core.modelsum$parameters[2,1])+core.modelsum$parameters[3,1],type="l",col="red")
dev.off();

shared<-read.csv(file="shared_genome.txt",header=F)
rownames(shared)<-shared[,1]
shared<-shared[2:length(shared[1,])]
shared<-t(shared)
shared.means<-NULL;shared.lb<-NULL;shared.ub<-NULL;for (f in 1:length(shared[1,])) {shared.means<-c(shared.means,mean(shared[,f]));shared.lb<-c(shared.lb,quantile(shared[,f],probs=0.025,names=FALSE));shared.ub<-c(shared.ub,quantile(shared[,f],probs=0.975,names=FALSE))}
pdf(file="shared_genes.pdf")
errbar(1:length(shared.means),shared.means,yplus=shared.ub,yminus=shared.lb,xlab="Number of sampled genomes",ylab="Number of shared CDSs")
y<-shared.means
x<-1:length(shared[1,])
shared.model<-nls(y~K*exp(-x/T)+O,start=c(K=2000,T=2,O=1700),control = list(maxiter = 100))
shared.modelsum<-summary(shared.model)
shared.modelsum
try(confint(shared.model))
points(x,shared.modelsum$parameters[1,1]*exp(-x/shared.modelsum$parameters[2,1])+shared.modelsum$parameters[3,1],type="l",col="green")
dev.off();
specific<-read.csv(file="strain_specific.txt",header=F)
rownames(specific)<-specific[,1]
specific<-specific[2:length(specific[1,])]
specific<-t(specific)
specific.means<-NULL;specific.lb<-NULL;specific.ub<-NULL;for (f in 1:length(specific[1,])) {specific.means<-c(specific.means,mean(specific[,f]));specific.lb<-c(specific.lb,quantile(specific[,f],probs=0.025,names=FALSE));specific.ub<-c(specific.ub,quantile(specific[,f],probs=0.975,names=FALSE))}
pdf(file="strain_specific_genes.pdf")
errbar(1:length(specific.means),specific.means,yplus=specific.ub,yminus=specific.lb,xlab="Number of sampled genomes",ylab="Number of strain specific CDSs")
y<-specific.means[2:length(specific.means)]
x<-2:length(specific[1,])
specific.model<-nls(y~K*exp(-x/T)+O,start=c(K=20,T=50,O=1),control = list(maxiter = 100))
specific.modelsum<-summary(specific.model)
specific.modelsum
try(confint(specific.model))
points(x,specific.modelsum$parameters[1,1]*exp(-x/specific.modelsum$parameters[2,1])+specific.modelsum$parameters[3,1],type="l",col="blue")
dev.off();
pan<-read.csv(file="pan_genome.txt",header=F)
rownames(pan)<-pan[,1]
pan<-pan[2:length(pan[1,])]
pan<-t(pan)
pan.means<-NULL;pan.lb<-NULL;pan.ub<-NULL;for (f in 1:length(pan[1,])) {pan.means<-c(pan.means,mean(pan[,f]));pan.lb<-c(pan.lb,quantile(pan[,f],probs=0.025,names=FALSE));pan.ub<-c(pan.ub,quantile(pan[,f],probs=0.975,names=FALSE))}
pdf(file="pan_genome.pdf")
errbar(1:length(pan.means),pan.means,yplus=pan.ub,yminus=pan.lb,xlab="Number of sampled genomes",ylab="Total number of observed CDSs")
y<-pan.means
x<-1:length(pan[1,])
points(x,pan.means[1]+specific.modelsum$parameters[3,1]*(x-1)+specific.modelsum$parameters[1,1]*exp(-2/specific.modelsum$parameters[2,1])*(1-exp(-(x-1)/specific.modelsum$parameters[2,1]))/(1-exp(-1/specific.modelsum$parameters[2,1])),col="red",type="l")
dev.off();
