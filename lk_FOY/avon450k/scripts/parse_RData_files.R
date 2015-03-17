setwd("~/projects/lk_FOY/avon450K/data/ARIES_data_B1361")

load("B1361.sampleData.releasev2_apr2014.Rdata")
write.csv(sampleDataForUser, "B1361.sampleData.releasev2_apr2014.csv")


load("B1361.pvalues.F7releasev2_apr2014.Rdata")
pvals.7 <- pvals
subset.sampleDataForUser <- sampleDataForUser[match(colnames(pvals.7),sampleDataForUser$collaboratorID),]
pvals.7.m <- pvals.7[,as.character(subset(subset.sampleDataForUser, gender=="m")$collaboratorID)]
pvals.7.f <- pvals.7[,as.character(subset(subset.sampleDataForUser, gender=="f")$collaboratorID)]
nums <- seq(1,nrow(pvals.7),150000)
nums <- c(nums, nrow(pvals.7))
for (x in 1:(length(nums)-1)){
  write.csv(pvals.7[nums[x]:nums[x+1],], file=paste("B1361.pvalues.F7releasev2_apr2014_", nums[x], "-", nums[x+1], ".csv", sep=""))
  write.csv(pvals.7.m[nums[x]:nums[x+1],], file=paste("B1361.pvalues.male.F7releasev2_apr2014_", nums[x], "-", nums[x+1], ".csv", sep=""))
  write.csv(pvals.7.f[nums[x]:nums[x+1],], file=paste("B1361.pvalues.female.F7releasev2_apr2014_", nums[x], "-", nums[x+1], ".csv", sep=""))
}


load("B1361.pvalues.15upreleasev2_apr2014.Rdata")
pvals.15 <- pvals
subset.sampleDataForUser <- sampleDataForUser[match(colnames(pvals.15),sampleDataForUser$collaboratorID),]
pvals.15.m <- pvals.15[,as.character(subset(subset.sampleDataForUser, gender=="m")$collaboratorID)]
pvals.15.f <- pvals.15[,as.character(subset(subset.sampleDataForUser, gender=="f")$collaboratorID)]
nums <- seq(1,nrow(pvals.15),150000)
nums <- c(nums, nrow(pvals.15))
for (x in 1:(length(nums)-1)){
  write.csv(pvals.15[nums[x]:nums[x+1],], file=paste("B1361.pvalues.15upreleasev2_apr2014_", nums[x], "-", nums[x+1], ".csv", sep=""))
  write.csv(pvals.15.m[nums[x]:nums[x+1],], file=paste("B1361.pvalues.male.15upreleasev2_apr2014_", nums[x], "-", nums[x+1], ".csv", sep=""))
  write.csv(pvals.15.f[nums[x]:nums[x+1],], file=paste("B1361.pvalues.female.15upreleasev2_apr2014_", nums[x], "-", nums[x+1], ".csv", sep=""))
}


load("B1361.rawBetas.F7releasev2_apr2014.Rdata")
betas.7 <- betas
subset.sampleDataForUser <- sampleDataForUser[match(colnames(betas.7),sampleDataForUser$collaboratorID),]
betas.7.m <- betas.7[,as.character(subset(subset.sampleDataForUser, gender=="m")$collaboratorID)]
betas.7.f <- betas.7[,as.character(subset(subset.sampleDataForUser, gender=="f")$collaboratorID)]
nums <- seq(1,nrow(betas.7),25000)
nums <- c(nums, nrow(betas.7))
for (x in 1:(length(nums)-1)){
  write.csv(betas.7[nums[x]:nums[x+1],], file=paste("B1361.rawBetas.F7releasev2_apr2014_", nums[x], "-", nums[x+1], ".csv", sep=""))
  write.csv(betas.7.m[nums[x]:nums[x+1],], file=paste("B1361.rawBetas.male.F7releasev2_apr2014_", nums[x], "-", nums[x+1], ".csv", sep=""))
  write.csv(betas.7.f[nums[x]:nums[x+1],], file=paste("B1361.rawBetas.female.F7releasev2_apr2014_", nums[x], "-", nums[x+1], ".csv", sep=""))
}


load("B1361.rawBetas.15upreleasev2_apr2014.Rdata")
betas.15 <- betas
subset.sampleDataForUser <- sampleDataForUser[match(colnames(betas.15),sampleDataForUser$collaboratorID),]
betas.15.m <- betas.15[,as.character(subset(subset.sampleDataForUser, gender=="m")$collaboratorID)]
betas.15.f <- betas.15[,as.character(subset(subset.sampleDataForUser, gender=="f")$collaboratorID)]
nums <- seq(1,nrow(betas.15),25000)
nums <- c(nums, nrow(betas.15))
for (x in 1:(length(nums)-1)){
  write.csv(betas.15[nums[x]:nums[x+1],], file=paste("B1361.rawBetas.15upreleasev2_apr2014_", nums[x], "-", nums[x+1], ".csv", sep=""))
  write.csv(betas.15.m[nums[x]:nums[x+1],], file=paste("B1361.rawBetas.male.15upreleasev2_apr2014_", nums[x], "-", nums[x+1], ".csv", sep=""))
  write.csv(betas.15.f[nums[x]:nums[x+1],], file=paste("B1361.rawBetas.female.15upreleasev2_apr2014_", nums[x], "-", nums[x+1], ".csv", sep=""))
}





