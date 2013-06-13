if (file.exists("/n/hsphS10/hsphfs1/chb/projects/lk_FOY")) {
  baseDir <- "/n/hsphS10/hsphfs1/chb/projects/lk_FOY"
}  else if (file.exists("/Volumes/home08/jhutchin/consults/lk_FOY/")) {
  baseDir <- "/Volumes/home08/jhutchin/consults/lk_FOY"
} else {
  baseDir <- "/Volumes/ody/consults/lk_FOY"
}
dataDir <- file.path(baseDir, "data")
resultsDir <- file.path(baseDir, "results")
metaDir <- file.path(baseDir, "meta")

load(file.path(resultsDir, "RDATA.raw_and_normalized_microarray.data.PBMC.U133Plus2.0"))
source("http://dl.dropboxusercontent.com/u/4253254/Resources/functions.r")

PCAplot(mic.norm.eset, categories="study", colorpalette=cbPalette, title="Uncorrected, studies", alpha=0.6)
PCAplot(mic.norm.eset, categories="stage", colorpalette=cbPalette, title="Uncorrected, stages", alpha=0.6)

edata <- exprs(mic.norm.eset)
pd <- pData(mic.norm.eset)
batch <- unlist(pd$study)
mod = model.matrix(~as.factor(stage), data=pd)
combat_edata = ComBat(dat=edata, batch=batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)
mic.norm.combat <- mic.norm
exprs(mic.norm.combat) <- combat_edata
PCAplot(mic.norm.combat, categories="study", colorpalette=cbPalette)