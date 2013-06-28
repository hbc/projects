#
# Analysis of Illumina arrays from the Cebon group
#
library(RankProd)
library(lumi)
library(lumiHumanAll.db)
library(lumiHumanIDMapping)
library(annotate)
library(pheatmap)

# The usual annotation issues. The early version of BeadStudio output 
# files use a numeric number as probe identifier, later on it uses the 
# new version of probe identifiers named as "ILMN_0000"
if (require(lumiHumanIDMapping))
  {lumiHumanIDMapping()}

# Read in array and control probes
fileName <- '../Data/iltp2565_Sample_Probe_Profile.txt'
covdesc <- '../Data/Sample_index_correct_drc.txt'
mel.lumi <- lumiR.batch(fileName, 
                        sampleInfoFile=covdesc,
                        lib.mapping='lumiHumanIDMapping',
                        convertNuID=T)

# Quick sanity check, what's the predicted chip type?
if (require(lumiHumanIDMapping))
  {getChipInfo(mel.lumi, species='Human')}
   
# Summary of data, QC Info
mel.lumi
summary(mel.lumi, 'QC')  
pData(mel.lumi)

# Preprocessing and quality control after normalization
mel.N.V <- lumiExpresso(mel.lumi, QC.evaluation=TRUE,
                        variance.stabilize=F,
                        normalize.param=list(method='vsn'))

# Summary of quality control info after processing
summary(mel.N.V, 'QC')

# Export preprocessed data for later
write.exprs(mel.N.V, file='processedMelanomaData_vsn.txt')

# Alternative probe / gene expression file using quantile normalization;
# use for now
mel.N.Q <- lumiExpresso(mel.lumi, QC.evaluation=TRUE, 
                        normalize.param=list(method='quantile'))        
summary(mel.N.Q, 'QC')
write.exprs(mel.N.Q, file='processedMelanomaData_quantiles.txt')

# QC plots
par(mfrow=c(1,1))
pdf('qc_density.pdf')
plot(mel.N.Q, what='density')
dev.off()

pdf('qc_densityCDF.pdf')
plotCDF(mel.N.Q, reverse=T)
dev.off()

pdf('qc_pairs.pdf')
pairs(mel.N.Q, smoothScatter=T, cex.labels=1)
dev.off()

pdf('qc_relations.pdf')
plot(mel.N.Q, what='sampleRelation')
dev.off()

# Check general sample relations; should be prettifed by
# using pData labels
pdf('qc_mdsRelations.pdf')
plot(mel.N.Q, what='sampleRelation', method='mds',
     color=c(1, 2, 1, 1, 2, 2, 1, 2))
dev.off()

#
# Get expression matrix and probe names, quantile version
#
dataMatrix <- exprs(mel.N.Q)
probes <- rownames(dataMatrix)
length(probes)

# Filtering based on probes being present
presentCount <- detectionCall(mel.lumi)
selDataMatrix <- dataMatrix[presentCount > 0, ]
selProbes <- rownames(selDataMatrix)
length(selProbes)

#
# Export for GSEA
#
symbols <- getSYMBOL(probes, 'lumiHumanAll.db')
export <- dataMatrix
rownames(export) <- symbols
export <- aggregate(export, list(rownames(export)), mean)
write.table(export, file='processedMelanomaData_quantiles_Symbol.txt',
            quote=F)

#
# Repeat (just the export) for VSN normalized data
#
dataMatrixVSN <- exprs(mel.N.V)
probesVSN <- rownames(dataMatrixVSN)
symbolsVSN <- getSYMBOL(probesVSN, 'lumiHumanAll.db')
exportVSN <- dataMatrixVSN
rownames(exportVSN) <- symbolsVSN
exportVSN <- aggregate(exportVSN, list(rownames(exportVSN)), mean)
write.table(exportVSN, file='processedMelanomaData_vsn_Symbol.txt',
            quote=F)

#                                        
# RankProd
#
# Comparison slow vs regular cells. Treat each cell line as coming
# from a different origin
stainCl <- c(1, 0, 1, 1, 0, 0, 1, 0)
stainOrigin <- c(1, 1, 2, 3, 3, 2, 4, 4)
stainExp.adv.out <- RPadvance(dataMatrix, stainCl,
                              stainOrigin, num.perm=100,
                              logged=T, rand=123);

# Scatter plot
pdf("stain_vs_unstained_rp.pdf");
plotRP(stainExp.adv.out, cutoff=0.1);
dev.off();
                                                          
# Limit by FDR 0.1
stainExp.genes <- topGene(stainExp.adv.out, cutoff=0.1,
                          method="pfp", logged=T,
                          logbase=2,
                          gene.names=unlist(lookUp(probes,
                          'lumiHumanAll.db', 'SYMBOL')))
stainExp.genes$Table1

write.table(stainExp.genes$Table1,
            file='stain_0.1_rp_upregulated.txt', sep="\t");
write.table(stainExp.genes$Table2,
            file='stain_0.1_rp_downregulated.txt', sep="\t");

# Again for the probe subset
stainExpSel.adv.out <- RPadvance(selDataMatrix, stainCl,
                              stainOrigin, num.perm=100,
                              logged=T, rand=123);                              

pdf("stain_vs_unstained_sel_rp.pdf");
plotRP(stainExpSel.adv.out, cutoff=0.1);
dev.off();                              

stainExpSel.genes <- topGene(stainExpSel.adv.out, cutoff=0.1,
                          method="pfp", logged=T,
                          logbase=2,
                          gene.names=unlist(lookUp(selProbes,
                          'lumiHumanAll.db', 'SYMBOL')))
                              
stainExpSel.genes$Table2

write.table(stainExpSel.genes$Table1,
            file='stain_0.1_sel_rp_upregulated.txt', sep="\t");
write.table(stainExpSel.genes$Table2,
            file='stain_0.1__sel_rp_downregulated.txt', sep="\t");                              

# Test visualization
library(pheatmap)                              
                            
# Merge up/down regulated genes
finalGenes <- topGene(stainExpSel.adv.out, cutoff=0.1,
                      method="pfp", logged=T,
                      logbase=2,
                      gene.names=selProbes)
finalGenesComb <- rbind(finalGenes$Table1, finalGenes$Table2)                              


m <- as.matrix(selDataMatrix[rownames(selDataMatrix) %in% rownames(finalGenesComb), ])
rownames(m) <- getSYMBOL(rownames(m), 'lumiHumanAll.db')

# Re-order samples in a meaningful away
colOrder <- c("LM34_Dil_unst", "LM42_Dil_unst", "LM44_Dil_unst", "LM28_Dil_unst",
              "LM34_Dil_stain", "LM42_Dil_stain", "LM44_Dil_stain", "LM28_Dil_stain")
mS <- m[, colOrder]       

# Difference in log intensities
mSF <- mS[, 1:4] - mS[, 5:8]                              
colnames(mSF) <- c('LM34', 'LM42', 'LM44', 'LM28')                               
                                                      
pdf('heatmap_0.1.pdf')                              
pheatmap(mSF, 
         cluster_cols=F,
         scale='column',
         fontsize=12)
dev.off()   

