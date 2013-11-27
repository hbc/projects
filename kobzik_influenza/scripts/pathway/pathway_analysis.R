# pathway enrichment, adapted from Emmanuel Dimont
# usage: pathEnrich(your_gene_list_of_interest, your_background_gene_list)
# gene lists should be a list of entrez ids
# background_gene_list should be tailored to your experiment, so for example if
# you are looking at RNA-seq data, your background should be all non-zero expressed genes,
# not all genes
pathEnrich = function (Genelist, geneset=pathways.Hs, bgGenelist)
{
    Nbackground = length(bgGenelist)
    genelist = unique(Genelist[!is.na(Genelist)])
    Nsig <- length(genelist)
    hyper <- as.data.frame(matrix(nrow = length(geneset), ncol = 1))
    colnames(hyper) <- c("p-value")
	hyper[,1] = as.numeric(lapply(geneset,function(x)
		{
			if(length(intersect(genelist,x))<1) return(1)
			else return(sum(dhyper(length(intersect(genelist,x)):Nsig,length(x), Nbackground - length(x), Nsig)))
		}))
    hyper[,2] <- p.adjust(hyper[, 1], method = "BH")
	overlap = lapply(geneset,function(x)
		{
			return(as.list(intersect(genelist,x)))
		})
    hyper[,3] = as.numeric(lapply(overlap,function(x) return(length(x))))
    hyper[,4] = as.numeric(lapply(geneset,function(x) return(length(x))))
    hyper[,5] <- names(geneset)
	genes = lapply(overlap, function(x) return(as.numeric(x)))
	hyper$genes = I(genes)
    colnames(hyper) <- c("p.value","FDR", "nGenes","nPathway","Name", "Genes")
    hyper = hyper[with(hyper, order(FDR)),]
    return(hyper)
	}
