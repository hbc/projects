install.packages(c("bookdown",
                   "cowplot",
                   "knitr",
                   "Matrix",
                   "mltools",
                   "pheatmap",
                   "RColorBrewer",
                   "readxl",
                   "rmarkdown",
                   "tidyverse",
                   "tools"))

source("https://bioconductor.org/biocLite.R")
biocLite(c("BiocStyle",
           "biomaRt",
           "scde",
           "steinbaugh/basejump"))

# [Seurat](http://satijalab.org/seurat/install.html)
devtools::install_url("https://github.com/satijalab/seurat/releases/download/v1.4.0/Seurat_1.4.0.9.tgz", binary = TRUE)
