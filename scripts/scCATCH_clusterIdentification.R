## Trying to use scCATCH for cluster identification in Karen's data

library(Seurat)
library(dplyr)
library(knitr)
library(cowplot)
library(DT)
library(scCATCH)


projectsDir <- "/camp/stp/babs/working/sopenam/projects/"
PI <- "tybulewiczv/"
scientist <- "micaela.sartoretti/"
project <- "SnRNAseq_of_human_adult_Down_Syndrome_hippocampi-ms358/"
resultsDir <- paste0( projectsDir, 
                      PI,
                      scientist,
                      project, 
                      "analysis_pilot/")


sceL3 <- readRDS( paste0(resultsDir,  "sceL3.rds"))
integrated  <- sceL3[[1]]
integrated  <- FindClusters(integrated , resolution = 0.4, verbose=FALSE)
# need normalized data (scale.data slot)
dat <- GetAssayData( integrated, slot="scale.data"  )

obj <- createscCATCH( data = dat, 
                       cluster = integrated@meta.data$seurat_clusters)
                       

obj <- findmarkergene(obj,
                               species = "Human",
                               marker = cellmatch,
                               tissue = "Brain")

obj <- findcelltype(object = obj)

write.table( obj@celltype, paste0( resultsDir, "scCATCH_cluster_annotation.txt"), sep="\t", row.names=FALSE, quote=FALSE)