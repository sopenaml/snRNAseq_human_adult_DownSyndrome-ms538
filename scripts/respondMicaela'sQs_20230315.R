library(Seurat)
library(dplyr)
library(knitr)
library(cowplot)
library(DT)
library(ggplot2)
library(openxlsx)
library( grid )
library(tidyverse)
library(patchwork)

projectsDir <- "/camp/stp/babs/working/sopenam/projects/"
PI <- "tybulewiczv/"
scientist <- "micaela.sartoretti/"
project <- "SnRNAseq_of_human_adult_Down_Syndrome_hippocampi-ms358/"


design <-read.xlsx( paste0( projectsDir, PI, scientist,project, "docs/", "SC22275.xlsx"), sheet=1, colNames =T, startRow=2)
rownames(design) <- design$Sample.limsid

#seq.dir <- "/camp/stp/babs/working/sopenam/projects/tybulewiczv/eva.lana-elola/03_scRNASeq_Dyrk1a_dosage/data/runs/"
seq.dir <-  paste0( projectsDir, PI, scientist,project, "cellranger/")
resultsDir <- paste0( projectsDir, PI, scientist,project, "analysis_pilot/" )

 sceL3 <- readRDS( paste0(resultsDir, "sceL3.rds"))

 sceL3[[1]] <- FindClusters(sceL3[[1]], resolution = 0.6, verbose =FALSE)
 markers <- c( "SLC17A7", "AQP4", "MOG", "PDGFRA", "GAD1", "RELN", "IGFBPL1","SOX11", "PLP1")
 FeaturePlot( sceL3[[1]], features=markers)
ggsave( paste0(resultsDir, "RequestedMarkers_20230315.jpeg"))

 DimPlot( sceL3[[1]], label=T)
ggsave( paste0(resultsDir, "UMAP_res0.6_20230315.jpeg"))

VlnPlot( sceL3[[1]], features=c("nFeature_RNA", "percent.mt") )
ggsave(paste0(resultsDir, "VlnPlot_nFeatures_20230315.jpeg") )

VlnPlot( sceL3[[1]], features=c("SLC17A7", "AQP4", "MOG", "PDGFRA", "GAD1", "RELN", "IGFBPL1", "PLP1","SOX11"))
ggsave(paste0(resultsDir, "VlnPlot_RequestedMarkers_20230315.jpeg") )

 FeaturePlot( sceL3[[1]], features=c("NOS1", "CXCL14", "CNR1" ))
 ggsave(paste0(resultsDir, "FeaturePlot_Nos1_Cxcl14_Cnr1_20230315.jpeg") )

VlnPlot( sceL3[[1]], features=c("NOS1", "CXCL14", "CNR1" ))
 ggsave(paste0(resultsDir, "VlnPlot_Nos1_Cxcl14_Cnr1_20230315.jpeg") )


## Micaela want's to remove cluster 7,8,17 because it's doublets
sce <- sceL3[[1]]
#remove cluster 7
sce <- subset( sce, idents=c(7,8,17), invert=TRUE)
## add cluster ids to the UMAP
DimPlot( sce, label=TRUE)
ggsave( "UMAP_20230320.jpeg")
library(ggrepel)

new.cluster.ids <- c( "Oligodendrocytes 1", 
                      "Astrocytes 1", 
                      "Oligodendrocytes 2", 
                      "Excitatory neurons 1",
                      "OPCs",
                      "Microglia",
                      "Astrocytes 2",
                      "Astrocytes 3",
                      "Inhibitory neurons 1",
                      "Inhibitory neurons 2",
                      "Excitatory neurons 2",
                      "Astrocytes 4",
                      "Oligodendrocytes 3",
                      "Excitatory neurons 3",
                      "Endothelial cells")

names(new.cluster.ids) <- levels(sce)
sce <- RenameIdents(sce, new.cluster.ids)
DimPlot(sce, reduction = "umap", label = TRUE, pt.size = 0.5, repel=TRUE) + NoLegend()
ggsave( "UMAP_with_labels_20230320.jpeg")
DimPlot(sce, reduction = "umap")
ggsave( "UMAP_no_labels_20230320.jpeg")

p1 <- DimPlot(sceL3[[1]], reduction="umap") 
p2 <- DimPlot(sce, reduction = "umap")
ggsave()

# more markers

markers2 <- c( "RBFOX3","FLT1","APBB1IP","VCAN","MOBP",  "GAD2")
FeaturePlot( sce, features= markers2, ncol=3)
ggsave( "FeaturePlot_newmarkers_20230318.jpeg")
VlnPlot( sce, features=markers2)
ggsave( "VlnPlot_newmarkers_20230318.jpeg")

# need to reorder clusters for a plot for Micaela

markers3 <- c("RBFOX3", "SLC17A7", "GAD2", "FLT1", "AQP4", "APBB1IP", "PDGFRA", "MOBP")

VlnPlot( sce, features=markers3,  stack=TRUE, fill.by="ident")  + 
   theme(legend.position = "none") 
   ggsave( "Stacked_VlnPlot_20230320.jpeg")


levels_order <- c( "Excitatory neurons 1",
"Excitatory neurons 2",
"Excitatory neurons 3",
"Inhibitory neurons 1",
"Inhibitory neurons 2",
"Endothelial cells",
"Astrocytes 1",
"Astrocytes 2",
"Astrocytes 3",
"Astrocytes 4",
"Microglia",
"OPCs",
"Oligodendrocytes 1",
"Oligodendrocytes 2",
"Oligodendrocytes 3"
)


#clusters <- levels(Idents(sce))
Idents(sce) <- factor( Idents(sce), levels = rev(levels_order))
VlnPlot( sce, features=markers3,  stack=TRUE, fill.by="ident")  + 
   theme(legend.position = "none") 
   ggsave( "Stacked_ordered_VlnPlot_20230320.jpeg")