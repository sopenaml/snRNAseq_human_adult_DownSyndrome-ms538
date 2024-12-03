# Run Clustifyr on cell markers


library( clustifyr)
library( clustifyrdatahub)
library( ExperimentHub)
library( Seurat)
library( dplyr)
# dir where to store results

projectsDir <- "/camp/stp/babs/working/sopenam/projects/"
PI <- "tybulewiczv/"
scientist <- "micaela.sartoretti/"
project <- "SnRNAseq_of_human_adult_Down_Syndrome_hippocampi-ms358/"
resultsDir <- paste0( projectsDir, 
                      PI,
                      scientist,
                      project, 
                      "analysis_pilot/")


#define species to use

species <- "homo_sapiens"
eh <- ExperimentHub::ExperimentHub(ask=FALSE)


if (species == "mus_musculus"){
    eh <- eh[eh$species=="Mus musculus"]
} else if (species == "homo_sapiens"){
    eh <- eh[eh$species=="Homo Sapiens"]
} else {
    stop("For species other than human or mouse the reference dataset needs to be manually curated.")
}


refs <- ExperimentHub::listResources(eh, "clustifyrdatahub")
if (species == "mus_musculus"){
    selVec <- c(
        "ref_MCA",
        "ref_hema_microarray",  
        "ref_tabula_muris_drop", 
        "ref_tabula_muris_facs", 
        "ref_mouse.rnaseq",
        "ref_moca_main",
        "ref_immgen" #,
         #"ref_mouse_atlas" This one does not work at present
    )
} else if (species == "homo_sapiens"){
    selVec <- c(
        "ref_hema_microarray",
        "ref_cortex_dev",
        "ref_pan_indrop",
        "ref_pan_smartseq2"
    )
    
} 

sceL3 <- readRDS( paste0(resultsDir, "sceL3.rds"))
integrated <- sceL3[[1]]
integrated  <- FindClusters(integrated , resolution = 0.4, verbose=FALSE)


for (i in 1:length(selVec)){
    print(paste0("Processing ", selVec[i]))
    ref_mat <- ExperimentHub::loadResources(
        eh, 
        "clustifyrdatahub",
        selVec[i]
    )[[1]]
    
    top2000 <- head(Seurat::VariableFeatures(integrated), 2000)
    
    ## Get count table ##
    counts <- integrated@assays$RNA@counts
    set.seed(123)
    select <- sample( colnames(counts), 5000)
    counts <- counts[, select]
    counts <- as.matrix(counts)
    metadata <- integrated@meta.data[select,]

    ## Create correlation matrix
    res <- clustifyr::clustify(
        input = counts , # matrix of normalized scRNA-seq counts (or SCE/Seurat object)
        metadata = metadata, # meta.data table containing cell clusters
        cluster_col = "seurat_clusters", # name of column in meta.data containing cell clusters
        ref_mat = ref_mat, # matrix of RNA-seq expression data for each cell type
        query_genes = top2000 # list of highly varible genes identified with Seurat
    )
    
    # Call cell types
    res2 <- clustifyr::cor_to_call(
      cor_mat = res,                  # matrix correlation coefficients
      cluster_col = "seurat_clusters" # name of column in meta.data containing cell clusters
    )
    
    res2 <- data.frame(res2)
    res2[["reference"]] <- selVec[i]
         
    if (i==1){
        dfRes <- res2
    } else {
        dfRes <- rbind(
            dfRes, 
            res2
        )
    }
}
## Remove low confidence hits
dfRes <- dfRes[dfRes$r > 0.5,]
dfRes <- dfRes[order(dfRes$seurat_clusters),]
dfRes$r <- round(dfRes$r, 4)
## Produce output table ##

#dfRes %>% group_by( seurat_clusters) %>% filter ( max(r))


  write.table( dfRes , paste0 ( resultsDir, "Celltype_Assignment.txt"), sep="\t",row.names=F, quote=FALSE )

