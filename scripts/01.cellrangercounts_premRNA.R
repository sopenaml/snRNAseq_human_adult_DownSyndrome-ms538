library(openxlsx)

projectsDir <- "/camp/stp/babs/working/sopenam/projects/"
PI <- "tybulewiczv/"
scientist <- "micaela.sartoretti/"
project <- "/SnRNAseq_of_human_adult_Down_Syndrome_hippocampi-ms358/"

design <- read.xlsx( paste0(projectsDir, PI, scientist, project, "docs/SC22275.xlsx"),
                     sheet = 1,  startRow = 2)

design <- design[,c(1:2)]

seq.dir1 <- "/camp/stp/babs/inputs/sequencing/fastq/230106_A01366_0336_AHKFGJDSX5/fastq/SC22275/"
seq.dir2 <- "/camp/stp/babs/inputs/sequencing/fastq/221222_A01366_0333_AHGMKHDMXY/fastq/SC22275/"
dataDir <- paste0(projectsDir, PI,scientist,project,"data/")
if(!file.exists( dataDir)){ dir.create( dataDir)}
docsDir <- paste0(projectsDir, PI,scientist,project,"docs/")
if(!file.exists( docsDir)){ dir.create( docsDir)}
cellrangerDir <- paste0( projectsDir, PI, scientist,project, "cellranger" )
if(!file.exists( cellrangerDir)){ dir.create( cellrangerDir)}
scriptsDir <-  paste(projectsDir, PI,scientist,project,"scripts/",sep="")
if(!file.exists( scriptsDir)){ dir.create( scriptsDir)}


samples1 <- dir(seq.dir1, pattern="_S[0-9]*_L00[1-9]*_R[1]_001.fastq.gz")
samples1 <- gsub( "_R[1]_001.fastq.gz", "", samples1)
samples1 <- unique(samples1)

samples2 <- dir(seq.dir2, pattern="_S[0-9]*_L00[1-9]*_R[1]_001.fastq.gz")
samples2 <- gsub( "_R[1]_001.fastq.gz", "", samples2)
samples2 <- unique(samples2)


shortsamplename <- gsub("_S[0-9]*_L00[1-9]*","", c(samples1, samples2)) # replace suffix with empty space
shortsamplename <- unique( shortsamplename )

ref <-"/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-gex-mm10-2020-A"


for (sample in shortsamplename) {
  

  cellranger <- paste("cellranger count",
                       paste0("--id=",  sample),
                       paste0("--transcriptome=", ref),
                       paste0("--fastqs=", paste( seq.dir1, seq.dir2, sep=",")),
                       paste0("--sample=", sample),
                       paste0("--localcores=", 24),
                       "--include-introns")
 

  ######## Write out to file
  
  sink(paste(scriptsDir,sample,"cellranger.sh",sep=""))
  cat("#!/bin/sh");cat("\n");
  cat("module purge") ;cat("\n");
  cat("module load CellRanger/6.1.2"); cat("\n");
  cat("cd ");cat( cellrangerDir); cat("\n");
  cat(cellranger); cat("\n");
  
  sink();
  
  sink(paste0( projectsDir, PI, scientist, project, "scripts/submit_cellranger.sh"), append = TRUE)
  cat("#!/bin/sh" ); cat ("\n\n");  
  cat("sbatch -c 24 --time=24:00:00"); cat(" "); 
  cat ( paste(scriptsDir,sample,"cellranger.sh",sep=""));cat("\n"); 
  sink()
  
}

system(paste( "bash", paste0( projectsDir, PI, scientist, project, "scripts/submit_cellranger.sh"), sep=" "))

