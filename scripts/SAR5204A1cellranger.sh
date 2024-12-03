#!/bin/sh
module purge
module load CellRanger/6.1.2
cd /camp/stp/babs/working/sopenam/projects/tybulewiczv/micaela.sartoretti//SnRNAseq_of_human_adult_Down_Syndrome_hippocampi-ms358/cellranger
cellranger count --id=SAR5204A1 --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-gex-GRCh38-2020-A --fastqs=/camp/stp/babs/inputs/sequencing/fastq/230106_A01366_0336_AHKFGJDSX5/fastq/SC22275/,/camp/stp/babs/inputs/sequencing/fastq/221222_A01366_0333_AHGMKHDMXY/fastq/SC22275/ --sample=SAR5204A1 --localcores=24 --include-introns
