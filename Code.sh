#! /bin/sh


####Premise 
# raw reads location on Premise 
/mnt/lz01/macmaneslab/shared/hypothalamus_seq/raw_reads

# location of STAR index 
/mnt/lz01/macmaneslab/shared/STAR_index

# location of MapCount.job
/mnt/lz01/macmaneslab/macmanes/MapCount.job


#### Purpose: To align reads to ref genome and count them 


#making directory for fastqc results
mkdir fastqc


#fastqc slurm script 

#!/bin/bash

#SBATCH --partition=macmanes
#SBATCH --cpus-per-task=40
#SBATCH --mem 310Gb
#SBATCH --open-mode=append
#SBATCH --exclude=node117,node118,node139
#SBATCH --output fastq_data.log

module load linuxbrew/colsa

fastqc /mnt/lz01/macmaneslab/shared/hypothalamus_seq/raw_reads/*fastq.gz -t 40 -o /mnt/lz01/macmaneslab/smc1079/fastq

#getting qc html files and downloading to local device 
get /mnt/lz01/macmaneslab/smc1079/fastq/*.html


#making directory for multiqc
mkdir multiqc

#multiqc script 
#!/bin/bash

#SBATCH --partition=macmanes
#SBATCH --cpus-per-task=40
#SBATCH --mem 310Gb
#SBATCH --open-mode=append
#SBATCH --exclude=node117,node118,node139
#SBATCH --output multiqc_data.log
 
module unload linuxbrew/colsa
module load anaconda/colsa
conda activate multiqc-1.10.1

multiqc /mnt/lz01/macmaneslab/sj1187/hypothalamus/mapping/fastq/*fastqc* -o /mnt/lz01/macmaneslab/smc1079/multiqc

#getting and downloading multiqc html files to local device 
get /mnt/lz01/macmaneslab/smc1079/fastq/*.html

#create sh file for mapping slurm script 
nano mapping_job.sh

#creating parameters and activating anaconda/colsa 
#!/bin/bash
#SBATCH --partition=macmanes
#SBATCH --cpus-per-task=40
#SBATCH --mem 310Gb
#SBATCH --exclude=node117,node118
#SBATCH --output hypo_data.log

ulimit -n 10000
module purge
module load anaconda/colsa
conda activate star-2.7.10b

#sbatch MapCount.job index_location reads_dir _results_dir

DIR=$(pwd)
INDEX=$1
READ_DIR=$2
RESULTS_DIR=$3

echo "index: $INDEX"
echo "Reads directory: $READ_DIR"
echo "RESULTS_DIR: $RESULTS_DIR"

# making a directory for star read mapping 
mkdir $RESULTS_DIR

#list all samples to map
samples=$(basename -a ${READ_DIR}/*1.fastq.gz | sed "s/R1.fastq.gz//g")

#
echo Mapping these samples: $samples

#mapping reads with star 
STAR --runMode alignReads \
--genomeDir $INDEX \
--readFilesIn ${READ_DIR}/${sample}R1.fastq.gz ${READ_DIR}/${sample}R2.fastq.gz \
--readFilesCommand zcat \
--runThreadN 40 \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 10 \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0.3 \
--quantMode GeneCounts \
--outFileNamePrefix $RESULTS_DIR/${sample} \
--outSAMtype BAM Unsorted \
--genomeLoad LoadAndKeep \
--limitBAMsortRAM 254748364800 \
--sjdbGTFtagExonParentTranscript Parent

#######################
##### count reads #####
#######################

conda deactivate
conda activate htseq-0.13.5

# count reads with htseq-count
echo Counting genes with htseq-count for these samples: $samples

for sample in $samples
do
htseq-count \
--stranded=reverse \
-f bam \
-t exon \
--idattr=Parent \
--additional-attr=gene \
--additional-attr=GeneID \
$RESULTS_DIR/${sample}Aligned.out.bam \
/mnt/lz01/macmaneslab/shared/genome/PerEre_H2_v1/PerEre_H2_v1_genomic.gff > $RESULTS_DIR/${sample}gene.counts
done

#mapping data analysis 
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "date" "sample" "reads" "unique" "multi" "toomany" "unmapped"> mappingdata.txt

for file in $(ls -lthr  *Log.final.out | awk '{print $9}')
do
sample=$(echo $file | sed "s/_Log.final.out//g" | sed "s/.Log.final.out//g")
date=$(grep "Started job" $file | cut -f 2 | cut -d ' ' -f1-2)
reads=$(grep "Number of input reads" $file | cut -f 2)
unique=$(grep "Uniquely mapped reads %" $file | cut -f 2)
multi=$(grep "% of reads mapped to multiple loci" $file | cut -f 2)
toomany=$(grep "% of reads mapped to too many loci" $file | cut -f 2)
unmapped=$(grep "Number of reads unmapped: too short" $file | cut -f 2)
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$date" "$sample" "$reads" "$unique" "$multi" "$toomany" "$unmapped" >> mappingdata.txt
done

#editing & downloading count files
mkdir count_data

for samples in *counts; 
do awk '{print $1, "\t", $NF}' $samples > count_data/$samples; 
done

#to make parsing in R easier - changing files names 
mv "1318FP_S28_L002_gene.counts" "1318_female_pit_no_S28_L002_gene.counts"
mv "1318FS_S27_L002_gene.counts" "1318_female_sonpvn_no_S27_L002_gene.counts"
mv "1327FP_S30_L002_gene.counts" "1327_female_pit_no_S30_L002_gene.counts"
mv "1327FS_S29_L002_gene.counts" "1327_female_sonpvn_no_S29_L002_gene.counts"
mv "1331MP_S32_L002_gene.counts" "1331_male_pit_no_S32_L002_gene.counts"
mv "1331MS_S31_L002_gene.counts" "1331_male_sonpvn_no_S31_L002_gene.counts"
mv "1365MP_S38_L002_gene.counts" "1365_male_pit_no_S38_L002_gene.counts"
mv "1365MS_S37_L002_gene.counts" "1365_male_sonpvn_no_S37_L002_gene.counts"
mv "1368MP_S36_L002_gene.counts" "1368_male_pit_no_S36_L002_gene.counts"
mv "1368MS_S35_L002_gene.counts" "1368_male_sonpvn_no_S35_L002_gene.counts"
mv "1369MP_S40_L002_gene.counts" "1369_male_pit_no_S40_L002_gene.counts"
mv "1369MS_S39_L002_gene.counts" "1369_male_sonpvn_no_S40_L002_gene.counts"
mv "1370FP_S26_L002_gene.counts" "1370_female_pit_no_S26_L002_gene.counts"
mv "1370FS_S25_L002_gene.counts" "1370_female_sonpvn_no_S25_L002_gene.counts"
mv "1372FP_S24_L002_gene.counts" "1372_female_pit_no_S24_L002_gene.counts"
mv "1372FS_S23_L002_gene.counts" "1372_female_sonpvn_no_S23_L002_gene.counts"
mv "1373FP_S22_L002_gene.counts" "1373_female_pit_no_S22_L002_gene.counts"
mv "1373FS_S21_L002_gene.counts" "1373_female_sonpvn_no_S21_L002_gene.counts"
mv "1374MP_S34_L002_gene.counts" "1374_male_pit_no_S34_L002_gene.counts"
mv "1374MS_S33_L002_gene.counts" "1374_male_sonpvn_no_S33_L002_gene.counts"
mv "1381FP_S2_L002_gene.counts" "1381_female_pit_yes_S2_L002_gene.counts"
mv "1381FS_S1_L002_gene.counts" "1381_female_sonpvn_yes_S1_L002_gene.counts"
mv "1384FP_S4_L002_gene.counts" "1384_female_pit_yes_S4_L002_gene.counts"
mv "1384FS_S3_L002_gene.counts" "1384_female_sonpvn_yes_S3_L002_gene.counts"
mv "1386MP_S14_L002_gene.counts" "1386_male_pit_yes_S14_L002_gene.counts"
mv "1386MS_S13_L002_gene.counts" "1386_male_sonpvn_yes_S13_L002_gene.counts"
mv "1619MP_S12_L002_gene.counts" "1619_male_pit_yes_S12_L002_gene.counts"
mv "1619MS_S11_L002_gene.counts" "1619_male_sonpvn_yes_S11_L002_gene.counts"
mv "1620MP_S16_L002_gene.counts" "1620_male_pit_yes_S16_L002_gene.counts"
mv "1620MS_S15_L002_gene.counts" "1620_male_sonpvn_yes_S15_L002_gene.counts"
mv "1638MP_S18_L002_gene.counts" "1638_male_pit_yes_S18_L002_gene.counts"
mv "1638MS_S17_L002_gene.counts" "1638_male_sonpvn_yes_S17_L002_gene.counts"
mv "1639MP_S20_L002_gene.counts" "1639_male_pit_yes_S20_L002_gene.counts"
mv "1639MS_S19_L002_gene.counts" "1639_male_sonpvn_yes_S19_L002_gene.counts"
mv "1640FP_S6_L002_gene.counts" "1640_female_pit_yes_S6_L002_gene.counts"
mv "1640FS_S5_L002_gene.counts" "1640_female_sonpvn_yes_S5_L002_gene.counts"
mv "1641FP_S8_L002_gene.counts" "1641_female_pit_yes_S8_L002_gene.counts"
mv "1641FS_S7_L002_gene.counts" "1641_female_sonpvn_S7_L002_gene.counts"
mv "1642FP_S10_L002_gene.counts" "1642_female_pit_yes_S10_L002_gene.counts"
mv "1642FS_S9_L002_gene.counts" "1642_female_sonpvn_S9_L002_gene.counts"

#running slurm script for mapping ~20 hours 
sbatch mapping_job.sh /mnt/lz01/macmaneslab/shared/STAR_index/ /mnt/lz01/macmaneslab/shared/hypothalamus_seq/raw_reads ~/RESULTS_DIR/results



#####RStudio 

#installing packages 
install.packages(c("tximport", "DESeq2", "dplyr", "foreach", "data.table", "splines", "ggthemes", "scales", "gridExtra", "tidyr", "pheatmap", "RColorBrewer", "ggplot2", "BiocParallel", "apeglm", "tidyverse", "topGO", "GO.db", "beanplot"))

#Loading packages required for analysis 
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

install(c("EnhancedVolcano", "pasilla", "topGO", "apeglm", "DESeq2"))
install.packages("gprofiler2")

library(DESeq2)
library(dplyr)
library(foreach)
library(data.table)
library(splines)
library(ggthemes)
library(scales)
library(gridExtra)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(apeglm)
library(tidyverse)
library(topGO)
library(GO.db)
library(readr)
library(readxl)
library(lubridate)
library(patchwork)
library(EnhancedVolcano)
library(Cairo)
library(biomaRt)
library(plotly)
library(reutils)
library(remotes)
library(GeneTonic)
library(gprofiler2)
library(beanplot)
library(reshape)   
library(gplots)
library(pheatmap)

#Loading samples into R
directory <- "/Users/sarahcouture/Desktop/hypo_seq_countdata"
sampleFiles <- grep("male",list.files(directory),value=TRUE)
sampleCondition <- sub("(.*yes).*","\\1",sampleFiles)
#setting up a df / values for sample table 
sample=c()
sample$filename <- as.data.frame(sampleFiles)
sample$sample <- sampleFiles  %>% gsub(pattern = "_L002.gene.counts", replacement = "")
sample$sample <- gsub("_S.*", "", sample$sample)
sample$tissue <- unlist(lapply(strsplit(sample$sample, split = "_"),"[[", 3))
sample$trt <- unlist(lapply(strsplit(sample$sample, split = "_"),"[[", 4))
sample$Animal_ID <- unlist(lapply(strsplit(sample$sample, split = "_"),"[[", 1))
sample$sex <- unlist(lapply(strsplit(sample$sample, split = "_"),"[[", 2))
#making a sample table 
sampleTable <- data.frame(sampleName = sample$sample,
                           fileName = sample$filename,
                           tissue = sample$tissue,
                           water = sample$trt,
                           sex = sample$sex,
                           Animal_ID = sample$Animal_ID)
#sampleTable$water<- factor(sampleTable$water)
sampleTable$sex<- factor(sampleTable$sex)
sampleTable$tissue<- factor(sampleTable$tissue)
sampleTable$group <- factor(paste0(sampleTable$tissue, sampleTable$water))

sampleTable
 


##DESeq
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ sex + group )
ddsHTSeq

#Prefiltering 
smallestGroupSize <- 6
keep <- rowSums(counts(ddsHTSeq) >= 10) >= smallestGroupSize
ddsHTSeq <- ddsHTSeq[keep,]
ddsHTSeq

#basic analysis

ddsHTSeq <- DESeq(ddsHTSeq)
resbasic <- results(ddsHTSeq)
summary(resbasic)
plotMA(resbasic, ylim=c(-4,4))


# tissue-wise analysis

#ddsHTSeq_tissue <- ddsHTSeq

#ddsHTSeq_tissue$group <- factor(paste0(ddsHTSeq_tissue$tissue, ddsHTSeq_tissue$water))
#design(ddsHTSeq_tissue) <- ~ group
#design(ddsHTSeq_tissue) <- formula(~ group + water)
#ddsHTSeq_tissue <- DESeq(ddsHTSeq_tissue)
#resultsNames(ddsHTSeq_tissue)


#pca plot
library(plotly)
vsd <- varianceStabilizingTransformation(ddsHTSeq, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("water", "tissue"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca <- ggplot(pcaData, aes(x = PC1, y = PC2, color = water, shape = tissue, name=name)) +
  geom_point(size = 3, show.legend = TRUE) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + 
  theme_bw() +
  guides(shape = guide_legend(order = 1),color = guide_legend(order = 2)) 
pca
ggplotly(pca)


#more MA plots for genes up/down regulated in both PVN/SON and pituitary 
res_pit <- results(ddsHTSeq, contrast=c("group", "pitno", "pityes"), lfcThreshold = 0.1)
summary(res_pit)
plotMA(res_pit, ylim=c(-3,3))

down_pit <- filter(as.data.frame(res_pit), log2FoldChange < 0, padj < 0.05)
up_pit <- filter(as.data.frame(res_pit), log2FoldChange > 0, padj < 0.05)


res_sonpvn <- results(ddsHTSeq, contrast=c("group", "sonpvnno", "sonpvnyes"), lfcThreshold = 0.1)
summary(res_sonpvn)
plotMA(res_sonpvn, ylim=c(-3,3))

down_sonpvn <- filter(as.data.frame(res_sonpvn), log2FoldChange < 0, padj < 0.05)
up_sonpvn <- filter(as.data.frame(res_sonpvn), log2FoldChange > 0, padj < 0.05)


#Print results
write.table(res_pit, "DE_pit.tsv", row.names = TRUE, sep = "\t")
write.table(res_sonpvn, "DE_hypo.tsv", row.names = TRUE, sep = "\t")


#Bean plot Pit 
hypo_counts_df <- as.data.frame(counts(ddsHTSeq, normalized=TRUE))
write.table(hypo_counts_df, "hypo_counts_df.tsv", row.names = TRUE, sep = "\t")

hypo_counts_df <- read_delim("hypo_counts_df.tsv", 
    delim = "\t", escape_double = FALSE, 
    col_types = cols(GeneID = col_character()), 
    trim_ws = TRUE)


beanplotPIT <- function(vals, genenames)
{
  gene <- hypo_counts_df[grep(paste("\\b", vals, "\\b", sep=""),hypo_counts_df$GeneID), ]
  PIT_NF <-  gene %>% select(contains("no")) %>%  select(contains("pit")) %>% select(contains("female")) %>% pivot_longer(everything()) 
  PIT_YF <-  gene %>% select(contains("yes")) %>%  select(contains("pit")) %>% select(contains("female")) %>% pivot_longer(everything()) 
  PIT_NM <-  gene %>% select(contains("no")) %>%  select(contains("pit")) %>% select(!contains("female")) %>% pivot_longer(everything())
  PIT_YM <-  gene %>% select(contains("yes")) %>%  select(contains("pit")) %>% select(!contains("female")) %>% pivot_longer(everything()) 
  beanplot(log(PIT_NF$value+1),log(PIT_YF$value+1), log(PIT_NM$value+1),log(PIT_YM$value+1), ll = 0, beanlinewd=0,
           side = "no", xlab="", ylab='logTPM',
           main = genenames, log="y",
           col = list("brown","dodgerblue2"),
           axes=F)
  points(rep(1, length(PIT_NF$value)), log(PIT_NF$value+1), col='black', pch=15)
  points(rep(2, length(PIT_YF$value)), log(PIT_YF$value+1), col='black', pch=15)
  points(rep(4, length(PIT_YM$value)), log(PIT_YM$value+1), col='black', pch=15)
  points(rep(3, length(PIT_NM$value)), log(PIT_NM$value+1), col='black', pch=15)
  axis(2)
  axis(1, at = 1:4, labels = c("F - Dehydrated", "F - Hydrated", "M - Dehydrated", "M - Hydrated"))
  #arrows(x0 = 1.5, y0 = median(log(NO$value+1)), x1=1.5, y1 = median(log(YES$value+1)), lwd=2, length=.1)
}
beanplotPIT("rna-XM_059276299.1" , "Pituitary: MUP4")


head(down_pit %>% arrange(log2FoldChange), n=20)
tail(up_pit %>% arrange(log2FoldChange), n=20)


#Bean Plot SON/PVN 
#hypo_counts_df2 <- as.data.frame(counts(ddsHTSeq, normalized=TRUE))
#write.table(hypo_counts_df2, "hypo_counts_df2.tsv", row.names = TRUE, sep = "\t")

hypo_counts_df2 <- read_delim("hypo_counts_df2.tsv", 
    delim = "\t", escape_double = FALSE, 
    col_types = cols(GeneID = col_character()), 
    trim_ws = TRUE)


beanplotPIT <- function(vals, genenames)
{
  gene <- hypo_counts_df2[grep(paste("\\b", vals, "\\b", sep=""),hypo_counts_df2$GeneID), ]
  SONPVN_NF <-  gene %>% select(contains("no")) %>%  select(contains("sonpvn")) %>% select(contains("female")) %>% pivot_longer(everything()) 
  SONPVN_YF <-  gene %>% select(contains("yes")) %>%  select(contains("sonpvn")) %>% select(contains("female")) %>% pivot_longer(everything()) 
  SONPVN_NM <-  gene %>% select(contains("no")) %>%  select(contains("sonpvn")) %>% select(!contains("female")) %>% pivot_longer(everything())
  SONPVN_YM <-  gene %>% select(contains("yes")) %>%  select(contains("sonpvn")) %>% select(!contains("female")) %>% pivot_longer(everything()) 
  beanplot(log(SONPVN_NF$value+1),log(SONPVN_YF$value+1), log(SONPVN_NM$value+1),log(SONPVN_YM$value+1), ll = 0, beanlinewd=0,
           side = "no", xlab="", ylab='logTPM',
           main = genenames, log="y",
           col = list("brown","dodgerblue2"),
           axes=F)
  points(rep(1, length(SONPVN_NF$value)), log(SONPVN_NF$value+1), col='black', pch=15)
  points(rep(2, length(SONPVN_YF$value)), log(SONPVN_YF$value+1), col='black', pch=15)
  points(rep(4, length(SONPVN_YM$value)), log(SONPVN_YM$value+1), col='black', pch=15)
  points(rep(3, length(SONPVN_NM$value)), log(SONPVN_NM$value+1), col='black', pch=15)
  axis(2)
  axis(1, at = 1:4, labels = c("F - Dehydrated", "F - Hydrated", "M - Dehydrated", "M - Hydrated"))
  #arrows(x0 = 1.5, y0 = median(log(NO$value+1)), x1=1.5, y1 = median(log(YES$value+1)), lwd=2, length=.1)
}
beanplotPIT("rna-XM_059276299.1" , "SonPvn: MUP4")


head(down_pit %>% arrange(log2FoldChange), n=20)
tail(up_pit %>% arrange(log2FoldChange), n=20)



#pheatmap
#if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
 
#BiocManager::install("pheatmap")
library("pheatmap")

coldata <- read.csv("sampletable.csv",header=TRUE,row.names=1)

coldata$tissue <- factor(coldata$tissue)
coldata$water <- factor(coldata$water)
coldata$sex <- factor(coldata$sex)
coldata$Animal_ID <- factor(coldata$Animal_ID)
coldata$group <- factor(coldata$group)
rownames(coldata) <- coldata$sampleName
 
coldata <- coldata[order(as.numeric(as.factor(coldata$group))),]

cts <- read.table("hypo_counts_df.tsv",sep="\t",row.names=1,header=TRUE, check.names = FALSE)
cts <- as.matrix(cts)
cts <- cts[ ,rownames(coldata)]

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ group)

dds <- DESeq(dds)
res <- results(dds)
 
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:40]

ntd <- normTransform(dds)
df <- data.frame(group = colData(dds)[,c("group")], row.names = rownames(colData(dds)))
 

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df)
 

rownames(cts)[select]

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ sex + water + tissue + sex:water + tissue:water + sex:tissue)

m1<- model.matrix(~ sex + water + tissue + sex:water + tissue:water + sex:tissue, colData(dds))


dds<- DESeq(dds, full=m1, betaPrior = F)

res <- results(dds)
 
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:40]

ntd <- normTransform(dds)

df <- data.frame(group = colData(dds)[,c("water" , "tissue", "sex")], row.names = rownames(colData(dds)))
 

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df)
