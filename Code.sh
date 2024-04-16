#! /bin/sh

# raw reads location on Premise 
/mnt/lz01/macmaneslab/shared/hypothalamus_seq/raw_reads

# location of STAR index 
/mnt/lz01/macmaneslab/shared/STAR_index

# location of MapCount.job
/mnt/lz01/macmaneslab/macmanes/MapCount.job


#### Purpose: To align reads to ref genome and count them 

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
cd count_data

for filename in *counts; do                  
    [ -f "$filename" ] || continue
    mv "$filename" "${filename/MP/_male_pit}"
    mv "$filename" "${filename/FP/_female_pit}"
    mv "$filename" "${filename/FS/_female_sonpvn}"
    mv "$filename" "${filename/MS/_male_sonpvn}"
done

#running slurm script for mapping ~20 hours 
sbatch mapping_job.sh /mnt/lz01/macmaneslab/shared/STAR_index/ /mnt/lz01/macmaneslab/shared/hypothalamus_seq/raw_reads ~/RESULTS_DIR/results


