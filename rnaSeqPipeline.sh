#!/bin/bash
#Assumes all necessary packages are already installed in env
#Configured for paired end reads
SECONDS=0

#Declare array containing srr accession numbers
srrArray=("SRR15852396" "SRR15852426")

# Step [0]: Retrieve SRA Data
# Use sra-tools 2.10, updated version does not work
#Set working directory
cd /Users/jeamin/Documents/bioinfo/bc_tissue_rnaSeq/sra/

for srr in ${srrArray[@]};
do 
  echo "starting data retrieval of" $srr
  prefetch $srr
  echo "begginning fasterq-dump, fastq files will be stored in rawData/"
  fasterq-dump $srr -O ../rawData/
done

# Delete files created from prefetch to free up disk space
rm -r SRR*

# Step [1]: Quality Control with fastqc
#Change working directory
cd /Users/jeamin/Documents/bioinfo/bc_tissue_rnaSeq/

for srr in ${srrArray[@]};
do
    #Set value for forward and reverse-read fastq files 
    fq_fwd=${srr}_1.fastq
    fq_rev=${srr}_2.fastq

    #FastQC on forward reads
    fastqc rawData/$fq_fwd -o rawData/

    #FastQC on reverse reads
    fastqc rawData/$fq_rev -o rawData/

    echo "FastQC done, reports stored in rawData/"

# Step [2]: Trimming with trimmomatic
#change adapter file accordingly

    trimmomatic PE \
    -threads 4 \
    -phred33 \
    -trimlog trimmedData/${srr}_trimmed.log \
    rawData/$fq_fwd rawData/$fq_rev \
    -baseout trimmedData/${srr}_trimmed.fastq \
    ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36  

    echo "Trimming completed for" $srr ", outputs stored in /trimmedData"
    #Have to delete the raw fastq after trimming to free disk space
    rm rawData/$fq_fwd
    rm rawData/$fq_rev
done

#Step [2.5]: quality checking the now trimmed data (ignoring unpaired reads)
fastqc trimmedData/*P.fastq -o trimmedData/
echo "QC on paired trimmed data complete, begginning alignment!"

# Step [3]: Aligning with HISAT2
# obtain genome index from:
# wget https://genome-idx.s3.amazon.aws.com/hisat/grch38_genome.tar.gz

for srr in ${srrArray[@]};
do
  #Set value for forward paired and reverse paired files
  fq_1p=${srr}_trimmed_1P.fastq
  fq_2p=${srr}_trimmed_2P.fastq

  echo begginning alignment of $srr to human genome
  hisat2 \
  -q --rna-strandness R \
  -x HISAT2/hisatGenomes/grch38/genome \
  -1 trimmedData/$fq_1p \
  -2 trimmedData/$fq_2p \
  | samtools sort -o mappedReads/${srr}.bam

done

# Step [4] Quantification with featureCounts
# download gtf file (variable to updates)
# wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz

# we will input both our bams in one command and get one output file
featureCounts \
-p -s 2 -a index/Homo_sapiens.GRCh38.109.gtf \
-o quants/bc_tissue_featureCounts.txt mappedReads/SRR15852396.bam mappedReads/SRR15852426.bam


# Step [5] Process counts file for DeSEQ2
# Removes columns and rows that are not needed for the next step and renames column headers to SRR#
(cat quants/bc_tissue_featureCounts.txt | cut -f1,7,8 | sed '1d' \
| sed -e "1s/mappedReads\/SRR15852396.bam/SRR15852396/g" \
-e "1s/mappedReads\/SRR15852426.bam/SRR15852426/g") > quants/bc_tissue_counts.txt


duration=SECONDS
echo "$((duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
