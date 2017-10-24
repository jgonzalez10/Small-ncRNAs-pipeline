#!/bin/sh

#This script analyzes Phytophthora infestans small ncRNAs at various stages of it life cycle: 24 hpi, 48 hpi, 72 hpi and in vitro mycelium. This data was originally published by Asman and collaborators (2014). Please make sure that you have the following tools installed and working in the computer or cluster where you are running this pipeline: sratoolkit/2.8.1, fastqc 0.11.2, cutadapt 1.12, bowtie 1.1.2 and htseq 0.6.1p1. If you have other versions available, please edit the module load command in the script.
  
printf "\nDownloading test data from NCBI\n"

module load sratoolkit/2.8.1

fastq-dump SRR1652586
fastq-dump SRR1652587
fastq-dump SRR1652588
fastq-dump SRR1652589

printf "\nChanging data names\n"

mv SRR1652586.fastq 88069_24H.fastq
mv SRR1652587.fastq 88069_48H.fastq
mv SRR1652588.fastq 88069_72H.fastq
mv SRR1652589.fastq 88069_mycelium.fastq

printf "\nTesting for data quality\n"

module load fastqc/0.11.2

fastqc 88069_24H.fastq
fastqc 88069_48H.fastq
fastqc 88069_72H.fastq
fastqc 88069_mycelium.fastq

printf "\nPre processing raw data\n"
printf "\nEliminating the adapter, trimming low quality bases and selecting by size between 18 and 33 bp\n"

module load cutadapt/1.12 

cutadapt --minimum-length 18 --maximum-length 33 --discard-untrimmed -a TGGAATTCTCGGGTGCCAAGG -u -1 -o 24.fastq 88069_24H.fastq
cutadapt --minimum-length 18 --maximum-length 33 --discard-untrimmed -a TGGAATTCTCGGGTGCCAAGG -u -1 -o 48.fastq 88069_48H.fastq 
cutadapt --minimum-length 18 --maximum-length 33 --discard-untrimmed -a TGGAATTCTCGGGTGCCAAGG -u -1 -o 72.fastq 88069_72H.fastq 
cutadapt --minimum-length 18 --maximum-length 33 --discard-untrimmed -a TGGAATTCTCGGGTGCCAAGG -u -1 -o mycelium.fastq 88069_mycelium.fastq 

printf "\nEliminating rRNA sequences\n"

module load bowtie/1.1.2

bowtie-build -f rRNA_Ensembl.txt rRNA

bowtie rRNA 24.fastq -S 24hrRNA.sam -t -q -v 0 -k 1 -p 8 --un 24h_no_rRNA.fq
bowtie rRNA 48.fastq -S 48hrRNA.sam -t -q -v 0 -k 1 -p 8 --un 48h_no_rRNA.fq
bowtie rRNA 72.fastq -S 72hrRNA.sam -t -q -v 0 -k 1 -p 8 --un 72h_no_rRNA.fq
bowtie rRNA mycelium.fastq -S myceliumrRNA.sam -t -q -v 0 -k 1 -p 8 --un mycelium_no_rRNA.fq

printf "\nEliminating tRNA sequences\n"

bowtie-build -f tRNA_Ensembl.txt tRNA

bowtie tRNA 24h_no_rRNA.fq -S 24htRNA.sam -t -q -v 0 -k 1 -p 8 --un 24h_no_rRNA_no_tRNA.fq
bowtie tRNA 48h_no_rRNA.fq -S 48htRNA.sam -t -q -v 0 -k 1 -p 8 --un 48h_no_rRNA_no_tRNA.fq
bowtie tRNA 72h_no_rRNA.fq -S 72htRNA.sam -t -q -v 0 -k 1 -p 8 --un 72h_no_rRNA_no_tRNA.fq
bowtie tRNA mycelium_no_rRNA.fq -S myceliumtRNA.sam -t -q -v 0 -k 1 -p 8 --un mycelium_no_rRNA_no_tRNA.fq

printf "\nEliminating small ncRNA sequences that align to the Solanum tuberosum genome\n"

bowtie-build -f S_tuberosum.fna STuber

bowtie STuber 24h_no_rRNA_no_tRNA.fq -S 24hSTuber.sam -t -q -v 0 -k 1 -p 8 --un 24h_no_rRNA_no_tRNA_no_Stuber.fq
bowtie STuber 48h_no_rRNA_no_tRNA.fq -S 48hSTuber.sam -t -q -v 0 -k 1 -p 8 --un 48h_no_rRNA_no_tRNA_no_Stuber.fq
bowtie STuber 72h_no_rRNA_no_tRNA.fq -S 72hSTuber.sam -t -q -v 0 -k 1 -p 8 --un 72h_no_rRNA_no_tRNA_no_Stuber.fq
bowtie STuber mycelium_no_rRNA_no_tRNA.fq -S myceliumSTuber.sam -t -q -v 0 -k 1 -p 8 --un mycelium_no_rRNA_no_tRNA_no_Stuber.fq

printf "\nAligning processed reads to the Phytophthora infestans  genome\n"

bowtie-build -f Phytophthora_infestans.fa Pinf
bowtie Pinf 24h_no_rRNA_no_tRNA_no_Stuber.fq -S 24hPinfestans.sam -t -q -v 0 -k 1 -p 8 
bowtie Pinf 48h_no_rRNA_no_tRNA_no_Stuber.fq -S 48hPinfestans.sam -t -q -v 0 -k 1 -p 8 
bowtie Pinf 72h_no_rRNA_no_tRNA_no_Stuber.fq -S 72hPinfestans.sam -t -q -v 0 -k 1 -p 8 
bowtie Pinf mycelium_no_rRNA_no_tRNA_no_Stuber.fq -S myceliumPinfestans.sam -t -q -v 0 -k 1 -p 8 
   
printf "\nGenerating count tables with HT-Seq\n"

module load htseq/0.6.1p1
htseq-count -s no -t exon -i gene_id -o 24h_htseqcount.sam 24hPinfestans.sam Pinfestans_definitivo_20172.gtf > 24h_htseqcount
htseq-count -s no -t exon -i gene_id -o 48h_htseqcount.sam 48hPinfestans.sam Pinfestans_definitivo_20172.gtf > 48h_htseqcount
htseq-count -s no -t exon -i gene_id -o 72h_htseqcount.sam 72hPinfestans.sam Pinfestans_definitivo_20172.gtf > 72h_htseqcount
htseq-count -s no -t exon -i gene_id -o mycelium_htseqcount.sam myceliumPinfestans.sam Pinfestans_definitivo_20172.gtf > mycelium_htseqcount
