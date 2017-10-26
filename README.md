# Small-ncRNAs-pipeline
This script analyzes Phytophthora infestans small ncRNAs at various stages of it life cycle: 24 hpi, 48 hpi, 72 hpi and in vitro mycelium. 
This data was originally published by Asman and collaborators (2014). 
1. Please make sure that you have the following tools installed and working in the computer or cluster where you are running this pipeline: 
sratoolkit/2.8.1
fastqc/0.11.2
cutadapt/1.12
bowtie/1.1.2
htseq/0.6.1p1. 

2. If you have other versions available, please edit the module load command in the script.

3. If these tools are installed and working, the script will function automatically. 
It is programmed to download the data sets for testing from NCBI using fastq-dump and to continue the analysis by itself.

4. If you want to use it for other data sets, please modify the identifiers indicated at the fastq-dump command.

5. As well, you will have to modify the parameters for cutadapt depending on your results from fastqc.

6. All required data sets are included in this GitHub page except for P. infestans and S. tuberosum genome sequences, because their file sizes exceeded the capacity of the repository. If required, these are available upon request. However, you can also download these genomes from Ensembl at http://protists.ensembl.org/Phytophthora_infestans/Info/Index and http://plants.ensembl.org/Solanum_tuberosum/Info/Index.
These genomes were downloaded in FASTA format and used without any alteration.

Please feel free to ask any questions! (j.gonzalez10@uniandes.edu.co)

