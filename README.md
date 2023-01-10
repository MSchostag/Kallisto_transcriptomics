# Kallisto_transcriptomics
This is a guide how to do basic transcriptomic (RNA seq) analysis from bacterial monoculture. This is ment as a guide for researchers at The Center for Microbial Secondary Metabolites (CeMiSt) - https://cemist.dtu.dk/.   

# Overview

- [Installing and setting up the environment in linux](#Installing-and-setting-up-the-environment-in-linux)
- [QC and filtering](#QC-and-filtering)
- [Running Kallisto](#Running-Kallisto)

# Installing and setting up the environment in linux

First you would need to install MiniConda. 
Go to https://docs.conda.io/en/latest/miniconda.html#linux-installers 

Find the newest version with 64-bit. Right-click and copy the hyperlink.

In terminal, download the file: 
```bash
wget LINK
```
Example: 
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-py310_22.11.1-1-Linux-x86_64.sh
```

Now run the bash file: 

```bash
./ Miniconda*.sh
```

Install Mamba using Miniconda:
```bash
conda install mamba -n base -c conda-forge
```
This is used to create an environment to work in. Have a look at this guide for commandlines when working in environments: 
https://www.imranabdullah.com/2021-08-21/Conda-and-Mamba-Commands-for-Managing-Virtual-Environments 

Create an environment for Kallisto analysis:
```bash
mamba create --name kallisto
```
Activate the Kallisto environment: 
```bash
mamba activate kallisto 
```
***From now on you would work in the kallisto environment.***

You would need the following programs for the transcriptomic analysis: 
- FastP (https://github.com/OpenGene/fastp) - Used to clean and filter your FastQ files 
- GFFread (https://github.com/gpertea/gffread) - Used to format gff files to transcript fasta file
- Kallisto (https://github.com/pachterlab/kallisto) - Program for quantifying abundances of transcripts from RNA-Seq data 

**Installing the programs:** 
```bash
mamba install -c bioconda fastp
mamba install -c bioconda gffread
mamba install -c bioconda kallisto
```
Say yes to the installation, when asked. 

# QC and filtering
Start of by create a folder with all you files for you analysis: 
```bash
mkdir PROJECT_NAME 
```
Example: 
```bash
mkdir phaeobacter_S26_rna_seq 
```
Transfere all raw RNA reads to server and store in folder. This could be called "raw": 
```bash
mkdir raw 
cd raw # transfere to here
```

Now filter the raw reads using default settings with FastP. This also includes adapter removal.  
Folder for filtered read files
```bash
mkdir filtered
```
Running FastP
```bash
fastp -i FORWARD_READ_1.fastq.gz -I REVERSE_READ_1.fastq.gz -o filtered/FORWARD_READ_1_f.fastq.gz -O filtered/REVERSE_READ_1_f.fastq.gz -h filtered/READS_1.html

example:
fastp -i RNAseq_of_Phaeobacter_piscinae_dtdaB_rep_1_72_hrs_fw_reads_1.fastq.gz -I RNAseq_of_Phaeobacter_piscinae_dtdaB_rep_1_72_hrs_fw_reads_2.fastq.gz -o filtered/RNAseq_of_Phaeobacter_piscinae_dtdaB_rep_1_72_hrs_fw_reads_1_f.fastq.gz -O filtered/RNAseq_of_Phaeobacter_piscinae_dtdaB_rep_1_72_hrs_fw_reads_2_f.fastq.gz -h filtered/dtdaB_rep_1_72_hrs_fw_reads_1.html
```

# Running Kallisto
For kallisto one would need a transcript fasta file to create a index file 

Example of transcript file: 
```bash
>gene1 CDS=1-10
ATGTCCGGAAACGAGCAGGCCCCCGCAGACTATGGCGCGGA.....
>gene2 CDS=1-20
ATGTCCGGAAACGAGCAGGCCCCCGCAGACTATGGCGCGGA.....
```
If one do not have a transcript fasta file, you could create it with a GFF file using GFFread. 

Generate a FASTA file with the DNA sequences for all transcripts in a GFF file. For this operation a fasta file with the genomic sequences has to be provided as well. This can be accomplished with a command line like this:
```bash
gffread -w transcripts.fa -g genome.fa annotation.gff
```

## setting up the files
