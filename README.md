# Kallisto_transcriptomics
This is a guide how to do basic transcriptomic (RNA seq) analysis from bacterial monoculture. This is ment as a guide for researchers at The Center for Microbial Secondary Metabolites (CeMiSt) - https://cemist.dtu.dk/.   

# Overview

- [Installing and setting up the environment in linux](#Installing-and-setting-up-the-environment-in-linux)
- [QC and filtering](#QC-and-filtering)
- [Running Kallisto](#Running-Kallisto)
  - [Input files](#Input-files)
  - [Build a Kallisto transcriptome index](#Build-a-Kallisto-transcriptome-index)
  - [Generate abundance estimates](#Generate-abundance-estimates)
- [Differential gene expression analysis](#Differential-gene-expression-analysis)
- [DeSeq2 analysis using RStudio](#DeSeq2-analysis-using-RStudio)
  - [Setting up files](#setting-up-the-files)

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
## Input files 
For kallisto one would need a transcript fasta file to create a index file.

Example of transcript file: 
```bash
>gene1 CDS=1-10
ATGTCCGGAAACGAGCAGGCCCCCGCAGACTATGGCGCGGA.....
>gene2 CDS=1-20
ATGTCCGGAAACGAGCAGGCCCCCGCAGACTATGGCGCGGA.....
```
If one do not have a transcript fasta file, you could create it with a GFF file and genome fasta file using GFFread. 

Generate a FASTA file with the DNA sequences for all transcripts in a GFF file. For this operation a fasta file with the genomic sequences has to be provided as well. This can be accomplished with a command line like this:
```bash
gffread -w transcripts.fa -g genome.fa annotation.gff

Example: 
gffread -w S26_cds.fna -g GCF_000826835.2_ASM82683v2_genomic.fna GCF_000826835.2_ASM82683v2_genomic.gff
```
## Build a Kallisto transcriptome index

Now create a Kallisto transcriptome index:
```bash
kallisto index -i transcripts.idx transcripts.fasta.gz

Example: 
kallisto index -i transcripts.idx S26_cds.fna
```

## Generate abundance estimates

First make a output folder: 
```bash
mkdir output
```

Create abundance estimates, for each sample: 
```bash
kallisto quant -i transcripts.idx -o output -b 100 FORWARD_READ_1_f.fastq.gz REVERSE_READ_1_f.fastq.gz -t 20

Example: 
kallisto quant -i transcripts.idx -o output/dtdaB_rep_1_72_hrs raw/filtered/RNAseq_of_Phaeobacter_piscinae_dtdaB_rep_1_72_hrs_fw_reads_1_f.fastq.gz raw/filtered/RNAseq_of_Phaeobacter_piscinae_dtdaB_rep_1_72_hrs_fw_reads_2_f.fastq.gz -t 20
```
In this example one would use 20 threads and bootstrap 100 times.

Do this for each sample in your analysis. 

For each sample you would now have created the following files: 
```bash
abundance.h5  
abundance.tsv
run_info.json
```

The abundance.h5 will be used in the differential gene expression analysis. 

# Differential gene expression analysis
# DeSeq2 analysis using RStudio
For differential gene expression analysis we would use DeSeq2 
http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

## setting up the files
First we need to setup the files and make additional files. 

Start by creating a folder were all the output from the analysis will be located. 

Create a folder called "kallisto". Within this folder one should transfere all the output folder from Kallisto (\output folder from linux) 

You would need to create two meta files: 
1) meta.txt
2) samples.txt

meta.txt should contain the meta data for your analysis. **Minimum: samplename and condition. OBS!!! use the same headers as in the example.** 

Example: 
```bash
sample	condition	rep
dtdaB_rep_1_24_hrs	m_24	1
dtdaB_rep_2_24_hrs	m_24	2
dtdaB_rep_3_24_hrs	m_24	3
wildtype_rep_1_24_hrs	wt_24	1
wildtype_rep_2_24_hrs	wt_24	2
wildtype_rep_3_24_hrs	wt_24	3
```
See also the file meta.txt

sample.txt should contain all the names of sample folders with in \kallisto. **See example in sample.txt, again the headers should be the same as in the file.**  
Example: 
```bash
run
dtdaB_rep_1_24_hrs
dtdaB_rep_2_24_hrs
dtdaB_rep_3_24_hrs
wildtype_rep_1_24_hrs
wildtype_rep_2_24_hrs
wildtype_rep_3_24_hrs
```

Now open you RStudio and create a new project (File->create project) and save it in the newly created analysis folder.  

Copy the kallisto_deseq.R and place it in the analysis folder. 

Now open the kallisto_deseq.R file with in the R project and run the analysis. 
