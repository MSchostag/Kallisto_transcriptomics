# Kallisto_transcriptomics
This is a guide how to do basic transcriptomic (RNA seq) analysis from bacterial monoculture. This is ment as a guide for researchers at The Center for Microbial Secondary Metabolites (CeMiSt) - https://cemist.dtu.dk/.   

## Setting up the environment in linux

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

## QC of raw sequencing reads
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
```



## setting up the files
