#!/bin/bash

## Script for downloading reference sequences from RefSeq database
## Author: Angela Crabtree

# Install wget if required (wget is not usually installed on MacOS)
wget || brew install wget || /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
wget || brew install wget

# Download Homo sapiens, GRCh38.p13, GCF_000001405.39
wget -O ./refseq/Hsap.fasta.gz https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz

# Download Saccharomyces paradoxus CBS432, GCF_002079055.1
wget -O ./refseq/Spara.fasta.gz https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_paradoxus/all_assembly_versions/GCF_002079055.1_ASM207905v1/GCF_002079055.1_ASM207905v1_genomic.fna.gz

# Download Saccharomyces cerevisiae S288C, R64, GCF_000146045.2
wget -O ./refseq/Scer.fasta.gz https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/all_assembly_versions/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
