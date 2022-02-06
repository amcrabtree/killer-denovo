# killer-denovo
Assembles de novo contigs from raw Illumina sequence reads. Optimized for dsRNA elements in yeast. 

### Dependencies

> Java, BBMap, SPAdes, R, fastp, Python3, ncbi-blast+, FastQC

This command line pipeline cleans and filters raw Illumina reads with fastp, filters genomic contaminants with bbsplit, assembles de novo contigs with SPAdes, and performs a local BLASTn of contigs using a curated database of notable yeast dsRNAs. The pipeline generates an html report containing graphs of contig length and coverage as well as numeric summary of filtering metrics. 