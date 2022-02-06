# killer-denovo
Assembles de novo contigs from raw Illumina sequence reads. Optimized for dsRNA elements in yeast. 

More specifically, this command line pipeline cleans and filters raw Illumina reads with fastp, filters genomic contaminants with bbsplit, assembles *de novo* contigs with SPAdes, and performs a local BLASTn of contigs using a curated database of notable yeast dsRNAs. The pipeline generates an html report containing graphs of contig length and coverage as well as numeric summary of filtering metrics. 

### **Dependencies**

* Java, BBMap, [SPAdes](https://github.com/ablab/spades), R, [fastp](https://github.com/OpenGene/fastp), Python3, ncbi-blast+, FastQC

<b>CLI Usage</b>
```
killer-denovo.sh -s sample_name -1 for_read.fastq.gz -2 rev_read.fastq.gz -o output_folder
```

<b>Options</b>
flag | description
------------ | -------------
-s	[ARG]	| (required) sample name, used for file naming purposes
-1	[ARG]	| (required) forward read file
-2	[ARG]	| (required) reverse read file
-o	[ARG]	| (required) destination of output folder
-t		| tests program with associated test file
-h		| help (print options)
<p>&nbsp;</p>


**Files Required for Program to Run:**

1. **Forward read file** – File contains raw NGS reads in forward orientation; “R1” should be in the filename somewhere (ex: “NCYC190_S3_L001_R1_001.fastq.gz”); can be .fastq or .fastq.gz

2. **Reverse read file** – File contains raw NGS reads in reverse orientation; “R2” should be in the filename somewhere (ex: “NCYC190_S3_L001_R2_001.fastq.gz”); can be .fastq or .fastq.gz

**Program Output**

|    Filename or Folder         |     Description                                                                                                |
|-------------------------------|----------------------------------------------------------------------------------------------------------------|
|     SAMPLE_contigs.fasta      |     This is a fasta-formatted file containing the de novo assembled contigs generated by SPAdes assembler.     |
|     SAMPLE_dsrna_blast.csv    |     A csv file containing local blastn of contigs against the Rowley lab virus blast database.                 |
|     SAMPLE_report.html        |     A summary of the de novo assembly results in html                                                          |
|     SAMPLE_report.md          |     A summary of the de novo assembly results in markdown                                                      |
|     spades_output             |     A folder containing all the SPAdes output. The contigs file is copied into the main “denovo” folder as “SAMPLE_contigs.fasta” so you don’t really need to open this spades folder unless you’re doing some weird behind-the-scenes stuff.                                      |
|     stats                     |     A folder containing FastQC reports about read quality before and after trimming/filtering with fastp (forward reads only). Also contains SAMPLE_fastp.html, a fastp report that contains information about the read quality before and after trimming/filtering with fastp. |
|     img                       |     A folder containing “SAMPLE_contigs.jpeg”, a dot plot of the de novo contig length versus coverage score   |

