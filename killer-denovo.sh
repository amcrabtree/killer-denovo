#!/bin/bash

## de novo assembly pipeline
## Author: Angela Crabtree

########## optionS ############
while getopts "s:1:2:o:th" opt; do
	case ${opt} in
		s) sample=$OPTARG ;;
		1) fread=$OPTARG ;;
		2) rread=$OPTARG ;;
		o) usroutdir=$OPTARG ;; 
		t) mode="test" ;;
		h)
			printf "\n\n-------------------------------------------------------\n"
			printf "\noptions:\n"
			printf "\n"
			printf "   -s [arg]	sample name (required)\n"
			printf "   -1 [arg]	forward read file (required)\n"
			printf "   -2 [arg]	reverse read file (required)\n"
			printf "   -o [arg]	destination of output folder (required)\n"
			printf "   -t		test (ensures required CL apps are working)\n"
			printf "   -h		help\n"
			printf "\n-------------------------------------------------------\n\n\n"
			exit 0
			;;
	esac
done

########## ESTABLISH BACKGROUND LOGISTICS ############

## LOAD MODULES (SERVER USE ONLY)
local=$(module load 2>&1 >/dev/null | grep 'command not found' | wc -l)
if [ $local != 1 ]; then
	module load java
	module load python
	module load spades
	module load R
	module load ncbi-blast
	module load fastqc
fi

## SET UP DIRECTORIES
	scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"   # location of this script file
	if [ -n "$mode" ]; then
		printf "\n\tOutput files will be saved to your current directory.\n"
		usroutdir=$(pwd)
		sample="TEST_STRAIN"
		fread=${scriptdir}/test/NCYC190_S3_L001_R1_001.fastq.gz
		rread=${scriptdir}/test/NCYC190_S3_L001_R2_001.fastq.gz
	fi	
	appdir=~/bin
	refseqdir="${scriptdir}"/refseq
	outdir="${usroutdir}"/"${sample}"/denovo
	mkdir -p ${outdir}/stats ${outdir}/img ${outdir}/polished-reads ${outdir}/contigs
	blast_db=${scriptdir}/refseq/dsrna.fasta

## ESTABLISH OUTPUT FILE NAMES & DIRECTORIES
	fastp_html=${outdir}/stats/${sample}_fastp.html
	contig_file=${outdir}/${sample}_contigs.fasta
	contig_img=${outdir}/img/${sample}_contigs.jpeg
	blast_output=${outdir}/${sample}_dsrna_blast.csv
	report_file=${outdir}/${sample}_denovo_report.md
	report_html=${outdir}/${sample}_denovo_report.html
	bbsplitR1=${outdir}/polished-reads/${sample}_R1_bbsplit.fastq.gz
	bbsplitR2=${outdir}/polished-reads/${sample}_R2_bbsplit.fastq.gz
	
########################## READ CLEANING/FILTERING ##############################

	printf "\n\n*** Cleaning reads with Fastp ***\n\n"
	# run fastp
	fastp \
		--trim_poly_x --trim_front1 6 --trim_tail1 6 \
		-i "${fread}" \
		-I "${rread}" \
		-o ${outdir}/${sample}_R1_filt.fastq.gz \
		-O ${outdir}/${sample}_R2_filt.fastq.gz \
		-h "$fastp_html"
	rm fastp.json
	# run FastQC on forward reads before and after read cleaneaninging
	fastqc -o ${outdir}/stats "${fread}"
	fastqc -o ${outdir}/stats ${outdir}/${sample}_R1_filt.fastq.gz
	
######################### GENOMIC CONTAMINANT FILTERING #########################

## SPECIFY REFERENCE GENOMES TO FILTER OUT
	if [[ ${mode} = "test" ]]; then
		refgen="${refseqdir}"/Scer.fasta.gz # this is the genome to filter out in the test condition
	else
		scer_ref="${refseqdir}"/Scer.fasta.gz
		spara_ref="${refseqdir}"/Spara.fasta.gz
		hsap_ref="${refseqdir}"/Hsap.fasta.gz
		refgen=${scer_ref},${spara_ref},${hsap_ref} # these are the genomes to filter out
	fi

## REMOVE READS THAT MAP TO REFERENCE GENOME(S)
	printf "\n\n*** Filtering S. cerevisiae, S. paradoxus, and H. sapiens reads with BBsplit ***\n\n"
	cd ${outdir} # you MUST be in the same working directory as the reads or bbsplit won't work
	bbsplit.sh -da \
		in1=${sample}_R1_filt.fastq.gz \
		in2=${sample}_R2_filt.fastq.gz \
		ref=${refgen} \
		path="${refseqdir}" \
		basename=dirty_reads_%.fastq.gz \
		outu1=$bbsplitR1 \
		outu2=$bbsplitR2
	
############################## ASSEMBLE CONTIGS ###################################

## de novo assembly with SPAdes assembler
	printf "\n\n*** Assembling contigs with SPAdes ***\n\n"
 	spades.py --cov-cutoff 10 --only-assembler \
 		-1 $bbsplitR1 \
 		-2 $bbsplitR2 \
 		-o ${outdir}/spades_output
 	cp ${outdir}/spades_output/contigs.fasta $contig_file
	# metaspades
	spades.py --meta --only-assembler \
 		-1 $bbsplitR1 \
 		-2 $bbsplitR2 \
 		-o ${outdir}/metaspades_output
 	cp ${outdir}/metaspades_output/contigs.fasta $contig_file
	
## Local blastn
	printf "\n\n*** Performing blastn Search using Rowley Lab Virus BLAST Database ***\n\n"
	# make necessary database (only needs to be done once, but this way no one has to check)
	makeblastdb -in $blast_db -dbtype nucl -parse_seqids
	blastn -db $blast_db \
		-query $contig_file \
		-max_target_seqs 1 \
		-outfmt "10 qacc sacc stitle pident length qstart qend sstart send" \
		-out $blast_output
	# Format the output file
	printf "\n\n*** Formatting Blast Output ***\n\n"
	# adds a column for the sample name
	sed -i'.og' -e "s/^/${sample}, /" $blast_output
	# adds column headers (doesn't work on linux for some reason)
	sed -i '1 i\sample,contig,match_acc,match_description,pident,length,qstart,qend,sstart,send' "$blast_output"
		
## Graph read depth using in-house R script
	printf "\n*** RScript - Contig Graph ***\n"
	cd ${scriptdir}/bin
	chmod a+x CL_contig_graph.R
	./CL_contig_graph.R $contig_file $contig_img
	# graph all contigs produced by various kmers 
	for dir in ${outdir}/spades_output/K*; do
		if [[ -f ${dir}/final_contigs.fasta ]]; then
			cntg=${dir}/final_contigs.fasta
		else
			cntg=${dir}/simplified_contigs.fasta
		fi
		cp $cntg ${outdir}/contigs/${sample}_$(basename $dir)_contigs.fasta
		./CL_contig_graph.R $cntg ${outdir}/img/${sample}_$(basename $dir)_contigs.jpeg
	done
	# graph all contigs produced by various kmers (metaspades)
	for dir in ${outdir}/metaspades_output/K*; do
		if [[ -f ${dir}/final_contigs.fasta ]]; then
			cntg=${dir}/final_contigs.fasta
		else
			cntg=${dir}/simplified_contigs.fasta
		fi
		cp $cntg ${outdir}/contigs/${sample}_$(basename $dir)_meta_contigs.fasta
		./CL_contig_graph.R $cntg ${outdir}/img/${sample}_$(basename $dir)_meta_contigs.jpeg
	done
		
#################################### STORE STATS ###################################

# calculate and store values
readwc1=$(echo "$(gzip -dc "${fread}" | wc -l | cut -d" " -f1)/4" | bc)
readwc2=$(echo "$(gzip -dc "${rread}" | wc -l | cut -d" " -f1)/4" | bc)
t_read_i=$(echo "$readwc1 + $readwc2"| bc)
readwc3=$(echo "$(gzip -dc ${bbsplitR1} | wc -l | cut -d" " -f1)/4" | bc)
readwc4=$(echo "$(gzip -dc ${bbsplitR2} | wc -l | cut -d" " -f1)/4" | bc)
t_read_f=$(echo "$readwc3 + $readwc4"| bc)
t_read_hs=$(echo "$(gzip -dc ${outdir}/dirty_reads_Hsap.fastq.gz | wc -l | cut -d" " -f1)/4" | bc)
t_read_sc=$(echo "$(gzip -dc ${outdir}/dirty_reads_Scer.fastq.gz | wc -l | cut -d" " -f1)/4" | bc)
t_read_sp=$(echo "$(gzip -dc ${outdir}/dirty_reads_Spara.fastq.gz | wc -l | cut -d" " -f1)/4" | bc)
contam1=$(printf %.2f $(echo "($t_read_hs/$t_read_i)*100"| bc -l))
contam2=$(printf %.2f $(echo "($t_read_sc/$t_read_i)*100"| bc -l))
contam3=$(printf %.2f $(echo "($t_read_sp/$t_read_i)*100"| bc -l))

spades_v=$(spades.py --version) 
java_v=$(java --version | head -n 3 | tail -n 1 )
python_v=$(python3 --version)
rscript_v=$(Rscript --version 2>&1 >/dev/null | head)
blast_v=$(blastn -version | head -n 1)
bbsplit_v=$(bbsplit.sh --version 2>&1 >/dev/null | head -n 4 | tail -n 1)
fastp_v=$(fastp --version 2>&1 >/dev/null | head)
fastqc_v=$(fastqc --version | head -n 3 | tail -n 1 )

################################ CREATE MARKDOWN FILE #####################################

touch $report_file
printf "*De novo* Contig Report\n" > $report_file
printf "==================================\n\n" >> $report_file
printf "## ${sample}\n\n" >> $report_file
printf "Forward reads file: ${fread}\n\n" >> $report_file
printf "Reverse reads file: ${rread}\n\n" >> $report_file
printf "* * * *\n\n" >> $report_file
printf "## Summary Table\n\n" >> $report_file
printf "| Metric                                       | Value               |\n" >> $report_file
printf "| :------------------------------------------- | :------------------ |\n" >> $report_file
printf "| Number of original reads:                    | $t_read_i               |\n" >> $report_file
printf "| Number of polished reads:                    | $t_read_f               |\n" >> $report_file
printf "| Human genome read contamination:             | $contam1 %%              |\n" >> $report_file
printf "| *S. cerevisiae* genome read contamination:   | $contam2 %%              |\n" >> $report_file
printf "| *S. paradoxus* genome read contamination:    | $contam3 %%              |\n\n" >> $report_file

# contig image
printf "![](img/${sample}_contigs.jpeg){ width=40%% }\n\n" >> $report_file
printf "* * * *\n\n" >> $report_file

# convert blastn results to markdown format and add top hits to report
# pip install csv2md (to install csv2mk)
csv2md $blast_output > ${outdir}/dsrna_blast.md
printf '%s\n\n' "### BLASTn Results for First 10 Contigs in $(basename ${contig_file}):" >> $report_file
head -n 2 ${outdir}/dsrna_blast.md >> $report_file
for i in `seq 10`; do
	grep "NODE_${i}_" ${outdir}/dsrna_blast.md | head -n 1 >> $report_file
done
printf "* * * *\n\n" >> $report_file

# kmer plots
printf '%s\n\n' "## All K-mer Contig Plots:" >> $report_file
for img in ${outdir}/img/${sample}_K*.jpeg
do
	printf "![$(basename $img)](img/$(basename $img)){ width=35%% }\n\n" >> $report_file
done

# print all the versions of the programs used in this script
printf "* * * *\n\n" >> $report_file
printf '%s\n\n' "### CL Applictions Used:" >> $report_file
printf '%s\n\n' "- ${spades_v}" >> $report_file
printf '\t\t%s\n\n' "$(head -n 1 ${outdir}/spades_output/spades.log | sed -n "s/^.\+\(spades.\+\)-1.\+/\1/p")" >> $report_file
printf '\t\t%s\n\n' "$(head -n 1 ${outdir}/metaspades_output/spades.log | sed -n "s/^.\+\(spades.\+\)-1.\+/\1/p")" >> $report_file
printf '%s\n\n' "- ${java_v}" >> $report_file
printf '%s\n\n' "- ${python_v}" >> $report_file
printf '%s\n\n' "- ${rscript_v}" >> $report_file
printf '%s\n\n' "- ${blast_v}" >> $report_file
printf '%s\n\n' "- ${bbsplit_v}" >> $report_file
printf '%s\n\n' "- ${fastp_v}" >> $report_file
printf '\t\t%s\n\n' "$(tail -n 1 ${fastp_html} | sed -n "s/^<div id='footer'> <p>\(.\+\)-i.\+/\1/p")" >> $report_file
printf '%s\n\n' "- ${fastqc_v}" >> $report_file

printf '\n%s\n\n' "*Script by Angela Crabtree*" >> $report_file

# convert markdown report to html using pandoc
pandoc --from markdown --to html $report_file > $report_html

# remove unnecessary files
rm ${outdir}/${sample}_R1_filt.fastq.gz
rm ${outdir}/${sample}_R2_filt.fastq.gz
rm ${outdir}/${sample}_dsrna_blast.csv.og
rm ${outdir}/dirty_reads_Hsap.fastq.gz
rm ${outdir}/dirty_reads_Scer.fastq.gz
rm ${outdir}/dirty_reads_Spara.fastq.gz
rm ${outdir}/dsrna_blast.md
rm -r ${outdir}/spades_output
rm -r ${outdir}/metaspades_output

printf "\n\tAll done!\n\n"