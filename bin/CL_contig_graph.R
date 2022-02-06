#!/usr/bin/env Rscript

main <- function() {
  ## open libraries
  packages <- c('ggplot2','stringr')
  # check.packages function: install and load multiple R packages.
  # Check to see if packages are installed. Install them if they are not, then load them into the R session.
  check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE, repo="http://cran.rstudio.com/")
    sapply(pkg, require, character.only = TRUE)
  }
  msg.trap <- capture.output( suppressMessages( check.packages(packages)) ) # silent mode
  
  ## parse arguments
  args <- commandArgs(trailingOnly = TRUE)
  ## parse arguments
  file <- args[1]  ## needed for importing data
  # option for output filename and type (argument #2)
  if (is.na(args[2])==FALSE){
    output_file <- args[2]
  } else {
    output_file <- paste(str_split(file, "[.]")[[1]][1],'.jpeg',sep="")
  }
  # graph height option (argument #3)
  if (is.na(args[3])==FALSE){
    grheight <- args[3]
  } else {
    grheight <- 5
  }
  # graph width option (argument #4)
  if (is.na(args[4])==FALSE){
    grwidth <- args[4]
  } else {
    grwidth <- 6
  }
  
  
  #****************** MAIN SCRIPT *******************
  
  ## Read data from pre-formatted dataset.
  outputname <- basename(output_file)
  samplename <- strsplit(outputname,"_")[[1]][[1]]
  contigs <- read.delim(file, header = FALSE, sep = '_') # parses file by underscore and newline
  contigs <- na.omit(contigs) # removes rows with 'NA'
  contigs$sample <- rep(samplename,nrow(contigs)) # name a new row with sample name
  contigs <- contigs[, c(7, 2, 4, 6)] # only include relevant columns
  row.names(contigs) <- 1:nrow(contigs) # renumber rows after reordering
  colnames(contigs) <- c('sample','node','len','cov') # rename column headers
  contigs <- contigs[which(contigs$len>150),] # delete contigs 150bp or less

  ############## AREA PLOT #################
  ggplot(contigs, aes(x=len, y=cov)) +
    coord_trans(y = "log10") +
    geom_point(size=3, alpha=0.5) +
    theme_classic() +
    labs(x="Contig Length", y="Contig Coverage", 
         title=paste(samplename," de novo Contigs")) +
    theme(axis.title = element_text(size = 16, face="bold"), 
          axis.text.y = element_text(size=14), 
          axis.text.x = element_text(size=14)) +
    geom_text(aes(label=ifelse(len>(max(len)*(1/2)),as.character(node),'')),
              hjust=1,vjust=1.5, color='red')
  
  ggsave(filename=output_file, height=as.numeric(grheight),
         width=as.numeric(grwidth))
  
}
main()