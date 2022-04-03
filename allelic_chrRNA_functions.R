#### Packages ####
library(reshape)
library(readr)
library(stringr)
library(ggplot2)
library(ggridges)
library(stringr)
library(plyr)
library(dplyr)
library(forcats)
library(minpack.lm)
library(mefa4)
library(DESeq2)
library(gridExtra)

#### Functions ####
#This function reads in .count matrices. 
#Upon reading individual count matrices, we create a "master" counts table appending the 3 counts columns (7:9) for each sample. 
LoadingCountFile <- function(file=""){
  Count <- read.table(file, header=T)
  return(Count)
}
#This function applies an "allelic filter" to only keep genes for which the number of allelically mapping reads exceeds a given cutoff. 
#We suggest a cut-off of 10 as reasonable
#Additionally, the 'fraction'. Set this to 1 to ensure only genes with 10 allelic reads in every sample analysed together, 
#or a smaller fraction (eg. 0.8) to reduce filtering stringency and enable analysis of more genes which may otherwise be filtered out by outlier samples (for example, from a sample sequenced at lower depth)
Allelic.Filter.OverFracSamples <- function(counts_table,cutoff,fraction){
  counts_table_g1 <- counts_table[,seq(8,ncol(counts_table),3)]
  counts_table_g2 <- counts_table[,seq(9,ncol(counts_table),3)]
  counts_g1g2 <- counts_table_g1 + counts_table_g2
  counts_table_filter <- data.frame()
  for (i in c(1:nrow(counts_g1g2))){
    pass <- length((which(counts_g1g2[i,] > cutoff) == "TRUE"))
    if (pass >= fraction*(ncol(counts_g1g2))){
      counts_table_filter <- rbind(counts_table_filter,counts_table[i,])
    }
  }
  return(counts_table_filter)
}
#Allelic ratio equation for iXist-ChrX-Dom lines, where g1_counts = Castaneous = Xa ad g2_counts = Domesticus/129 = Xi
Allelic.Ratio <- function(g1_counts,g2_counts){
  Xa <- as.numeric(g1_counts)
  Xi <- as.numeric(g2_counts)
  AllelicRatio <- ( Xi / (Xi + Xa))
  return(AllelicRatio)
}
#Allelic ratio equation for iXist-ChrX-Cast lines, where g1_counts = Castaneous = Xi ad g2_counts = Domesticus/129 = Xa
Allelic.Ratio.Reverse <- function(g1_counts,g2_counts){
  Xi <- as.numeric(g1_counts)
  Xa <- as.numeric(g2_counts)
  AllelicRatio <- ( Xi / (Xi + Xa))
  return(AllelicRatio)
}
#Restricts analysis to just genes on ChrX proximal to the Xist locus (0-103Mb), for iXist-ChrX-Dom lines
chrX1.filter <- function(counts_table){
  chrRNA_counts_table_chrX <- filter(counts_table, !grepl("chrX",counts_table$Chr))
  chrRNA_counts_table_chrX <-  counts_table[which(grepl("chrX",counts_table$Chr) == TRUE), ]
  numextract <- function(string){
    str_extract(string, "\\-*\\d+\\.*\\d*")
  }
  chrRNA_chrX1_filter_End <- cbind(chrRNA_counts_table_chrX[1], as.numeric(numextract(chrRNA_counts_table_chrX$End)))
  chrRNA_chrX1_filter_End <- chrRNA_chrX1_filter_End[which(chrRNA_chrX1_filter_End$`as.numeric(numextract(chrRNA_counts_table_chrX$End))` < 103483233), ]
  chrRNA_filter_chrX1 <- counts_table[
    which(counts_table$Geneid %in% chrRNA_chrX1_filter_End$Geneid ), ]
  return(chrRNA_filter_chrX1)
}
#For iXist-ChrX-Dom lines, applies the allelic ratio equation to a counts table of the form specified above to create an 'allelic ratio' table.
#Columns 7 onwards in this table are allelic ratios for each sample
AR.table <- function(counts_table){
  AllelicRatio_table <- counts_table[1:6]
  for (j in seq(8,ncol(counts_table),3)){
    AllelicRatio_table <- cbind(AllelicRatio_table,
                                Allelic.Ratio(as.matrix(counts_table[j]),
                                              as.matrix(counts_table[j+1]) ) )
  }
  return(AllelicRatio_table)
}
#For iXist-ChrX-Cast lines, applies the allelic ratio equation to a counts table of the form specified above.
AR.table.reverse <- function(counts_table){
  AllelicRatio_table <- counts_table[1:6]
  for (j in seq(8,ncol(counts_table),3)){
    AllelicRatio_table <- cbind(AllelicRatio_table,
                                Allelic.Ratio.Reverse(as.matrix(counts_table[j]),
                                              as.matrix(counts_table[j+1]) ) )
  }
  return(AllelicRatio_table)
}
#Filters out genes from allelic ratio tables which are heavily allelically skewed in untreated mESCs (where biallelic expression is expected). 
#Upper and lower allelic ratio thresholds can be set, we recommend 0.1 and 0.9, with the reasoning that extreme allelic skewing of certain genes in untreated samples are artifacts from misannotated/incongruous SNPs.
initialAR.filter <- function(AR_table,lower,upper){
  AR_table_filter1 <- subset(AR_table,
                                       lower <= AR_table[7])
  AR_table_filter2 <- subset(AR_table_filter1,
                                      upper >= AR_table_filter1[7])
  return(AR_table_filter2)
}
#Calculates the "reads per million" for each gene from a counts table 
chrRNA.RPM <- function(counts_table){
  #Calculate sum of counts for each sample and divide through (per million) to make RPM
  countsSums <- colSums(counts_table[7:ncol(counts_table)])
  GeneNormFactors <- (rep(countsSums[seq(1,length(countsSums),3)], each=3))/1000000
  chrRNA_RPM <- cbind(counts_table[1:6], sweep(counts_table[7:ncol(counts_table)], 2, GeneNormFactors, "/"))
  return(chrRNA_RPM)
}

