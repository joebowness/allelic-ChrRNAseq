#### Functions ####
LoadingCountFile <- function(file=""){
  Count <- read.table(file, header=T)
  return(Count)
}
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
Allelic.Ratio <- function(g1_counts,g2_counts){
  Xa <- as.numeric(g1_counts)
  Xi <- as.numeric(g2_counts)
  AllelicRatio <- ( Xi / (Xi + Xa))
  return(AllelicRatio)
}
Allelic.Ratio.Reverse <- function(g1_counts,g2_counts){
  Xi <- as.numeric(g1_counts)
  Xa <- as.numeric(g2_counts)
  AllelicRatio <- ( Xi / (Xi + Xa))
  return(AllelicRatio)
}
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
AR.table <- function(counts_table){
  AllelicRatio_table <- counts_table[1:6]
  for (j in seq(8,ncol(counts_table),3)){
    AllelicRatio_table <- cbind(AllelicRatio_table,
                                Allelic.Ratio(as.matrix(counts_table[j]),
                                              as.matrix(counts_table[j+1]) ) )
  }
  return(AllelicRatio_table)
}
AR.table.reverse <- function(counts_table){
  AllelicRatio_table <- counts_table[1:6]
  for (j in seq(8,ncol(counts_table),3)){
    AllelicRatio_table <- cbind(AllelicRatio_table,
                                Allelic.Ratio.Reverse(as.matrix(counts_table[j]),
                                              as.matrix(counts_table[j+1]) ) )
  }
  return(AllelicRatio_table)
}
initialAR.filter <- function(AR_table,lower,upper){
  AR_table_filter1 <- subset(AR_table,
                                       lower <= AR_table[7])
  AR_table_filter2 <- subset(AR_table_filter1,
                                      upper >= AR_table_filter1[7])
  return(AR_table_filter2)
}
chrRNA.RPM <- function(counts_table){
  #Calculate sum of counts for each sample and divide through (per million) to make RPM
  countsSums <- colSums(counts_table[7:ncol(counts_table)])
  GeneNormFactors <- (rep(countsSums[seq(1,length(countsSums),3)], each=3))/1000000
  chrRNA_RPM <- cbind(counts_table[1:6], sweep(counts_table[7:ncol(counts_table)], 2, GeneNormFactors, "/"))
  return(chrRNA_RPM)
}

