# allelic-ChrRNAseq
**R functions and scripts for allelic analysis of ChrRNA-seq data sets for investigating X chromosome silencing**


This repository provides numerous functions written in R which streamline the allele-specific analysis of Chromatin RNA-seq data sets in order to quantify Xist-mediated silencing.  


These functions and scripts are designed to take as inputs '.count' matrices. These are created from our UNIX command-line pipeline, which does the following:
1. Individual RNA-seq libraries are mapped by [STAR](https://github.com/alexdobin/STAR).
2. Libraries are allelically split by [SNPsplit](https://github.com/FelixKrueger/SNPsplit).  
3. Unsplit, genome1, and genome2 .bam files are sorted using samtools
4. We then run [featureCounts](http://bioconductor.org/packages/release/bioc/html/Rsubread.html) using a .gtf file of non-redundant genes. This generates the counts table in which columns 1-6 contain gene information and columns 7-9 specify the read/fragment counts from unsplit, genome1, and genome2 .bam files respectively. 
