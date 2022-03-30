# allelic-ChrRNAseq
**R functions and scripts for allelic analysis of ChrRNA-seq data sets for investigating X chromosome silencing**


This repository provides numerous functions written in R which streamline the allele-specific analysis of Chromatin RNA-seq data sets in order to quantify Xist-mediated silencing.  


These functions and scripts are designed to take as inputs '.count' matrices produced by [featureCounts](http://bioconductor.org/packages/release/bioc/html/Rsubread.html).  In our UNIX command-line pipeline, RNA-seq libraries are mapped by [STAR](https://github.com/alexdobin/STAR), filtered and sorted, then allelically split by [SNPsplit](https://github.com/FelixKrueger/SNPsplit).  We then run featureCounts using a .gtf file of a list of non-redundant genes and generate counts files in which columns 7-9 file specify the read/fragment counts from unsplit, genome1, and genome2 .bam files respectively. 
