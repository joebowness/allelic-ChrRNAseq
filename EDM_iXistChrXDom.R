### REQUIRED FILES FOR THE EXPONENTIAL DECAY MODEL ####
# 1. An allelic ratio table, with columns from 7 onwards samples of the ES-to-NPC differentiation timecourse to fit to an exponential decay model
# 2. A list of escapee genes, as defined from the final NPC timepoints 

library(gplots)

### SET-UP OF MODEL (for iXist-ChrX-Dom) ####
AR_table_forEDM <- A11B2_AR_table #this is an allelic ratio table from two replicates of ES-to-NPC differentiation timecourse in iXist-ChrX-Dom cells
NPC_escapees <- A11B2_NPC_escapees #this is a list of escapee genes in iXist-ChrX-Dom cells, defined as genes with an allelic ratio > 0.1 in an average of two NPC chrRNA-seq samples
AR_table_forEDM_merged <- A11B2_AR_table_merged #this the allelic ratio table with replicates of the same timepoint averaged together

options(scipen=999)
# From column names, collect timepoints of samples in the allelic ratio table
colnames(AR_table_forEDM)
timepoints <- c(0,0,24,24,48,48,72,72,120,120,168,168,360,504)
colnames(AR_table_forEDM)[7:ncol(AR_table_forEDM)] <- timepoints

# Set some approximate initial parameter estimates for running nls on whole dataset
y_f <- 0.1
y_0 <- 0.4
theta <- 2
TimeFrac <- 0.5 # set to 0.5 to calculate to halftime

### RUN MODEL ####
#estimate initial parameters from running whole dataset with model
AR_table_forEDM_melt <- melt(AR_table_forEDM, id=c("Geneid","Chr","Start","End","Strand","Length"), stringsAsFactors = FALSE)
colnames(AR_table_forEDM_melt)[7:8] <- c("time", "ratio")
AR_table_forEDM_melt$time <- as.numeric(as.character(AR_table_forEDM_melt$time))
plot(AR_table_forEDM_melt$time,AR_table_forEDM_melt$ratio, ylim=c(0,0.8))
overall_lm <- nlsLM(ratio ~ yf + y0*exp(-time/th),
                    data = AR_table_forEDM_melt[,7:8], start = list(y0 = y_0, yf = y_f, th = theta))
coef(overall_lm)
y_0 <- coef(overall_lm)[1]
y_f <- coef(overall_lm)[2]
theta <- coef(overall_lm)[3]
curve((y_f + y_0*exp(-x/theta)), add=TRUE, col="red" )
#Calculate overall halftime
-theta*log2((TimeFrac*(y_0+y_f)-y_f)/y_0)

#simple boxplot of all samples to input into model, check agree with defined timpeoints
AR_matrix_forEDM <- as.matrix(AR_table_forEDM[7:ncol(AR_table_forEDM)])
par(mfrow=c(1,1))
boxplot.matrix(AR_matrix_forEDM, use.cols=TRUE,
               ylab = 'chrRNA Allelic Ratio (Xi/(Xi+Xa))',
               ylim=c(0,1),
               names = colnames(AR_matrix_forEDM))
abline(h=0.5)

#model
par(mfrow=c(3,3))
OutTableEDM <- data.frame()
data_2D <- data.frame() # for plotting overall data in 2D
res_table <- data.frame() #for testing output of overall fit
EDM_FitValues <- data.frame() #for testing output of overall fit
for (row in 1:nrow(AR_matrix_forEDM)){
  #set up for single genes, plot data points 
  Geneid <- AR_table_forEDM[row,1]
  Start <- AR_table_forEDM[row,3]
  End <- AR_table_forEDM[row,3]
  gene <- as.data.frame(cbind(timepoints,AR_matrix_forEDM[row,]))
  colnames(gene) <- c('time','ratio')
  ratio <- gene$ratio
  time <- gene$time
  #data_2D <- rbind(data_2D,gene)
  #plot(log(time),ratio) #plot log-log
  plot(time,ratio, ylim=c(0,0.8), xlab = Geneid)
  #remove problematic timepoints for special genes
  if (Geneid == "Mbnl3"){gene <- subset(gene,timepoints != 24)}
  if (Geneid == "2610018G03Rik"){gene <- subset(gene,timepoints < 360)}
  #do model to non-zero asymptote if an escapee 
  if (Geneid %in% NPC_escapees$Geneid){
    cc <- try(nlsLM(ratio ~ yf + y0*exp(-time/th),
                    data = gene, start = list(y0 = y_0, yf = y_f, th = theta)))
    if (is(cc,"try-error")) {
      EDM_start <- NA
      EDM_final <- NA
      EDM_rate <- NA
      ResStError <- NA
      EDM_halftime <- NA
    } 
    else {
      #perform model fitting
      EDMfit <- nlsLM(ratio ~ (yf + y0*exp(-time/th)),
                      data = gene, start = list(y0 = y_0, yf = y_f, th = theta))
      gene_y0 <- coef(EDMfit)[1]
      gene_yf <- coef(EDMfit)[2]
      gene_th <- coef(EDMfit)[3]
      curve((gene_yf + gene_y0*exp(-x/gene_th)), add=TRUE, col="red" )
      #Compute model information
      residuals <- residuals(EDMfit)
      fitval <- fitted.values(EDMfit)
      res_add <- cbind(fitval,residuals)
      res_table <- rbind(res_table,res_add)
      #points(fitval,residuals)
      ResStError <- sqrt(sum((residuals(EDMfit)*(residuals(EDMfit)))))/2
      #Rate parameter
      EDM_start <- as.numeric(gene_y0) + as.numeric(gene_yf)   #initial AR
      abline(h=EDM_start)
      abline(h=EDM_start/2)
      EDM_final <- as.numeric(gene_yf)
      abline(h=EDM_final)
      EDM_halftime <- -gene_th*log((TimeFrac*(gene_y0+gene_yf)-gene_yf)/gene_y0)
      EDM_halftime[is.nan(EDM_halftime)] <- NA
      abline(v = EDM_halftime)
      if (EDM_final < 0){EDM_final <- 0} #set final AR to 0 if the asymptote is below zero
      EDM_rate <- as.numeric(gene_th) #rate constant
      #make NA if slope reversed or AR increases
      if (EDM_rate < 0 | EDM_final > EDM_start){EDM_start <- EDM_final <- EDM_rate <- ResStError <- EDM_halftime <- NA}
      #make everything NA if poor model fit
      #if (ResStError > 0.1){ EDM_start <- EDM_final <- EDM_rate <- ResStError <- EDM_halftime <- NA} 
      #make NA if can't produce halftime
      if (is.na(EDM_halftime) == TRUE | EDM_halftime > 500) {EDM_start <- EDM_final <- EDM_rate <- ResStError <- EDM_halftime <- NA} 
    }
  }
  #do model to completion if a silenced gene
  else{
    cc <- try(nlsLM(ratio ~ y0*exp(-time/th),
                    data = gene, start = list(y0 = y_0, th = theta)))
    if (is(cc,"try-error")) {
      EDM_start <- NA
      EDM_final <- NA
      EDM_rate <- NA
      ResStError <- NA
      EDM_halftime <- NA
    } 
    else {
      #perform model fitting
      EDMfit <- nlsLM(ratio ~ (y0*exp(-time/th)),
                      data = gene, start = list(y0 = y_0, th = theta))
      gene_y0 <- coef(EDMfit)[1]
      gene_th <- coef(EDMfit)[2]
      curve((gene_y0*exp(-x/gene_th)), add=TRUE, col="red" )
      #Compute model information
      residuals <- residuals(EDMfit)
      fitval <- fitted.values(EDMfit)
      res_add <- cbind(fitval,residuals)
      res_table <- rbind(res_table,res_add)
      #points(fitval,residuals)
      ResStError <- sqrt(sum((residuals(EDMfit)*(residuals(EDMfit)))))/2
      #Show relevent measures on graph
      EDM_start <- as.numeric(gene_y0)   #initial AR
      abline(h=EDM_start)
      abline(h=EDM_start/2)
      EDM_final <- 0
      abline(h=EDM_final)
      EDM_halftime <- -gene_th*log((TimeFrac*(gene_y0))/gene_y0)
      EDM_halftime[is.nan(EDM_halftime)] <- NA
      abline(v = EDM_halftime)
      EDM_rate <- as.numeric(gene_th) #rate constant 
      #make NA if poor model fit
      #if (ResStError > 0.1){EDM_start <- EDM_final <- EDM_rate <- ResStError <- EDM_halftime <- NA} 
    }
  }
  
  EDMinfo <- cbind(Geneid, Start, EDM_start, EDM_final, EDM_rate, EDM_halftime, ResStError)
  #EDMpeaklist <- rbind(EDMpeaklist,peakname,stringsAsFactors=FALSE)
  OutTableEDM <- rbind(OutTableEDM,EDMinfo)
  OutTableEDM[,3:7] <- as.numeric(unlist(OutTableEDM[,3:7]))
  EDM_FitValues <- rbind.fill(EDM_FitValues, as.data.frame(t(fitval)))
}

### MODEL OUTPUT ####
#band plot of residuals
par(mfrow=c(1,1))
bandplot(res_table$fitval,res_table$residuals,cex=0.6)
abline(h=0,col="green") 
abline(v=median(AR_table_forEDM_merged$'0'),col="green")
abline(v=0.5*(median(AR_table_forEDM_merged$'0')),col="green")

#compare boxplots of real data to model data
par(mfrow=c(1,2))
AR_matrix_forEDM_merged <- as.matrix(AR_table_forEDM_merged[7:ncol(AR_table_forEDM_merged)])
boxplot.matrix(AR_matrix_forEDM_merged, use.cols=TRUE,
               ylab = 'chrRNA Allelic Ratio (Xi/(Xi+Xa))',
               xlab = 'real data',
               ylim=c(0,1))
abline(h=0.5)
EDM_Fit_Matrix <- as.matrix(EDM_FitValues)
colnames(EDM_Fit_Matrix) <- timepoints
boxplot.matrix(EDM_Fit_Matrix[,c(1,3,5,7,9,11,13,14)], use.cols=TRUE,
               ylab = 'chrRNA Allelic Ratio (Xi/(Xi+Xa))',
               xlab = 'model data',
               ylim=c(0,1))
abline(h=0.5)
#summary values and density of halftimes
summary(OutTableEDM$EDM_halftime)
par(mfrow=c(1,1))
plot(density(OutTableEDM$EDM_halftime, na.rm = TRUE,adjust = 0.8), xlim=c(0,150), ylim=c(0,0.03), xaxs="i",yaxs="i") 
#define slow, medium and fast-silencing genes
fast_threshold <- quantile(OutTableEDM$EDM_halftime, 0.33333, na.rm=TRUE)
slow_threshold <- quantile(OutTableEDM$EDM_halftime, 0.66666, na.rm=TRUE)
abline(v=(slow_threshold))
abline(v=(fast_threshold))

slow_silencing_genes <- subset(OutTableEDM,
                                     OutTableEDM$EDM_halftime > slow_threshold)
medium_silencing_genes <- subset(OutTableEDM,
                                       OutTableEDM$EDM_halftime < slow_threshold & OutTableEDM$EDM_halftime > fast_threshold)
fast_silencing_genes <- subset(OutTableEDM,
                                     OutTableEDM$EDM_halftime < fast_threshold)
