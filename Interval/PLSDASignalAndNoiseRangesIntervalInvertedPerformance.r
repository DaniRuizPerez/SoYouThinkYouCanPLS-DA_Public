#########################################################################################################################
#########################################################################################################################
###### PAPER:          So you think you can PLS-DA?
###### NAME:           PLSDASignalAndNoiseRangesIntervalInvertedPerformance.r
###### AUTHOR:         Daniel Ruiz-Perez, PhD Student
###### AFFILIATION:    Florida International University
#########################################################################################################################
#########################################################################################################################

rm(list=ls())

library(mixOmics)
library(plotly)


ncomp = 1 #Number of components for plsda
perturbation = 0 #The perturbtion to add to the signal
#Values to iterate over the features lenght and signal length
noiseSeq = seq(1, 1000, by = 75)
samplesSeq = seq(10, 150, by = 10)
nRepetitions = 100
repetitionsSeq = seq(1, nRepetitions, by = 1)
nSignal = 10

outputMatrixPLSDA = array(0, dim=c(length(noiseSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 

#We execute with different matrixes to minimize random effects
for (repetition in repetitionsSeq) {
  print(as.numeric(repetition/nRepetitions)*100)
  
  #We execute with a random seed
  set.seed(sample(1:100, 1))
  set.seed(runif(1,0,10))
  
  outputMatrixPLSDAAux = array(0, dim=c(length(noiseSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 
  
  
  # limitsAlsoLeft = c(0,1/3,2/3,1)
  limitsAlsoLeft = c(0,0.2,0.4,0.6,0.8,1)
  #limitsAlsoLeft = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
  
  
  OriginalM = c()
  for (j in seq(1:(noiseSeq[length(noiseSeq)]))) {
    OriginalM = cbind(runif(samplesSeq[length(samplesSeq)],limitsAlsoLeft[j%%10+1],limitsAlsoLeft[j%%10+2]),OriginalM)
  }
  
  
  #We iterate over different number of features
  noiseCount = 0
  for (nNoise in noiseSeq) {
    noiseCount = noiseCount +1
    
    samplesCount = 0
    #We iterate over different number of signals
    for (nSamples in samplesSeq) {
      samplesCount = samplesCount +1
      
      #We create the matrix, sampling from the original
      if (nNoise == 0){
        m = matrix(nrow = nSamples, ncol = 0 )
        
      }
      else{
        m = OriginalM[1:nSamples,1:nNoise]
      }
      
      #Array that contains the limits that can take the variables to be classified as class 1
 
      
      Y <- rep(1:0, each=nSamples/2)
      
      
      
      for (j in seq(1:(nSignal))) {
        m = cbind(c(runif(nSamples/2,0,1),runif(nSamples/2,limitsAlsoLeft[j%%10+1],limitsAlsoLeft[j%%10+2])),m)
      }
      
      
      X <- as.data.frame(m)
      #As the first and second half of the colums have different signal, the classes are different for them
      
      plsda <- splsda(X, Y, ncomp = ncomp, scale = TRUE)
      
      perf <- perf(plsda, validation = "Mfold",dist = "all", folds = 5, progressBar = FALSE, nrepeat = 10,auc = TRUE)
      
      outputMatrixPLSDAAux[noiseCount,samplesCount] = as.numeric(perf$auc$`comp 1`[1])
      
      
      
    }
    # setwd("..")
  }
  outputMatrixPLSDA = outputMatrixPLSDA + outputMatrixPLSDAAux
}

outputMatrixPLSDA = outputMatrixPLSDA/nRepetitions




plsdaP <- plot_ly(y = noiseSeq, x = samplesSeq, z = outputMatrixPLSDA,showlegend=TRUE, name="PLSDA") %>% 
  layout(title = "Performance of PLSDA vs PCA",
         annotations = list(x = 1.1,y = 0.30,text = 'PLSDA',showarrow = FALSE),
         scene = list(
           xaxis = list(title = "# Signal"), 
           yaxis = list(title = "# Noise"), 
           zaxis = list(title = "Performance")
         )) %>% add_surface()
plsdaP




title = paste("IntervalDataIntervalsCVPerformanceInverted","nRep",nRepetitions,".html",sep="")
Sys.setenv("plotly_username" = "druiz072")
Sys.setenv("plotly_api_key" = "Gnv8DXsx4cKYFt8uVD90")
setwd("D:/Google Drive/FIU/RESEARCH/PLS-DA/graphs")
htmlwidgets::saveWidget(plsdaP,title)


