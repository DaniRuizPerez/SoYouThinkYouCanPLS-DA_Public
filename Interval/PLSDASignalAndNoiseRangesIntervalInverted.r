#########################################################################################################################
#########################################################################################################################
###### PAPER:          So you think you can PLS-DA?
###### NAME:           PLSDASignalAndNoiseRangesIntervalInverted.r
###### AUTHOR:         Daniel Ruiz-Perez, PhD Student
###### AFFILIATION:    Florida International University
#########################################################################################################################
#########################################################################################################################

rm(list=ls())

library(mixOmics)
library(plotly)


nSamples = 200 #Number of samples
ncomp = 1 #Number of components for plsda
perturbation = 0 #The perturbtion to add to the signal
#Values to iterate over the features lenght and signal length
noiseSeq = seq(100, 400, by = 25)
singalSeq = seq(1, 30, by = 2)
nRepetitions = 100
repetitionsSeq = seq(1, nRepetitions, by = 1)

outputMatrixPCA = array(0, dim=c(length(noiseSeq),length(singalSeq))) #First dimension number of features, second number of signal, 
outputMatrixPLSDA = array(0, dim=c(length(noiseSeq),length(singalSeq))) #First dimension number of features, second number of signal, 

#We execute with different matrixes to minimize random effects
for (repetition in repetitionsSeq) {
  print(as.numeric(repetition/nRepetitions)*100)
  
  #We execute with a random seed
  set.seed(sample(1:100, 1))
  set.seed(runif(1,0,10))
  
  outputMatrixPCAAux = array(0, dim=c(length(noiseSeq),length(singalSeq))) #First dimension number of features, second number of signal, 
  outputMatrixPLSDAAux = array(0, dim=c(length(noiseSeq),length(singalSeq))) #First dimension number of features, second number of signal, 
  
  limitsAlsoLeft = c(0,1/3,2/3,1)
  # limitsAlsoLeft = c(0,0.2,0.4,0.6,0.8,1)
  # limitsAlsoLeft = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
  

  OriginalM = c()
  for (j in seq(1:(noiseSeq[length(noiseSeq)]))) {
    OriginalM = cbind(runif(nSamples,limitsAlsoLeft[j%%3+1],limitsAlsoLeft[j%%3+2]),OriginalM)
  }
  
  
  #We iterate over different number of features
  noiseCount = 0
  for (nNoise in noiseSeq) {
    noiseCount = noiseCount +1
   
    signalCount = 0
    #We iterate over different number of signals
    for (nSignal in singalSeq) {
      signalCount = signalCount +1
  
      #We create the matrix, sampling from the original
      if (nNoise == 0){
        m = matrix(nrow = nSamples, ncol = 0 )

      }
      else{
      m = OriginalM[1:nSamples,1:nNoise]
      }
      
      Y <- rep(1:0, each=nSamples/2)
      
      
      for (j in seq(1:(nSignal))) {
        m = cbind(c(runif(nSamples/2,0,1),runif(nSamples/2,limitsAlsoLeft[j%%3+1],limitsAlsoLeft[j%%3+2])),m)
      }
      
      X <- as.data.frame(m)
      #As the first and second half of the colums have different signal, the classes are different for them
      
      plsda <- splsda(X, Y, ncomp = ncomp, scale = TRUE)
        
      # #PCA
      pca <- pca(X, ncomp = ncomp)
 

      #We compute the performance of the models as the percentage of signal features in the important features
      importantVariablesPCA = selectVar(pca, comp = 1)$name[1:nSignal]
      importantVariablesPLSDA = selectVar(plsda, comp = 1)$name[1:nSignal]
      for (n in seq(1:(nSignal))) {
        if (paste('V',n,sep = "") %in% importantVariablesPCA){
          outputMatrixPCAAux[noiseCount,signalCount] = outputMatrixPCAAux[noiseCount,signalCount] + 1/nSignal
        }
        if (paste('V',n,sep = "") %in% importantVariablesPLSDA){
          outputMatrixPLSDAAux[noiseCount,signalCount] = outputMatrixPLSDAAux[noiseCount,signalCount] + 1/nSignal
        }
      }

     
      
      
      
    }
    # setwd("..")
  }
  outputMatrixPLSDA = outputMatrixPLSDA + outputMatrixPLSDAAux
  outputMatrixPCA = outputMatrixPCA + outputMatrixPCAAux
}

outputMatrixPLSDA = outputMatrixPLSDA/nRepetitions
outputMatrixPCA = outputMatrixPCA/nRepetitions





pcaP <- plot_ly(colors=c("dark blue","orange") ,showlegend=TRUE, name = "PCA", y = noiseSeq,x = singalSeq, z = outputMatrixPCA)%>%  
  layout(title = "Performance of PCA",
         annotations = list(x = 2.17,y = 1.03,text = 'PCA',showarrow = FALSE),
         scene = list(
           xaxis = list(title = "# Signal"), 
           yaxis = list(title = "# Noise"),
           zaxis = list(title = "Performance"))) %>% add_surface()
#pcaP

plsdaP <- plot_ly(y = noiseSeq, x = singalSeq, z = outputMatrixPLSDA,showlegend=TRUE, name="PLSDA") %>% 
  layout(title = "Performance of PLSDA vs PCA",
         annotations = list(x = 1.1,y = 0.30,text = 'PLSDA',showarrow = FALSE),
         scene = list(
           xaxis = list(title = "# Signal"), 
           yaxis = list(title = "# Noise"), 
           zaxis = list(title = "Performance")
         )) %>% add_surface()
#plsdaP


p <- subplot(pcaP,plsdaP)
p


