#########################################################################################################################
#########################################################################################################################
###### PAPER:          So you think you can PLS-DA?
###### NAME:           PLSDASignalAndNoiseCLUSTERClusterSizeSeparation.r
###### AUTHOR:         Daniel Ruiz-Perez, PhD Student
###### AFFILIATION:    Florida International University
#########################################################################################################################
#########################################################################################################################


rm(list=ls())

library(mixOmics)
library(plotly)
library(clusterGeneration)


ncomp = 1 #Number of components for plsda
separationSeq = seq(-0.9,0.9 , by = 0.1) #Never more than the samples
samplesSeq = seq(10, 100010, by = 5000)
nRepetitions = 100
repetitionsSeq = seq(1, nRepetitions, by = 1)
nSignal = 10
nNoise = 200

outputMatrixPCA = array(0, dim=c(length(separationSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 
outputMatrixPLSDA = array(0, dim=c(length(separationSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 

#We execute with different matrixes to minimize random effects
for (repetition in repetitionsSeq) {
  print(as.numeric(repetition/nRepetitions)*100)
  
  outputMatrixPCAAux = array(0, dim=c(length(separationSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 
  outputMatrixPLSDAAux = array(0, dim=c(length(separationSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 

 
  start.time <- Sys.time()

  
  #We iterate over different number of features
  separatioCount = 0
  for (nSeparation in separationSeq) {
    separatioCount = separatioCount +1
    # dir.create(paste(nSamples,"x",nSeparation))
    # setwd(paste(nSamples,"x",nSeparation))
    print (nSeparation)
    
    samplesCount = 0
    #We iterate over different number of signals
    for (nSamples in samplesSeq) {
      samplesCount = samplesCount +1

      #We create the data
      c <- genRandomClust(numClust=2, sepVal=nSeparation, numNonNoisy=nSignal, numNoisy=nNoise, 
                          numOutlier=0, numReplicate=1,clustszind= 1,clustSizeEq=(nSamples+250))

      X <- as.data.frame(c$datList$test_1[1:(nSamples*2),1:(nSignal+nNoise)])
      Y <- c$memList$test_1[1:(nSamples*2)]
      
      
      
      plsda <- splsda(X, Y, ncomp = ncomp, scale = TRUE)

      pca <- pca(X, ncomp = ncomp)
   
      #We obain the variables that carry the signal
      signalVariables = setdiff(1:ncol(X),c$noisyList$test_1)
      
      #We compute the performance of the models as the percentage of signal features in the important features
      importantVariablesPCA = selectVar(pca, comp = 1)$name[1:nSignal]
      importantVariablesPLSDA = selectVar(plsda, comp = 1)$name[1:nSignal]
      for (n in signalVariables) {
        if (paste('x',n,sep = "") %in% importantVariablesPCA){
          outputMatrixPCAAux[separatioCount,samplesCount] = outputMatrixPCAAux[separatioCount,samplesCount] + 1/nSignal
        }
        if (paste('x',n,sep = "") %in% importantVariablesPLSDA){
          outputMatrixPLSDAAux[separatioCount,samplesCount] = outputMatrixPLSDAAux[separatioCount,samplesCount] + 1/nSignal
        }
      }
    }
  }
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
  outputMatrixPLSDA = outputMatrixPLSDA + outputMatrixPLSDAAux
  outputMatrixPCA = outputMatrixPCA + outputMatrixPCAAux
}

outputMatrixPLSDA = outputMatrixPLSDA/nRepetitions
outputMatrixPCA = outputMatrixPCA/nRepetitions


pcaP <- plot_ly(colors=c("dark blue","orange") ,showlegend=TRUE, name = "PCA", y = separationSeq,x = samplesSeq, z = outputMatrixPCA)%>%  
  layout(title = "Performance of PCA",
         annotations = list(x = 2.17,y = 1.03,text = 'PCA',showarrow = FALSE),
         scene = list(
           xaxis = list(title = "Samples"), 
           yaxis = list(title = "# Separation"), 
           zaxis = list(title = "Performance"))) %>% add_surface()
#pcaP

plsdaP <- plot_ly(y = separationSeq, x = samplesSeq, z = outputMatrixPLSDA,showlegend=TRUE, name="PLSDA") %>% 
  layout(title = "Performance of PLSDA vs PCA",
         annotations = list(x = 1.1,y = 0.30,text = 'PLSDA',showarrow = FALSE),
         scene = list(
           xaxis = list(title = "Samples"), 
           yaxis = list(title = "# Separation"), 
           zaxis = list(title = "Performance")
         )) %>% add_surface()
#plsdaP


p <- subplot(pcaP,plsdaP)
p

