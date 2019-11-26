#########################################################################################################################
#########################################################################################################################
###### PAPER:          So you think you can PLS-DA?
###### NAME:           PLSDASignalAndNoiseCLUSTERPerformance.r
###### AUTHOR:         Daniel Ruiz-Perez, PhD Student
###### AFFILIATION:    Florida International University
#########################################################################################################################
#########################################################################################################################



#Now the data is generated with the package clustergeneration, where plsda does very good

rm(list=ls())

library(mixOmics)
library(plotly)
library(clusterGeneration)


ncomp = 1 #Number of components for plsda
nonSignalSeq = seq(0,1000 , by = 75) 
samplesSeq = seq(10, 150, by = 10)
nRepetitions = 100
repetitionsSeq = seq(1, nRepetitions, by = 1)
nSignal = 10

outputMatrixPLSDA = array(0, dim=c(length(nonSignalSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 

start.time <- Sys.time()


#We execute with different matrixes to minimize random effects
for (repetition in repetitionsSeq) {
  print(as.numeric(repetition/nRepetitions)*100)
  
  outputMatrixPLSDAAux = array(0, dim=c(length(nonSignalSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 

 
  #We iterate over different number of features
  noiseCount = 0
  for (nNoise in nonSignalSeq) {
    noiseCount = noiseCount +1
    print(nNoise)

    
    samplesCount = 0
    #We iterate over different number of signals
    for (nSamples in samplesSeq) {
      samplesCount = samplesCount +1
      print(nSamples)
      #We create the data
      
      
      c <- genRandomClust(numClust=2, sepVal=0.5, numNonNoisy=nSignal, numNoisy=nNoise, 
                          numOutlier=0, numReplicate=1,clustszind= 1,clustSizeEq=(nSamples+nNoise*(1.5)))
      
      
      X <- as.data.frame(c$datList$test_1[1:(nSamples*2),1:(nSignal+nNoise)])
      Y <- c$memList$test_1[1:(nSamples*2)]

      
      plsda <- splsda(X, Y, ncomp = ncomp, scale = TRUE)
      
      perf <- perf(plsda, validation = "Mfold",dist = "all", folds = 5, progressBar = FALSE, nrepeat = 10,auc = TRUE)

      outputMatrixPLSDAAux[noiseCount,samplesCount] = as.numeric(perf$auc$`comp 1`[1])
    }
  }
  outputMatrixPLSDA = outputMatrixPLSDA + outputMatrixPLSDAAux
}

outputMatrixPLSDA = outputMatrixPLSDA/nRepetitions




plsdaP <- plot_ly(y = nonSignalSeq, x = samplesSeq, z = outputMatrixPLSDA) %>%   layout(title = "Performance of PLSDA",
  scene = list(
  xaxis = list(title = "# Samples"), 
  yaxis = list(title = "# Noise"), 
  zaxis = list(title = "Performance"))) %>% add_surface()
plsdaP


end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

