#########################################################################################################################
#########################################################################################################################
###### PAPER:          So you think you can PLS-DA?
###### NAME:           PLSDASignalAndNoisePerformance.r
###### AUTHOR:         Daniel Ruiz-Perez, PhD Student
###### AFFILIATION:    Florida International University
#########################################################################################################################
#########################################################################################################################

rm(list=ls())

library(mixOmics)
library(plotly)


nSamples = 200 #Number of samples
constantFirstHalf = 0.5 #value that the sum of the signal of the first half of samples add up to
constantSecondHalf = -0.5 #value that the sum of the signal of the second half of samples add up to
ncomp = 1 #Number of components for plsda
perturbation = 0 #The perturbtion to add to the signal
#Values to iterate over the features lenght and signal length
noiseSeq = seq(0, 1000, by = 75)
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

  #values are %of good yfeatures in the first #signal values
  OriginalM = matrix(rnorm(nSamples*(noiseSeq[length(noiseSeq)]+nSignal)),
                     ncol=(noiseSeq[length(noiseSeq)]+nSignal)) #To avoid having a new matrix with every execution, we just 
  #samples this original one

  #We iterate over different number of features
  noiseCount = 0
  for (nNoise in noiseSeq) {
    noiseCount = noiseCount +1
    # dir.create(paste(nSamples,"x",nNoise))
    # setwd(paste(nSamples,"x",nNoise))
    
    samplesCount = 0
    #We iterate over different number of signals
    for (nSamples in samplesSeq) {
      samplesCount = samplesCount +1
      
      # name = paste(nSamples,"x",nNoise,"(sxf)", "sig", nSignal,"comp",ncomp,"pert",perturbation, sep = "") #Name of the output files
      # nametxt = paste(nSamples,"x",nNoise,"(sxf)","comp",ncomp,"pert",perturbation, sep = "") #Name of the output files
      
      #We create the matrix, sampling from the original
      m = OriginalM[1:nSamples,1:(nNoise+nSignal)]

      
      #We make the signal features to add up to a constant iff class 1 and another iff class 0
      for (i in seq(1:nSamples)) {
        rSig <- rnorm(nSignal)
        rSig <- rSig/sum(rSig)
        for (j in seq(1:nSignal)) {
          if (i <= (nSamples/2)){
            m[i,j] <- rSig[j]*constantFirstHalf + runif(1,-perturbation,perturbation)
          }else {
            m[i,j] <- rSig[j]*constantSecondHalf+ runif(1,-perturbation,perturbation)
          }
        }
      }
      
      X <- as.data.frame(m)
      #As the first and second half of the colums have different signal, the classes are different for them
      Y <- rep(1:0, each=nSamples/2)
      
      plsda <- splsda(X, Y, ncomp = ncomp, scale = TRUE)
   
      
      perf <- perf(plsda, validation = "Mfold",dist = "all", folds = 5, progressBar = FALSE, nrepeat = 10,auc = TRUE)
      
      outputMatrixPLSDAAux[noiseCount,samplesCount] = as.numeric(perf$auc$`comp 1`[1])
      
    }
    # setwd("..")
  }
  outputMatrixPLSDA = outputMatrixPLSDA + outputMatrixPLSDAAux
}

outputMatrixPLSDA = outputMatrixPLSDA/nRepetitions





plsdaP <- plot_ly(y = noiseSeq, x = samplesSeq, z = outputMatrixPLSDA) %>%   layout(title = "Performance of PLSDA",
scene = list(
xaxis = list(title = "# Samples"), 
yaxis = list(title = "# Noise"), 
zaxis = list(title = "Performance"))) %>% add_surface()
plsdaP

