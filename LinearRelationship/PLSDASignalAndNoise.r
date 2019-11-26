#########################################################################################################################
#########################################################################################################################
###### PAPER:          So you think you can PLS-DA?
###### NAME:           PLSDASignalAndNoise.r
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
noiseSeq = seq(100, 800, by = 30)
singalSeq = seq(1, 20, by = 1)
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

  #values are %of good yfeatures in the first #signal values
  OriginalM = matrix(rnorm(nSamples*(noiseSeq[length(noiseSeq)]+singalSeq[length(singalSeq)])),
                     ncol=(noiseSeq[length(noiseSeq)]+singalSeq[length(singalSeq)])) #To avoid having a new matrix with every execution, we just 
  #samples this original one

  #We iterate over different number of features
  noiseCount = 0
  for (nNoise in noiseSeq) {
    noiseCount = noiseCount +1
    # dir.create(paste(nSamples,"x",nNoise))
    # setwd(paste(nSamples,"x",nNoise))
    
    signalCount = 0
    #We iterate over different number of signals
    for (nSignal in singalSeq) {
      signalCount = signalCount +1
      
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
      
      # We print the scores for the most important features of each component
      # aux = cbind(selectVar(plsda, comp = 1)$value,selectVar(plsda, comp = 2)$value)
      # for (c in seq(3:(ncomp))) {
      #  aux = cbind(aux,selectVar(plsda, comp = c)$value)
      # }
      # print(head(aux))
      # 
      # png(file=file.path(paste(name,"PLSDA.png", sep = "")))
      # plotIndiv(plsda , comp = c(1,2),group = Y, ind.names = TRUE,ellipse = TRUE, legend = TRUE, title = paste('PLSDA comp 1-2',name))
      # dev.off()
       
        
      # #PCA
      pca <- pca(X, ncomp = ncomp)
      # print(head(selectVar(pca)['value']$value))
      #plotIndiv(pca ,group = Y, comp = c(1,2), ind.names = TRUE,ellipse = TRUE, legend = TRUE, title = 'PCA comp 1 - 2')
      
      #   
      # cat(paste("\n\nSIGNAL:",nSignal,"\n"),file= file.path(paste(nametxt,".txt", sep = "")),append=TRUE)
      # #We save in the files the most important features for each PLSDA component
      # cat("\nPLSDA",file= file.path(paste(nametxt,".txt", sep = "")),append=TRUE)
      # for (c in seq(1:(ncomp))) {
      #   cat(paste("\n","comp",c,": ", sep = ""),file= file.path(paste(nametxt,".txt", sep = "")),append=TRUE)
      #   cat(head(selectVar(plsda, comp = c)$name),file= file.path(paste(nametxt,".txt", sep = "")),append=TRUE)
      # }
      # 
      # #Then for each PCA component
      # cat("\n\nPCA",file= file.path(paste(nametxt,".txt", sep = "")),append=TRUE)
      # for (c in seq(1:(ncomp))) {
      #   cat(paste("\n","comp",c,": ", sep = ""),file= file.path(paste(nametxt,".txt", sep = "")),append=TRUE)
      #   cat(head(selectVar(pca, comp = c)$name),file= file.path(paste(nametxt,".txt", sep = "")),append=TRUE)
      # }

      
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

