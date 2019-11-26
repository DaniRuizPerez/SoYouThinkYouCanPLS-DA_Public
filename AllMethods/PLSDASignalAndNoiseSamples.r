#########################################################################################################################
#########################################################################################################################
###### PAPER:          So you think you can PLS-DA?
###### NAME:           PLSDASignalAndNoiseCLUSTERClusterSizeSeparation.r
###### AUTHOR:         Daniel Ruiz-Perez, PhD Student. and Haibin Guan
###### AFFILIATION:    Florida International University
#########################################################################################################################
#########################################################################################################################



#Original one, plots features vs signal vs performance


rm(list=ls())
packages_required <- c("mixOmics", "plotly", "clusterGeneration","Rdimtools","doParallel","stringi","caret")
packages_required <- packages_required[!(packages_required %in% installed.packages()[,"Package"])]
if (length(packages_required)) {
  install.packages(packages_required,repos="http://cran.us.r-project.org")
}

getwd()
library(mixOmics)
library(plotly)
library(clusterGeneration)
library(Rdimtools)
library(doParallel)
library(caret)
cl<-makeCluster(6)
registerDoParallel(cl)

load(file="NoiseSamples_30.RData")
print(ls())
### load outputMatix from previous .RDATA
outputMatrixPLSDA = outputMatrixPLSDA*repetition
outputMatrixPCA = outputMatrixPCA*repetition
outputMatrixICA = outputMatrixICA*repetition
outputMatrixRLDA = outputMatrixRLDA*repetition
outputMatrixSPCA = outputMatrixSPCA*repetition



constantFirstHalf = 0.5 #value that the sum of the signal of the first half of samples add up to
constantSecondHalf = -0.5 #value that the sum of the signal of the second half of samples add up to
ncomp = 1 #Number of components for plsda
perturbation = 0 #The perturbtion to add to the signal
#Values to iterate over the features lenght and signal length
noiseSeq = seq(150, 1000, by = 20)
samplesSeq = seq(10, 400, by = 10)
nRepetitions =50
repetitionsSeq = seq(repetition+1, nRepetitions, by = 1)
nSignal= 10

#outputMatrixPCA = array(0, dim=c(length(noiseSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 
#outputMatrixPLSDA = array(0, dim=c(length(noiseSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 
#outputMatrixICA = array(0, dim=c(length(noiseSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 
#outputMatrixRLDA = array(0, dim=c(length(noiseSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 
#outputMatrixSPCA = array(0, dim=c(length(noiseSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 



#We execute with different matrixes to minimize random effects
for (repetition in repetitionsSeq) {
  print(as.numeric(repetition/nRepetitions)*100)
  
  #We execute with a random seed
  set.seed(sample(1:100, 1))
  set.seed(runif(1,0,10))
  
  outputMatrixPCAAux = array(0, dim=c(length(noiseSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 
  outputMatrixPLSDAAux = array(0, dim=c(length(noiseSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 
  outputMatrixICAAux = array(0, dim=c(length(noiseSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 
  outputMatrixRLDAAux = array(0, dim=c(length(noiseSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 
  outputMatrixSPCAAux = array(0, dim=c(length(noiseSeq),length(samplesSeq))) #First dimension number of features, second number of signal, 
  
  
  #values are %of good yfeatures in the first #signal values
  OriginalM = matrix(rnorm(samplesSeq[length(samplesSeq)]*(noiseSeq[length(noiseSeq)] +nSignal)),
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
      
      plsda <- mixOmics::splsda(X, Y, ncomp = ncomp, scale = TRUE)
   
        
      # #PCA
      pca <- mixOmics::pca(X, ncomp = ncomp)

      ica<-do.ica(X,ndim=ncomp,type='logcosh',tpar=1,maxiter=100)
      
      rlda<-do.rlda(X,Y,ndim=ncomp)
      
      spca<-mixOmics::spca(X,ncomp=ncomp)
      
 
      #We compute the performance of the models as the percentage of signal features in the important features
      importantVariablesPCA = selectVar(pca, comp = 1)$name[1:nSignal]
      importantVariablesPLSDA = selectVar(plsda, comp = 1)$name[1:nSignal]
      importantVariablesICA=head(order(abs(ica$projection),decreasing=TRUE),n=nSignal)
      importantVariablesRLDA=head(order(abs(rlda$projection),decreasing=TRUE),n=nSignal)
      importantVariablesSPCA = mixOmics::selectVar(spca, comp = 1)$name[1:nSignal]

      for (n in seq(1:(nSignal))) {
        if (paste('V',n,sep = "") %in% importantVariablesPCA){
          outputMatrixPCAAux[noiseCount,samplesCount] = outputMatrixPCAAux[noiseCount,samplesCount] + 1/nSignal
        }
        if (paste('V',n,sep = "") %in% importantVariablesPLSDA){
          outputMatrixPLSDAAux[noiseCount,samplesCount] = outputMatrixPLSDAAux[noiseCount,samplesCount] + 1/nSignal
        }
        if (paste(n,sep = "") %in% importantVariablesICA){
          outputMatrixICAAux[noiseCount,samplesCount] = outputMatrixICAAux[noiseCount,samplesCount] + 1/nSignal
        }
        if (paste(n,sep = "") %in% importantVariablesRLDA){
          outputMatrixRLDAAux[noiseCount,samplesCount] = outputMatrixRLDAAux[noiseCount,samplesCount] + 1/nSignal
        }
        if (paste('V',n,sep = "") %in% importantVariablesSPCA){
          outputMatrixSPCAAux[noiseCount,samplesCount] = outputMatrixSPCAAux[noiseCount,samplesCount] + 1/nSignal
        }
        
      }
    }
    # setwd("..")
  }
  outputMatrixPLSDA = outputMatrixPLSDA + outputMatrixPLSDAAux
  outputMatrixPCA = outputMatrixPCA + outputMatrixPCAAux
  outputMatrixICA = outputMatrixICA + outputMatrixICAAux
  outputMatrixRLDA = outputMatrixRLDA + outputMatrixRLDAAux
  outputMatrixSPCA = outputMatrixSPCA + outputMatrixSPCAAux
save(list=ls(all=TRUE),file="NoiseSamples_50.RData")

}

outputMatrixPLSDA = outputMatrixPLSDA/repetition
outputMatrixPCA = outputMatrixPCA/repetition
outputMatrixICA = outputMatrixICA/repetition
outputMatrixRLDA = outputMatrixRLDA/repetition
outputMatrixSPCA = outputMatrixSPCA/repetition

plsdaP <- plot_ly(y = noiseSeq,x = samplesSeq, z = outputMatrixPLSDA,showlegend=TRUE, name="PLSDA") %>% 
  layout(title = "Performance of PLSDA ",
         annotations = list(x = 6.,y = 1.01,text = 'PLSDA',showarrow = FALSE),
         scene = list(
           xaxis = list(title = "# Samples"), 
           yaxis = list(title = "# Noise"), 
           zaxis = list(title = "Performance")
         )) %>% add_surface()
#plsdaP

pcaP <- plot_ly(colors=c("dark blue","orange"),showlegend=TRUE, name = "PCA",y = noiseSeq,x = samplesSeq, z = outputMatrixPCA)%>%  
  layout(title = "Performance of PLSDA vs PCA",
         annotations = list(x = 5.3,y = 0.84,text = 'PCA',showarrow = FALSE),
         scene = list(
           xaxis = list(title = "# Samples"), 
           yaxis = list(title = "# Noise"), 
           zaxis = list(title = "Performance"))) %>% add_surface()
#pcaP


icaP <- plot_ly(colors=c("grey50","navajowhite3"),y = noiseSeq,x = samplesSeq, z = outputMatrixICA,showlegend=TRUE, name="ICA") %>% 
  layout(title = "Performance of PLSDA vs PCA vs ICA ",
         annotations = list(x = 4.1,y = 0.66,text = 'ICA',showarrow = FALSE),
         scene = list(
           xaxis = list(title = "# Samples"), 
           yaxis = list(title = "# Noise"), 
           zaxis = list(title = "Performance")
         )) %>% add_surface()
#ICA


rldaP <- plot_ly(colors=c("dark red","salmon"),y = noiseSeq,x = samplesSeq, z = outputMatrixRLDA,showlegend=TRUE, name="RLDA") %>% 
  layout(title = "Performance of PLSDA vs PCA vs ICA vs RLDA ",
         annotations = list(x = 2.8,y = 0.495,text = 'RLDA',showarrow = FALSE),
         scene = list(
           xaxis = list(title = "# Samples"), 
           yaxis = list(title = "# Noise"), 
           zaxis = list(title = "Performance")
         )) %>% add_surface()
#RLDA

spcaP <- plot_ly(colors=c("slategray1","slategray4"),y = noiseSeq,x = samplesSeq, z = outputMatrixSPCA,showlegend=TRUE, name="SPCA") %>% 
  layout(title = "Performance of PLSDA vs PCA vs ICA vs RLDA vs SPCA",
         annotations = list(x = 1.4,y = 0.32,text = 'SPCA',showarrow = FALSE),
         scene = list(
           xaxis = list(title = "# Samples"), 
           yaxis = list(title = "# Noise"), 
           zaxis = list(title = "Performance")
         )) %>% add_surface()
#SPCA


p<-subplot(plsdaP,pcaP,icaP,rldaP,spcaP)
p

htmlwidgets::saveWidget(as_widget(p), "NoiseSamples_50.html")


save(list=ls(all=TRUE),file="NoiseSamples_50.RData")

