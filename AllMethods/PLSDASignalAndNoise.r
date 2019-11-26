#########################################################################################################################
#########################################################################################################################
###### PAPER:          So you think you can PLS-DA?
###### NAME:           PLSDASignalAndNoiseCLUSTERClusterSizeSeparation.r
###### AUTHOR:         Daniel Ruiz-Perez, PhD Student. and Haibin Guan
###### AFFILIATION:    Florida International University
#########################################################################################################################
#########################################################################################################################



#rm(list=ls())

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
outputMatrixICA = array(0, dim=c(length(noiseSeq),length(singalSeq))) #First dimension number of features, second number of signal, 
outputMatrixRLDA = array(0, dim=c(length(noiseSeq),length(singalSeq))) #First dimension number of features, second number of signal, 
outputMatrixSPCA = array(0, dim=c(length(noiseSeq),length(singalSeq))) #First dimension number of features, second number of signal, 

#We execute with different matrixes to minimize random effects
for (repetition in repetitionsSeq) {
  print(as.numeric(repetition/nRepetitions)*100)
  
  #We execute with a random seed
  set.seed(sample(1:100, 1))
  set.seed(runif(1,0,10))
  
  outputMatrixPCAAux = array(0, dim=c(length(noiseSeq),length(singalSeq))) #First dimension number of features, second number of signal, 
  outputMatrixPLSDAAux = array(0, dim=c(length(noiseSeq),length(singalSeq))) #First dimension number of features, second number of signal, 
  outputMatrixICAAux = array(0, dim=c(length(noiseSeq),length(singalSeq))) #First dimension number of features, second number of signal, 
  outputMatrixRLDAAux = array(0, dim=c(length(noiseSeq),length(singalSeq))) #First dimension number of features, second number of signal, 
  outputMatrixSPCAAux = array(0, dim=c(length(noiseSeq),length(singalSeq))) #First dimension number of features, second number of signal,   
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
      
      plsda <- mixOmics::splsda(X, Y, ncomp = ncomp, scale = TRUE)
      
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
      pca <- mixOmics::pca(X, ncomp = ncomp)
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

      
      ica<-do.ica(X,ndim=1,type='logcosh',tpar=1,maxiter=100)
      
      
      rlda<-do.rlda(X,Y,ndim=1)
      
      spca<-mixOmics::spca(X,ncomp=ncomp)
      
      #We compute the performance of the models as the percentage of signal features in the important features
      importantVariablesPCA = mixOmics::selectVar(pca, comp = 1)$name[1:nSignal]
      importantVariablesPLSDA = mixOmics::selectVar(plsda, comp = 1)$name[1:nSignal]
      importantVariablesSPCA = mixOmics::selectVar(spca, comp = 1)$name[1:nSignal]
      importantVariablesICA=head(order(abs(ica$projection),decreasing=TRUE),n=nSignal)
      importantVariablesRLDA=head(order(abs(rlda$projection),decreasing=TRUE),n=nSignal)

      for (n in seq(1:(nSignal))) {
        if (paste('V',n,sep = "") %in% importantVariablesPCA){
          outputMatrixPCAAux[noiseCount,signalCount] = outputMatrixPCAAux[noiseCount,signalCount] + 1/nSignal
        }
        if (paste('V',n,sep = "") %in% importantVariablesPLSDA){
          outputMatrixPLSDAAux[noiseCount,signalCount] = outputMatrixPLSDAAux[noiseCount,signalCount] + 1/nSignal
        }
        if (paste(n,sep = "") %in% importantVariablesICA){
          outputMatrixICAAux[noiseCount,signalCount] = outputMatrixICAAux[noiseCount,signalCount] + 1/nSignal
        }
        if (paste(n,sep = "") %in% importantVariablesRLDA){
          outputMatrixRLDAAux[noiseCount,signalCount] = outputMatrixRLDAAux[noiseCount,signalCount] + 1/nSignal
        }
        if (paste('V',n,sep = "") %in% importantVariablesSPCA){
          outputMatrixSPCAAux[noiseCount,signalCount] = outputMatrixSPCAAux[noiseCount,signalCount] + 1/nSignal
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
save(list=ls(all=TRUE),file="SignalAndNoise_100.RData")
}

outputMatrixPLSDA = outputMatrixPLSDA/repetition
outputMatrixPCA = outputMatrixPCA/repetition
outputMatrixICA = outputMatrixICA/repetition
outputMatrixRLDA = outputMatrixRLDA/repetition
outputMatrixSPCA = outputMatrixSPCA/repetition


plsdaP <- plot_ly(colors=c("palegreen","palegreen4"),y = noiseSeq,x = singalSeq, z = outputMatrixPLSDA,showlegend=TRUE, name="PLSDA") %>% 
  layout(title = "Performance of PLSDA ",
         annotations = list(x = 8.71,y = 1.01,text = 'PLSDA',showarrow = FALSE),
         scene = list(
           xaxis = list(title = "# Signal"), 
           yaxis = list(title = "# Noise"), 
           zaxis = list(title = "Performance")
         )) %>% add_surface()
#plsdaP

pcaP <- plot_ly(colors=c("royalblue","royalblue4") ,showlegend=TRUE, name = "PCA",y = noiseSeq,x = singalSeq, z = outputMatrixPCA)%>%  
  layout(title = "Performance of PLSDA vs PCA",
         annotations = list(x = 6.99,y = 0.87,text = 'PCA',showarrow = FALSE),
         scene = list(
           xaxis = list(title = "# Signal"), 
           yaxis = list(title = "# Noise"), 
           zaxis = list(title = "Performance"))) %>% add_surface()
#pcaP


icaP <- plot_ly(colors=c("pink","pink4"),y = noiseSeq,x = singalSeq, z = outputMatrixICA,showlegend=TRUE, name="ICA") %>% 
  layout(title = "Performance of PLSDA vs PCA vs ICA ",
         annotations = list(x = 5.69,y = 0.73,text = 'ICA',showarrow = FALSE),
         scene = list(
           xaxis = list(title = "# Signal"), 
           yaxis = list(title = "# Noise"), 
           zaxis = list(title = "Performance")
         )) %>% add_surface()
#ICA


rldaP <- plot_ly(colors=c("grey50","dark grey"),y = noiseSeq,x = singalSeq,z = outputMatrixRLDA,showlegend=TRUE, name="RLDA") %>% 
  layout(title = "Performance of PLSDA vs PCA vs ICA vs RLDA ",
         annotations = list(x = 4.368,y = 0.57,text = 'RLDA',showarrow = FALSE),
         scene = list(
           xaxis = list(title = "# Signal"), 
           yaxis = list(title = "# Noise"), 
           zaxis = list(title = "Performance")
         )) %>% add_surface()
#RLDA

spcaP <- plot_ly(colors=c("slategray1","slategray4"),y = noiseSeq,x = singalSeq,z = outputMatrixSPCA,showlegend=TRUE, name="SPCA") %>% 
  layout(title = "Performance of PLSDA vs PCA vs ICA vs RLDA vs SPCA",
         annotations = list(x = 3.1,y = 0.425,text = 'SPCA',showarrow = FALSE),
         scene = list(
           xaxis = list(title = "# Signal"), 
           yaxis = list(title = "# Noise"), 
           zaxis = list(title = "Performance")
         )) %>% add_surface()
#SPCA

p<-subplot(plsdaP,pcaP,icaP,rldaP,spcaP)
p


htmlwidgets::saveWidget(as_widget(p), "SignalAndNoise_100.html")

save(list=ls(all=TRUE),file="SignalAndNoise_100.RData")


