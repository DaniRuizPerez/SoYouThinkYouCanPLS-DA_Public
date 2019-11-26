#########################################################################################################################
#########################################################################################################################
###### PAPER:          So you think you can PLS-DA?
###### NAME:           PLSDABVDataset.r.r
###### AUTHOR:         Daniel Ruiz-Perez, PhD Student
###### AFFILIATION:    Florida International University
#########################################################################################################################
#########################################################################################################################



rm(list=ls())

library(mixOmics)
library(plotly)


abundances=read.csv("abund.csv",header=T,sep=",")
lab=read.csv("specificity.csv",header=T)



Y = abundances[,2]
abundances = abundances[,-2]
rownames(abundances) = (abundances[,1])
abundances = abundances[,-1]
X = abundances


for (i in 1:5){
  lab[,2*i] = lab[,2*i]/(sum(lab[,2*i]))
}


AUWC = c()
#We do the one vs all approach. For every cst, we try to differenciate from all the others
for (i in 1:5){
    newY = as.numeric(Y) == i
    
    ######### calculate cumulative specificity for plsda
    plsda <- splsda(X, newY, ncomp = 2, scale = TRUE)
    #plotIndiv(plsda ,group = YIteration, comp = c(1,2), ind.names = TRUE,ellipse = TRUE, legend = TRUE, title = 'PLSDA comp 1 - 2')
    loadingsPLSDA <- (cbind(selectVar(plsda, comp = 1)$value,selectVar(plsda, comp = 2)$value))
    acumPLSDA = c(0)
    for (j in 2:(dim(abundances)[2]+1)){
      variableFound = grep(paste("\\b",rownames(loadingsPLSDA)[j-1],"\\b",sep=""),as.character(lab[,2*i-1]),value = F,ignore.case = T)
      valueOfJFeature = 0
      if (length(variableFound) == 1)
        valueOfJFeature = lab[,2*i][variableFound]
      
      acumPLSDA[j] = acumPLSDA[j-1]+valueOfJFeature
    }
    acumPLSDA = acumPLSDA[-1]    

    rownames(loadingsPLSDA)[1:10]
    
    ######### calculate cumulative specificity for pca
    pca <- pca(X, ncomp = 2)
    #plotIndiv(pca ,group = newY, comp = c(1,2), ind.names = TRUE,ellipse = TRUE, legend = TRUE, title = 'PCA comp 1 - 2')
    loadingsPCA <- cbind(selectVar(pca, comp = 1)$value,selectVar(pca, comp = 2)$value)
    acumPCA = c(0)
    for (j in 2:(dim(abundances)[2]+1)){
      variableFound = grep(paste("\\b",rownames(loadingsPCA)[j-1],"\\b",sep=""),as.character(lab[,2*i-1]),value = F,ignore.case = T, )
      valueOfJFeature = 0
      if (length(variableFound) == 1){
        valueOfJFeature = lab[,2*i][variableFound]
      }
      acumPCA[j] = acumPCA[j-1]+valueOfJFeature
    }
    acumPCA = acumPCA[-1]    
    
    

    ######### We plot the results, cumulatively    
    if (i == 1)
    plot(1:ncol(X)/length(acumPCA)*100,acumPCA,type = "l", col = i, xlab = "% features selected",
         ylab = "Specificity",main = "",lty=2,lwd=2,xlim=c(0, 10))
    else
      lines(1:ncol(X)/length(acumPCA)*100,acumPCA,type = "l", col = i,lty=2,lwd=2)
    
    # lines(1:ncol(X)/length(acumPCA)*100,acumPLSDA,type = "l", col = i,lwd=2)

    
    
    AUWC = cbind(AUWC,c(sum(acumPCA)/length(acumPCA)))
    
    
}

# legend(90, 0.20, legend=c("PCA", "PLS-DA"),col=c(1,1), lty=2:1, cex=0.8,lwd=2)
# legend(92.5, 0.55, legend=c("CST-I", "CST-II", "CST-III","CST-IV","CST-V"),col=1:5, cex=0.8,pch = 19)

# legend(0.4, 0.20, legend=c("PCA", "PLS-DA"),col=c(1,1), lty=2:1, cex=0.8,lwd=2)
# legend(0.43, 0.55, legend=c("CST-I", "CST-II", "CST-III","CST-IV","CST-V"),col=1:5, cex=0.8,pch = 19)

# legend(40, 0.20, legend=c("PCA", "PLS-DA"),col=c(1,1), lty=2:1, cex=0.8,lwd=2)
# legend(43, 0.55, legend=c("CST-I", "CST-II", "CST-III","CST-IV","CST-V"),col=1:5, cex=0.8,pch = 19)

# rownames(AUWC) = c("PCA","PLSDA")
AUWC

