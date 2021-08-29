##############data generation function#########
genSimpleData <- function(conErrLev, catErrLev) {
  tmpDat <- genMixedData(sampSize = 1000, nConVar = 2, nCatVar = 2,
                         nCatLevels = 4, nConWithErr = 2, nCatWithErr = 2,
                         popProportions = c(.5, .5), conErrLev = conErrLev,
                         catErrLev = catErrLev)
  tmpDat <- within(tmpDat, conVars <- as.data.frame(scale(conVars)))
  tmpDat <- within(tmpDat, catVarsFac <- as.data.frame(lapply(
    as.data.frame(catVars), factor)))
  within(tmpDat, catVarsDum <- as.data.frame(dummyCodeFactorDf(
    catVarsFac)))
}


###########clustering fucntion#############
clust_ord<-function(sim,Nnorms){
  result<-vector(length = 2)
  cluster<-sim[["trueID"]]
  data<-cbind(scale(sim[["conVars"]]), sim[["catVars"]])
  size<-nrow(data)
  model_types<-c("EII", "VII", "EEI", "VEI", "EVI", "VVI")
  model_list<-lapply(model_types,function(x) tryCatch(clustMD(X=data, G=2, CnsIndx=2, OrdIndx = 4,Nnorms = Nnorms,
                                                              MaxIter = 1e3, model = x, store.params = FALSE,
                                                              scale = FALSE, startCL = "kmeans", autoStop = TRUE, ma.band = 50,
                                                              stop.tol = 0.0001),error=function(e) NULL))
  if (is.null(unlist(model_list))){
    result[1]=NA
    result[2]=NA
  }
  else{
    bics <- vector()
    index <- vector()
    for (i in c(1:6)){
      if (!is.null(model_list[[i]])){
        bics<- c(bics,model_list[[i]][["BIChat"]])
        index<-c(index, i)
      }
    }
    clustmd<-model_list[[index[bics==min(bics)]]]
    clustmd_clust<-clustmd[["cl"]]
    result[1]<-adjustedRandIndex(cluster,clustmd_clust)
    result[2]<-ifelse(sum(cluster!=clustmd_clust)<size-sum(cluster!=clustmd_clust),
                      sum(cluster!=clustmd_clust)/size,1-sum(cluster!=clustmd_clust)/size)
  }
  return(result)
}

clust_ord_multiclust<-function(sim,Nnorms,nClust){
  result<-vector(length = 2)
  cluster<-sim[["trueID"]]
  data<-cbind(scale(sim[["conVars"]]), sim[["catVars"]])
  size<-nrow(data)
  model_types<-c("EII", "VII", "EEI", "VEI", "EVI", "VVI")
  model_list<-lapply(model_types,function(x) tryCatch(clustMD(X=data, G=nClust, CnsIndx=2, OrdIndx = 4,Nnorms = Nnorms,
                                                              MaxIter = 1e3, model = x, store.params = FALSE,
                                                              scale = FALSE, startCL = "kmeans", autoStop = TRUE, ma.band = 50,
                                                              stop.tol = 0.0001),error=function(e) NULL))
  if (is.null(unlist(bic_list))){
    result[1]=NA
    result[2]=NA
  }
  else{
    bics <- vector()
    index <- vector()
    for (i in c(1:6)){
      if (!is.null(model_list[[i]])){
        bics<- c(bics,model_list[[i]][["BIChat"]])
        index<-c(index, i)
      }
    }
    clustmd<-model_list[[index[bics==min(bics)]]]
    clustmd_clust<-clustmd[["cl"]]
    result[1]<-adjustedRandIndex(cluster,clustmd_clust)
    agree=0
    for (i in (1:nClust)){
      agree = agree+ max(sapply(1:nClust,function(x) sum(clustmd_clust[cluster==i]==x)))
    }
    result[2]<-(size-agree)/size
    
  }
  return(result)
}


############Simulation study#################
######normal same overlap
overlap_prop<-c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6)
n<-length(overlap_prop)
set.seed(1)
datasets<-lapply(1:n,function(x)
  genSimpleData(conErrLev = overlap_prop[x], catErrLev = overlap_prop[x])
)
clust_ord_same<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clust_ord(datasets[[x]],1e4),error=function(e) NULL))
})


######normal diff overlap
overlap_prop<-c(0.05,0.1,0.3,0.5)
overlap_levels<-expand.grid(overlap_prop,overlap_prop)
n<-nrow(overlap_levels)
set.seed(1)
datasets<-lapply(1:n,function(x)
  genSimpleData(conErrLev = overlap_levels[x,1], catErrLev = overlap_levels[x,2])
)
clust_ord_diff<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clust_ord(datasets[[x]],1e4),error=function(e) NULL))
})




######pgNorm diff overlap
overlap_prop<-c(0.05,0.2,0.35,0.5)
overlap_levels<-expand.grid(overlap_prop,overlap_prop)
n<-nrow(overlap_levels)
########kurtosis8
set.seed(1)
datasets<-lapply(1:n,function(x)
  genPGNormSimpleData(size=1e3,conErrLev = overlap_levels[x,1], 
                      catErrLev = overlap_levels[x,2],p=0.7019)
)
clustx3_8<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clust_ord(datasets[[x]],1e4),error=function(e) NULL))
})



########kurtosis6
set.seed(1)
datasets<-lapply(1:n,function(x)
  genPGNormSimpleData(size=1e3,conErrLev = overlap_levels[x,1], 
                      catErrLev = overlap_levels[x,2],p=0.7785)
)
clustx3_6<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clust_ord(datasets[[x]],1e4),error=function(e) NULL))
})



########kurtosis4
set.seed(1)
datasets<-lapply(1:n,function(x)
  genPGNormSimpleData(size=1e3,conErrLev = overlap_levels[x,1], 
                      catErrLev = overlap_levels[x,2],p=0.9012)
)
clustx3_4<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clust_ord(datasets[[x]],1e4),error=function(e) NULL))
})





######LogNorm diff overlap

########skewness1
set.seed(1)
datasets<-lapply(1:n,function(x)
  genLogNormSimpleData(size=1e3,conErrLev = overlap_levels[x,1], 
                       catErrLev = overlap_levels[x,2],sdlog = 0.3143)
)
clustx3_1<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clust_ord(datasets[[x]],1e4),error=function(e) NULL))
})

########skewness2.5
set.seed(1)
datasets<-lapply(1:n,function(x)
  genLogNormSimpleData(size=1e3,conErrLev = overlap_levels[x,1], 
                       catErrLev = overlap_levels[x,2],sdlog = 0.6409)
)
clustx3_2<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clust_ord(datasets[[x]],1e4),error=function(e) NULL))
})


########skewness9
set.seed(1)
datasets<-lapply(1:n,function(x)
  genLogNormSimpleData(size=1e3,conErrLev = overlap_levels[x,1], 
                       catErrLev = overlap_levels[x,2],sdlog = 1.1310)
)
clustx3_9<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clust_ord(datasets[[x]],1e4),error=function(e) NULL))
})


########diff number of variables
clust_ord_diffVar<-function(sim,Nnorms,nConVar,nOrdVar){
  result<-vector(length = 2)
  cluster<-sim[["trueID"]]
  data<-cbind(scale(sim[["conVars"]]), sim[["catVars"]])
  size<-nrow(data)
  model_types<-c("EII", "VII", "EEI", "VEI", "EVI", "VVI")
  bic_list<-lapply(model_types,function(x) tryCatch(clustMD(X=data, G=2, CnsIndx=nConVar, OrdIndx = nOrdVar+nConVar,Nnorms = Nnorms,
                                                            MaxIter = 1e3, model = x, store.params = FALSE,
                                                            scale = FALSE, startCL = "kmeans", autoStop = TRUE, ma.band = 50,
                                                            stop.tol = 0.0001)[["BIChat"]],error=function(e) NULL))
  num<-sapply(1:6,function(x) if(is.null(bic_list[[x]])){bic_list[[x]]=min(unlist(bic_list)+1)}else{bic_list[[x]]=bic_list[[x]]})
  model_ind<-model_types[which.min(num)]#####choose by BIC
  clustmd<-clustMD(X=data, G=2, CnsIndx=nConVar, OrdIndx = nOrdVar+nConVar,Nnorms = Nnorms,
                   MaxIter = 1e3, model = model_ind, store.params = FALSE,
                   scale = FALSE, startCL = "kmeans", autoStop = TRUE, ma.band = 50,
                   stop.tol = 0.0001)
  clustmd_clust<-clustmd[["cl"]]
  result[1]<-adjustedRandIndex(cluster,clustmd_clust)
  result[2]<-ifelse(sum(cluster!=clustmd_clust)<size-sum(cluster!=clustmd_clust),
                    sum(cluster!=clustmd_clust)/size,1-sum(cluster!=clustmd_clust)/size)
  return(result)
}
####################################################################
#######################10 variables case############################
####################################################################
set.seed(1)
overlap_prop<-c(0.05,0.2,0.35,0.5)
nvar<-matrix(c(1,9,9,1,2,8,8,2,3,7,7,3,4,6,6,4,5,5),byrow = T,ncol = 2)
n<-nrow(nvar)
########0.05
set.seed(1)
datasets<-lapply(1:n,function(x)
  genDiffNVarData(size=1e3,ErrLev = overlap_prop[1], 
                  nConVar = nvar[x,1],nCatVar = nvar[x,2])
)
clustx3_005<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clust_ord_diffVar(datasets[[x]],1e4,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})


########0.2
set.seed(1)
datasets<-lapply(1:n,function(x)
  genDiffNVarData(size=1e3,ErrLev = overlap_prop[2], 
                  nConVar = nvar[x,1],nCatVar = nvar[x,2])
)
clustx3_02<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clust_ord_diffVar(datasets[[x]],1e4,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})


######0.35
set.seed(1)
datasets<-lapply(1:n,function(x)
  genDiffNVarData(size=1e3,ErrLev = overlap_prop[3], 
                  nConVar = nvar[x,1],nCatVar = nvar[x,2])
)
clustx3_035<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clust_ord_diffVar(datasets[[x]],1e4,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})



######0.5
set.seed(1)
datasets<-lapply(1:n,function(x)
  genDiffNVarData(size=1e3,ErrLev = overlap_prop[4], 
                  nConVar = nvar[x,1],nCatVar = nvar[x,2])
)
clustx3_05<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clust_ord_diffVar(datasets[[x]],1e4,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})


#############multi clust##########
set.seed(1)
overlap_prop<-c(0.05,0.1,0.2,0.3,0.4,0.5)
props = c(2.5,2,1.6,1.4,1.2,1.2)
n<-length(overlap_prop)
####################################################################
##########################3 clusters###############################
####################################################################


set.seed(1)
datasets<-lapply(1:n,function(x)
  genMultiClustData(1e3,4,overlap_prop[x]*props[x],overlap_prop[x],3,2,2)
)

clust_ord_3<-lapply(1:n,function(x){
  set.seed(1)
  return(tryCatch(clust_ord_multiclust(datasets[[x]],1e4,3),error = function(e) NULL))
})
####################################################################
#######################4 variables case############################
####################################################################
set.seed(1)
props = c(2.5,2,1.6,1.4,1.2,1.2)
datasets<-lapply(1:n,function(x)
  genMultiClustData(1e3,4,overlap_prop[x]*props[x],overlap_prop[x],4,2,2)
)
clust_ord_4<-lapply(1:n,function(x){
  set.seed(1)
  return(tryCatch(clust_ord_multiclust(datasets[[x]],1e4,4),error = function(e) NULL))
})

####################################################################
#######################5 variables case############################
####################################################################

set.seed(1)
props = c(2.5,2,1.6,1.4,1.2,1.2)
datasets<-lapply(1:n,function(x)
  genMultiClustData(1e3,4,overlap_prop[x]*props[x],overlap_prop[x],5,2,2)
)

clust_ord_5<-lapply(1:n,function(x){
  set.seed(1)
  return(tryCatch(clust_ord_multiclust(datasets[[x]],1e4,5),error = function(e) NULL))
})


####################################################################
#######################5 variables case############################
####################################################################

set.seed(1)
props = c(2.5,2,1.6,1.4,1.2,1.2)
datasets<-lapply(1:n,function(x)
  genMultiClustData(1e3,4,overlap_prop[x]*props[x],overlap_prop[x],6,2,2)
)

clust_ord_6<-lapply(1:n,function(x){
  set.seed(1)
  return(tryCatch(clust_ord_multiclust(datasets[[x]],1e4,6),error = function(e) NULL))
})





######################plot################
library(readxl)
plot_data<-read_xlsx('ord_vs_nom.xlsx')
plot_data<-as.data.frame(plot_data)
plot(x=plot_data[,1],y=plot_data[,2],xlab = 'clustering ARI for ordinal data',
     ylab='clustering ARI for nominal data',main = 'Pairwise comparison for two types of categorical data',
     cex.main=0.9,cex.lab=0.9)
lines(x=seq(0,1,0.001),y=seq(0,1,0.001),col='red')
