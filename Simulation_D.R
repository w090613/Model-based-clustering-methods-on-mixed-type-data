library(MixSim)
library(arules)

##############data generation function#########
genMultiClustData<-function(size,nLevels,MaxOverlap, MeanOverlap, nClust, nConVar,nCatVar){
  Q <- MixSim(MaxOmega = MaxOverlap, BarOmega = MeanOverlap, K = nClust, p = nConVar+nCatVar)
  A <- simdataset(n = size, Pi = Q$Pi, Mu = Q$Mu, S = Q$S)
  cat_mat<-A$X
  catVars = matrix(nrow=size,ncol = nCatVar)
  for (i in c(1:nCatVar)){
    catVars[,i] = discretize(cat_mat[,nConVar+i],breaks = nLevels,labels=FALSE)
  }
  conVars = A$X[,c(1:nConVar)]
  trueID= A$id
  simdata<- list(trueID = trueID, conVars= conVars,catVars=catVars)
  return(simdata)
}
#############plot example simulated data##########
set.seed(2)
simdata1<-genMultiClustData(300,4,0.125,0.05,3,2,2)
with(data = simdata1,
     expr = pairs(cbind(conVars, jitter(catVars)), pch = c(1, 2,3)[trueID],
                  col = c(1,2,4)[trueID], main =  "Pair plot of pairwise overlaping 0.05 and maximum overlapping 0.125 \nfor both types of variables",cex.main=0.8))
set.seed(1)
simdata1<-genMultiClustData(300,4,0.4,0.3,3,2,2)
with(data = simdata1,
     expr = pairs(cbind(conVars, jitter(catVars)), pch = c(1, 3,4)[trueID],
                  col = c(1,2,4)[trueID],main =  "Pair plot of pairwise overlaping 0.3 and maximum overlapping 0.4 \nfor both types of variables",cex.main=0.8))

###########clustering fucntion#############
clustx3_MultiClust<-function(sim,Nnorms, nClust){
  result<-vector(length = 6)
  cluster<-sim[["trueID"]]
  data<-cbind(scale(sim[["conVars"]]), sim[["catVars"]])
  size<-nrow(data)
  ######finite mixture
  fmm<-flexmixedruns(x = data, xvarsorted = TRUE,continuous = 2,
                     discrete = 2,ppdim = c(4,4),n.cluster = nClust,allout = FALSE)
  fmm_clust<-fmm$flexout@cluster
  result[1]<-adjustedRandIndex(cluster,fmm_clust)
  
  agree=0
  for (i in (1:nClust)){
    agree = agree+ max(sapply(1:nClust,function(x) sum(fmm_clust[cluster==i]==x)))
  }
  result[2]<-(size-agree)/size
  
  ######clustMD
  model_types<-c("EII", "VII", "EEI", "VEI", "EVI", "VVI")
  model_list<-lapply(model_types,function(x) tryCatch(clustMD(X=data, G=nClust, CnsIndx=2, OrdIndx = 2,Nnorms = Nnorms,
                                                              MaxIter = 1e3, model = x, store.params = FALSE,
                                                              scale = FALSE, startCL = "kmeans", autoStop = TRUE, ma.band = 50,
                                                              stop.tol = 0.0001),error=function(e) NULL))
  if (is.null(unlist(model_list))){
    result[3]=NA
    result[4]=NA
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
    result[3]<-adjustedRandIndex(cluster,clustmd_clust)
    agree=0
    for (i in (1:nClust)){
      agree = agree+ max(sapply(1:nClust,function(x) sum(clustmd_clust[cluster==i]==x)))
    }
    result[4]<-(size-agree)/size
    
  }
  #####KAMILA
  catDf <- data.frame(apply(sim$catVars, 2, factor), stringsAsFactors = TRUE)
  conDf <- data.frame(scale(sim$conVars), stringsAsFactors = TRUE)
  kam<-kamila(conDf, catDf, numClust = nClust, numInit = 10)
  kam_clust<-kam$finalMemb
  result[5]<-adjustedRandIndex(cluster,kam_clust)
  agree=0
  for (i in (1:nClust)){
    agree = agree+ max(sapply(1:nClust,function(x) sum(kam_clust[cluster==i]==x)))
  }
  result[6]<-(size-agree)/size
  return(result)
  
}


############Simulation study#################
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


clustx3_3<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_MultiClust(datasets[[x]],1e4,3),error=function(e) NULL))
})

tf_null = sapply(1:n,function(x) is.null(clustx3_3[[x]]))
id_more = c(1:n)[tf_null]
clustx3_3more<-lapply(id_more,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_MultiClust(datasets[[x]],1e5,3),error=function(e) NULL))
})


####################################################################
#######################4 variables case############################
####################################################################
set.seed(1)
props = c(2.5,2,1.6,1.4,1.2,1.2)
datasets<-lapply(1:n,function(x)
  genMultiClustData(1e3,4,overlap_prop[x]*props[x],overlap_prop[x],4,2,2)
)
clustx3_4<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_MultiClust(datasets[[x]],1e4,4),error=function(e) NULL))
})

tf_null = sapply(1:n,function(x) is.null(clustx3_4[[x]]))
id_more = c(1:n)[tf_null]
clustx3_4more<-lapply(id_more,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_MultiClust(datasets[[x]],1e5,4),error=function(e) NULL))
})


####################################################################
#######################5 variables case############################
####################################################################

set.seed(1)
props = c(2.5,2,1.6,1.4,1.2,1.2)
datasets<-lapply(1:n,function(x)
  genMultiClustData(1e3,4,overlap_prop[x]*props[x],overlap_prop[x],5,2,2)
)


clustx3_5<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_MultiClust(datasets[[x]],1e4,5),error=function(e) NULL))
})

tf_null = sapply(1:n,function(x) is.null(clustx3_5[[x]]))
id_more = c(1:n)[tf_null]
clustx3_5more<-lapply(id_more,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_MultiClust(datasets[[x]],1e5,5),error=function(e) NULL))
})



####################################################################
#######################5 variables case############################
####################################################################

set.seed(1)
props = c(2.5,2,1.6,1.4,1.2,1.2)
datasets<-lapply(1:n,function(x)
  genMultiClustData(1e3,4,overlap_prop[x]*props[x],overlap_prop[x],6,2,2)
)


clustx3_6<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_MultiClust(datasets[[x]],1e4,6),error=function(e) NULL))
})

tf_null = sapply(1:n,function(x) is.null(clustx3_6[[x]]))
id_more = c(1:n)[tf_null]
clustx3_6more<-lapply(id_more,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_MultiClust(datasets[[x]],1e5,6),error=function(e) NULL))
})




