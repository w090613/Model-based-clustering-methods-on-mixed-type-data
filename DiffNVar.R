##############data generation function#########
genDiffNVarData <- function(size,nConVar,nCatVar,ErrLev) {
  tmpDat <- genMixedData(sampSize = size, nConVar = nConVar, nCatVar = nCatVar,
                         nCatLevels = 4, nConWithErr =nConVar, nCatWithErr = nCatVar,
                         popProportions = c(.5, .5), conErrLev = ErrLev,
                         catErrLev = ErrLev)
  tmpDat <- within(tmpDat, conVars <- as.data.frame(scale(conVars)))
  tmpDat <- within(tmpDat, catVarsFac <- as.data.frame(lapply(
    as.data.frame(catVars), factor)))
  within(tmpDat, catVarsDum <- as.data.frame(dummyCodeFactorDf(
    catVarsFac)))
}
###########clustering fucntion#############
clustx3_diffVar<-function(sim,Nnorms, nConVar,nCatVar){
  result<-vector(length = 6)
  cluster<-sim[["trueID"]]
  data<-cbind(scale(sim[["conVars"]]), sim[["catVars"]])
  size<-nrow(data)
  ######finite mixture
  fmm<-flexmixedruns(x = data, xvarsorted = TRUE,continuous = nConVar,
                     discrete = nCatVar,ppdim = c(4,4),n.cluster = 2,allout = FALSE)
  fmm_clust<-fmm$flexout@cluster
  result[1]<-adjustedRandIndex(cluster,fmm_clust)
  result[2]<-ifelse(sum(cluster!=fmm_clust)<size-sum(cluster!=fmm_clust),
                    sum(cluster!=fmm_clust)/size,1-sum(cluster!=fmm_clust)/size)
  
  ######clustMD
  model_types<-c("EII", "VII", "EEI", "VEI", "EVI", "VVI")
  model_list<-lapply(model_types,function(x) tryCatch(clustMD(X=data, G=2, CnsIndx=nConVar, OrdIndx = nConVar,Nnorms = Nnorms,
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
    result[4]<-ifelse(sum(cluster!=clustmd_clust)<size-sum(cluster!=clustmd_clust),
                      sum(cluster!=clustmd_clust)/size,1-sum(cluster!=clustmd_clust)/size)
  }
  
  #####KAMILA
  catDf <- data.frame(apply(sim$catVars, 2, factor), stringsAsFactors = TRUE)
  conDf <- data.frame(scale(sim$conVars), stringsAsFactors = TRUE)
  kam<-kamila(conDf, catDf, numClust = 2, numInit = 10)
  kam_clust<-kam$finalMemb
  result[5]<-adjustedRandIndex(cluster,kam_clust)
  result[6]<-ifelse(sum(cluster!=kam_clust)<size-sum(cluster!=kam_clust),
                    sum(cluster!=kam_clust)/size,1-sum(cluster!=kam_clust)/size)
  return(result)
  
}
############Simulation study#################

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
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e4,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})

tf_null = sapply(1:n,function(x) is.null(clustx3_005[[x]]))
id_more = c(1:n)[tf_null]
clustx3_005more<-lapply(id_more,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e5,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})

########0.2
set.seed(1)
datasets<-lapply(1:n,function(x)
  genDiffNVarData(size=1e3,ErrLev = overlap_prop[2], 
                  nConVar = nvar[x,1],nCatVar = nvar[x,2])
)
clustx3_02<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e4,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})

tf_null = sapply(1:n,function(x) is.null(clustx3_02[[x]]))
id_more = c(1:n)[tf_null]
clustx3_02more<-lapply(id_more,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e5,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})
######0.35
set.seed(1)
datasets<-lapply(1:n,function(x)
  genDiffNVarData(size=1e3,ErrLev = overlap_prop[3], 
                  nConVar = nvar[x,1],nCatVar = nvar[x,2])
)
clustx3_035<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e4,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})

tf_null = sapply(1:n,function(x) is.null(clustx3_035[[x]]))
id_more = c(1:n)[tf_null]
clustx3_035more<-lapply(id_more,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e5,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})

######0.5
set.seed(1)
datasets<-lapply(1:n,function(x)
  genDiffNVarData(size=1e3,ErrLev = overlap_prop[4], 
                  nConVar = nvar[x,1],nCatVar = nvar[x,2])
)
clustx3_05<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e4,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})

tf_null = sapply(1:n,function(x) is.null(clustx3_05[[x]]))
id_more = c(1:n)[tf_null]
clustx3_05more<-lapply(id_more[1],function(x) {
  set.seed(1)
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e5,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})






####################################################################
#######################5 variables case############################
####################################################################
set.seed(1)
overlap_prop<-c(0.05,0.2,0.35,0.5)
overlap_levels<-expand.grid(overlap_prop,overlap_prop)
n<-nrow(overlap_levels)
########0.05
set.seed(1)
datasets<-lapply(1:n,function(x)
  genDiffNVarData(size=1e3,ErrLev = overlap_prop[1], 
                  nConVar = nvar[x,1],nCatVar = nvar[x,2])
)
clustx3_005<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e4,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})

tf_null = sapply(1:n,function(x) is.null(clustx3_005[[x]]))
id_more = c(1:n)[tf_null]
clustx3_005more<-lapply(id_more,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e5,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})

########0.2
set.seed(1)
datasets<-lapply(1:n,function(x)
  genDiffNVarData(size=1e3,ErrLev = overlap_prop[2], 
                  nConVar = nvar[x,1],nCatVar = nvar[x,2])
)
clustx3_02<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e4,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})

tf_null = sapply(1:n,function(x) is.null(clustx3_02[[x]]))
id_more = c(1:n)[tf_null]
clustx3_02more<-lapply(id_more,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e5,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})
######0.35
set.seed(1)
datasets<-lapply(1:n,function(x)
  genDiffNVarData(size=1e3,ErrLev = overlap_prop[3], 
                  nConVar = nvar[x,1],nCatVar = nvar[x,2])
)
clustx3_035<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e4,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})

tf_null = sapply(1:n,function(x) is.null(clustx3_035[[x]]))
id_more = c(1:n)[tf_null]
clustx3_035more<-lapply(id_more,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e5,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})

######0.5
set.seed(1)
datasets<-lapply(1:n,function(x)
  genDiffNVarData(size=1e3,ErrLev = overlap_prop[4], 
                  nConVar = nvar[x,1],nCatVar = nvar[x,2])
)
clustx3_05<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e4,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})

tf_null = sapply(1:n,function(x) is.null(clustx3_05[[x]]))
id_more = c(1:n)[tf_null]
clustx3_05more<-lapply(id_more[2],function(x) {
  set.seed(1)
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e5,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})

####################################################################
#######################15 variables case############################
####################################################################
set.seed(1)
overlap_prop<-c(0.05,0.2,0.35,0.5)
nvar<-matrix(c(3,12,12,3,5,10,10,5,7,8,8,7),byrow = T,ncol = 2)
n<-nrow(nvar)
########0.05
set.seed(1)
datasets<-lapply(1:n,function(x)
  genDiffNVarData(size=1e3,ErrLev = overlap_prop[1], 
                  nConVar = nvar[x,1],nCatVar = nvar[x,2])
)
clustx3_005<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e4,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})



########0.2
set.seed(1)
datasets<-lapply(1:n,function(x)
  genDiffNVarData(size=1e3,ErrLev = overlap_prop[2], 
                  nConVar = nvar[x,1],nCatVar = nvar[x,2])
)
clustx3_02<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e4,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})


######0.35
set.seed(1)
datasets<-lapply(1:n,function(x)
  genDiffNVarData(size=1e3,ErrLev = overlap_prop[3], 
                  nConVar = nvar[x,1],nCatVar = nvar[x,2])
)
clustx3_035<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e4,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})



######0.5
set.seed(1)
datasets<-lapply(1:n,function(x)
  genDiffNVarData(size=1e3,ErrLev = overlap_prop[4], 
                  nConVar = nvar[x,1],nCatVar = nvar[x,2])
)
clustx3_05<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3_diffVar(datasets[[x]],1e4,nvar[x,1],nvar[x,2]),error=function(e) NULL))
})





