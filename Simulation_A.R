##############data generation function#########
genSimpleData <- function(conErrLev, catErrLev) {
  tmpDat <- genMixedData(sampSize = 200, nConVar = 2, nCatVar = 2,
                         nCatLevels = 4, nConWithErr = 2, nCatWithErr = 2,
                         popProportions = c(.5, .5), conErrLev = conErrLev,
                         catErrLev = catErrLev)
  tmpDat <- within(tmpDat, conVars <- as.data.frame(scale(conVars)))
  tmpDat <- within(tmpDat, catVarsFac <- as.data.frame(lapply(
    as.data.frame(catVars), factor)))
  within(tmpDat, catVarsDum <- as.data.frame(dummyCodeFactorDf(
    catVarsFac)))
}


#############plot example simulated data##########

set.seed(1)
with(data = genSimpleData(conErrLev = 0.3, catErrLev = 0.3),
     expr = pairs(cbind(conVars, jitter(catVars)), pch = c(1, 3)[trueID],
                  col = trueID, main = "Pair plot of overlaping 0.3 \nfor both types of variables",cex.main=0.9))
set.seed(1)
with(data = genSimpleData(conErrLev = 0.1, catErrLev = 0.5),
     expr = pairs(cbind(conVars, jitter(catVars)), pch = c(1, 3)[trueID],
                  col = trueID, main = "Pair plot of overlaping 0.1 for continuous variables \nand 0.5 for categorical variables",cex.main=0.9))




###########clustering fucntion#############
clustx3<-function(sim,Nnorms){
  result<-vector(length = 6)
  cluster<-sim[["trueID"]]
  data<-cbind(scale(sim[["conVars"]]), sim[["catVars"]])
  size<-nrow(data)
  ######finite mixture
  fmm<-flexmixedruns(x = data, xvarsorted = TRUE,continuous = 2,
                     discrete = 2,ppdim = c(4,4),n.cluster = 2,allout = FALSE)
  fmm_clust<-fmm$flexout@cluster
  result[1]<-adjustedRandIndex(cluster,fmm_clust)
  result[2]<-ifelse(sum(cluster!=fmm_clust)<size-sum(cluster!=fmm_clust),
         sum(cluster!=fmm_clust)/size,1-sum(cluster!=fmm_clust)/size)
  
  ######clustMD
  model_types<-c("EII", "VII", "EEI", "VEI", "EVI", "VVI")
  model_list<-lapply(model_types,function(x) tryCatch(clustMD(X=data, G=2, CnsIndx=2, OrdIndx = 2,Nnorms = Nnorms,
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

############Simulatin study#################
overlap_prop<-c(0.05,0.1,0.2,0.3,0.4,0.5,0.6)
n<-length(overlap_prop)
set.seed(1)
datasets<-lapply(1:n,function(x)
  genSimpleData(conErrLev = overlap_prop[x], catErrLev = overlap_prop[x])
)
clustx3_same<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3(datasets[[x]],1e4),error=function(e) NULL))
})
tf_null = sapply(1:n,function(x) is.null(clustx3_same[[x]]))
id_more = c(1:n)[tf_null]
clustx3_same_more<-lapply(id_more,function(x) {
  set.seed(1)
  return(tryCatch(clustx3(datasets[[x]],1e5),error=function(e) NULL))
})


overlap_prop<-c(0.05,0.1,0.3,0.5)
overlap_levels<-expand.grid(overlap_prop,overlap_prop)
n<-nrow(overlap_levels)
set.seed(1)
datasets<-lapply(1:n,function(x)
  genSimpleData(conErrLev = overlap_levels[x,1], catErrLev = overlap_levels[x,2])
)
clustx3_diff<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3(datasets[[x]],1e4),error=function(e) NULL))
})
tf_null = sapply(1:n,function(x) is.null(clustx3_diff[[x]]))
id_more = c(1:n)[tf_null]
clustx3_diff_more<-lapply(id_more,function(x) {
  set.seed(1)
  return(tryCatch(clustx3(datasets[[x]],1e5),error=function(e) NULL))
})


############plot default overlap#####################
mu2 = -2 * qnorm(OVERLAP_DEFAULT/2)
a=seq(-5,10,length.out=1000)
b=seq(-5,10,length.out = 1000)
plot(a,dnorm(a),type='l',ylab='density',ylim=c(-0.03,0.4),
     main='Overlapping level 0.01 \nbetween two normal densities',xlab='')
axis(1, at=mu2/2,labels='x',cex=0.5)
lines(b,dnorm(b,mean = mu2),col='red')
legend('topright',c('density 1','density 2'),col=c('black','red'),lty=c(1,1),cex=0.5)




############The ARI is imported form the result#####################

overlap_prop<-c(0.05,0.1,0.2,0.3,0.4,0.5,0.6)
fmm_ari = c(1,0.996,0.953,0.843,0.736,0.53,0.371)
clustmd_ari = c(1,1,0.956,0.828,0.719,0.541,0.333)
kamila_ari =c(1,0.996,0.956,0.854,0.729,0.538,0.343)
plot(overlap_prop,fmm_ari,pch=19,type ='b',ylab='ARI',xlab='Overlapping level',
     main ='Plot of ARI against Overlapping Level \nfor clustering methods',ylim=c(0.3,1))
points(overlap_prop,clustmd_ari,pch=19,col='red',type ='b')
points(overlap_prop,kamila_ari,pch=19,col='blue',type ='b')
legend('topright',c('Latent Class Model','ClustMD',"KAMILA"),pch=c(19,19,19),cex=0.8,col=c("black",'red','blue'))
