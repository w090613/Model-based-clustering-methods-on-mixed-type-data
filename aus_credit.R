#############clustering function ##############

clustx3_real1<-function(data,Nnorms){
  result<-vector(length = 6)
  cluster<-data[["trueID"]]
  ppdim<- as.vector(apply(data[['catVars']],2,function(x) length(unique(x))))
  size<-length(cluster)
  OrdVar <- data[['catVars']][,ppdim==2]
  NomVar<- data[['catVars']][,ppdim!=2]
  dataset<-cbind(scale(data[["conVars"]]), OrdVar,NomVar)
  nConVar<-ncol(data[['conVars']])
  nCatVar<-ncol(data[['catVars']])
  nOrdVar<-ncol(OrdVar)+nConVar
  ######finite mixture
  fmm<-flexmixedruns(x = dataset, xvarsorted = TRUE,continuous = nConVar,
                     discrete = nCatVar,ppdim = ppdim,n.cluster = 2,allout = FALSE,recode = TRUE)
  fmm_clust<-fmm$flexout@cluster
  result[1]<-adjustedRandIndex(cluster,fmm_clust)
  result[2]<-ifelse(sum(cluster!=fmm_clust)<size-sum(cluster!=fmm_clust),
                    sum(cluster!=fmm_clust)/size,1-sum(cluster!=fmm_clust)/size)
  
  ######clustMD
  model_types<-c("EII", "VII", "EEI", "VEI", "EVI", "VVI")
  model_list<-lapply(model_types,function(x) tryCatch(clustMD(X=data, G=2, CnsIndx=nConVar, OrdIndx = nOrdVar,Nnorms = Nnorms,
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
  catDf <- data.frame(apply(data$catVars, 2, factor), stringsAsFactors = TRUE)
  conDf <- data.frame(scale(data$conVars), stringsAsFactors = TRUE)
  kam<-kamila(conDf, catDf, numClust = 2, numInit = 10)
  kam_clust<-kam$finalMemb
  result[5]<-adjustedRandIndex(cluster,kam_clust)
  result[6]<-ifelse(sum(cluster!=kam_clust)<size-sum(cluster!=kam_clust),
                    sum(cluster!=kam_clust)/size,1-sum(cluster!=kam_clust)/size)
  return(result)
  
}
aus_credit<-read.table(file = 'australian.dat')
View(aus_credit)
aus_credit<-na.omit(aus_credit)
trueID = aus_credit[,15]+1
sum(trueID==1)
sum(trueID==2)
conVars = aus_credit[,-c(1,4,5,6,8,9,11,12,15)]
catVars = aus_credit[,c(1,4,5,6,8,9,11,12)]
catVars = apply(catVars,2,function(x) {
  if (min(x)==0) x+1
  else x})
aus_list<-list(trueID = trueID, conVars= conVars,catVars=catVars)
set.seed(1)
aus_result<-clustx3_real1(aus_list,1e5)

catVars_less <- catVars[,-c(3:4)]
aus_list_alt<-list(trueID = trueID, conVars= conVars,catVars=catVars_less)
set.seed(1)
aus_result_alt<-clustx3_real1(aus_list_alt,1e4)


with(data = aus_list,
     expr = pairs(scale(conVars), pch = c(1, 3)[trueID],
                  col= trueID, main = "Australia Credit Data \ncontinuous variable pair plot",
                  cex.main=0.9))
par(mfrow=c(2,3))
plot(density(conVars[,1]),main='Density of variable 2',cex.main=0.8)
plot(density(conVars[,2]),main='Density of variable 3',cex.main=0.8)
plot(density(conVars[,3]),main='Density of variable 7',cex.main=0.8)
plot(density(conVars[,4]),main='Density of variable 10',cex.main=0.8)
plot(density(conVars[,5]),main='Density of variable 13',cex.main=0.8)
plot(density(conVars[,6]),main='Density of variable 14',cex.main=0.8)



par(mfrow=c(1,1))
with(data = aus_list,
     expr = pairs( jitter(catVars), pch = c(1, 3)[trueID],
                  col= trueID, main = "Australia Credit Data \ncategorical variable pair plot",
                  cex.main=0.9))
