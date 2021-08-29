#############clustering function ##############
clustx3_real2<-function(data,Nnorms,nOrd){
  result<-vector(length = 6)
  cluster<-data[["trueID"]]
  ppdim<- as.vector(apply(data[['catVars']],2,function(x) length(unique(x))))
  size<-nrow(dataset)
  OrdVar<-data[['catVars']][,1:nOrd]
  NomVar<-data[['catVars']][,-c(1:nOrd)]
  dataset<-cbind(scale(data[["conVars"]]), OrdVar,NomVar)
  nConVar<-ncol(data[['conVars']])
  nCatVar<-ncol(data[['catVars']])
  nOrdVar<-nOrd+nConVar
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







############alternative dataset
german_credit<-read.csv(file = '~/Desktop/Thesis/code/data/german_credit_data.csv')
View(german_credit)
german_credit$trueID = read.table(file = '~/Desktop/Thesis/code/data/german.data')[,21]
german_credit<-na.omit(german_credit)
german_credit <- german_credit[,-1]
trueID = german_credit$trueID

conVars = german_credit[,-c(10,2,3,4,5,6,9)]
catVars = german_credit[,c(2,3,4,5,6,9)]

catVars = apply(catVars,2,function(x) 
  as.numeric(factor(x)))

catVars_reordered = catVars[,c(1,2,4,5,3,6)]
View(catVars_reordered)
german_list<-list(trueID = trueID, conVars= conVars,catVars=catVars_reordered)
set.seed(1)
ger_result<-clustx3_real2(german_list,1e4,1)


par(mfrow=c(2,2))
plot(density(conVars[,1]),main='Age',cex.main=0.9)
plot(density(conVars[,2]),main='Credit_amount',cex.main=0.9)
plot(density(conVars[,3]),main='Duration',cex.main=0.9)



par(mfrow=c(1,1))
with(data = german_list,
     expr = pairs(scale(conVars), pch = c(1, 3)[trueID],
                  col= trueID, main = "German Credit Data \ncontinuous variable pair plot",
                  cex.main=0.9))

with(data = german_list,
     expr = pairs(jitter(catVars), pch = c(1, 3)[trueID],
                  col= trueID, main = "Australia Credit Data \ncategorical variable pair plot",
                  cex.main=0.9))



