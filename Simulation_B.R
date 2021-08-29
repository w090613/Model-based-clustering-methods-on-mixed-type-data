################################################
#######################p-gnorm #################
################################################

##############data generation function#########
genPGNormMixedData<-function (sampSize, nConVar, nCatVar, nCatLevels, nConWithErr, 
          nCatWithErr, popProportions, conErrLev, catErrLev,p) 
  
{
  require("pgnorm")
  if (length(popProportions) > 2) {
    stop("Error in genMixedData: More than two populations not currently supported")
  }
  if (is.logical(nConWithErr)) {
    if (length(nConWithErr) != nConVar) {
      stop("Error in genMixedData: length of nConWithErr must equal nConVar")
    }
  }
  else if (length(nConWithErr) != 1 | !is.numeric(nConWithErr)) {
    stop("Error in genMixedData: improper input to nConWithErr")
  }
  else if (nConWithErr > nConVar) {
    stop("Error in genMixedData: nConWithErr must be less than nConVar")
  }
  else {
    nConWithErr = rep(c(T, F), c(nConWithErr, nConVar - nConWithErr))
  }
  if (nCatLevels%%length(popProportions) != 0) {
    stop("Error in genMixedData: number of categorical levels must\n        be a multiple of number of populations")
  }
  OVERLAP_DEFAULT = 0.01
  trueID = sample(x = 1:length(popProportions), size = sampSize, 
                  replace = T, prob = popProportions)
  trueMus = c(0, -2 * qnorm(OVERLAP_DEFAULT/2))
  muVec = trueMus[trueID]
  conVars = matrix(rpgnorm(sampSize * nConVar,p, mean = muVec, 
                         sigma = 1), nrow = sampSize, ncol = nConVar)
  errVariance = (abs(diff(trueMus))/(2 * qnorm(conErrLev/2)))^2 - 
    1
  for (i in 1:nConVar) {
    if (nConWithErr[i]) {
      ithError = rnorm(sampSize, mean = 0, sd = sqrt(errVariance))
      conVars[, i] = conVars[, i] + ithError
    }
  }
  rightCatProb = (1 - OVERLAP_DEFAULT/2)/(nCatLevels/2)
  wrongCatProb = (OVERLAP_DEFAULT/2)/(nCatLevels/2)
  popProb1 = rep(c(rightCatProb, wrongCatProb), each = nCatLevels/2)
  popProb2 = rep(c(wrongCatProb, rightCatProb), each = nCatLevels/2)
  numInPop1 = sum(trueID == 1)
  if (nCatWithErr < nCatVar) {
    pop1CatNoErr = sample(x = 1:nCatLevels, size = numInPop1 * 
                            (nCatVar - nCatWithErr), replace = T, prob = popProb1)
    pop2CatNoErr = sample(x = 1:nCatLevels, size = (sampSize - 
                                                      numInPop1) * (nCatVar - nCatWithErr), replace = T, 
                          prob = popProb2)
    catNoErr = rep(NA, sampSize * (nCatVar - nCatWithErr))
    vectIdNoErr = rep(trueID, nCatVar - nCatWithErr)
    catNoErr[vectIdNoErr == 1] = pop1CatNoErr
    catNoErr[vectIdNoErr == 2] = pop2CatNoErr
    catNoErr = matrix(catNoErr, ncol = nCatVar - nCatWithErr)
  }
  else {
    catNoErr = c()
  }
  if (nCatWithErr > 0) {
    rightCatProbErr = (1 - catErrLev/2)/(nCatLevels/2)
    wrongCatProbErr = (catErrLev/2)/(nCatLevels/2)
    popProb1Err = rep(c(rightCatProbErr, wrongCatProbErr), 
                      each = nCatLevels/2)
    popProb2Err = rep(c(wrongCatProbErr, rightCatProbErr), 
                      each = nCatLevels/2)
    pop1CatWithErr = sample(x = 1:nCatLevels, size = numInPop1 * 
                              nCatWithErr, replace = T, prob = popProb1Err)
    pop2CatWithErr = sample(x = 1:nCatLevels, size = (sampSize - 
                                                        numInPop1) * nCatWithErr, replace = T, prob = popProb2Err)
    catWithErr = rep(NA, sampSize * nCatWithErr)
    vectIdWithErr = rep(trueID, nCatWithErr)
    catWithErr[vectIdWithErr == 1] = pop1CatWithErr
    catWithErr[vectIdWithErr == 2] = pop2CatWithErr
    catWithErr = matrix(catWithErr, ncol = nCatWithErr)
  }
  else {
    catWithErr = c()
  }
  catVars = cbind(catNoErr, catWithErr)
  return(list(trueID = trueID, trueMus = trueMus, conVars = conVars, 
              errVariance = errVariance, popProbsNoErr = data.frame(popProb1 = popProb1, 
                                                                    popProb2 = popProb2, stringsAsFactors = TRUE), popProbsWithErr = data.frame(popProb1Err = popProb1Err, 
                                                                                                                                                popProb2Err = popProb2Err, stringsAsFactors = TRUE), 
              catVars = catVars))
  }
genPGNormSimpleData <- function(size,conErrLev, catErrLev,p) {
  tmpDat <- genPGNormMixedData(sampSize = size, nConVar = 2, nCatVar = 2,
                         nCatLevels = 4, nConWithErr = 2, nCatWithErr = 2,
                         popProportions = c(.5, .5), conErrLev = conErrLev,
                         catErrLev = catErrLev,p=p)
  tmpDat <- within(tmpDat, conVars <- as.data.frame(scale(conVars)))
  tmpDat <- within(tmpDat, catVarsFac <- as.data.frame(lapply(
    as.data.frame(catVars), factor)))
  within(tmpDat, catVarsDum <- as.data.frame(dummyCodeFactorDf(
    catVarsFac)))
}


#############plot example simulated data##########
set.seed(1)
with(data = genPGNormSimpleData(size=200,conErrLev = 0.3, catErrLev = 0.3,p=0.7019),
     expr = pairs(cbind(conVars, jitter(catVars)), pch = c(1, 3)[trueID],
                  col = trueID, main = "Pair plot of overlaping 0.3 \nfor both types of variables",cex.main=0.9))
set.seed(1)
with(data = genPGNormSimpleData(size=200,conErrLev = 0.1, catErrLev = 0.5,p=0.7019),
     expr = pairs(cbind(conVars, jitter(catVars)), pch = c(1, 3)[trueID],
                  col = trueID, main = "Pair plot of overlaping 0.1 for continuous variables \nand 0.5 for categorical variables",cex.main=0.9))
set.seed(1)
sim1<-genPGNormSimpleData(size = 1e4,conErrLev = 0.05, catErrLev = 0.05,p=0.7019)
cluster1<-sim1[["trueID"]]
data1<-cbind(scale(sim1[["conVars"]]), sim1[["catVars"]])
View(data1)
plot(density(data1[,1]))



############Simulation study#################
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
  return(tryCatch(clustx3(datasets[[x]],1e4),error=function(e) NULL))
})

tf_null = sapply(1:n,function(x) is.null(clustx3_8[[x]]))
id_more = c(1:n)[tf_null]
clustx3_8more<-lapply(id_more,function(x) {
  set.seed(1)
  return(tryCatch(clustx3(datasets[[x]],1e5),error=function(e) NULL))
})

########kurtosis6
set.seed(1)
datasets<-lapply(1:n,function(x)
  genPGNormSimpleData(size=1e3,conErrLev = overlap_levels[x,1], 
                      catErrLev = overlap_levels[x,2],p=0.7785)
)
clustx3_6<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3(datasets[[x]],1e4),error=function(e) NULL))
})

tf_null = sapply(1:n,function(x) is.null(clustx3_6[[x]]))
id_more = c(1:n)[tf_null]
clustx3_6more<-lapply(id_more,function(x) {
  set.seed(1)
  return(tryCatch(clustx3(datasets[[x]],1e5),error=function(e) NULL))
})

########kurtosis4
set.seed(1)
datasets<-lapply(1:n,function(x)
  genPGNormSimpleData(size=1e3,conErrLev = overlap_levels[x,1], 
                      catErrLev = overlap_levels[x,2],p=0.9012)
)
clustx3_4<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3(datasets[[x]],1e4),error=function(e) NULL))
})

tf_null = sapply(1:n,function(x) is.null(clustx3_4[[x]]))
id_more = c(1:n)[tf_null]
clustx3_4more<-lapply(id_more[3],function(x) {
  set.seed(1)
  return(tryCatch(clustx3(datasets[[x]],1e5),error=function(e) NULL))
})


################################################
#######################log-norm ################
################################################

##############data generation function#########
genLogNormMixedData<-function (sampSize, nConVar, nCatVar, nCatLevels, nConWithErr, 
                               nCatWithErr, popProportions, conErrLev, catErrLev,sdlog) 
  
{
  require("pgnorm")
  if (length(popProportions) > 2) {
    stop("Error in genMixedData: More than two populations not currently supported")
  }
  if (is.logical(nConWithErr)) {
    if (length(nConWithErr) != nConVar) {
      stop("Error in genMixedData: length of nConWithErr must equal nConVar")
    }
  }
  else if (length(nConWithErr) != 1 | !is.numeric(nConWithErr)) {
    stop("Error in genMixedData: improper input to nConWithErr")
  }
  else if (nConWithErr > nConVar) {
    stop("Error in genMixedData: nConWithErr must be less than nConVar")
  }
  else {
    nConWithErr = rep(c(T, F), c(nConWithErr, nConVar - nConWithErr))
  }
  if (nCatLevels%%length(popProportions) != 0) {
    stop("Error in genMixedData: number of categorical levels must\n        be a multiple of number of populations")
  }
  #OVERLAP_DEFAULT = 0.01
  trueID = sample(x = 1:length(popProportions), size = sampSize, 
                  replace = T, prob = popProportions)
  #trueMus = c(0, -2 * qnorm(OVERLAP_DEFAULT/2))
  #muVec = trueMus[trueID]
  conVars = matrix(rlnorm(sampSize * nConVar, meanlog = 1, 
                          sdlog = sdlog), nrow = sampSize, ncol = nConVar)
  f = function(x){
    x1=(1-sdlog^2)-sqrt((1-sdlog^2)^2-(1+2*sdlog^2 * log(x*sdlog*sqrt(2*pi))))
    x2=(1-sdlog^2)+sqrt((1-sdlog^2)^2-(1+2*sdlog^2 * log(x*sdlog*sqrt(2*pi))))
    return(abs(plnorm(exp(x2),meanlog = 1, sdlog = sdlog)-plnorm(exp(x1),meanlog = 1, sdlog = sdlog)-1+conErrLev))
  }
  result = optim(par=0.2,fn=f,method = "SANN")$par
  x1=(1-sdlog^2)-sqrt((1-sdlog^2)^2-(1+2*sdlog^2 * log(result*sdlog*sqrt(2*pi))))
  x2=(1-sdlog^2)+sqrt((1-sdlog^2)^2-(1+2*sdlog^2 * log(result*sdlog*sqrt(2*pi))))
  a = abs(exp(x2)-exp(x1))
  conVars[trueID==1,]<-conVars[trueID==1,]+a
  
  rightCatProb = (1 - OVERLAP_DEFAULT/2)/(nCatLevels/2)
  wrongCatProb = (OVERLAP_DEFAULT/2)/(nCatLevels/2)
  popProb1 = rep(c(rightCatProb, wrongCatProb), each = nCatLevels/2)
  popProb2 = rep(c(wrongCatProb, rightCatProb), each = nCatLevels/2)
  numInPop1 = sum(trueID == 1)
  if (nCatWithErr < nCatVar) {
    pop1CatNoErr = sample(x = 1:nCatLevels, size = numInPop1 * 
                            (nCatVar - nCatWithErr), replace = T, prob = popProb1)
    pop2CatNoErr = sample(x = 1:nCatLevels, size = (sampSize - 
                                                      numInPop1) * (nCatVar - nCatWithErr), replace = T, 
                          prob = popProb2)
    catNoErr = rep(NA, sampSize * (nCatVar - nCatWithErr))
    vectIdNoErr = rep(trueID, nCatVar - nCatWithErr)
    catNoErr[vectIdNoErr == 1] = pop1CatNoErr
    catNoErr[vectIdNoErr == 2] = pop2CatNoErr
    catNoErr = matrix(catNoErr, ncol = nCatVar - nCatWithErr)
  }
  else {
    catNoErr = c()
  }
  if (nCatWithErr > 0) {
    rightCatProbErr = (1 - catErrLev/2)/(nCatLevels/2)
    wrongCatProbErr = (catErrLev/2)/(nCatLevels/2)
    popProb1Err = rep(c(rightCatProbErr, wrongCatProbErr), 
                      each = nCatLevels/2)
    popProb2Err = rep(c(wrongCatProbErr, rightCatProbErr), 
                      each = nCatLevels/2)
    pop1CatWithErr = sample(x = 1:nCatLevels, size = numInPop1 * 
                              nCatWithErr, replace = T, prob = popProb1Err)
    pop2CatWithErr = sample(x = 1:nCatLevels, size = (sampSize - 
                                                        numInPop1) * nCatWithErr, replace = T, prob = popProb2Err)
    catWithErr = rep(NA, sampSize * nCatWithErr)
    vectIdWithErr = rep(trueID, nCatWithErr)
    catWithErr[vectIdWithErr == 1] = pop1CatWithErr
    catWithErr[vectIdWithErr == 2] = pop2CatWithErr
    catWithErr = matrix(catWithErr, ncol = nCatWithErr)
  }
  else {
    catWithErr = c()
  }
  catVars = cbind(catNoErr, catWithErr)
  return(list(trueID = trueID, trueMus = trueMus, conVars = conVars, 
              errVariance = errVariance, popProbsNoErr = data.frame(popProb1 = popProb1, 
                                                                    popProb2 = popProb2, stringsAsFactors = TRUE), popProbsWithErr = data.frame(popProb1Err = popProb1Err, 
                                                                                                                                                popProb2Err = popProb2Err, stringsAsFactors = TRUE), 
              catVars = catVars))
}

genLogNormSimpleData <- function(size,conErrLev, catErrLev,sdlog) {
  tmpDat <- genLogNormMixedData(sampSize = size, nConVar = 2, nCatVar = 2,
                                nCatLevels = 4, nConWithErr = 2, nCatWithErr = 2,
                                popProportions = c(.5, .5), conErrLev = conErrLev,
                                catErrLev = catErrLev,sdlog=sdlog)
  tmpDat <- within(tmpDat, conVars <- as.data.frame(scale(conVars)))
  tmpDat <- within(tmpDat, catVarsFac <- as.data.frame(lapply(
    as.data.frame(catVars), factor)))
  within(tmpDat, catVarsDum <- as.data.frame(dummyCodeFactorDf(
    catVarsFac)))
}

#############plot example simulated data##########
set.seed(1)
with(data = genLogNormSimpleData(size =200,conErrLev = 0.3, catErrLev = 0.3,sdlog=0.3143),
     expr = pairs(cbind(conVars, jitter(catVars)), pch = c(1, 3)[trueID],
                  col = trueID, main="Pair plot of overlaping 0.3 \nfor both types of variables",cex.main=0.9))
set.seed(1)
with(data = genLogNormSimpleData(size =200,conErrLev = 0.1, catErrLev = 0.5,sdlog=0.3143),
     expr = pairs(cbind(conVars, jitter(catVars)), pch = c(1, 3)[trueID],
                  col = trueID, main = "Pair plot of overlaping 0.1 for continuous variables \nand 0.5 for categorical variables",cex.main=0.9))

############Simulation study#################
overlap_prop<-c(0.05,0.2,0.35,0.5)
overlap_levels<-expand.grid(overlap_prop,overlap_prop)
n<-nrow(overlap_levels)
########skewness1
set.seed(1)
datasets<-lapply(1:n,function(x)
  genLogNormSimpleData(size=1e3,conErrLev = overlap_levels[x,1], 
                       catErrLev = overlap_levels[x,2],sdlog = 0.3143)
)
clustx3_1<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3(datasets[[x]],1e4),error=function(e) NULL))
})

tf_null = sapply(1:n,function(x) anyNA(clustx3_1[[x]]))
id_more = c(1:n)[tf_null]
clustx3_1more<-lapply(id_more,function(x) {
  set.seed(1)
  return(tryCatch(clustx3(datasets[[x]],1e5),error=function(e) NULL))
})

########skewness2.5
set.seed(1)
datasets<-lapply(1:n,function(x)
  genLogNormSimpleData(size=1e3,conErrLev = overlap_levels[x,1], 
                       catErrLev = overlap_levels[x,2],sdlog = 0.6409)
)
clustx3_2<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3(datasets[[x]],1e4),error=function(e) NULL))
})




########skewness9
set.seed(1)
datasets<-lapply(1:n,function(x)
  genLogNormSimpleData(size=1e3,conErrLev = overlap_levels[x,1], 
                       catErrLev = overlap_levels[x,2],sdlog = 1.1310)
)
clustx3_9<-lapply(1:n,function(x) {
  set.seed(1)
  return(tryCatch(clustx3(datasets[[x]],1e4),error=function(e) NULL))
})






############plot default overlap#####################

b=seq(0,10,length.out=1000)
sdlog=0.3143
conErrLev=0.01
f = function(x){
  x1=(1-sdlog^2)-sqrt((1-sdlog^2)^2-(1+2*sdlog^2 * log(x*sdlog*sqrt(2*pi))))
  x2=(1-sdlog^2)+sqrt((1-sdlog^2)^2-(1+2*sdlog^2 * log(x*sdlog*sqrt(2*pi))))
  return(abs(plnorm(exp(x2),meanlog = 1, sdlog = sdlog)-plnorm(exp(x1),meanlog = 1, sdlog = sdlog)-1+conErrLev))
}
result = optim(par=0.2,fn=f,method = "SANN")$par
x1=(1-sdlog^2)-sqrt((1-sdlog^2)^2-(1+2*sdlog^2 * log(x*sdlog*sqrt(2*pi))))
x2=(1-sdlog^2)+sqrt((1-sdlog^2)^2-(1+2*sdlog^2 * log(x*sdlog*sqrt(2*pi))))
a = abs(exp(x2)-exp(x1))

plot(b,dlnorm(b,meanlog=1,sdlog = 0.3143),type='l',ylab='density',
     main='Overlapping level 0.01 \nbetween two log-normal densities',xlab='',cex.axis=0.7)
x1=(1-sdlog^2)-sqrt((1-sdlog^2)^2-(1+2*sdlog^2 * log(x*sdlog*sqrt(2*pi))))
axis(2, at=dlnorm(exp(x1),meanlog = 1,sdlog=0.3143),labels='y',cex.axis=0.9)
axis(1, at=exp(x1),labels='x1',cex.axis=0.5)
axis(1, at=exp(x2),labels='x2',cex.axis=0.5)
lines(b+4.690593,dlnorm(b,meanlog=1,sdlog=0.3143),col='red')
lines(b,rep(dlnorm(exp(x1),meanlog = 1,sdlog=0.3143),length(b)),lty=2)
legend('topright',c('density 1','density 2'),col=c('black','red'),lty=c(1,1),cex=0.5)


