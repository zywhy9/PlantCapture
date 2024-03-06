library(MCMCvis)

## Read Data

nset <- 1000 ## Number of datasets
niter <- 15000 ## Final number of iteration per chain (after burnin)
nchain <- 3 ## Number of MCMC chains
nitert <- niter * nchain ## Total iterations
npar <- 3 ## Number of parameters of interest
results <- matrix(NA, nrow=nitert*nset, ncol=npar)
colnames(results) <- c("H","p.c","p.maybe")

true.value <- c(1500,0.7,0.2) ## True value of parameters


cp <- rep(0,npar)
psd <- matrix(0,nrow=nset,ncol=npar)
pmean <- matrix(0,nrow=nset,ncol=npar)
pmed <- matrix(0,nrow=nset,ncol=npar)
ci.l <- matrix(0,nrow=nset,ncol=npar)
ci.u <- matrix(0,nrow=nset,ncol=npar)

for(i in 1:nset){
  temp <- readRDS(paste0("res15/res",i,".rds")) ## Read dataset for small cities. Change the address for large cities.
  tempH <- as.vector(temp[,,"H"]) ## Transform matrix to vector for each parameter
  tempps <- as.vector(temp[,,"p.c"])
  temppmay <- as.vector(temp[,,"p.maybe"])
  temp <- cbind(tempH,tempps,temppmay) ## Save data to a matrix with columns as variables and rows as iterations
  results[(nitert*(i-1)+1):(nitert*i),] <- temp ## Save this dataset to the final results matrix
  mcmc.output <- list(chain1=temp[1:niter,],chain2=temp[(niter+1):(niter*2),],chain3=temp[(niter*2+1):nitert,])
  res <- MCMCsummary(object = mcmc.output, func=median, func_name = "median")
  ci.l[i,] <- res$`2.5%` ## 2.5 percentile of MCMC samples
  ci.u[i,] <- res$`97.5%` ## 97.5 percentile of MCMC samples
  psd[i,] <- res$sd ## Standard deviation of MCMC samples
  pmean[i,] <- res$mean ## Mean of MCMC samples
  pmed[i,] <- res$median ## Median of MCMC samples
  for(j in 1:npar){
    if(true.value[j]<=res$`97.5%`[j] && true.value[j]>=res$`2.5%`[j]){
      cp[j] <- cp[j] + 1 ## Calculate coverage probability
    }
  }
}

rm(temp,tempH,tempps,temppmay) ## Remove temp variables to save space

## Evaluation

tpmed <- apply(pmed,2,mean) ## Expected posterior median
tpmn <- apply(results,2,mean) ## Expected posterior mean
psemean <- apply(pmean,2,sd) ## Standard error of posterior mean
psemed <- apply(pmed,2,sd) ## Standard error of posterior median
epsd <- apply(psd,2,mean) ## Expected posterior standard deviation

bias <- tpmed - true.value ## Bias
rbias <- bias / true.value ## Relative bias

mse <- rep(0,npar)
for(i in 1:npar){
  mse[i] <- mean((pmed[,i]-true.value[i])^2) ## Mean squared error
}
rmse <- sqrt(mse) ## Root mean squared error
rrmse <- rmse/true.value ## Relative RMSE

cp <- cp/nset ## Coverage probability
lci <- apply(ci.u-ci.l,2,mean) ## Average length of credible interval
list(mean=round(tpmn,4),
     median=round(tpmed,4),
     epsd=round(epsd,4),
     se=round(psemed,4),
     rbias=round(rbias,4),
     rrmse=round(rrmse,4),
     cp=cp,
     lci=round(lci,4))

# saveRDS(list(mean=tpmn,
#              median=tpmed,
#              mode=tpmode,
#              epsd=epsd,
#              se=psemed,
#              semean=psemean,
#              rbias=rbias,
#              rrmse=rrmse,
#              cp=cp,lci=lci),file="eva100.rds")
