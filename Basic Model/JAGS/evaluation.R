library(MCMCvis)

## Read Data
nset <- 1000
niter <- 15000
nchain <- 3
nitert <- niter * nchain
npar <- 3
results <- matrix(NA, nrow=nitert*nset, ncol=npar)
colnames(results) <- c("H","p.s","p.maybe")

true.value <- c(1500,0.7,0.2)

cp <- rep(0,npar)
psd <- matrix(0,nrow=nset,ncol=npar)
ci.l <- matrix(0,nrow=nset,ncol=npar)
ci.u <- matrix(0,nrow=nset,ncol=npar)

for(i in 1:nset){
  temp <- readRDS(paste0("res25/res",i,".rds"))
  tempH <- as.vector(temp[,,"H"])
  tempps <- as.vector(temp[,,"p.s"])
  temppmay <- as.vector(temp[,,"p.maybe"])
  temp <- cbind(tempH,tempps,temppmay)
  results[(nitert*(i-1)+1):(nitert*i),] <- temp
  mcmc.output <- list(chain1=temp[1:niter,],chain2=temp[(niter+1):(niter*2),],chain3=temp[(niter*2+1):nitert,])
  res <- MCMCsummary(object = mcmc.output)
  ci.l[i,] <- res$`2.5%`
  ci.u[i,] <- res$`97.5%`
  psd[i,] <- res$sd
  for(j in 1:npar){
    if(true.value[j]<=res$`97.5%`[j] && true.value[j]>=res$`2.5%`[j]){
      cp[j] <- cp[j] + 1
    }
  }
}

rm(temp,tempH,tempps,temppmay)

## Evaluation
pmn <- apply(results,2,mean)
pse <- apply(results,2,sd)

bias <- pmn - true.value
rbias <- bias / true.value

mse <- rep(0,npar)
for(i in 1:npar){
  mse[i] <- mean((results[,i]-true.value[i])^2)
}
rmse <- sqrt(mse)
rrmse <- rmse/true.value

cp <- cp/nset
lcp <- apply(ci.u-ci.l,2,mean)
epsd <- apply(psd,2,mean)
list(mean=round(pmn,4),epsd=round(epsd,4),se=round(pse,4),rbias=round(rbias,4),rrmse=round(rrmse,4),cp=cp,lcp=round(lcp,4))

saveRDS(list(mean=pmn,epsd=epsd,se=pse,rbias=rbias,rrmse=rrmse,cp=cp,lcp=lcp),file="eva.rds")
