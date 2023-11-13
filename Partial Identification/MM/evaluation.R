library(MCMCvis)

## Read Data
nset <- 1000
niter <- 15000
nchain <- 3
nitert <- niter * nchain
npar <- 6
results <- matrix(NA, nrow=nitert*nset, ncol=npar)
colnames(results) <- c("H1.t", "H2.t", "p.s[1]", "p.s[2]", "p.ni.maybe", "p.s.i")

true.value <- c(1500, 1500, 0.9, 0.4, 0.2, 0.8)

cp <- rep(0,npar)
psd <- matrix(0,nrow=nset,ncol=npar)
ci.l <- matrix(0,nrow=nset,ncol=npar)
ci.u <- matrix(0,nrow=nset,ncol=npar)

for(i in 1:nset){
  temp <- readRDS(paste0("res25/res",i,".rds"))
  tempH1 <- as.vector(temp[,,"H1.t"])
  tempH2 <- as.vector(temp[,,"H2.t"])
  tempps1 <- as.vector(temp[,,"p.s[1]"])
  tempps2 <- as.vector(temp[,,"p.s[2]"])
  temppmay <- as.vector(temp[,,"p.ni.maybe"])
  temppsi <- as.vector(temp[,,"p.s.i"])
  temp <- cbind(tempH1,tempH2,tempps1,tempps2,temppmay,temppsi)
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

rm(temp,tempH,tempps,temppmay,temppsi)

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

saveRDS(list(mean=pmn,epsd=epsd,se=pse,rbias=rbias,rrmse=rrmse,cp=cp,lcp=lcp),file="eva25.rds")
