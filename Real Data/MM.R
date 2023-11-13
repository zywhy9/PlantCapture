library(rjags)
library(R2jags)
library(MCMCvis)

## Extract Data from Literature
cityname <- c("Chicago", "New Orleans", "Phoenix", "New York", "Los Angeles")
M <- c(13, 58, 26, 94, 25) ### Number of plants in the sites where enumerators were seen for each city from page 18 of Martin et al. (1997)
Y <- c(11, 109, 104, 1240, 217) ### Census count in the sites where enumerators were seen for each city from page 20 of Martin et al. (1997)
p.interview <- c(0.07, 0.67, 0.44, 0.37, 0.33) ### Percent of interviewed plants for each city from page 15 of Martin et al. (1997)
p.yes <- c(0, 0.1, 0.08, 0.17, 0.02) ### Percent of plants answering yes for each city from page 15 of Martin et al. (1997)
p.maybe <- c(0.18, 0.07, 0.03, 0.12, 0.04) ### Percent of plants answering maybe for each city from page 15 of Martin et al. (1997)
p.no <- c(0.25, 0.10, 0.10, 0.20, 0.13) ### Percent of plants answering no for each city from page 15 of Martin et al. (1997)
p.invalid <- c(0.5, 0.05, 0.36, 0.14, 0.48) ### Percent of plants did not see enumerators for each city from page 15 of Martin et al. (1997)

## Final Data
M.i <- round(M * p.interview / (1 - p.invalid))
M.yes <- round(M * p.yes / (1 - p.invalid))
M.no <- round(M * p.no / (1 - p.invalid))
M.maybe <- M - M.i - M.yes - M.no
census.unkn <- Y - M.yes - M.i


## Model
ModelCode <- "model{
  #Priors
  p.s ~ dunif(0,1) #probability to be seen
  p.ni.maybe ~ dunif(0,1) #probability a plant that says maybe they were seen and not interviewed
  p.s.i ~ dunif(0,1) #probability a plant was seen that he got interviewed
  
  #Transformation
  p.i <- p.s * p.s.i
  p.s.maybe <- p.s * (1 - p.s.i) * p.ni.maybe
  p.yes <- p.s * (1 - p.s.i) * (1 - p.ni.maybe)
  p.ns.maybe <- (1 - p.s) * p.ni.maybe
  p.maybe <- p.s.maybe + p.ns.maybe
  
  #Model
  M.i ~ dbin(p.i, M)
  M.yes ~ dbin(p.yes / (1 - p.i), M - M.i)
  M.maybe ~ dbin(p.maybe / (1 - p.yes - p.i), M - M.i - M.yes)
  M.s.maybe ~ dbin(p.s.maybe / p.maybe, M.maybe)
  H.s <- census.unkn - M.s.maybe
  H0 <- H.s / p.s
  H ~ dnorm(H0, 1/(H0*(1-p.s)/p.s))
}
"

## Analysis
for(set in 1:5){
  
  set.seed(set)
  
  M.s.maybe.inits <- round(M.maybe[set]*0.3)
  
  initial.values <- list(list("p.s" = 0.5,
                              "p.s.i" = 0.5,
                              "p.ni.maybe" = 0.5,
                              "M.s.maybe" = M.s.maybe.inits),
                         list("p.s" = 0.6,
                              "p.s.i" = 0.3,
                              "p.ni.maybe" = 0.8,
                              "M.s.maybe" = M.s.maybe.inits),
                         list("p.s" = 0.4,
                              "p.s.i" = 0.8,
                              "p.ni.maybe" = 0.1,
                              "M.s.maybe" = M.s.maybe.inits))
  
  vars.monitor <- c("p.s", "p.ni.maybe", "p.s.i", "H0", "H")
  data <- list("M" = M[set],
               "M.i" = M.i[set],
               "M.yes" = M.yes[set],
               "M.maybe" = M.maybe[set],
               "census.unkn" = census.unkn[set])
  
  jagsfit <- jags(data=data, n.chains=3, inits=initial.values,
                  parameters.to.save=vars.monitor, n.iter=30000, n.burnin=15000,n.thin=1,
                  DIC=TRUE, model.file=textConnection(ModelCode))
  
  saveRDS(jagsfit$BUGSoutput$sims.array,paste0("MM/",cityname[set],".rds"))
}

## Evaluation
nset <- 5
niter <- 15000
nchain <- 3
nitert <- niter * nchain
npar <- 5
results <- matrix(NA, nrow=nitert*nset, ncol=npar)
colnames(results) <- c("H", "H0", "p.s", "p.ni.maybe", "p.s.i")

psd <- matrix(0,nrow=nset,ncol=npar)
pmean <- matrix(0,nrow=nset,ncol=npar)
pmed <- matrix(0,nrow=nset,ncol=npar)
ci.l <- matrix(0,nrow=nset,ncol=npar)
ci.u <- matrix(0,nrow=nset,ncol=npar)

for(i in 1:nset){
  temp <- readRDS(paste0("MM/",cityname[i],".rds"))
  tempH <- as.vector(temp[,,"H"])
  tempH0 <- as.vector(temp[,,"H0"])
  tempps <- as.vector(temp[,,"p.s"])
  temppmay <- as.vector(temp[,,"p.ni.maybe"])
  temppsi <- as.vector(temp[,,"p.s.i"])
  temp <- cbind(tempH,tempH0,tempps,temppmay,temppsi)
  results[(nitert*(i-1)+1):(nitert*i),] <- temp
  mcmc.output <- list(chain1=temp[1:niter,],chain2=temp[(niter+1):(niter*2),],chain3=temp[(niter*2+1):nitert,])
  res <- MCMCsummary(object = mcmc.output, func="median", func_name = "median")
  ci.l[i,] <- res$`2.5%`
  ci.u[i,] <- res$`97.5%`
  psd[i,] <- res$sd
  pmed[i,] <- res$median
  pmean[i,] <- res$mean
}

rm(temp,tempH,tempH0,tempps,temppmay,temppsi)

eva <- cbind(as.vector(t(pmean)),as.vector(t(pmed)),as.vector(t(psd)),as.vector(t(ci.l)),as.vector(t(ci.u)))
colnames(eva) <- c("Mean","Median", "SD", "Lower", "Upper")
rownames(eva) <- rep(c("H", "H0", "p.s", "p.ni.maybe", "p.s.i"),5)
round(eva,2)