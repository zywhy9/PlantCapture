library(rjags)
library(R2jags)

ModelCode <- "model{
  #Priors
  p.s ~ dunif(0,1) #probability to be seen
  p.ni.maybe ~ dunif(0,1) #probability a plant that says maybe they were seen and not interviewed
  p.s.i ~ dunif(0,1) #probability a plant was seen that he got interviewed
  H.prior ~ dlnorm(0, 1/10^2) 
  H <- round(H.prior)
  
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
  H.s ~ dbin(p.s, H)
  H.i ~ dbin(p.s.i, H.s)
  census.unkn ~ dsum(M.s.maybe, H.s)
}
"

for(set in 1:1000){
  
  set.seed(set)
  
  # simdata <- readRDS(paste0("sim15/sim",set,".rds")) ## Read simulated dataset for small cities
  simdata <- readRDS(paste0("sim100/sim",set,".rds")) ## Read simulated dataset for large cities
  
  M.s.maybe.inits <- round(simdata$M.maybe*0.7) ## Set initial values
  H.s.inits <- simdata$census.unkn - M.s.maybe.inits ## Set initial values
  H.inits <- H.s.inits/0.7 ## Set initial values
  
  initial.values <- list(list("p.s" = 0.5,
                              "p.s.i" = 0.5,
                              "p.ni.maybe" = 0.5,
                              "M.s.maybe" = M.s.maybe.inits,
                              "H.s" = H.s.inits,
                              "H.prior" = H.inits),
                         list("p.s" = 0.6,
                              "p.s.i" = 0.3,
                              "p.ni.maybe" = 0.8,
                              "M.s.maybe" = M.s.maybe.inits,
                              "H.s" = H.s.inits,
                              "H.prior" = H.inits),
                         list("p.s" = 0.4,
                              "p.s.i" = 0.8,
                              "p.ni.maybe" = 0.1,
                              "M.s.maybe" = M.s.maybe.inits,
                              "H.s" = H.s.inits,
                              "H.prior" = H.inits))
  
  vars.monitor <- c("p.s", "p.ni.maybe", "p.s.i", "H")
  data <- list("M" = simdata$M,
               "M.i" = simdata$M.i,
               "M.yes" = simdata$M.yes,
               "M.maybe" = simdata$M.maybe,
               "census.unkn" = simdata$census.unkn,
               "H.i" = simdata$H.i)
  
  jagsfit <- jags(data=data, n.chains=3, inits=initial.values,
                  parameters.to.save=vars.monitor, n.iter=30000, n.burnin=15000,n.thin=1,
                  DIC=TRUE, model.file=textConnection(ModelCode))
  
  # saveRDS(jagsfit$BUGSoutput$sims.array,paste0("res15/res",set,".rds")) ## Save results for small cities
  saveRDS(jagsfit$BUGSoutput$sims.array,paste0("res100/res",set,".rds")) ## Save results for large cities
}
