library(rjags)
library(R2jags)

## Model for Analysis
ModelCode <- "model{
  ## Priors
  
  p.s ~ dunif(0,1) #probability to be seen
  p.maybe ~ dunif(0,1) #probability a plant that says maybe they were seen
  H ~ dlnorm(0, 1/10^2)
  
  ## Model
  
  M.maybe ~ dbin(p.maybe, M)
  M.yes ~ dbin(p.s, M-M.maybe)
  M.no <- M - M.maybe
  census.unkn ~ dbin(p.s,total)
  total <- round(H) + M.maybe #number of real homeless
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
                              "p.maybe" = 0.5,
                              "H" = H.inits),
                         list("p.s" = 0.6,
                              "p.maybe" = 0.3,
                              "H" = H.inits),
                         list("p.s" = 0.4,
                              "p.maybe" = 0.1,
                              "H" = H.inits))
  vars.monitor <- c("p.s", "p.maybe", "H")
  data <- list("M" = simdata$M,
               "M.yes" = simdata$M.yes,
               "M.maybe" = simdata$M.maybe,
               "census.unkn" = simdata$census.unkn)
  
  jagsfit <- jags(data=data, n.chains=3, inits=initial.values,
                  parameters.to.save=vars.monitor, n.iter=30000, n.burnin=15000,n.thin=1,
                  DIC=TRUE, model.file=textConnection(ModelCode))
  
  # saveRDS(jagsfit$BUGSoutput$sims.array,paste0("res15/res",set,".rds")) ## Save results for small cities
  saveRDS(jagsfit$BUGSoutput$sims.array,paste0("res100/res",set,".rds")) ## Save results for large cities
}
