library(rjags)
library(R2jags)
load("simmodel.RData")

for(set in 1:1000){
  
  set.seed(set)
  
  simdata <- readRDS(paste0("../simdata100/sim",set,".rds"))
  
  M.s.maybe.inits <- round(simdata$M.maybe*0.7)
  H.s.inits <- simdata$census.unkn - M.s.maybe.inits
  H.inits <- H.s.inits/0.7
  
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
  
  saveRDS(jagsfit$BUGSoutput$sims.array,paste0("res100/res",set,".rds"))
}
