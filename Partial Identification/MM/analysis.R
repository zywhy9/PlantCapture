library(rjags)
library(R2jags)
load("simmodel.RData")

for(set in 1:1000){
  
  set.seed(set)
  
  simdata <- readRDS(paste0("../sim100/sim",set,".rds"))
  
  M.s.maybe.inits <- round(simdata$M.maybe*0.3)
  # H.s.inits <- simdata$temp - M.s.maybe.inits
  # H.inits <- H.s.inits/0.7
  
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
  
  vars.monitor <- c("p.s", "p.ni.maybe", "p.s.i", "H1", "H2")
  data <- list("M" = simdata$M,
               "M.i" = simdata$M.i,
               "M.yes" = simdata$M.yes,
               "M.maybe" = simdata$M.maybe,
               "census.unkn" = simdata$census.unkn)
  
  jagsfit <- jags(data=data, n.chains=3, inits=initial.values,
                  parameters.to.save=vars.monitor, n.iter=30000, n.burnin=15000,n.thin=1,
                  DIC=TRUE, model.file=textConnection(ModelCode))
  
  # saveRDS(jagsfit$BUGSoutput$sims.array,paste0("res100/res",set,".rds"))
}
