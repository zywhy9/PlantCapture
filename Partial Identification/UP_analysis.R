library(rjags)
library(R2jags)

ModelCode <- "model{
  #Priors
  p.c ~ dunif(0,1) #probability to be seen
  p.ni.maybe ~ dunif(0,1) #probability a plant that says maybe they were seen and not interviewed
  p.c.i ~ dunif(0,1) #probability a plant was seen that he got interviewed
  
  #Transformation
  p.i <- p.c * p.c.i
  p.c.maybe <- p.c * (1 - p.c.i) * p.ni.maybe
  p.yes <- p.c * (1 - p.c.i) * (1 - p.ni.maybe)
  p.ns.maybe <- (1 - p.c) * p.ni.maybe
  p.maybe <- p.c.maybe + p.ns.maybe
  
  #Model
  M.i ~ dbin(p.i, M)
  M.yes ~ dbin(p.yes / (1 - p.i), M - M.i)
  M.maybe ~ dbin(p.maybe / (1 - p.yes - p.i), M - M.i - M.yes)
  M.s.maybe ~ dbin(p.c.maybe / p.maybe, M.maybe)
  H.s <- census.unkn - M.s.maybe
  H.i ~ dbin(p.c.i, H.s)
  H1 <- H.s / p.c
  H2 ~ dnorm(H1, 1/(H1*(1-p.c)/p.c))
}
"

for(set in 1:1000){
  
  set.seed(set)
  
  # simdata <- readRDS(paste0("sim15/sim",set,".rds")) ## Read simulated dataset for small cities
  simdata <- readRDS(paste0("sim100/sim",set,".rds")) ## Read simulated dataset for large cities
  
  M.s.maybe.inits <- round(simdata$M.maybe*0.3)
  
  initial.values <- list(list("p.c" = 0.5,
                              "p.c.i" = 0.5,
                              "p.ni.maybe" = 0.5,
                              "M.s.maybe" = M.s.maybe.inits),
                         list("p.c" = 0.6,
                              "p.c.i" = 0.3,
                              "p.ni.maybe" = 0.8,
                              "M.s.maybe" = M.s.maybe.inits),
                         list("p.c" = 0.4,
                              "p.c.i" = 0.8,
                              "p.ni.maybe" = 0.1,
                              "M.s.maybe" = M.s.maybe.inits))
  
  vars.monitor <- c("p.c", "p.ni.maybe", "p.c.i", "H1", "H2")
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
