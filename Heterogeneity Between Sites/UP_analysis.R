library(rjags)
library(R2jags)

## Inference
ModelCode <- "model{
  #Priors
  p.c[1] ~ dunif(0,1)
  p.c[2] ~ dunif(0,1)
  p.ni.maybe ~ dunif(0,1) #probability a plant that says maybe they were seen and not interviewed
  p.c.i ~ dunif(0,1) #probability a plant was seen that he got interviewed
  
  #Tranformation
  for(i in 1:2){
    p.i[i] <- p.c[i] * p.c.i
    p.c.maybe[i] <- p.c[i] * (1 - p.c.i) * p.ni.maybe
    p.yes[i] <- p.c[i] * (1 - p.c.i) * (1 - p.ni.maybe)
    p.ns.maybe[i] <- (1 - p.c[i]) * p.ni.maybe
    p.maybe[i] <- p.c.maybe[i] + p.ns.maybe[i]
  }
  
  #Model
  for(i in 1:2){
    M.i[i] ~ dbin(p.i[i], M[i])
    M.yes[i] ~ dbin(p.yes[i] / (1 - p.i[i]), M[i] - M.i[i])
    M.maybe[i] ~ dbin(p.maybe[i] / (1 - p.yes[i] - p.i[i]), M[i] - M.i[i] - M.yes[i])
    M.s.maybe[i] ~ dbin(p.c.maybe[i] / p.maybe[i], M.maybe[i])
    H.s[i] <-  temp[i] - M.s.maybe[i]
    H.i[i] ~ dbin(p.c.i, H.s[i])
    H1[i] <- H.s[i] / p.c[i]
    H2[i] ~ dnorm(H1[i], 1/(H1[i]*(1-p.c[i])/p.c[i]))
  }
  H1.t <- sum(H1[1:2])
  H2.t <- sum(H2[1:2])
}
"

for(set in 1:1000){
  
  set.seed(set)
  
  simdata <- readRDS(paste0("sim30/sim",set,".rds"))
  
  M.s.maybe.inits <- round(simdata$M.maybe*0.3)
  H.s.inits <- simdata$temp - M.s.maybe.inits
  H.inits <- H.s.inits/0.7
  
  initial.values <- list(list("p.c" = c(0.8,0.2),
                              "p.c.i" = 0.5,
                              "p.ni.maybe" = 0.5,
                              "M.s.maybe" = M.s.maybe.inits),
                         list("p.c" = c(0.6, 0.4),
                              "p.c.i" = 0.3,
                              "p.ni.maybe" = 0.8,
                              "M.s.maybe" = M.s.maybe.inits),
                         list("p.c" = c(0.7, 0.3),
                              "p.c.i" = 0.8,
                              "p.ni.maybe" = 0.1,
                              "M.s.maybe" = M.s.maybe.inits))
  
  vars.monitor <- c("p.c", "p.ni.maybe", "p.c.i", "H1.t", "H2.t","H1","H2")
  data <- list("M" = simdata$M,
               "M.i" = simdata$M.i,
               "M.yes" = simdata$M.yes,
               "M.maybe" = simdata$M.maybe,
               "temp" = simdata$temp,
               "H.i" = simdata$H.i)
  print(set)
  jagsfit <- jags(data=data, n.chains=3, inits=initial.values,
                  parameters.to.save=vars.monitor, n.iter=30000, n.burnin=15000,n.thin=1,
                  DIC=TRUE, model.file=textConnection(ModelCode))
  
  saveRDS(jagsfit$BUGSoutput$sims.array,paste0("res30/res",set,".rds"))
}


