## Model codes
ModelCode <- "model{
  #Priors
  p.s[1] ~ dunif(0,1)
  p.s[2] ~ dunif(0,1)
  p.ni.maybe ~ dunif(0,1) #probability a plant that says maybe they were seen and not interviewed
  p.s.i ~ dunif(0,1) #probability a plant was seen that he got interviewed
  H.prior[1] ~ dnorm(1000, 1/10000^2) T(0,)
  H.prior[2] ~ dnorm(1000, 1/10000^2) T(0,)
  H <- round(H.prior)
  
  #Tranformation
  for(i in 1:2){
    p.i[i] <- p.s[i] * p.s.i
    p.s.maybe[i] <- p.s[i] * (1 - p.s.i) * p.ni.maybe
    p.yes[i] <- p.s[i] * (1 - p.s.i) * (1 - p.ni.maybe)
    p.ns.maybe[i] <- (1 - p.s[i]) * p.ni.maybe
    p.maybe[i] <- p.s.maybe[i] + p.ns.maybe[i]
  }
  
  #Model
  for(i in 1:2){
    M.i[i] ~ dbin(p.i[i], M[i])
    M.yes[i] ~ dbin(p.yes[i] / (1 - p.i[i]), M[i] - M.i[i])
    M.maybe[i] ~ dbin(p.maybe[i] / (1 - p.yes[i] - p.i[i]), M[i] - M.i[i] - M.yes[i])
    M.s.maybe[i] ~ dbin(p.s.maybe[i] / p.maybe[i], M.maybe[i])
    H.s[i] ~ dbin(p.s[i], H[i])
    H.i[i] ~ dbin(p.s.i, H.s[i])
    temp[i] ~ dsum(M.s.maybe[i], H.s[i])
  }
}
"
save(ModelCode,file="simmodel.RData")
