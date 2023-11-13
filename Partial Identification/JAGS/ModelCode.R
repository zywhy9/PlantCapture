## Model codes
ModelCode <- "model{
  #Priors
  p.s ~ dunif(0,1) #probability to be seen
  p.ni.maybe ~ dunif(0,1) #probability a plant that says maybe they were seen and not interviewed
  p.s.i ~ dunif(0,1) #probability a plant was seen that he got interviewed
  H.prior ~ dnorm(1000, 1/10000^2) T(0,)
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
save(ModelCode,file="simmodel.RData")
