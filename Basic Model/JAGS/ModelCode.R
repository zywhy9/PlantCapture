## Model codes
ModelCode <- "model{
  #Priors
  p.s ~ dunif(0,1) #probability to be seen
  p.maybe ~ dunif(0,1) #probability a plant that says maybe they were seen
  H ~ dnorm(1000, 1/10000^2) T(0,)
  #Model
  M.maybe ~ dbin(p.maybe, M)
  M.yes ~ dbin(p.s, M-M.maybe)
  M.no <- M - M.maybe
  census.unkn ~ dbin(p.s,total)
  total <- round(H) + M.maybe #number of real homeless
}
"
save(ModelCode,file="simmodel.RData")
