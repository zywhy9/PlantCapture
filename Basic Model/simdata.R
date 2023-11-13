## Small City
set.seed(1234)
M <- 15
H <- 150
p.s <- 0.7
p.maybe <- 0.2


for(set in 1:1000){
  H.s <- rbinom(1, H, p.s)
  flag <- T
  while(flag){
    M.yes <- rbinom(1, M, p.s*(1-p.maybe))
    if(M.yes==0){next} ## Avoid zeros which cause error
    M.s.maybe <- rbinom(1, M - M.yes, p.s*p.maybe / (1-p.s*(1-p.maybe)))
    M.ns.maybe <- rbinom(1, M - M.yes - M.s.maybe, (1-p.s)*p.maybe/(1-p.s*(1-p.maybe)-p.s*p.maybe))
    M.maybe <- M.s.maybe + M.ns.maybe
    M.no <- M - M.yes - M.maybe
    if(M.no==0){next} ## Avoid zeros which cause error
    flag <- F
  }
  Y <- M.yes + H.s + M.s.maybe
  
  simdata <- list(
    "M.yes" = M.yes,
    "M.maybe" = M.maybe,
    "M.no" = M.no,
    "M" = M,
    "census.unkn" = Y-M.yes)
  saveRDS(simdata, file=paste0("sim15nz/sim",set,".rds"))
}

## Large City
set.seed(1234)
M <- 100
H <- 1500
p.s <- 0.7
p.maybe <- 0.2


for(set in 1:1000){
  H.s <- rbinom(1, H, p.s)
  flag <- T
  while(flag){
    M.yes <- rbinom(1, M, p.s*(1-p.maybe))
    if(M.yes==0){next} ## Avoid zeros which cause error
    M.s.maybe <- rbinom(1, M - M.yes, p.s*p.maybe / (1-p.s*(1-p.maybe)))
    M.ns.maybe <- rbinom(1, M - M.yes - M.s.maybe, (1-p.s)*p.maybe/(1-p.s*(1-p.maybe)-p.s*p.maybe))
    M.maybe <- M.s.maybe + M.ns.maybe
    M.no <- M - M.yes - M.maybe
    if(M.no==0){next} ## Avoid zeros which cause error
    flag <- F
  }
  Y <- M.yes + H.s + M.s.maybe
  
  simdata <- list(
    "M.yes" = M.yes,
    "M.maybe" = M.maybe,
    "M.no" = M.no,
    "M" = M,
    "census.unkn" = Y-M.yes)
  saveRDS(simdata, file=paste0("sim100/sim",set,".rds"))
}
