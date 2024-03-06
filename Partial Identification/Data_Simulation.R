set.seed(1234)

p.c <- 0.7
p.ni.maybe <- 0.2
p.c.i <- 0.8

p.i <- p.c * p.c.i
p.yes <- p.c * (1 - p.c.i) * (1 - p.ni.maybe)
p.c.maybe <- p.c * (1 - p.c.i) * p.ni.maybe
p.ns.maybe <- (1 - p.c) * p.ni.maybe
p.no <- (1 - p.c) * (1 - p.ni.maybe)

## Small city
H <- 1500
M <- 100

for(set in 1:1000){
  H.s <- rbinom(1, H, p.c)
  H.i <- rbinom(1, H.s, p.c.i)
  flag <- T
  while(flag){
    M.i <- rbinom(1, M, p.i)
    M.yes <- rbinom(1, M - M.i, p.yes/(1 - p.i))
    if(M.i==0 & M.yes==0){next}
    M.s.maybe <- rbinom(1, M - M.i - M.yes, p.c.maybe / (1 - p.i - p.yes))
    M.ns.maybe <- rbinom(1, M - M.i - M.yes - M.s.maybe, p.ns.maybe / (1 - p.i - p.yes - p.c.maybe))
    M.maybe <- M.s.maybe + M.ns.maybe
    M.no <- M - M.i - M.yes - M.maybe
    if(M.no==0){next}
    flag <- F
  }
  census.unkn <- H.s + M.s.maybe
  
  simdata <- list(
    "M.i" = M.i,
    "M.yes" = M.yes,
    "M.maybe" = M.maybe,
    "M" = M,
    "census.unkn" = census.unkn,
    "H.i" = H.i)
  saveRDS(simdata, file=paste0("sim100/sim",set,".rds"))
}

## Small city
H <- 150
M <- 15

for(set in 1:1000){
  H.s <- rbinom(1, H, p.c)
  H.i <- rbinom(1, H.s, p.c.i)
  flag <- T
  while(flag){
    M.i <- rbinom(1, M, p.i)
    M.yes <- rbinom(1, M - M.i, p.yes/(1 - p.i))
    if(M.i==0 & M.yes==0){next}
    M.s.maybe <- rbinom(1, M - M.i - M.yes, p.c.maybe / (1 - p.i - p.yes))
    M.ns.maybe <- rbinom(1, M - M.i - M.yes - M.s.maybe, p.ns.maybe / (1 - p.i - p.yes - p.c.maybe))
    M.maybe <- M.s.maybe + M.ns.maybe
    M.no <- M - M.i - M.yes - M.maybe
    if(M.no==0){next}
    flag <- F
  }
  census.unkn <- H.s + M.s.maybe
  
  simdata <- list(
    "M.i" = M.i,
    "M.yes" = M.yes,
    "M.maybe" = M.maybe,
    "M" = M,
    "census.unkn" = census.unkn,
    "H.i" = H.i)
  saveRDS(simdata, file=paste0("sim15/sim",set,".rds"))
}
