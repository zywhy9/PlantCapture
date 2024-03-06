set.seed(1234)

w.e <- 0.6
p.c.e <- 0.9
p.c.h <- 0.4
p.c <- c(p.c.e, p.c.h)
p.ni.maybe <- 0.2
p.c.i <- 0.8

p.i <- p.c * p.c.i
p.yes <- p.c * (1 - p.c.i) * (1 - p.ni.maybe)
p.c.maybe <- p.c * (1 - p.c.i) * p.ni.maybe
p.ns.maybe <- (1 - p.c) * p.ni.maybe
p.no <- (1 - p.c) * (1 - p.ni.maybe)


## Small City
H <- 300
M.t <- 30

for(set in 1:1000){
  H.s <- c(rbinom(1, H * w.e, p.c[1]), rbinom(1, H * (1 - w.e), p.c[2]))
  H.i <- c(rbinom(1, H.s[1], p.c.i), rbinom(1, H.s[2], p.c.i))
  M.e <- M.t * w.e
  M.h <- M.t - M.e
  M <- c(M.e, M.h)
  
  flag <- T
  while(flag){
    M.i <- c(rbinom(1, M[1], p.i[1]), rbinom(1, M[2], p.i[2]))
    M.yes <- c(rbinom(1, M[1] - M.i[1], p.yes[1]/(1 - p.i[1])), rbinom(1, M[2] - M.i[2], p.yes[2]/(1 - p.i[2])))
    if(M.i[1]==0 & M.yes[1]==0){next}
    if(M.i[2]==0 & M.yes[2]==0){next}
    M.s.maybe <- c(rbinom(1, M[1] - M.i[1] - M.yes[1], p.c.maybe[1] / (1 - p.i[1] - p.yes[1])), 
                   rbinom(1, M[2] - M.i[2] - M.yes[2], p.c.maybe[2] / (1 - p.i[2] - p.yes[2])))
    M.ns.maybe <- c(rbinom(1, M[1] - M.i[1] - M.yes[1] - M.s.maybe[1], p.ns.maybe[1] / (1 - p.i[1] - p.yes[1] - p.c.maybe[1])),
                    rbinom(1, M[2] - M.i[2] - M.yes[2] - M.s.maybe[2], p.ns.maybe[2] / (1 - p.i[2] - p.yes[2] - p.c.maybe[2])))
    M.maybe <- M.s.maybe + M.ns.maybe
    M.no <- M - M.i - M.yes - M.maybe
    if(0 %in% M.no){next}
    flag <- F
  }
  census.unkn <- H.s - H.i + M.s.maybe
  
  simdata <- list(
    "M.i" = M.i,
    "M.yes" = M.yes,
    "M.maybe" = M.maybe,
    "M" = M,
    "census.unkn" = census.unkn,
    "temp" = H.s + M.s.maybe,
    "H.i" = H.i)
  saveRDS(simdata, file=paste0("sim30/sim",set,".rds"))
}

## Small City
H <- 1500
M.t <- 100

for(set in 1:1000){
  H.s <- c(rbinom(1, H * w.e, p.c[1]), rbinom(1, H * (1 - w.e), p.c[2]))
  H.i <- c(rbinom(1, H.s[1], p.c.i), rbinom(1, H.s[2], p.c.i))
  M.e <- M.t * w.e
  M.h <- M.t - M.e
  M <- c(M.e, M.h)
  
  flag <- T
  while(flag){
    M.i <- c(rbinom(1, M[1], p.i[1]), rbinom(1, M[2], p.i[2]))
    M.yes <- c(rbinom(1, M[1] - M.i[1], p.yes[1]/(1 - p.i[1])), rbinom(1, M[2] - M.i[2], p.yes[2]/(1 - p.i[2])))
    if(M.i[1]==0 & M.yes[1]==0){next}
    if(M.i[2]==0 & M.yes[2]==0){next}
    M.s.maybe <- c(rbinom(1, M[1] - M.i[1] - M.yes[1], p.c.maybe[1] / (1 - p.i[1] - p.yes[1])), 
                   rbinom(1, M[2] - M.i[2] - M.yes[2], p.c.maybe[2] / (1 - p.i[2] - p.yes[2])))
    M.ns.maybe <- c(rbinom(1, M[1] - M.i[1] - M.yes[1] - M.s.maybe[1], p.ns.maybe[1] / (1 - p.i[1] - p.yes[1] - p.c.maybe[1])),
                    rbinom(1, M[2] - M.i[2] - M.yes[2] - M.s.maybe[2], p.ns.maybe[2] / (1 - p.i[2] - p.yes[2] - p.c.maybe[2])))
    M.maybe <- M.s.maybe + M.ns.maybe
    M.no <- M - M.i - M.yes - M.maybe
    if(0 %in% M.no){next}
    flag <- F
  }
  census.unkn <- H.s - H.i + M.s.maybe
  
  simdata <- list(
    "M.i" = M.i,
    "M.yes" = M.yes,
    "M.maybe" = M.maybe,
    "M" = M,
    "census.unkn" = census.unkn,
    "temp" = H.s + M.s.maybe,
    "H.i" = H.i)
  saveRDS(simdata, file=paste0("sim100/sim",set,".rds"))
}
