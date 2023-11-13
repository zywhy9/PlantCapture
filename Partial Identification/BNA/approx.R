## Calculate Likelihood
ld_complete <- function(M.i,M.yes,M.maybe,M.no,Y,H.i,param){
  ## Target
  H <- exp(param[1])
  p.s <- boot::inv.logit(param[2])
  p.s.i <- boot::inv.logit(param[3])
  p.ni.maybe <- boot::inv.logit(param[4])
  ## Transformation
  p.i <- p.s * p.s.i
  p.yes <- p.s * (1 - p.s.i) * (1 - p.ni.maybe)
  p.s.maybe <- p.s * (1 - p.s.i) * p.ni.maybe
  p.ns.maybe <- (1 - p.s) * p.ni.maybe
  p.no <- (1 - p.s) * (1 - p.ni.maybe)
  ## Density
  res <- c()
  uppb <- min(M.maybe,Y-M.i-M.yes)
  lowb <- min(max(0,round(Y-M.i-M.yes-H)),uppb)
  for(M.maybe.s in lowb:uppb){
    M.maybe.ns <- M.maybe - M.maybe.s
    H.s <- Y - M.i - M.yes - M.maybe.s
    M.comp <- c(M.i,M.yes,M.maybe.s,M.maybe.ns,M.no)
    ld1 <- dmultinom(M.comp, prob=c(p.i,p.yes,p.s.maybe,p.ns.maybe,p.no), log=T) ## Multinomial(M,theta)
    ld2 <- dnorm(H.s, H*p.s, sqrt(H*p.s*(1-p.s)), log=T) ## Binomial(H,p.s)
    ld3 <- dbinom(H.i, H.s, p.s.i, log=T) ## Binomial(H.s,p.s.i)
    logitpar <- param[2:4]
    res <- c(res, ld1+ld2+ld3+sum(-logitpar-2*log(1+exp(-logitpar))))
  }
  res <- matrixStats::logSumExp(res)
  if(res==-Inf){
    return(-1e10) ## Avoid error with Inf
  }else{
    return(res)
  }
}

## Delta Method for log
var_log <- function(e,v){
  ex <- exp(e)
  varx <- v*ex^2
  return(varx)
}

## Delta Method for logit
var_logit <- function(e,v){
  ds <- boot::inv.logit(e) * (1 - boot::inv.logit(e))
  varx <- v * ds^2
  return(varx)
}

set.seed(1234)
nset <- 1000
npar <- 4
meanmt <- matrix(NA, nrow = nset, ncol = npar)
hmt <- matrix(NA, nrow = nset, ncol = npar)
fsdmt <- matrix(NA, nrow = nset, ncol = npar)
cilmt <- matrix(NA, nrow = nset, ncol = npar)
ciumt <- matrix(NA, nrow = nset, ncol = npar)
start_time <- Sys.time()
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = nset, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar

for(i in 1:nset){
  data <- readRDS(paste0("../sim15/sim",i,".rds"))
  M.i <- data$M.i
  M.yes <- data$M.yes
  M.maybe <- data$M.maybe
  M <- data$M
  H.i <- data$H.i
  census.unkn <- data$census.unkn
  M.no <- M - M.i - M.yes - M.maybe
  Y <- census.unkn + M.i + M.yes
  res <- optim(par=c(log(200),boot::logit(0.5),boot::logit(0.5),boot::logit(0.5)), fn=ld_complete,
              M.i=M.i,M.yes=M.yes,M.maybe=M.maybe,M.no=M.no,Y=Y,H.i=H.i,
              control = list(fnscale=-1), hessian = T)
  meanmt[i,] <- c(exp(res$par[1]),boot::inv.logit(res$par[2:4]))
  hmt[i,] <- sqrt(diag(solve(-numDeriv::hessian(ld_complete,M.i=M.i,M.yes=M.yes,M.maybe=M.maybe,M.no=M.no,Y=Y,H.i=H.i,x=res$par))))
  fsdmt[i,] <- sqrt(c(var_log(res$par[1],hmt[i,1]^2), var_logit(e=res$par[2:4],v=hmt[i,2:4]^2)))
  cilmt[i,] <- c(exp(res$par[1] - 1.96 * hmt[i,1]),boot::inv.logit(res$par[2:4] - 1.96 * hmt[i,2:4]))
  ciumt[i,] <- c(exp(res$par[1] + 1.96 * hmt[i,1]),boot::inv.logit(res$par[2:4] + 1.96 * hmt[i,2:4]))
  setTxtProgressBar(pb, i)
}
close(pb)
end_time <- Sys.time()
(end_time - start_time)/nset

true.value <- c(150,0.7,0.8,0.2)

pmn <- apply(meanmt,2,mean)
psd <- apply(meanmt,2,sd)
efpsd <- apply(fsdmt,2,mean)

bias <- pmn - true.value
rbias <- bias / true.value

mse <- rep(0,npar)
for(i in 1:npar){
  mse[i] <- mean((meanmt[,i]-true.value[i])^2)
}
rmse <- sqrt(mse)
rrmse <- rmse/true.value

cp <- rep(0,npar)
for(set in 1:nset){
  for(i in 1:npar){
    if(true.value[i]<=ciumt[set,i] && true.value[i]>=cilmt[set,i]){
      cp[i] <- cp[i] + 1
    }
  }
}
cp <- cp/nset
lcp <- apply(ciumt-cilmt,2,mean)
list(mean=round(pmn,4),epsd=round(efpsd,4),sd=round(psd,4),rbias=round(rbias,4),rrmse=round(rrmse,4),cp=cp,lcp=round(lcp,4))
