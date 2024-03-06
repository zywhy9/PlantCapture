## Calculate Likelihood
ld_complete <- function(M.i,M.yes,M.maybe,M.no,Y,H.i,param){
  ## Target
  H <- exp(param[1])
  p.c <- boot::inv.logit(param[2])
  p.c.i <- boot::inv.logit(param[3])
  p.ni.maybe <- boot::inv.logit(param[4])
  ## Transformation
  p.i <- p.c * p.c.i
  p.yes <- p.c * (1 - p.c.i) * (1 - p.ni.maybe)
  p.c.maybe <- p.c * (1 - p.c.i) * p.ni.maybe
  p.ns.maybe <- (1 - p.c) * p.ni.maybe
  p.no <- (1 - p.c) * (1 - p.ni.maybe)
  ## Density
  res <- c()
  uppb <- min(M.maybe,Y-M.i-M.yes)
  lowb <- min(max(0,round(Y-M.i-M.yes-H)),uppb)
  for(M.maybe.s in lowb:uppb){
    M.maybe.ns <- M.maybe - M.maybe.s
    H.s <- Y - M.i - M.yes - M.maybe.s
    M.comp <- c(M.i,M.yes,M.maybe.s,M.maybe.ns,M.no)
    ld1 <- dmultinom(M.comp, prob=c(p.i,p.yes,p.c.maybe,p.ns.maybe,p.no), log=T) ## Multinomial(M,theta)
    ld2 <- dnorm(H.s, H*p.c, sqrt(H*p.c*(1-p.c)), log=T) ## Binomial(H,p.c)
    ld3 <- dbinom(H.i, H.s, p.c.i, log=T) ## Binomial(H.s,p.c.i)
    res <- c(res, ld1+ld2+ld3)
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
sdmt <- matrix(NA, nrow = nset, ncol = npar)
cilmt <- matrix(NA, nrow = nset, ncol = npar)
ciumt <- matrix(NA, nrow = nset, ncol = npar)
start_time <- Sys.time()
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = nset, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar

for(i in 1:nset){
  data <- readRDS(paste0("sim15/sim",i,".rds")) ## Small city example
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
  temp <- diag(solve(-numDeriv::hessian(ld_complete,M.i=M.i,M.yes=M.yes,M.maybe=M.maybe,M.no=M.no,Y=Y,H.i=H.i,x=res$par),tol=1e-50))
  hmt[i,] <- ifelse(temp<0,0,sqrt(temp)) ## Standard deviation on transformed scale
  sdmt[i,] <- sqrt(c(var_log(res$par[1],hmt[i,1]^2), var_logit(e=res$par[2:4],v=hmt[i,2:4]^2))) ## Standard deviation on original scale
  cilmt[i,] <- c(exp(res$par[1] - 1.96 * hmt[i,1]),boot::inv.logit(res$par[2:4] - 1.96 * hmt[i,2:4])) ## Lower bound of 95% CrI
  ciumt[i,] <- c(exp(res$par[1] + 1.96 * hmt[i,1]),boot::inv.logit(res$par[2:4] + 1.96 * hmt[i,2:4])) ## Upper bound of 95% CrI
  setTxtProgressBar(pb, i)
}
close(pb)
end_time <- Sys.time()
(end_time - start_time)/nset

true.value <- c(150,0.7,0.8,0.2) ## True value for small cities

pmn <- apply(meanmt,2,mean) ## Expected posterior mean
psd <- apply(meanmt,2,sd) ## Standard error of posterior mean
epsd <- apply(sdmt,2,mean) ## Expected posterior standard deviation

bias <- pmn - true.value ## Bias
rbias <- bias / true.value ## Relative bias

mse <- rep(0,npar)
for(i in 1:npar){
  mse[i] <- mean((meanmt[,i]-true.value[i])^2) ## Mean squared error
}
rmse <- sqrt(mse) ## Root mean squared error
rrmse <- rmse/true.value ## Relative RMSE

cp <- rep(0,npar)
for(set in 1:nset){
  for(i in 1:npar){
    if(true.value[i]<=ciumt[set,i] && true.value[i]>=cilmt[set,i]){
      cp[i] <- cp[i] + 1
    }
  }
}
cp <- cp/nset ## Coverage probability
lci <- apply(ciumt-cilmt,2,mean) ## Average length of credible interval
list(mean=round(pmn,4),epsd=round(efpsd,4),se=round(psd,4),rbias=round(rbias,4),rrmse=round(rrmse,4),cp=cp,lci=round(lci,4))

