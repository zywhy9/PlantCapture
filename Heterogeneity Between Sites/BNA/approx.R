## Calculate Likelihood
ld_complete <- function(M.i,M.yes,M.maybe,M.no,Y,H.i,param){
  ## Target
  H.e <- exp(param[1])
  H.h <- exp(param[2])
  p.s.e <- boot::inv.logit(param[3])
  p.s.h <- boot::inv.logit(param[4])
  p.s.i <- boot::inv.logit(param[5])
  p.ni.maybe <- boot::inv.logit(param[6])
  ## Transformation
  p.s <- c(p.s.e, p.s.h)
  H <- c(H.e, H.h)
  
  p.i <- p.yes <- p.s.maybe <- p.ns.maybe <- p.no <- H.s <- rep(NA, 2)
  for(i in 1:2){
    p.i[i] <- p.s[i] * p.s.i
    p.yes[i] <- p.s[i] * (1 - p.s.i) * (1 - p.ni.maybe)
    p.s.maybe[i] <- p.s[i] * (1 - p.s.i) * p.ni.maybe
    p.ns.maybe[i] <- (1 - p.s[i]) * p.ni.maybe
    p.no[i] <- (1 - p.s[i]) * (1 - p.ni.maybe)
  }
  
  ## Density
  res <- c()
  uppb.e <- min(M.maybe[1], Y[1] - M.i[1] - M.yes[1] - H.i[1])
  lowb.e <- min(max(0, round(Y[1] - M.i[1] - M.yes[1] - H[1])), uppb.e)
  for(M.maybe.s.e in lowb.e:uppb.e){
    uppb.h <- min(M.maybe[2], Y[2] - M.i[2] - M.yes[2] - H.i[2])
    lowb.h <- min(max(0, round(Y[2] - M.i[2] - M.yes[2] - H[2])), uppb.h)
    for(M.maybe.s.h in lowb.h:uppb.h){
      # print(c(M.maybe.s.e,H[1]))
      M.maybe.s <- c(M.maybe.s.e, M.maybe.s.h)
      M.maybe.ns <- M.maybe - M.maybe.s
      M.comp <- cbind(M.i,M.yes,M.maybe.s,M.maybe.ns,M.no)
      for(i in 1:2){
        H.s[i] <- Y[i] - M.i[i] - M.yes[i] - M.maybe.s[i]
      }
      ld <- matrix(0,nrow=2,ncol=3)
      for(i in 1:2){
        ld[i,1] <- dmultinom(M.comp[i,], prob=c(p.i[i],p.yes[i],p.s.maybe[i],p.ns.maybe[i],p.no[i]), log=T) ## Multinomial(M,theta)
        ld[i,2] <- dnorm(H.s[i], H[i]*p.s[i], sqrt(H[i]*p.s[i]*(1-p.s[i])), log=T) ## Binomial(H,p.s)
        ld[i,3] <- dbinom(H.i[i], H.s[i], p.s.i, log=T) ## Binomial(H.s,p.s.i)
      }
      logitpar <- param[3:6]
      res <- c(res, sum(ld)+sum(-logitpar-2*log(1+exp(-logitpar)))+sum(param[1:2])) ## data likelihood + flat prior for logit + flat prior for log
    }
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
npar <- 6
# meanmt <- matrix(NA, nrow = nset, ncol = (npar-1))
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
  data <- readRDS(paste0("../sim100/sim",i,".rds"))
  M.i <- data$M.i
  M.yes <- data$M.yes
  M.maybe <- data$M.maybe
  M <- data$M
  H.i <- data$H.i
  census.unkn <- data$census.unkn
  M.no <- M - M.i - M.yes - M.maybe
  Y <- census.unkn + M.i + M.yes + H.i
  res <- optim(par=c(log(1000),log(1000),boot::logit(0.8),boot::logit(0.3),boot::logit(0.9),boot::logit(0.1)), fn=ld_complete,
              M.i=M.i,M.yes=M.yes,M.maybe=M.maybe,M.no=M.no,Y=Y,H.i=H.i,
              control = list(fnscale=-1,maxit=10000), hessian = T)
  # meanmt[i,] <- c(sum(exp(res$par[1:2])),boot::inv.logit(res$par[3:npar]))
  meanmt[i,] <- c(exp(res$par[1:2]),boot::inv.logit(res$par[3:npar]))
  # hmt[i,] <- sqrt(diag(solve(-numDeriv::hessian(ld_complete,M.i=M.i,M.yes=M.yes,M.maybe=M.maybe,M.no=M.no,Y=Y,H.i=H.i,x=res$par))))
  temp <- diag(solve(-numDeriv::hessian(ld_complete,M.i=M.i,M.yes=M.yes,M.maybe=M.maybe,M.no=M.no,Y=Y,H.i=H.i,x=res$par),tol=1e-50))
  hmt[i,] <- ifelse(is.na(sqrt(temp)),0,sqrt(temp))
  fsdmt[i,] <- sqrt(c(var_log(res$par[1:2],hmt[i,1:2]^2), var_logit(e=res$par[3:npar],v=hmt[i,3:npar]^2)))
  cilmt[i,] <- c(exp(res$par[1:2] - 1.96 * hmt[i,1:2]),boot::inv.logit(res$par[3:npar] - 1.96 * hmt[i,3:npar]))
  ciumt[i,] <- c(exp(res$par[1:2] + 1.96 * hmt[i,1:2]),boot::inv.logit(res$par[3:npar] + 1.96 * hmt[i,3:npar]))
  setTxtProgressBar(pb, i)
}
close(pb)
end_time <- Sys.time()
(end_time - start_time)/nset

# saveRDS(list(meanmt=meanmt,hmt=hmt,fsdmt=fsdmt,cilmt=cilmt,ciumt=ciumt),"res100.rds")

resmeanmt <- cbind(meanmt[,1]+meanmt[,2],meanmt[,3:npar])
resfsdmt <- cbind(sqrt(fsdmt[,1]^2+fsdmt[,2]^2),fsdmt[,3:npar])
rescilmt <- cbind(cilmt[,1]+cilmt[,2],cilmt[,3:npar])
resciumt <- cbind(ciumt[,1]+ciumt[,2],ciumt[,3:npar])

true.value <- c(1500,0.9,0.4,0.8,0.2)

pmn <- apply(resmeanmt,2,mean)
psd <- apply(resmeanmt,2,sd)
efpsd <- apply(resfsdmt,2,mean)

bias <- pmn - true.value
rbias <- bias / true.value

mse <- rep(0,npar-1)
for(i in 1:(npar-1)){
  mse[i] <- mean((resmeanmt[,i]-true.value[i])^2)
}
rmse <- sqrt(mse)
rrmse <- rmse/true.value

cp <- rep(0,npar-1)
for(set in 1:nset){
  for(i in 1:(npar-1)){
    if(true.value[i]<=resciumt[set,i] && true.value[i]>=rescilmt[set,i]){
      cp[i] <- cp[i] + 1
    }
  }
}
cp <- cp/nset
lcp <- apply(resciumt-rescilmt,2,mean)
list(mean=round(pmn,4),epsd=round(efpsd,4),sd=round(psd,4),rbias=round(rbias,4),rrmse=round(rrmse,4),cp=cp,lcp=round(lcp,4))
