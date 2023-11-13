## Extract Data from Literature
cityname <- c("Chicago", "New Orleans", "Phoenix", "New York", "Los Angeles")
M <- c(13, 58, 26, 94, 25) ### Number of plants in the sites where enumerators were seen for each city from page 18 of Martin et al. (1997)
Y <- c(11, 109, 104, 1240, 217) ### Census count in the sites where enumerators were seen for each city from page 20 of Martin et al. (1997)
p.interview <- c(0.07, 0.67, 0.44, 0.37, 0.33) ### Percent of interviewed plants for each city from page 15 of Martin et al. (1997)
p.yes <- c(0, 0.1, 0.08, 0.17, 0.02) ### Percent of plants answering yes for each city from page 15 of Martin et al. (1997)
p.maybe <- c(0.18, 0.07, 0.03, 0.12, 0.04) ### Percent of plants answering maybe for each city from page 15 of Martin et al. (1997)
p.no <- c(0.25, 0.10, 0.10, 0.20, 0.13) ### Percent of plants answering no for each city from page 15 of Martin et al. (1997)
p.invalid <- c(0.5, 0.05, 0.36, 0.14, 0.48) ### Percent of plants did not see enumerators for each city from page 15 of Martin et al. (1997)

## Final Data
M.i <- round(M * p.interview / (1 - p.invalid))
M.yes <- round(M * p.yes / (1 - p.invalid))
M.no <- round(M * p.no / (1 - p.invalid))
M.maybe <- M - M.i - M.yes - M.no
census.unkn <- Y - M.yes - M.i

# Calculate Likelihood
ld_complete <- function(M.i,M.yes,M.maybe,M.no,Y,param){
  ## Target
  H <- exp(param[1])
  p.s <- boot::inv.logit(param[2])
  p.ni.maybe <- boot::inv.logit(param[3])
  p.s.i <- boot::inv.logit(param[4])
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
    logitpar <- param[2:4]
    res <- c(res, ld1+ld2+sum(-logitpar-2*log(1+exp(-logitpar))))
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
nset <- 5
npar <- 4
pmean <- matrix(NA, nrow = nset, ncol = npar)
hmt <- matrix(NA, nrow = nset, ncol = npar)
psd <- matrix(NA, nrow = nset, ncol = npar)
ci.l <- matrix(NA, nrow = nset, ncol = npar)
ci.u <- matrix(NA, nrow = nset, ncol = npar)

for(i in 1:nset){
  res <- optim(par=c(log(200),boot::logit(0.5),boot::logit(0.5),boot::logit(0.5)), fn=ld_complete,
               M.i=M.i[i],M.yes=M.yes[i],M.maybe=M.maybe[i],M.no=M.no[i],Y=Y[i],
               control = list(fnscale=-1), hessian = T)
  pmean[i,] <- c(exp(res$par[1]),boot::inv.logit(res$par[2:4]))
  hmt[i,] <- sqrt(diag(solve(-numDeriv::hessian(ld_complete,M.i=M.i[i],M.yes=M.yes[i],M.maybe=M.maybe[i],M.no=M.no[i],Y=Y[i],x=res$par))))
  psd[i,] <- sqrt(c(var_log(res$par[1],hmt[i,1]^2), var_logit(e=res$par[2:4],v=hmt[i,2:4]^2)))
  ci.l[i,] <- c(exp(res$par[1] - 1.96 * hmt[i,1]),boot::inv.logit(res$par[2:4] - 1.96 * hmt[i,2:4]))
  ci.u[i,] <- c(exp(res$par[1] + 1.96 * hmt[i,1]),boot::inv.logit(res$par[2:4] + 1.96 * hmt[i,2:4]))
}

eva <- cbind(as.vector(t(pmean)),as.vector(t(psd)),as.vector(t(ci.l)),as.vector(t(ci.u)))
colnames(eva) <- c("Median", "SD", "Lower", "Upper")
rownames(eva) <- rep(c("H", "p.s", "p.ni.maybe", "p.s.i"),5)
round(eva,2)
