M.i <- c(2, 41, 18, 40, 16)
M.yes <- c(0, 6, 3, 19, 1)
M.maybe <- c(5, 5, 1, 13, 2)
M.no <- c(6, 6, 4, 22, 6)
Y <- c(11, 109, 104, 1240, 217)


for(i in 1:5){
  res <- optim(par=c(log(2000),boot::logit(0.5),boot::logit(0.5),boot::logit(0.5)), fn=ld_complete,
               M.i=M.i[i],M.yes=M.yes[i],M.maybe=M.maybe[i],M.no=M.no[i],Y=Y[i],
               control = list(fnscale=-1), hessian = T)
  meanmt[i,] <- c(exp(res$par[1]),boot::inv.logit(res$par[2:4]))
  hmt[i,] <- sqrt(diag(solve(-numDeriv::hessian(ld_complete,M.i=M.i[i],M.yes=M.yes[i],M.maybe=M.maybe[i],M.no=M.no[i],Y=Y[i],x=res$par))))
  fsdmt[i,] <- sqrt(c(var_log(res$par[1],hmt[i,1]^2), var_logit(e=res$par[2:4],v=hmt[i,2:4]^2)))
  cilmt[i,] <- c(exp(res$par[1] - 1.96 * hmt[i,1]),boot::inv.logit(res$par[2:4] - 1.96 * hmt[i,2:4]))
  ciumt[i,] <- c(exp(res$par[1] + 1.96 * hmt[i,1]),boot::inv.logit(res$par[2:4] + 1.96 * hmt[i,2:4]))
}
