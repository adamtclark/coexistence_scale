pout<-log(c(0.1, 0.3, 0.4, 0.02, 0.02, 0.02))

levmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dN1<-State[1]*Pars[1]*(1-State[1])-Pars[4]*State[1]
    dN2<-State[2]*Pars[2]*(1-State[1]-State[2])-Pars[5]*State[2]-
      State[2]*State[1]*Pars[1]
    dN3<-State[3]*Pars[3]*(1-State[1]-State[2]-State[3])-Pars[6]*State[3]-
      State[3]*State[1]*Pars[1]-State[3]*State[2]*Pars[2]
    
    return(list(c(dN1, dN2, dN3)))
  })
}

optfun<-function(pout) {
  clst_meta = exp(c(pout[1], pout[2], pout[3]))
  mlst_meta = exp(c(pout[4], pout[5], pout[6]))#rep(0.02, 3)
  
  #simtime<-79
  #population<-populate(gridout, nlst = s0_ind,
  #                     clst = clst_meta, radlst = Inf, mlst = mlst_meta)
  #grid_sub<-grid_subset_BLOCKSIZE(gridout, size = 1)
  
  #out<-run_metapopulation(tmax=simtime, nsteps = simtime, gridout, population, talktime = 0, sites_sub = grid_sub$sites)
  #tmp<-out$output[,c(4:2)]/prod(gridout$lng)
  
  tmp<-ode(s0_ind/100^2, timebands_small, levmod, c(clst_meta, mlst_meta))
  
  sum((tmp[,-1]-abunds_obs[,c(3,2,1)])^2, na.rm=T)
}

optout<-optim(pout, optfun, control=list(trace=4))



pout<-optout$par
#[1]  -2.1361578   0.5237058  -2.1162532
#[4]  -2.5475349   0.3291012 -10.2173169

matplot(timebands_small, tmp[,-1], lty=1, type="l", col=2:4, lwd=2)
matlines(timebands_small, abunds_obs[,c(3,2,1)], lty=2, lwd=2, col=2:4)
