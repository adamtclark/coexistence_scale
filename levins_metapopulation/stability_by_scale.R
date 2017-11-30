error
rm(list=ls())
setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")


#TODO:
#1. Set up wrapper functions for stability
#2. Apply to both models
#3. Set up wrapper functions for spatial subsetting
#4. Think about PSF model to use

#load functions
require(rEDM)
source("run_metapopulation_wrapper.R")

#par(mfcol=c(4,2), mar=c(4,4,2,2))
set.seed(1152)

##### Try Tilman metapopulation model
gridout<-makegrid(xlng = 20, ylng = 20)
clst = c(0.15, 0.3, 0.8)
#getceq(clst)
population<-populate(gridout, nlst = floor(getceq(clst)*prod(gridout$lng)), clst = clst, radlst = Inf)


#tmax=1000; nsteps = 1000; talktime = 0; runtype="metapopulation"
out_meta<-run_metapopulation(tmax=1000, nsteps = 1000, gridout, population, talktime = 0)
#plot_metapop(out_meta)
#plot_map(out_meta, gridout)



#Estimate Eigenvalue
simtime<-100

out_meta_1<-NULL
for(i in 1:length(out_meta$plotdata$ceq)) {
  pt<-rep(0, length(out_meta$plotdata$ceq))
  pt[i]<-0.1
  out_meta_1[[i]]<-rerunrun_metapopulation(out=out_meta, tmax=simtime, talktime = 0, runtype = "metapopulation", perturb=pt, perturbsites=1:out_meta$plotdata$ngrid)
}
  
  
#out_meta_1[[1]]<-rerunrun_metapopulation(out=out_meta, tmax=simtime, talktime = 0, runtype = "metapopulation", perturb=c(0.1, 0, 0), perturbsites=1:out_meta$plotdata$ngrid)
#out_meta_1[[2]]<-rerunrun_metapopulation(out=out_meta, tmax=simtime, talktime = 0, runtype = "metapopulation", perturb=c(0, 0.1, 0), perturbsites=1:out_meta$plotdata$ngrid)
#out_meta_1[[3]]<-rerunrun_metapopulation(out=out_meta, tmax=simtime, talktime = 0, runtype = "metapopulation", perturb=c(0, 0, 0.1), perturbsites=1:out_meta$plotdata$ngrid)

Euse<-6
#sppos<-1
eigenlst<-matrix(ncol=length(out_meta_1), nrow=(simtime-1))

for(sppos in 1:length(out_meta_1)) {
  tmp<-simplex(time_series = c(out_meta$output[,sppos+1], out_meta_1[[sppos]]$output[,sppos+1]),
               E=Euse,
               lib=c(1,length(out_meta$output[,2])),
               pred=c(nrow(out_meta$output)-Euse+1, nrow(out_meta$output)+nrow(out_meta_1[[sppos]]$output)), stats_only = FALSE)
  
  pred_diff<-abs(tmp$model_output[[1]]$pred-tmp$model_output[[1]]$obs)/(out_meta$plotdata$ngrid)
  
  eigenlst[,sppos]<-cumsum(log(pred_diff[-1]/pred_diff[-length(pred_diff)]))/(1:(length(pred_diff)-1))

  #tm<-1:length(pred_diff)
  #b<-pred_diff[1]
  
  
  #mod<-nls(pred_diff~b0*exp(-exp(lmd)*tm)+b1,
  #         start=c(b0=b, lmd=log(abs(mean(log(pred_diff/b)/tm))), b1=pred_diff[length(pred_diff)]))
  #summary(mod)
  #plot(tm, pred_diff, type="l"); abline(h=0, v=0, lty=3)
  #lines(tm, coef(mod)[1]*exp(-exp(coef(mod)[2])*tm)+coef(mod)[3], col=2)
  
  #eigenlst[sppos]<-exp(coef(mod)[2])*(-sign(coef(mod)[1]-coef(mod)[3]))
}

matplot(1:99, eigenlst, type="l", col=1:3, lty=1); abline(h=0, lty=3)


#Estimate increase when rare
out_meta_2<-NULL
burnin<-100

for(i in 1:length(out_meta$plotdata$ceq)) {
  pt<-rep(0, length(out_meta$plotdata$ceq))
  pt[i]<-1
  tmp<-rerunrun_metapopulation(out=out_meta, tmax=burnin, talktime = 0, runtype = "metapopulation", perturb=pt, perturbsites=1:out_meta$plotdata$ngrid)
  
  newpop<-tmp$output[nrow(tmp$output),-1]
  newpop[i]<-ceiling((out_meta$plotdata$ngrid-sum(newpop[-i]))*0.05)
  
  population<-populate(gridout, nlst = newpop, clst = clst, radlst = Inf)
  
  out_meta_2[[i]]<-run_metapopulation(tmax=simtime, gridout=gridout, population=population, talktime = 0)
}


#out_meta_2[[1]]<-rerunrun_metapopulation(out=out_meta, tmax=simtime, talktime = 0, runtype = "metapopulation", perturb=c(0.95, 0, 0), perturbsites=1:out_meta$plotdata$ngrid)
#out_meta_2[[2]]<-rerunrun_metapopulation(out=out_meta, tmax=simtime, talktime = 0, runtype = "metapopulation", perturb=c(0, 0.95, 0), perturbsites=1:out_meta$plotdata$ngrid)
#out_meta_2[[3]]<-rerunrun_metapopulation(out=out_meta, tmax=simtime, talktime = 0, runtype = "metapopulation", perturb=c(0, 0, 0.95), perturbsites=1:out_meta$plotdata$ngrid)

grwrare<-matrix(ncol=length(out_meta_2), nrow=(simtime-1))

for(sppos in 1:length(out_meta_1)) {
  pred_grw<-out_meta_2[[sppos]]$output[,sppos+1]/out_meta$plotdata$ngrid
  
  grwrare[,sppos]<-cumsum(log(pred_grw[-1]/pred_grw[-length(pred_grw)]))/(1:(length(pred_grw)-1))
}

matplot(1:99, grwrare, type="l", col=1:3, lty=1); abline(h=0, lty=3)















#Test rEDM
Elst<-2:10
outcol_meta<-out_meta$output[,2]
simplout_meta<-suppressWarnings(simplex(outcol_meta, E=Elst))
E_meta<-Elst[min(which((max(simplout_meta$rho)-simplout_meta$rho)/diff(range(simplout_meta$rho))<0.1))]
plot(simplout_meta$E, simplout_meta$rho, type="l", xlab="E", ylab="rho"); abline(v=E_meta, lty=3)
predL_meta<-predict_vs_L(outcol = outcol_meta, E=E_meta)
predlag_meta<-test_predict_tlag(outcol_meta, Luse=min(c(floor(length(outcol_meta)/5), predL_meta$Lmin)), E=E_meta)

##### Try Hubbell NZNS neutral model
clst<-rep(0.5, 3)
population<-populate(gridout, nlst = round(rep(unique(abs(getceq(clst)))/length(clst), length(clst))*prod(gridout$lng)),
                     clst = clst, radlst = Inf)
out_neut<-run_metapopulation(tmax=1000, nsteps = 1000, gridout, population, talktime = 0, runtype = "neutral")
plot_metapop(out_neut)

#Test rEDM
outcol_neut<-out_neut$output[,2]
simplout_neut<-suppressWarnings(simplex(outcol_neut, E=Elst))
E_neu<-Elst[min(which((max(simplout_neut$rho)-simplout_neut$rho)/diff(range(simplout_neut$rho))<0.1))]
plot(simplout_neut$E, simplout_neut$rho, type="l", xlab="E", ylab="rho"); abline(v=E_neu, lty=3)
predL_neut<-predict_vs_L(outcol = outcol_neut, E=E_neu)
predlag_neut<-test_predict_tlag(outcol=outcol_neut, Luse=min(c(floor(length(outcol_neut)/5), predL_neut$Lmin)), E=E_neu)





