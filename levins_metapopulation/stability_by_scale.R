error
rm(list=ls())
setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

#TODO:
#1. Will have to find out how tore-package event times.
#2. Need to calculate eigenvalues/Chesson's metrics for comparison.


#load functions
require(rEDM)
source("run_metapopulation_wrapper.R")

par(mfcol=c(4,2), mar=c(4,4,2,2))
set.seed(1152)

##### Try Tilman metapopulation model
gridout<-makegrid(xlng = 20, ylng = 20)
clst = c(0.15, 0.3, 0.8)
#getceq(clst)
population<-populate(gridout, nlst = floor(getceq(clst)*prod(gridout$lng)), clst = clst, radlst = Inf)

out_meta<-run_metapopulation(tmax=1000, nsteps = 1000, gridout, population, talktime = 0)
plot_metapop(out_meta)
#plot_map(out_meta, gridout)

#Test rEDM
Elst<-2:10
outcol_meta<-out_meta$output[,2]
simplout_meta<-suppressWarnings(simplex(outcol_meta, E=Elst))
E_meta<-Elst[min(which((max(simplout_meta$rho)-simplout_meta$rho)/diff(range(simplout_meta$rho))<0.1))]
plot(simplout_meta$E, simplout_meta$rho, type="l", xlab="E", ylab="rho"); abline(v=E_meta, lty=3)
predL_meta<-predict_vs_L(outcol = outcol_meta, E=E_meta)
predlag_meta<-test_predict_tlag(outcol_meta, Luse=min(c(floor(length(outcol_meta)/5), predL_meta$Lmin)), E=E_meta)

##### Try Hubbell NZNS neutral model
population<-populate(gridout, nlst = rep(floor(prod(gridout$lng)/3), 3), clst = rep(0.5, 3), radlst = Inf)
out_neut<-run_metapopulation(tmax=1000, nsteps = 1000, gridout, population, talktime = 0, runtype = "neutral")
plot_metapop(out_neut)

#Test rEDM
outcol_neut<-out_neut$output[,2]
simplout_neut<-suppressWarnings(simplex(outcol_neut, E=Elst))
E_neu<-Elst[min(which((max(simplout_neut$rho)-simplout_neut$rho)/diff(range(simplout_neut$rho))<0.1))]
plot(simplout_neut$E, simplout_neut$rho, type="l", xlab="E", ylab="rho"); abline(v=E_neu, lty=3)
predL_neut<-predict_vs_L(outcol = outcol_neut, E=E_neu)
predlag_neut<-test_predict_tlag(outcol=outcol_neut, Luse=min(c(floor(length(outcol_neut)/5), predL_neut$Lmin)), E=E_neu)




#TODO:
#1. Update L function, and find E
#2. Update functions for spatial (re-package)
#3. Find ways to calculate eigens and Chesson




#re-wrapping example:
#population1<-rewrap_pop(out1, population)
#out1<-run_metapopulation(tmax=100, nsteps = 100, gridout, population1, talktime = 0)
#plot_metapop(out1)
#plot_map(out1, gridout)

