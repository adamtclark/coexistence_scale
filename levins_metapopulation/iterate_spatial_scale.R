#!/usr/bin/env Rscript
#error
rm(list=ls())
#setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

#load packages
require(rEDM)
require(parallel)

#load scripts
source("run_metapopulation_wrapper.R")

#set up for runs
niterations<-1000   #CHANGE TO ALTER NUMBER OF ITERATIONS
scalelst<-c(0.01, 0.05, 0.1, 0.5, 0.75, 1)
radlst<-Inf

#set up simulations
gridout<-makegrid(xlng = 100, ylng = 100)
xfac<-3
xfac_fast<-7
ptb<-0.2

tmax<-300 #timeseries length
tmax_long<-1000
burnin<-200 #burning for growth rate when rare method
simtime<-200 #time spans for equilibria dectection
invarburn<-200

lglst<-round(seq(0, tmax*0.8-10, length=10)) #lags for invar test

clst_meta = c(0.15, 0.35)*xfac
mlst_meta = rep(0.1, length(clst_meta))*xfac
population_meta<-populate(gridout, nlst = floor(getceq(clst_meta, mlst_meta)*prod(gridout$lng)),
                          clst = clst_meta, radlst = Inf, mlst = mlst_meta)

clst_neut<-rep(0.5, 2)*xfac
mlst_neut<-rep(0.1, length(clst_neut))*xfac_fast
tmp<-abs(getceq(clst_neut, mlst_neut))
population_neut<-populate(gridout, nlst = round(rep(unique(tmp[tmp>0])/length(clst_neut), length(clst_neut))*prod(gridout$lng)),
                          clst = clst_neut, radlst = Inf, mlst = mlst_neut)

clst_dist= c(0.145, 0.2)*xfac
mlst_dist = rep(0.1, length(clst_dist))*xfac
distlst<-c(0.95, 0)
distfrq<-50
population_dist<-populate(gridout, nlst = rep(floor(prod(gridout$lng)/length(clst_dist)*0.8), length(clst_dist)),
                          clst = clst_dist, radlst = Inf, mlst = mlst_dist)

clst_psf= c(1.5, 1)*xfac
mlst_psf = rep(0.1, length(clst_psf))*xfac
population_psf<-populate(gridout, nlst = rep(floor(prod(gridout$lng)/length(clst_psf)*0.8), length(clst_psf)),
                         clst = clst_psf, radlst = Inf, mlst = mlst_psf)

intmat_rps<-rbind(c(1,0,0,1),
                  c(1,1,0,0),
                  c(0,1,1,0),
                  c(0,0,1,1))
clst_rps<-rep(0.4, 4)*xfac*0.4
mlst_rps<-rep(0.08, length(clst_rps))*xfac_fast*0.4
tmp<-abs(getceq(clst_rps, mlst_rps))
population_rps<-populate(gridout, nlst = round(rep(unique(tmp[tmp>0])/length(clst_rps), length(clst_rps))*prod(gridout$lng)),
                         clst = clst_rps, radlst = Inf, mlst = mlst_rps)

#open cluster
if(!exists("cl") & niterations>1) {
  cl <- makeCluster(mc <- getOption("cl.cores", min(c(niterations, detectCores())))) #cluters for simulations
}

#explor needed variables
clusterExport(cl, c("invarburn",
                    "gridout",
                    "population_meta", "population_neut", "population_dist", "population_psf", "population_rps",
                    "ptb", "intmat_rps",
                    "tmax", "tmax_long", "burnin", "simtime", "lglst",
                    "clst_meta", "mlst_meta", "clst_neut", "mlst_neut", "clst_dist", "mlst_dist", "clst_psf", "mlst_psf", "clst_rps", "mlst_rps",
                    "distlst", "distfrq",
                    "psfwrapper",
                    "run_metapopulation", "rerunrun_metapopulation", "getceq", "getE", "loadrun", "unloadrun", "getrunname",
                    "estimate_eqreturn", "estimate_rarereturn", "estimate_invar", "beta_estimate", "predict_vs_L", "test_predict_tlag",
                    "makemat", "makemat_inv"))

#run simulations
matout_tot<-NULL

for(i in 1:length(scalelst)) {
  grid_sub<-grid_subset(gridout, size = scalelst[i])
  
  clusterExport(cl, c("grid_sub"))
  
  #run parallel program for predicting community biomass
  clusterout<-try(parLapply(cl=cl, 1:niterations, fun=runpar))
  
  if(!is.character(clusterout)) {
    tmp<-t(matrix(nrow=nrow(clusterout[[1]]), unlist(clusterout)))
    
    matout_tot<-rbind(matout_tot, cbind(scale=scalelst[i], tscale=1:ncol(clusterout[[1]]), iter=rep(1:niterations, each=ncol(clusterout[[1]])), tmp))
    
    #save outputs to csv
    write.csv(matout_tot, "output/matout_tot.csv", row.names=F)
  }
  
  print(round(i/length(scalelst),2))
}

#close cluster
if(exists("cl")) {
  stopCluster(cl)
}

