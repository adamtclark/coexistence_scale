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
niterations<-2
scalelst<-c(0.01, 0.02, 0.05, 0.1, 0.3, 0.5, 1)

#set up simulations
gridout<-makegrid(xlng = 100, ylng = 100) #grid for simulation
xfac<-5 #scaling factor for c and m paramters
ptb<-0.2 #size of perturbations for equilibrium method (in fraction of initial value)

tmax<-1000 #timeseries length
burnin<-100 #burning for growth rate when rare method
simtime<-100 #time spans for equilibria dectection

lglst<-round(seq(0, tmax, length=20)) #lags for invar test

minE<-4 #minimum E for simplex algorithm

clst_meta = c(0.15, 0.3, 0.8, 3)*xfac
mlst_meta = rep(0.1, length(clst_meta))*xfac
population_meta<-populate(gridout, nlst = floor(getceq(clst_meta, mlst_meta)*prod(gridout$lng)),
                          clst = clst_meta, radlst = Inf, mlst = mlst_meta)

clst_neut<-rep(0.5, 4)*xfac
mlst_neut<-rep(0.1, length(clst_neut))*xfac
population_neut<-populate(gridout, nlst = round(rep(unique(abs(getceq(clst_neut, mlst_neut)))/length(clst_neut), length(clst_neut))*prod(gridout$lng)),
                          clst = clst_neut, radlst = Inf, mlst = mlst_neut)

#open cluster
if(!exists("cl") & niterations>1) {
  cl <- makeCluster(mc <- getOption("cl.cores", min(c(niterations, detectCores())))) #cluters for simulations
}

#explor needed variables
clusterExport(cl, c("gridout", "population_meta", "population_neut", "ptb",
                    "tmax", "burnin", "simtime", "lglst",
                    "clst_meta", "mlst_meta", "clst_neut", "mlst_neut",
                    "run_metapopulation", "rerunrun_metapopulation", "getceq", "getE", "minE",
                    "estimate_eqreturn", "estimate_rarereturn", "estimate_invar", "predict_vs_L", "test_predict_tlag")) #need to add in all functions and sub-functions here

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

