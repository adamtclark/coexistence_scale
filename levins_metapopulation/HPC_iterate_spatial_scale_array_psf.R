#!/usr/bin/env Rscript
#error
rm(list=ls())
#setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

#load packages
require(parallel)

#load scripts
source("run_metapopulation_wrapper.R")

########################################
# NEW functions
########################################
makemat<-function(eig1out, eig2out, r0out, nsp) {
  tmpmat<-matrix(ncol=nsp*6, nrow=max(c(nrow(eig1out$eigenlst),
                                           nrow(eig2out$eigenlst),
                                           nrow(r0out$grwrare),
                                           nrow(eig1out$eigenlst_tot),
                                           nrow(eig2out$eigenlst_tot),
                                           nrow(r0out$grwrare_tot),na.rm=T)))
  
  nn<-0
  tmpmat[1:nrow(eig1out$eigenlst),(1:nsp)+nsp*(nn)]<-eig1out$eigenlst; nn<-nn+1
  tmpmat[1:nrow(eig2out$eigenlst),(1:nsp)+nsp*(nn)]<-eig2out$eigenlst; nn<-nn+1
  tmpmat[1:nrow(r0out$grwrare),(1:nsp)+nsp*(nn)]<-r0out$grwrare; nn<-nn+1
  
  tmpmat[1:nrow(eig1out$eigenlst_tot),(1:nsp)+nsp*(nn)]<-eig1out$eigenlst_tot; nn<-nn+1
  tmpmat[1:nrow(eig2out$eigenlst_tot),(1:nsp)+nsp*(nn)]<-eig2out$eigenlst_tot; nn<-nn+1
  tmpmat[1:nrow(r0out$grwrare_tot),(1:nsp)+nsp*(nn)]<-r0out$grwrare_tot; nn<-nn+1
  
  return(tmpmat)
}

runpar_psf<-function(...) {
  #Function to automate stability tests across models
  #Is embarrassingly parallel, and can be used to iterate runs
  
  tst<-try({
    ##### run simulations
    #run psf
    sd<-round(runif(1)*1e9)
    
    set.seed(sd)
    out_psf<-run_metapopulation(tmax=tmax, gridout = gridout, population = population_psf, talktime = 0, runtype = "psf", sites_sub = grid_sub$sites)
    set.seed(sd)
    out_psf_long<-run_metapopulation(tmax=tmax_long, gridout = gridout, population = population_psf, talktime = 0, runtype = "psf", sites_sub = grid_sub$sites)
    
    eig_psf1<-estimate_eqreturn(out_psf, simtime=simtime, runtype="psf", replace_perturb = 1, talktime=0, prtb=ptb, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites, doplot = FALSE)
    eig_psf2<-estimate_eqreturn(out_psf, simtime=simtime, runtype="psf", replace_perturb = 0, talktime=0, prtb=ptb, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites, doplot = FALSE)
    r0_psf<-estimate_rarereturn(out_psf, simtime=simtime, burnin=burnin, runtype="psf", perturbsites = grid_sub$sites, sites_sub = grid_sub$sites, doplot = FALSE)
    
    ##### save output
    #make output matrix
    nsp<-length(clst_psf)
    
    #collect and output data
    matout<-cbind(makemat(eig1out=eig_psf1, eig2out=eig_psf2, r0out=r0_psf, nsp=nsp))
    
    matout<-t(matout)
  })
  
  if(is.character(tst)) {
    matout<-NULL
  }
  
  return(matout)
}

psfwrapper<-function(population, gridout, sites_sub, abundances, gridsize, c_sptraits, tmax, speciesid, xylim, nsteps, destroyed, tmpspdestroy, soilstart=0) {
  ##### Prepare run
  
  #adjust for extra time step in this function
  #tmax<-(tmax-1)
  #nsteps<-(nsteps-1)
  
  if(is.infinite(population$radlst)) {
    sp1dis=1; sp2dis=1        #global dispersal
  } else {
    sp1dis=0; sp2dis=0        #local dispersal
  }
  initabund=(abundances/gridsize)*100   #initial native abundance (sp1)
  seed=c_sptraits[1]/c_sptraits[2]      #ratio of seed production (exotic:native)
  
  #extract parameters from simulation run
  dim=ceiling(sqrt(gridsize)) #grid edge size
  tmax=tmax                   #maximum time
  outmat=numeric((tmax+1)*3)  #matrix for storing species abundances
  outmat_sub=outmat           #matrix for storing species abundances in subset
  outmap0=numeric((dim)^2)    #matrix for storing initial conditions
  outmap=outmap0              #matrix for storing end conditions
  
  #Fixed model parameters
  sp1fb=80; sp2fb=80        #feedback
  sp1m=5; sp2m=5              #stochastic mortality
  scenario=3                  #begin with neutral soils
  pr_nocol=0.5                #reduction in probability of colonization event - allows for empty cells
  stepsize=1                  #step size in by which changes in soil occur (in percent)
  edge=0                      #absorbing edge conditions
  
  #Set up starting conditions
  #locations of species
  speciesiduse<-speciesid #initially, largest number means "empty"
  speciesiduse[speciesiduse>2]<-2
  speciesiduse<-speciesiduse+1
  speciesiduse[speciesiduse==3]<-0 #now, zero means "empty"
  
  #Soil conditions for each species
  if(sum(abs(soilstart))==0) {
    soilstate<-matrix(nrow=(dim+2), ncol=(dim+2), data=0)
    if(scenario==0) {
      soilstate[]<-0
    } else if(scenario==1) {
      soilstate[]<-(-abs(sp1fb))
    } else if(scenario==2) {
      soilstate[]<-abs(sp2fb)
    } else if(scenario==3) {
      #start with state of occupied species
      soilstate[cbind(gridout$xpos+1, gridout$ypos+1)]<-c(0, (-abs(sp1fb)), abs(sp2fb))[speciesid+1]
    }
  } else {
    soilstate<-soilstart
  }
  
  #expand vectors to include edges
  speciesiduse_m<-matrix(speciesiduse, nrow=dim)
  speciesiduse_m<-rbind(speciesiduse_m, rep(0, dim))
  speciesiduse_m<-cbind(speciesiduse_m, rep(0, dim+1))
  speciesiduse_m<-cbind(rep(0, dim+1), speciesiduse_m)
  speciesiduse_m<-rbind(rep(0, dim+2), speciesiduse_m)
  
  ##### Run C code
  out<-.C("sis_metapopulation",
          psp1dis=as.integer(sp1dis), psp2dis=as.integer(sp2dis),
          ps1fb=as.integer(sp1fb), psp2fb=as.integer(sp2fb), 
          psp1m=as.integer(sp1m), psp2m=as.integer(sp2m),
          pinitabund=as.integer(initabund), pseed=as.double(seed),
          pedge=as.integer(edge), pdim=as.integer(dim), ptmax=as.integer(tmax),
          ppr_nocol=as.double(pr_nocol), pstepsize=as.integer(stepsize),
          outmat=as.integer(outmat), outmat_sub=as.integer(outmat_sub),
          outmap0=as.integer(outmap0), outmap=as.integer(outmap),
          speciesid=as.integer(c(speciesiduse_m)), 
          c_sites_sub=as.integer(sites_sub-1), plngsub=as.integer(length(sites_sub)),
          soilstate=as.integer(c(soilstate)),
          destroyed=as.integer(destroyed), spdestroy=as.integer(tmpspdestroy))
  
  ##### repack output
  #total abundance
  m<-matrix(out$outmat, ncol=3)
  #Exclude empty cells or invaded cells from the periphery
  if(edge==0) {
    m[,1]<-m[,1]-(4*(dim+2)-4)
  } else {
    m[,1]<-m[,1]-(3*(dim+2)-2)
    m[,2]<-m[,2]-(dim)
  }
  
  #subset abundance
  m_sub<-matrix(out$outmat_sub, ncol=3)
  
  #species locations
  speciesid<-c(out$outmap)
  speciesid<-speciesid-1
  speciesid[speciesid==-1]<-2 #now 2 is "empty"
  
  cout<-list(ptmax=tmax, pgridsize=gridsize, pnsp=length(population$nlst), xylim=xylim, destroyed=destroyed, spdestroy=tmpspdestroy, #grid
             c_sptraits=population$clst, m_sptraits=population$mlst, abundances=m[nrow(m),-1], colsites=NA, pncolsites=NA, #traits
             eventtimes_c=rep(NA, gridsize), eventtimes_m=rep(NA, gridsize), #events
             speciesid=speciesid, #species
             output=c(cbind(1:(tmax+1), m[,-1])), pnsteps=nsteps,
             ptalktime=NA,
             abundances_sub=m_sub[nrow(m_sub),-1], sites_sub=(sites_sub-1), pnsites_sub=length(sites_sub), output_sub=c(cbind(1:(tmax+1), m_sub[,-1])),
             population=population, gridout=gridout, soilstate=out$soilstate)
  
  return(cout)
}


#set up for runs
niterations<-20000   #CHANGE TO ALTER NUMBER OF ITERATIONS
scalelst<-c(0.005, 0.01, 0.05, 0.1, 0.5, 0.75, 1)
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

lglst<-round(seq(0, tmax_long*0.8-10, length=10)) #lags for invar test

clst_psf= c(1.5, 1)*xfac
mlst_psf = rep(0.1, length(clst_psf))*xfac
population_psf<-populate(gridout, nlst = rep(floor(prod(gridout$lng)/length(clst_psf)*0.8), length(clst_psf)),
                         clst = clst_psf, radlst = Inf, mlst = mlst_psf)

#open cluster
if(!exists("cl") & niterations>1) {
  cl <- makeForkCluster(mc <- getOption("cl.cores", min(c(niterations, detectCores())))) #cluters for simulations
}

#explor needed variables
clusterExport(cl, c("invarburn",
                    "gridout",
                    "population_psf",
                    "ptb",
                    "tmax", "tmax_long", "burnin", "simtime", "lglst",
                    "clst_psf", "mlst_psf",
                    "psfwrapper",
                    "run_metapopulation", "rerunrun_metapopulation", "getceq", "getE", "loadrun", "unloadrun", "getrunname",
                    "estimate_eqreturn", "estimate_rarereturn", "estimate_invar", "beta_estimate", "predict_vs_L", "test_predict_tlag",
                    "makemat", "makemat_inv"))

#run simulations
#for(i in 1:length(scalelst)) {
i<-as.integer(Sys.getenv("SGE_TASK_ID", "1"))

grid_sub<-grid_subset(gridout, size = scalelst[i])

clusterExport(cl, c("grid_sub"))

#run parallel program for predicting community biomass
clusterout<-try(parLapplyLB(cl=cl, 1:niterations, fun=runpar_psf))

if(!is.character(clusterout)) {
  tmp<-t(matrix(nrow=nrow(clusterout[[1]]), unlist(clusterout)))
  
  matout_tot<-cbind(scale=scalelst[i], tscale=1:ncol(clusterout[[1]]), iter=rep(1:niterations, each=ncol(clusterout[[1]])), tmp)
  
  #separate
  matout_dyn<-matout_tot
  
  #save outputs to csv
  save(list=c("matout_dyn"), file = paste("output/matout_dyn_psf_", scalelst[i], ".rda", sep=""))
}

print(round(i/length(scalelst),2))
#}

#close cluster
if(exists("cl")) {
  stopCluster(cl)
}

