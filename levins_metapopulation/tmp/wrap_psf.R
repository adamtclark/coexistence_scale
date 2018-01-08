error
rm(list=ls())
setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

#load functions
require(rEDM)
source("run_metapopulation_wrapper.R")
require(RColorBrewer)
source("~/Dropbox/Rfunctions/figure_functions.R")

##### Try PSF model
gridout<-makegrid(xlng = 100, ylng = 100)
grid_sub<-grid_subset(gridout, size = 0.05)

xfac<-5
ptb<-0.1
collst<-c("black", brewer.pal(4, "Dark2"))

set.seed(171205)
clst_psf = c(1.5, 1)*xfac
mlst_psf = rep(0.1, length(clst_psf))*xfac

population_psf<-populate(gridout, nlst = c(0.4, 0.4),
                          clst = clst_psf, radlst = Inf, mlst = mlst_psf)


tmax=200; gridout = gridout; population = population_psf; runtype = "psf"; talktime = 0
nsteps=tmax; sites_sub=grid_sub$sites; prt=0; prtfrq=0

#out_psf<-run_metapopulation(tmax=200, gridout = gridout, population = population_psf, runtype = "psf", talktime = 0)






psfwrapper<-function(population, gridout, grid_sub, abundances, gridsize, c_sptraits, tmax, speciesid) {
  ##### Prepare run
  
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
  outmap0=numeric((dim)^2)  #matrix for storing initial conditions
  outmap=outmap0              #matrix for storing end conditions
  
  #Fixed model parameters
  sp1fb=-80; sp2fb=-80        #feedback
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
          c_sites_sub=as.integer(grid_sub$sites), plngsub=as.integer(length(grid_sub$sites)),
          soilstate=as.integer(c(soilstate)))
  
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
  
  cout<-list(ptmax=tmax, pgridsize=gridsize, pnsp=nsp, xylim=xylim, destroyed=destroyed, spdestroy=rep(0, nsp), #grid
            c_sptraits=c_sptraits, m_sptraits=m_sptraits, abundances=m[nrow(m),-1], colsites=NA, pncolsites=NA, #traits
            eventtimes_c=NA, eventtimes_m=NA, #events
            speciesid=speciesid, #species
            output=c(cbind(1:(tmax+1), m[,-1])), pnsteps=nsteps,
            ptalktime=talktime,
            abundances_sub=m_sub[nrow(m_sub),-1], sites_sub=c_sites_sub, pnsites_sub=pnsites_sub, output_sub=c(cbind(1:(tmax+1), m_sub[,-1])))
  
  return(cout)
}

str(out$full)

system("R CMD SHLIB sis_metapopulation.c")
dyn.unload("sis_metapopulation.so")
dyn.load("sis_metapopulation.so")











cout<-.C(runname,
         ptmax= as.double(tmax), pgridsize=as.integer(gridsize), pnsp=as.integer(nsp), xylim=as.integer(xylim), destroyed=as.integer(destroyed), spdestroy=as.integer(rep(0, nsp)), #grid
         c_sptraits=as.double(c_sptraits), m_sptraits=as.double(m_sptraits), abundances=as.integer(abundances), colsites=as.integer(colsites), pncolsites=as.integer(ncolsites), #traits
         eventtimes_c=as.double(eventtimes_c), eventtimes_m=as.double(eventtimes_m), #events
         speciesid=as.integer(speciesid), #species
         output=as.double(output), pnsteps=as.integer(nsteps),
         ptalktime=as.integer(talktime),
         abundances_sub=as.integer(abundances_sub), sites_sub=as.integer(c_sites_sub), pnsites_sub=as.integer(pnsites_sub), output_sub=as.double(output_sub))

#save output
out_spatial<-matrix(cout$output_sub, nrow=nsteps+1)
out_spatial<-out_spatial[out_spatial[,1]!=0 | (1:nrow(out_spatial))==1,]





