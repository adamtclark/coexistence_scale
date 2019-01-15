#error
rm(list=ls())
setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

#load scripts
require(RColorBrewer)
source("run_metapopulation_wrapper.R")
source("util/figure_functions.R")

########################################
# NEW functions
########################################
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
  sp1fb=80; sp2fb=80          #feedback
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
      soilstate[cbind(gridout$xpos+1, gridout$ypos+1)]<-c(0, (abs(sp1fb)), -abs(sp2fb))[speciesid+1]
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



modplotfun_psf<-function(out, r0out, collst, burnin=0, doaxis1=T, plotpos=1, ...) {
  #original fxn
  abunds<-out$output
  if(burnin>0) {
    abunds<-abunds[-c(1:burnin),]
    abunds[,1]<-abunds[,1]-burnin
  }
  
  mxt<-ceiling(max(abunds[,1]))
  mxt_r0_0<-ceiling(max(r0out$out0_lst[[plotpos]]$output[,1]))
  mxt_r0<-ceiling(max(r0out$out_lst[[plotpos]]$output[,1]))
  
  #eig & r0
  tmp_r0_0<-r0out$out0_lst[[plotpos]]$output
  tmp_r0_0[,1]<-tmp_r0_0[,1]+mxt#+mxt_eig
  
  tmp_r0<-r0out$out_lst[[plotpos]]$output
  tmp_r0[,1]<-tmp_r0[,1]+mxt+mxt_r0_0#+mxt_eig
  
  #combine
  abunds2<-rbind(abunds, tmp_r0_0, tmp_r0)
  
  pabunds2<-abunds2/out$plotdata$ngrid
  pabunds2[,1]<-abunds2[,1]
  
  #plot r0
  suppressWarnings(matplot(pabunds2[,1], pabunds2[,-1], type="l", lty=1, col=collst, lwd=c(1.5, rep(1.5, ncol(pabunds2)-2)), xlab="", ylab="", xaxs="i", axes=F, xlim=range(c(0, pabunds2[,1])), ...))
  
  abline(v=c(mxt+mxt_r0_0), lty=2); abline(h=c(0,1), lty=3, lwd=1)
  
  
  #Add in disturbance lines
  ap1<-tmp_r0[1,plotpos+1]/out$plotdata$ngrid
  ap2<-0
  ap1<-max(c(ap1, ap2+diff(range(pabunds2[,-1],na.rm=T))*0.1))
  arrows(mxt+mxt_r0_0, ap2, mxt+mxt_r0_0, ap1, lend=2, length = 0.06, col=1, lwd=1.5)
  
  if(doaxis1) {
    axis(1, cex.axis=1.6)
  }
  axis(2, las=2, cex.axis=1.6); box()
  
  return(list(pabunds2=pabunds2))
}


###########################
#set up for runs
scalelst<-c(1)
radlst<-Inf

#set up simulations
gridout<-makegrid(xlng = 100, ylng = 100)
xfac<-3
xfac_fast<-7
ptb<-0.2

tmax<-300 #timeseries length
tmax_long<-1000
burnin<-100 #burning for growth rate when rare method
simtime<-200 #time spans for equilibria dectection
invarburn<-200

lglst<-round(seq(0, tmax_long*0.8-10, length=10)) #lags for invar test

clst_psf= c(1, 1.5)*xfac
mlst_psf = rep(0.1, length(clst_psf))*xfac


#Species 1
population_psf1<-populate(gridout, nlst = c(0, 0.1),
                         clst = clst_psf, radlst = Inf, mlst = mlst_psf)

grid_sub<-grid_subset(gridout, size = scalelst[1])

out_psf_long1<-run_metapopulation(tmax=tmax_long, gridout = gridout, population = population_psf1, talktime = 0, runtype = "psf", sites_sub = grid_sub$sites)

out_psf1<-rerunrun_metapopulation(out_psf_long1, tmax=100, runtype="psf", talktime = 0)
r0_psf1<-estimate_rarereturn(out_psf1, simtime=simtime, burnin=burnin, runtype="psf", perturbsites = grid_sub$sites, sites_sub = grid_sub$sites, doplot = FALSE, add_amount = 1)


#Species 2
population_psf2<-populate(gridout, nlst = c(0.1, 0),
                         clst = clst_psf, radlst = Inf, mlst = mlst_psf)

out_psf_long2<-run_metapopulation(tmax=tmax_long, gridout = gridout, population = population_psf2, talktime = 0, runtype = "psf", sites_sub = grid_sub$sites)

out_psf2<-rerunrun_metapopulation(out_psf_long2, tmax=100, runtype="psf", talktime = 0)
r0_psf2<-estimate_rarereturn(out_psf2, simtime=simtime, burnin=burnin, runtype="psf", perturbsites = grid_sub$sites, sites_sub = grid_sub$sites, doplot = FALSE, add_amount = 1)


#Plot results
collst<-adjustcolor(c(brewer.pal(4, "Set1")), alpha.f = 0.7)

par(mfrow=c(2,1), mar=c(4,4,2,2))
plot1<-modplotfun_psf(out_psf1, r0_psf1, collst, plotpos=1)
plot2<-modplotfun_psf(out_psf2, r0_psf2, collst, plotpos=2)


