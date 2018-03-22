#!/usr/bin/env Rscript
#error
rm(list=ls())

setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

#Load functions
source("run_metapopulation_wrapper.R")
source("util/filled.contour3.R")
source("util/figure_functions.R")
source("util/plot_grid_functions.R")
require(data.table)


if(FALSE) {
  ##############################
  # Data prep
  ##############################
  #Load data
  d.bm<-read.csv("data/e014/TableS3a_abundance_data.csv")
  d.pc<-read.csv("data/e014/TableS3b_percentcover_data.csv")
  d.fld<-read.csv("data/e014/TableS1_field_information.csv")
  d.spc<-read.csv("data/e014/TableS2_species_information.csv")
  
  #pc wide to long
  d.pc<-dcast(d.pc[!is.na(d.pc$Group),], Year+Field+Transect+Plot+Cover_Bareground+Cover_Litter~Group, value.var = "Cover_Group")
  d.pc[,c("Annual", "C3_Perennial", "C4S_Perennial", "Forb_Perennial", "Leg_Perennial")][is.na(d.pc[,c("Annual", "C3_Perennial", "C4S_Perennial", "Forb_Perennial", "Leg_Perennial")])]<-0
  
  #Add field info to data
  index.bm<-paste(d.bm$Field, d.bm$Transect)
  index.pc<-paste(d.pc$Field, d.pc$Transect)
  index.fld<-paste(d.fld$Field, d.fld$Transect)
  
  addcol<-c("ForestStatus", "Burned", "FirstYearAfterBurning", "YearAbandoned")
  d.bm<-data.frame(d.bm, d.fld[match(index.bm, index.fld),addcol])
  d.pc<-data.frame(d.pc, d.fld[match(index.pc, index.fld),addcol])
  
  d.bm$Age<-d.bm$Year-d.bm$YearAbandoned
  d.pc$Age<-d.pc$Year-d.pc$YearAbandoned
  
  d.bm<-d.bm[,c("Year", "Field","Transect",
                "ForestStatus", "Burned", "FirstYearAfterBurning",
                "YearAbandoned","Age",
                "Biomass_Litter", "Biomass_TotalVegetation",
                "Biomass_Annual", "Biomass_C3_Perennial", "Biomass_C4S_Perennial",
                "Biomass_Forb_Perennial", "Biomass_Leg_Perennial")]
  names(d.bm)[names(d.bm)%in%c("Biomass_Annual", "Biomass_C3_Perennial", "Biomass_C4S_Perennial",
                               "Biomass_Forb_Perennial", "Biomass_Leg_Perennial")]<-
    c("Annual", "C3_Perennial", "C4S_Perennial", "Forb_Perennial", "Leg_Perennial")
  
  d.pc<-d.pc[,c("Year", "Field","Transect", "Plot",
                "ForestStatus", "Burned", "FirstYearAfterBurning",
                "YearAbandoned","Age",
                "Cover_Bareground", "Cover_Litter",
                "Annual", "C3_Perennial", "C4S_Perennial",
                "Forb_Perennial", "Leg_Perennial")]
  
  #Get subset of data
  sd.pc<-d.pc[d.pc$ForestStatus!="HF" & d.pc$Burned==0,]
  sd.bm<-d.bm[d.bm$ForestStatus!="HF" & d.bm$Burned==0,]
  
  ##############################
  # Get rates
  ##############################
  # Get estimate of "start state"
  tmp<-sd.pc[sd.pc$Age<=5,c("Field", "Annual", "C3_Perennial", "C4S_Perennial","Forb_Perennial", "Leg_Perennial")]
  tmp2<-as.matrix(log(tmp[,-1]))
  tmp2[!is.finite(tmp2)]<-NA
  dmin<-0.0001
  state0<-exp(colMeans(tmp2, na.rm=T))*colMeans(tmp[,-1]>0)+dmin
  
  sd.pc$Annual_r0<-log((sd.pc$Annual+dmin)/state0[1])#*(1/sd.pc$Age)
  sd.pc$C3_Perennial_r0<-log((sd.pc$C3_Perennial+dmin)/state0[2])#*(1/sd.pc$Age)
  sd.pc$C4S_Perennial_r0<-log((sd.pc$C4S_Perennial+dmin)/state0[3])#*(1/sd.pc$Age)
  sd.pc$Forb_Perennial_r0<-log((sd.pc$Forb_Perennial+dmin)/state0[4])#*(1/sd.pc$Age)
  sd.pc$Leg_Perennial_r0<-log((sd.pc$Leg_Perennial+dmin)/state0[5])#*(1/sd.pc$Age)
  
  #Get means by scales
  bks<-seq(0, 90, by=5)
  dtmax<-5
  sd.pc$Age.cat<-cut(sd.pc$Age, breaks = bks)
  timebands<-sort(unique(sd.pc$Age))#levels(sd.pc$Age.cat)
  tmp<-table(sd.pc$Age.cat, sd.pc$Field)
  mxplot<-max(rowSums(tmp))
  mnplot<-min(rowSums(tmp))
  mxfld<-max(rowSums(tmp>0))
  niter<-1000
  
  arrayout<-array(dim=c(length(timebands), 16*mnplot, 3, niter))
  
  r0fun<-function(vec, s0) {
    #get mean for each scale
    tmp<-unname(log((cumsum(vec)/(1:length(vec))+0.0001)/s0))
    
    tmp[!is.finite(tmp)]<-NA
    tmp
  }
  
  for(iter in 1:niter) {
    for(i in 1:length(timebands)) {
      #ps<-which(sd.pc$Age.cat==timebands[i])
      ps<-which(abs(sd.pc$Age-timebands[i])<=dtmax)
      fldlst_tmp<-sort(unique(sd.pc$Field[ps]))
      fldord<-sample(1:length(fldlst_tmp), length(fldlst_tmp), rep=F)
      
      psinclude<-matrix(nrow=mnplot, ncol=length(fldlst_tmp))
      
      for(j in 1:length(fldlst_tmp)) {
        ps2<-which(sd.pc$Field[ps]==fldlst_tmp[fldord[j]])
        
        if(length(ps2)>50) {
          #select year/transect
          tmp<-unique(sd.pc$Transect[ps[ps2]])
          if(length(tmp)>1) {
            trns_smp<-sample(tmp, 2)
          } else {
            trns_smp<-tmp
          }
          
          tmp<-unique(sd.pc$Year[ps[ps2]])
          if(length(tmp)>1) {
            yr_smp<-sample(tmp, 1)
          } else {
            yr_smp<-tmp
          }
          
          #save in psinclude
          tmp<-ps[ps2[which((sd.pc$Transect[ps[ps2]] %in% trns_smp)&
                              (sd.pc$Year[ps[ps2]] %in% yr_smp))]]
          psinclude[1:length(tmp),j]<-tmp
          
        } else {
          #save in psinclude
          psinclude[c(1:length(ps2)),j]<-ps[ps2]
        }
      }
      
      arrayout[i,1:length(psinclude),1,iter]<-r0fun(sd.pc$Annual[c(psinclude[,1:j])], state0[1])/timebands[i]
      arrayout[i,1:length(psinclude),2,iter]<-r0fun(sd.pc$C3_Perennial[c(psinclude[,1:j])], state0[2])/timebands[i]
      arrayout[i,1:length(psinclude),3,iter]<-r0fun(sd.pc$C4S_Perennial[c(psinclude[,1:j])], state0[3])/timebands[i]
    }
    if(iter/10==floor(iter/10)) {
      print(round(iter/niter,2))
    }
  }
  
  #Get quantiles for each...
  array_quant<-apply(arrayout, 1:3, function(x) quantile(x, c(0.5), na.rm=T))
  xscl<-(1:dim(array_quant)[2])/2
  
  
  #Shrink matrix
  #spatial
  xscl_small<-sort(unique(round(c(1, exp(seq(log(1), log(dim(array_quant)[2]), length=50)), dim(array_quant)[2]))))
  xscl_small<-xscl_small[xscl_small<=dim(array_quant)[2]]
  
  fuzzmat<-matrix(ncol=length(xscl_small), nrow=dim(array_quant)[2], data=0)
  for(i in 1:ncol(fuzzmat)) {
    fuzzmat[(c(0, xscl_small)[i]):(c(0, xscl_small)[i+1]),i]<-1/length((c(0, xscl_small)[i]):(c(0, xscl_small)[i+1]))
  }
  array_quant_small<-array(dim=c(dim(array_quant)[1], length(xscl_small), dim(array_quant)[3]))
  for(i in 1:dim(array_quant)[3]) {
    tmp<-array_quant[,,i]
    tmp[is.na(tmp)]<-0
    array_quant_small[,,i]<-tmp%*%fuzzmat
    array_quant_small[,,i]<-array_quant_small[,,i]*timebands
    
    array_quant_small[,,i][array_quant_small[,,i]==0]<-NA
  }
  
  #temporal
  timebands_small<-sort(unique(round(c(1, exp(seq(log(1), log(max(timebands)), length=40)), max(timebands)))))
  timebands_small<-timebands_small[timebands_small<=max(timebands)]
  
  fuzzmat<-matrix(ncol=length(timebands_small), nrow=dim(array_quant_small)[1], data=0)
  for(i in 1:ncol(fuzzmat)) {
    tmp<-which((timebands>c(0,timebands_small)[i]) & (timebands<=timebands_small[i]))
    fuzzmat[tmp,i]<-1/length(tmp)
  }
  tmp<-which(colSums(fuzzmat)>0)
  fuzzmat<-fuzzmat[,tmp]
  timebands_small<-timebands_small[tmp]
  
  array_quant_small2<-array(dim=c(length(timebands_small), dim(array_quant_small)[2], dim(array_quant_small)[3]))
  for(i in 1:dim(array_quant)[3]) {
    tmp<-array_quant_small[,,i]
    tmp[is.na(tmp)]<-0
    array_quant_small2[,,i]<-t(t(tmp)%*%fuzzmat)
    array_quant_small2[,,i][array_quant_small2[,,i]==0]<-NA
  }
  
    
  ##############################
  # Run simulation...
  ##############################
  #set up for runs
  niterations<-niter   #CHANGE TO ALTER NUMBER OF ITERATIONS
  scalelst<-c(1,2,3,5,7,10,20,30,50,70,100,200,300,400)/c(100^2)
  radlst<-Inf
  
  #set up simulations
  gridout<-makegrid(xlng = 100, ylng = 100)
  
  simtime<-81 #time spans for equilibria dectection
  
  clst_meta = c(0.1,0.3,0.4)
  mlst_meta = c(0.02, 0.02, 0.02)
  
  array_sim<-array(dim=c(simtime-1, length(scalelst), 3, niter))
  s0<-unname(state0[3:1])
  s0_ind<-round(s0*prod(gridout$lng))
  s0_ind_back<-s0_ind/prod(gridout$lng)
  
  for(i in 1:niter) {
    population<-populate(gridout, nlst = s0_ind,
                            clst = clst_meta, radlst = Inf, mlst = mlst_meta)
    for(j in 1:length(scalelst)) {
      grid_sub<-grid_subset(gridout, size = scalelst[j])
      
      out<-run_metapopulation(tmax=simtime, nsteps = simtime, gridout, population, talktime = 0, sites_sub = grid_sub$sites)
      ngrid<-length(grid_sub$sites)
      r0<-t(log(t(out$output_spatial[-1,-1]/ngrid+dmin)/(s0_ind_back)))[,c(3:1)]
      
      array_sim[,j,,i]<-r0
    }
    if(i/10==floor(i/10)) {
      print(round(i/niter,2))
    }
  }
  
  #Get quantiles
  array_sim_quant<-apply(array_sim, 1:3, function(x) quantile(x, c(0.5), na.rm=T))
  
  #Shrink matrix temporally
  timebands_small<-sort(unique(round(c(1, exp(seq(log(1), log(max(timebands)), length=40)), max(timebands)))))
  timebands_small<-timebands_small[timebands_small<=max(timebands)]
  
  fuzzmat<-matrix(ncol=length(timebands_small), nrow=dim(array_sim_quant)[1], data=0)
  for(i in 1:ncol(fuzzmat)) {
    tmp<-which((timebands>c(0,timebands_small)[i]) & (timebands<=timebands_small[i]))
    fuzzmat[tmp,i]<-1/length(tmp)
  }
  tmp<-which(colSums(fuzzmat)>0)
  fuzzmat<-fuzzmat[,tmp]
  timebands_small<-timebands_small[tmp]
  
  array_sim_quant_small2<-array(dim=c(length(timebands_small), dim(array_sim_quant)[2], dim(array_sim_quant)[3]))
  for(i in 1:dim(array_sim_quant)[3]) {
    tmp<-array_sim_quant[,,i]
    tmp[is.na(tmp)]<-0
    array_sim_quant_small2[,,i]<-t(t(tmp)%*%fuzzmat)
    array_sim_quant_small2[,,i][array_sim_quant_small2[,,i]==0]<-NA
  }
  
  
  
  save.image("output/save_e014_comp.RData")
} else {
  load("output/save_e014_comp.RData")
}





##############################
#Plot output
##############################
m<-(matrix(nrow=3, 1:6))
m<-cbind(m, 7)
layout(m, widths=c(1,1,0.3))
par(mar=c(2,3,1,0), oma=c(2.5,2,2,3))
#plotout<-plot_cont_emp(arrayout=array_quant_small, xscalslst=log10(xscl_small), xlst=log10(timebands), nlevels=10, logx=TRUE, logy=TRUE, logz=FALSE, nstart=1, ofs1=c(0, 0))
squse<-c(-8,-4,-2,-1,0,2,4,8)#c(-100, seq(-4,0, by=1), seq(1,4,by=1), 100)
plotout<-plot_cont_emp(arrayout=array_quant_small2, xscalslst=log10(xscl_small), xlst=log10(timebands_small), nlevels=10, logx=TRUE, logy=TRUE, logz=FALSE, nstart=1, sqlst=squse, ofs1=c(0, 0))
#plotout<-plot_cont_emp(arrayout=array_quant*timebands, xscalslst=log10(xscl*2), xlst=log10(timebands), nlevels=10, logx=TRUE, logy=TRUE, logz=FALSE, nstart=1, sqlst=squse, ofs1=c(0, 0))
plotout<-plot_cont_emp(arrayout=array_sim_quant_small2, xscalslst=log10(scalelst*prod(gridout$lng)), xlst=log10(timebands_small), nlevels=10, logx=TRUE, logy=TRUE, logz=FALSE, nstart=4, sqlst=squse, ofs1=c(0, 0))

par(mar=c(2,3.5,1,1))
sqplot<-c(-8,-4,-2,-1,0,2,4,8)#c(-5, seq(-4,0, by=1), seq(1,4,by=1), 5)
filled.legend(z=matrix(sqplot), levels=sqplot, col=plotout$collst2, key.axes = axis(4, at = sqplot, labels = sqplot, las=2))



##############################
#Plot temporal trend
##############################
#Get average trajectory
tmtmp<-timebands
tmtmp[timebands>=79]<-timebands[timebands>=79]+1
abunds_obs<-exp(array_quant_small[,which.min(abs(xscl_small-100)),])*rep(state0[1:3], each=dim(array_quant_small)[1])-0.0001
#matplot(tmtmp, abunds_obs, type="l", col=c("blue", "green", "red"), lty=1, lwd=2, xlab="Age", ylab="pcover", xlim=c(1,80), xaxs="i"); abline(h=0, lty=3)

par(mfrow=c(2,1), mar=c(4,4,2,2))
matplot(out$output[,1], out$output[,-1]/out$plotdata$ngrid, type="l", col=c("green", "blue", "red"), lty=1, lwd=2, xlab="Age", ylab="pcover", xlim=c(1,80), xaxs="i"); abline(h=0, lty=3)
matplot(tmtmp, abunds_obs, type="l", col=c("red", "blue", "green"), lty=1, lwd=2, xlab="Age", ylab="pcover", xlim=c(1,80), xaxs="i"); abline(h=0, lty=3)


#Can do same for simulation output...

grid_sub<-grid_subset(gridout, size = scalelst[1])
out_sub<-run_metapopulation(tmax=simtime, gridout = gridout, population = population, talktime = 0, runtype = "metapopulation", sites_sub = grid_sub$sites)
matplot(out_sub$output_spatial[,1], out_sub$output_spatial[,-1]/length(grid_sub$sites), type="l", col=c("green", "blue", "red"), lty=1, lwd=2, xlab="Age", ylab="pcover", xlim=c(1,80), xaxs="i", ylim=c(0, 0.85)); abline(h=0, lty=3)


