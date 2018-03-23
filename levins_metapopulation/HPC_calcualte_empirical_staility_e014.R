#!/usr/bin/env Rscript
#error
rm(list=ls())
tstarg<-commandArgs(1)
#setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

#Load functions
source("run_metapopulation_wrapper.R")
source("util/filled.contour3.R")
source("util/figure_functions.R")
source("util/plot_grid_functions.R")

#main input variables
niter<-100

if(FALSE) {
  require(data.table)
  
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
  
  
  save(list = c("sd.pc"), file = "output/sd.pc.rd")
}
#Get empirical data
load("output/sd.pc.rd")

##############################
# Get rates
##############################
# Get estimate of "start state"
tmp<-sd.pc[sd.pc$Age<=5,c("Field", "Annual", "C3_Perennial", "C4S_Perennial","Forb_Perennial", "Leg_Perennial")]
tmp2<-as.matrix(log(tmp[,-1]))
tmp2[!is.finite(tmp2)]<-NA
dmin<-0.0001
state0<-exp(colMeans(tmp2, na.rm=T))*colMeans(tmp[,-1]>0)+dmin

#Get means by scales
bks<-seq(0, 90, by=5)
dtmax<-5
timebands<-sort(unique(sd.pc$Age))
sd.pc$Age.cat<-cut(sd.pc$Age, breaks = bks)
tmp<-table(sd.pc$Age.cat, sd.pc$Field)
mxplot<-max(rowSums(tmp))
mnplot<-min(rowSums(tmp))
mxfld<-max(rowSums(tmp>0))
mxfldplt<-16


#Shrink matrix
#spatial
dmdsp<-mxfldplt*mnplot

xscl_small<-sort(unique(round(c(1, exp(seq(log(1), log(dmdsp), length=50)), dmdsp))))
xscl_small<-xscl_small[xscl_small<=dmdsp]

fuzzmat_space<-matrix(ncol=length(xscl_small), nrow=dmdsp, data=0)
for(i in 1:ncol(fuzzmat_space)) {
  fuzzmat_space[(c(0, xscl_small)[i]):(c(0, xscl_small)[i+1]),i]<-1/length((c(0, xscl_small)[i]):(c(0, xscl_small)[i+1]))
}

#temporal
timebands_small<-sort(unique(round(c(1, exp(seq(log(1), log(max(timebands)), length=40)), max(timebands)))))
timebands_small<-timebands_small[timebands_small<=max(timebands)]

fuzzmat_time<-matrix(ncol=length(timebands_small), nrow=length(timebands), data=0)
for(i in 1:ncol(fuzzmat_time)) {
  tmp<-which((timebands>c(0,timebands_small)[i]) & (timebands<=timebands_small[i]))
  fuzzmat_time[tmp,i]<-1/length(tmp)
}
tmp<-which(colSums(fuzzmat_time)>0)
fuzzmat_time<-fuzzmat_time[,tmp]
timebands_small<-timebands_small[tmp]


arrayout<-array(dim=c(length(timebands_small), length(xscl_small), 3, niter))
tmparray<-array(dim=c(length(timebands), 16*mnplot, 3))
  
r0fun<-function(vec, s0) {
  #get mean for each scale
  tmp<-unname(log((cumsum(vec)/(1:length(vec))+0.0001)/s0))
  
  tmp[!is.finite(tmp)]<-NA
  tmp
}

if(FALSE) {
  save(list = "xscl_small", file = "output/tmp/xscl_small.rda")
  save(list = "timebands_small", file = "output/tmp/timebands_small.rda")
  save(list = "state0", file = "output/tmp/state0.rda")
}


for(iter in 1:niter) {
  tmparray[]<-NA
  
  for(i in 1:length(timebands)) {
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
    
    tmparray[i,1:length(psinclude),1]<-r0fun(sd.pc$Annual[c(psinclude[,1:j])], state0[1])/timebands[i]
    tmparray[i,1:length(psinclude),2]<-r0fun(sd.pc$C3_Perennial[c(psinclude[,1:j])], state0[2])/timebands[i]
    tmparray[i,1:length(psinclude),3]<-r0fun(sd.pc$C4S_Perennial[c(psinclude[,1:j])], state0[3])/timebands[i]
  }
  #fuzz and save
  for(i in 1:dim(tmparray)[3]) {
    tmp<-tmparray[,,i]
    tmp[is.na(tmp)]<-0
    tmp2<-tmp%*%fuzzmat_space
    tmp2<-tmp2*timebands
   
    tmp3<-t(t(tmp2)%*%fuzzmat_time)
    tmp3[tmp3==0]<-NA

    arrayout[,,i,iter]<-tmp3
  }
  
  
  #if(iter/10==floor(iter/10)) {
  #  print(round(iter/niter,2))
  #}
}

#save output
if(length(tstarg)==0) {tstarg<-runif(1)}
fname<-paste("/work/clarka/e014_emp_r0out_", tstarg[1], ".rda", sep="")
save(list=c("arrayout"), file = fname)


##############################
# Run simulation...
##############################
#set up for runs
niterations<-niter
scalelst<-c(1,2,3,5,7,10,20,30,50,70,100,200,300,400)/c(100^2)
radlst<-Inf

#set up simulations
gridout<-makegrid(xlng = 100, ylng = 100)

simtime<-81 #time spans for equilibria dectection

clst_meta = c(0.1,0.3,0.4)
mlst_meta = c(0.02, 0.02, 0.02)

array_sim<-array(dim=c(length(timebands_small), length(scalelst), 3, niter))
array_sim_tmp<-array(dim=c(simtime-1, length(scalelst), 3))

s0<-unname(state0[3:1])
s0_ind<-round(s0*prod(gridout$lng))
s0_ind_back<-s0_ind/prod(gridout$lng)

#averaging matrix
fuzzmat_time2<-matrix(ncol=length(timebands_small), nrow=simtime-1, data=0)
for(i in 1:ncol(fuzzmat_time2)) {
  tmp<-which((timebands>c(0,timebands_small)[i]) & (timebands<=timebands_small[i]))
  fuzzmat_time2[tmp,i]<-1/length(tmp)
}
tmp<-which(colSums(fuzzmat_time2)>0)
fuzzmat_time2<-fuzzmat_time2[,tmp]

for(i in 1:niter) {
  population<-populate(gridout, nlst = s0_ind,
                       clst = clst_meta, radlst = Inf, mlst = mlst_meta)
  for(j in 1:length(scalelst)) {
    grid_sub<-grid_subset(gridout, size = scalelst[j])
    
    out<-run_metapopulation(tmax=simtime, nsteps = simtime, gridout, population, talktime = 0, sites_sub = grid_sub$sites)
    ngrid<-length(grid_sub$sites)
    r0<-t(log(t(out$output_spatial[-1,-1]/ngrid+dmin)/(s0_ind_back)))[,c(3:1)]
    
    array_sim_tmp[,j,]<-r0
  }
  #if(i/10==floor(i/10)) {
  #  print(round(i/niter,2))
  #}
  
  for(j in 1:dim(array_sim_tmp)[3]) {
    tmp<-array_sim_tmp[,,j]
    tmp[is.na(tmp)]<-0
    array_sim[,,j,i]<-t(t(tmp)%*%fuzzmat_time2)
    array_sim[,,j,i][array_sim[,,j,i]==0]<-NA
  }
}

#save output
fname<-paste("/work/clarka/e014_sim_r0out_", tstarg[1], ".rda", sep="")
save(list=c("array_sim"), file = fname)


