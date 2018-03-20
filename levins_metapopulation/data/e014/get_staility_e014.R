#!/usr/bin/env Rscript
#error
rm(list=ls())

setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

#Load functions
source("run_metapopulation_wrapper.R")
source("util/filled.contour3.R")
source("util/figure_functions.R")
require(data.table)

##############################
# Run simulation...
##############################



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

#Add tm+1 to data
#e14_yrlst<-sort(unique(d.pc$Year))
#e54_yrlst<-sort(unique(d.bm$Year))
#index.bm.yr<-paste(d.bm$Field, d.bm$Transect, d.bm$Year)
#index.pc.yr<-paste(d.pc$Field, d.pc$Transect, d.pc$Plot, d.pc$Year)
#index.bm.yrp1<-paste(d.bm$Field, d.bm$Transect, d.bm$Year+1)
#index.pc.yrp1<-paste(d.pc$Field, d.pc$Transect, d.pc$Plot, e14_yrlst[match(d.pc$Year, e14_yrlst)+1])

#d.pc$Annual_tp1<-d.pc$Annual[match(index.pc.yrp1, index.pc.yr)]
#d.pc$C3_Perennial_tp1<-d.pc$C3_Perennial[match(index.pc.yrp1, index.pc.yr)]
#d.pc$C4S_Perennial_tp1<-d.pc$C4S_Perennial[match(index.pc.yrp1, index.pc.yr)]
#d.pc$Forb_Perennial_tp1<-d.pc$Forb_Perennial[match(index.pc.yrp1, index.pc.yr)]
#d.pc$Leg_Perennial_tp1<-d.pc$Leg_Perennial[match(index.pc.yrp1, index.pc.yr)]

#d.bm$Annual_tp1<-d.bm$Annual[match(index.bm.yrp1, index.bm.yr)]
#d.bm$C3_Perennial_tp1<-d.bm$C3_Perennial[match(index.bm.yrp1, index.bm.yr)]
#d.bm$C4S_Perennial_tp1<-d.bm$C4S_Perennial[match(index.bm.yrp1, index.bm.yr)]
#d.bm$Forb_Perennial_tp1<-d.bm$Forb_Perennial[match(index.bm.yrp1, index.bm.yr)]
#d.bm$Leg_Perennial_tp1<-d.bm$Leg_Perennial[match(index.bm.yrp1, index.bm.yr)]

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
niter<-100

arrayout<-array(dim=c(length(timebands), 16*mnplot, 3, niter))

r0fun<-function(vec, s0) {
  #tmp<-unname(log((mean(vec,na.rm=T)+0.0001)/s0))
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

#Shrink matrix

#Get quantiles for each...
array_quant<-apply(arrayout, 1:3, function(x) quantile(x, c(0.5), na.rm=T))
xscl<-(1:dim(array_quant)[2])/2


#Plot output
m<-(matrix(nrow=3, 1:6))
m<-cbind(m, 7)
layout(m, widths=c(1,0.5))
par(mar=c(2,3,1,0), oma=c(2.5,2,2,3))

plotout<-plot_cont_emp(arrayout=array_quant, xscalslst=log10(xscl), xlst=log10(timebands), nlevels=10, logx=TRUE, logy=TRUE, logz=FALSE, nstart=1, ofs1=c(0, 0))
  

