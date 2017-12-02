error
rm(list=ls())
setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

#TODO:
#2. Neet to think upbetter "return to equilibrium" method?

#3. Set up wrapper functions for spatial subsetting
#-. FIRST - update wrapper for rerun
#-. THEN - update neutral c code
#-. THEN - build spatial output into other functions

#4. Think about PSF model to use
#-. Needs both intransitivity etc. and cycles that drop off-cycle with disturbance.
#-. Maybe one example can be rock-paper-scissors?

#load functions
require(rEDM)
source("run_metapopulation_wrapper.R")



############################################################
# "Global" run
############################################################
par(mfcol=c(4,2), mar=c(4,4,2,2))
set.seed(171201)

##### Try Tilman metapopulation model
gridout<-makegrid(xlng = 50, ylng = 50)
clst_meta = c(0.15, 0.3, 0.8)
population<-populate(gridout, nlst = floor(getceq(clst_meta)*prod(gridout$lng)), clst = clst_meta, radlst = Inf)
out_meta<-run_metapopulation(tmax=1000, gridout = gridout, population = population, talktime = 0)

#getceq(clst_meta)
plot_metapop(out_meta)
#plot_map(out_meta, gridout)

getEmeta<-getE(out_meta, Elst = 2:10)
E_meta<-getEmeta$Eout

#Estimate Eigenvalue
#eig_meta<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation", useeq=getceq(clst_meta), replace_perturb = 1)
#eig_meta<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation", replace_perturb = 1)
#eig_meta<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation", E=E_meta, replace_perturb = 1)

#Estimate increase when rare
r0_meta<-estimate_rarereturn(out_meta, simtime=20, burnin=100, runtype="metapopulation")

#Estimate invariability
invar_out<-estimate_invar(out_meta, E=E_meta, burnin=0, doplot=TRUE)



##### Try Hubbell NZNS neutral model
clst_neut<-rep(0.5, 3)
population<-populate(gridout, nlst = round(rep(unique(abs(getceq(clst_neut)))/length(clst_neut), length(clst_neut))*prod(gridout$lng)),
                     clst = clst_neut, radlst = Inf)
out_neut<-run_metapopulation(tmax=1000, gridout = gridout, population = population, talktime = 0, runtype = "neutral")
plot_metapop(out_neut)

getEneut<-getE(out_neut, Elst = 2:10)
E_neut<-getEneut$Eout

#Estimate Eigenvalue
#eig_neut<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral", useeq = rep(unique(abs(getceq(clst_neut)))/length(clst_neut), length(clst_neut)), replace_perturb = 1)
eig_neut<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral", replace_perturb = 1)
#eig_neut<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral", E=E_neut, replace_perturb = 1)

#Estimate increase when rare
r0_neut<-estimate_rarereturn(out_neut, simtime=20, burnin=100, runtype="neutral")

#Estimate invariability
invar_neut<-estimate_invar(out_neut, E=E_neut, burnin=0, doplot=TRUE)



############################################################
# Spatial subsetting
############################################################
out_meta<-run_metapopulation(tmax=1000, gridout = gridout, population = population, talktime = 0, runtype = "metapopulation_spatial", sites_sub = 1:100)



population<-populate(gridout, nlst = floor(getceq(clst_meta)*prod(gridout$lng)), clst = clst_meta, radlst = Inf)

out_meta_sp<-run_metapopulation(tmax=1, gridout = gridout, population = population, talktime = 0)

for(i in 1:100) {
  out_meta_sp<-rerunrun_metapopulation(out=out_meta_sp, tmax=1, talktime = 0)
}






