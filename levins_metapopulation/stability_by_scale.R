error
rm(list=ls())
setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

#TODO:
#2. Neet to think upbetter "return to equilibrium" method?
#-. Think about increasing size to decrease variability - show improvement in equilibrium model
#-. Note that perturbation needs to be larger than stochastic variation - when large, this becomes similar to "increase when rare".

#3. Set up wrapper functions for spatial subsetting
#-. FIRST - update wrapper for rerun
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
par(mfcol=c(5,2), mar=c(4,4,2,2))
set.seed(171201)

##### Try Tilman metapopulation model
gridout<-makegrid(xlng = 50, ylng = 50)
xfac<-5
ptb<-0.2

clst_meta = c(0.15, 0.3, 0.8)*xfac
mlst_meta = rep(0.1, 3)*xfac
#getceq(clst_meta, mlst_meta)
population_meta<-populate(gridout, nlst = floor(getceq(clst_meta, mlst_meta)*prod(gridout$lng)), clst = clst_meta, radlst = 5, mlst = mlst_meta)
out_meta<-run_metapopulation(tmax=100, gridout = gridout, population = population_meta, talktime = 0)

#getceq(clst_meta)
plot_metapop(out_meta)
#plot_map(out_meta, gridout)

getEmeta<-getE(out_meta, Elst = 2:10)
E_meta<-getEmeta$Eout
E_meta[E_meta<5]<-5

#Estimate Eigenvalue
#eig_meta<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation", useeq=getceq(clst_meta), replace_perturb = 1)

eig_meta<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation", replace_perturb = 1, talktime=0, prtb=ptb)
mtext("replace", side=3)

eig_meta<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation", replace_perturb = 0, talktime=0, prtb=ptb)
mtext("no replace", side=3)

#eig_meta<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation", E=E_meta, replace_perturb = 1)

#Estimate increase when rare
r0_meta<-estimate_rarereturn(out_meta, simtime=100, burnin=100, runtype="metapopulation")

#Estimate invariability
invar_out<-estimate_invar(out_meta, E=E_meta, burnin=0, doplot=TRUE)



##### Try Hubbell NZNS neutral model
clst_neut<-rep(0.5, 3)*xfac
mlst_neut<-rep(0.1, 3)*xfac
population_neut<-populate(gridout, nlst = round(rep(unique(abs(getceq(clst_neut, mlst_neut)))/length(clst_neut), length(clst_neut))*prod(gridout$lng)),
                     clst = clst_neut, radlst = 5, mlst = mlst_neut)
out_neut<-run_metapopulation(tmax=100, gridout = gridout, population = population_neut, talktime = 0, runtype = "neutral")
plot_metapop(out_neut)

getEneut<-getE(out_neut, Elst = 2:10)
E_neut<-getEneut$Eout
E_neut[E_neut<5]<-5

#Estimate Eigenvalue
#eig_neut<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral", useeq = rep(unique(abs(getceq(clst_neut)))/length(clst_neut), length(clst_neut)), replace_perturb = 1)

eig_neut<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral", replace_perturb = 1, talktime=0, prtb=ptb)
mtext("replace", side=3)
eig_neut<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral", replace_perturb = 0, talktime=0, prtb=ptb)
mtext("no replace", side=3)

#eig_neut<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral", E=E_neut, replace_perturb = 1)

#Estimate increase when rare
r0_neut<-estimate_rarereturn(out_neut, simtime=100, burnin=100, runtype="neutral")

#Estimate invariability
invar_neut<-estimate_invar(out_neut, E=E_neut, burnin=0, doplot=TRUE)



############################################################
# Spatial subsetting
############################################################
grid_sub<-grid_subset(gridout, size = 0.05)

out_meta<-run_metapopulation(tmax=1000, gridout = gridout, population = population_meta, talktime = 0, runtype = "metapopulation_spatial", sites_sub = grid_sub$sites)
plot_metapop(out_meta)
plot_metapop(out_meta, sites=1)
plot_map(out_meta, gridout, grid_sub=grid_sub)

out_neut<-run_metapopulation(tmax=1000, gridout = gridout, population = population_neut, talktime = 0, runtype = "neutral_spatial", sites_sub = grid_sub$sites)
plot_metapop(out_neut)
plot_metapop(out_neut, sites=1)
plot_map(out_neut, gridout, grid_sub=grid_sub)






