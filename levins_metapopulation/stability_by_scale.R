error
rm(list=ls())
setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

#TODO:
#1. Update destruction algorithm to allow limits for individual species

#2. Neet to think up better "return to equilibrium" method?
#-. Think about increasing size to decrease variability - show improvement in equilibrium model
#-. Note that perturbation needs to be larger than stochastic variation - when large, this becomes similar to "increase when rare".

#4. Add Stan's PSF model
#-. Needs both intransitivity etc. and cycles that drop off-cycle with disturbance.

#5. One example of rock-paper-scissors



#load functions
require(rEDM)
source("run_metapopulation_wrapper.R")



############################################################
# "Global" run
############################################################
pdf("figures/stability_by_scale.pdf", width=6, height=8, colormodel = "cmyk")

par(mfcol=c(5,2), mar=c(4,4,2,2), oma=c(0.1,0.1,2,0))
set.seed(171205)

##### Try Tilman metapopulation model
gridout<-makegrid(xlng = 100, ylng = 100)
xfac<-5
ptb<-0.2

clst_meta = c(0.15, 0.3, 0.8, 3)*xfac
mlst_meta = rep(0.1, length(clst_meta))*xfac
#getceq(clst_meta, mlst_meta)
population_meta<-populate(gridout, nlst = floor(getceq(clst_meta, mlst_meta)*prod(gridout$lng)),
                          clst = clst_meta, radlst = Inf, mlst = mlst_meta)
out_meta<-run_metapopulation(tmax=1000, gridout = gridout, population = population_meta, talktime = 0)

#getceq(clst_meta)

plot_metapop(out_meta)
mtext("Levins Model", side=3)

#plot_map(out_meta, gridout)

getEmeta<-getE(out_meta, Elst = 2:10)
E_meta<-getEmeta$Eout
E_meta[E_meta<4]<-4

#Estimate Eigenvalue
#eig_meta<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation", useeq=getceq(clst_meta), replace_perturb = 1)

eig_meta1<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation", replace_perturb = 1, talktime=0, prtb=ptb)
mtext("replace", side=3, cex=0.8)

eig_meta2<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation", replace_perturb = 0, talktime=0, prtb=ptb)
mtext("no replace", side=3, cex=0.8)

#eig_meta<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation", E=E_meta, replace_perturb = 1)

#Estimate increase when rare
r0_meta<-estimate_rarereturn(out_meta, simtime=100, burnin=100, runtype="metapopulation")

#Estimate invariability
invar_out<-estimate_invar(out_meta, E=E_meta, burnin=0, doplot=TRUE)



##### Try Hubbell N-ZNS neutral model
clst_neut<-rep(0.5, 4)*xfac
mlst_neut<-rep(0.1, length(clst_neut))*xfac
population_neut<-populate(gridout, nlst = round(rep(unique(abs(getceq(clst_neut, mlst_neut)))/length(clst_neut), length(clst_neut))*prod(gridout$lng)),
                     clst = clst_neut, radlst = Inf, mlst = mlst_neut)
out_neut<-run_metapopulation(tmax=1000, gridout = gridout, population = population_neut, talktime = 0, runtype = "neutral")

plot_metapop(out_neut)
mtext("Neutral Model", side=3)

getEneut<-getE(out_neut, Elst = 2:10)
E_neut<-getEneut$Eout
E_neut[E_neut<4]<-4

#Estimate Eigenvalue
#eig_neut<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral", useeq = rep(unique(abs(getceq(clst_neut)))/length(clst_neut), length(clst_neut)), replace_perturb = 1)

eig_neut1<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral", replace_perturb = 1, talktime=0, prtb=ptb)
mtext("replace", side=3, cex=0.8)
eig_neut2<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral", replace_perturb = 0, talktime=0, prtb=ptb)
mtext("no replace", side=3, cex=0.8)

#eig_neut<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral", E=E_neut, replace_perturb = 1)

#Estimate increase when rare
r0_neut<-estimate_rarereturn(out_neut, simtime=100, burnin=100, runtype="neutral")

#Estimate invariability
invar_neut<-estimate_invar(out_neut, E=E_neut, burnin=0, doplot=TRUE)

mtext("Global Spatial Scale", side = 3, outer = TRUE)

############################################################
# Spatial subsetting
############################################################
#META
grid_sub<-grid_subset(gridout, size = 0.2)
ptb<-0.2

out_meta<-run_metapopulation(tmax=1000, gridout = gridout, population = population_meta, talktime = 0, runtype = "metapopulation_spatial", sites_sub = grid_sub$sites)
#plot_metapop(out_meta)

plot_metapop(out_meta, sites=1)
mtext("Levins Model", side=3)
#plot_map(out_meta, gridout, grid_sub=grid_sub)

getEmeta<-getE(out_meta, Elst = 2:10, sites_sub = grid_sub$sites)
E_meta<-getEmeta$Eout
E_meta[E_meta<4]<-4

#eig_meta<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation", replace_perturb = 0, talktime=0, prtb=ptb)
eig_meta1<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation_spatial", replace_perturb = 1, talktime=0, prtb=ptb, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)
mtext("replace", side=3, cex=0.8)

eig_meta2<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation_spatial", replace_perturb = 0, talktime=0, prtb=ptb, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)
mtext("no replace", side=3, cex=0.8)

#r0_meta<-estimate_rarereturn(out_meta, simtime=100, burnin=100, runtype="metapopulation")
r0_meta<-estimate_rarereturn(out_meta, simtime=100, burnin=100, runtype="metapopulation_spatial", perturbsites = grid_sub$sites, sites_sub = grid_sub$sites)

invar_meta<-estimate_invar(out_meta, E=E_meta, burnin=0, doplot=TRUE, sites_sub = grid_sub$sites)

#NEUTRAL
out_neut<-run_metapopulation(tmax=1000, gridout = gridout, population = population_neut, talktime = 0, runtype = "neutral_spatial", sites_sub = grid_sub$sites)
#plot_metapop(out_neut)

plot_metapop(out_neut, sites=1)
mtext("Neutral Model", side=3)
#plot_map(out_neut, gridout, grid_sub=grid_sub)

getEneut<-getE(out_neut, Elst = 2:10, sites_sub = grid_sub$sites)
E_neut<-getEneut$Eout
E_neut[E_neut<4]<-4

#eig_neut<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral", replace_perturb = 1, talktime=0, prtb=ptb)
eig_neut1<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral_spatial", replace_perturb = 1, talktime=0, prtb=ptb, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)
mtext("replace", side=3, cex=0.8)
eig_neut2<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral_spatial", replace_perturb = 0, talktime=0, prtb=ptb, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)
mtext("no replace", side=3, cex=0.8)

#r0_neut<-estimate_rarereturn(out_neut, simtime=100, burnin=100, runtype="neutral")
r0_neut<-estimate_rarereturn(out_neut, simtime=100, burnin=100, runtype="neutral_spatial", perturbsites = grid_sub$sites, sites_sub = grid_sub$sites)

invar_neut<-estimate_invar(out_neut, E=E_neut, burnin=0, doplot=TRUE, sites_sub = grid_sub$sites)

mtext("20% Spatial Subset", side = 3, outer = TRUE)


dev.off()