error
rm(list=ls())
setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

require(rEDM)
source("run_metapopulation_wrapper.R")

##### Try disturbance model
gridout<-makegrid(xlng = 100, ylng = 100)
xfac<-5
ptb<-0.2

clst_dist= c(0.145, 0.16, 0.18, 0.21)*xfac
mlst_dist = rep(0.1, length(clst_dist))*xfac
distlst<-c(0.95, 0.7, 0.2, 0)
getceq(clst_dist, mlst_dist)


population_dist<-populate(gridout, nlst = rep(floor(prod(gridout$lng)/length(clst_dist)*0.8), length(clst_dist)),
                          clst = clst_dist, radlst = Inf, mlst = mlst_dist)
getceq(clst_dist, mlst_dist)
grid_sub<-grid_subset(gridout, size = 0.2)

tmax=1000; nsteps=tmax; gridout = gridout; population = population_dist; talktime = 0; runtype = "disturbance";
prt=c(0.95, 0.5, 0.2, 0); prtfrq=20; sites_sub = grid_sub$sites




out_dist<-run_metapopulation(tmax=200, gridout = gridout, population = population_dist, talktime = 0, runtype = "disturbance",
                             prt=distlst, prtfrq=20)
plot_metapop(out_dist, dotot = FALSE, ylim=c(0, 0.4))

out_rep<-rerunrun_metapopulation(out=out_dist, tmax=200, talktime = 0, runtype = "disturbance", perturb = distlst, prt=distlst, prtfrq=20)
plot_metapop(out_rep, dotot = FALSE, ylim=c(0, 0.4))


out_meta<-run_metapopulation(tmax=1000, gridout = gridout, population = population_dist, talktime = 0, runtype = "metapopulation")
plot_metapop(out_meta, dotot = FALSE, ylim=c(0, 0.4))

out_rep<-rerunrun_metapopulation(out=out_meta, tmax=200, talktime = 0, runtype = "metapopulation", prt=distlst, prtfrq=20)
plot_metapop(out_rep, dotot = FALSE, ylim=c(0, 0.4))





runtype = "disturbance_spatial"

out_dist<-run_metapopulation(tmax=1000, gridout = gridout, population = population_dist, talktime = 0, runtype = "disturbance_spatial",
                             prt=distlst, prtfrq=20, sites_sub = grid_sub$sites)
plot_metapop(out_dist, dotot = FALSE, ylim=c(0, 0.4))
plot_metapop(out_dist, dotot = FALSE, ylim=c(0, 0.4), sites = 1)

out_rep<-rerunrun_metapopulation(out=out_dist, tmax=200, talktime = 0, runtype = "disturbance_spatial", perturb = distlst, prt=distlst, prtfrq=20, sites_sub = grid_sub$sites, perturbsites = 1:out_dist$plotdata$ngrid)
plot_metapop(out_rep, dotot = FALSE, ylim=c(0, 0.4))
plot_metapop(out_rep, dotot = FALSE, ylim=c(0, 0.4), sites = 1)

#out<-out_dist
#perturb=rep(0, out$full$pnsp); perturbsites=1:out$plotdata$ngrid; addn=0; addsites=perturbsites; replace_perturb=0




out_meta<-run_metapopulation(tmax=1000, gridout = gridout, population = population_dist, talktime = 0, runtype = "metapopulation_spatial", sites_sub = grid_sub$sites)

plot_metapop(out_meta, dotot = FALSE, ylim=c(0, 0.4))
plot_metapop(out_meta, dotot = FALSE, ylim=c(0, 0.4), sites = 1)

out_rep<-rerunrun_metapopulation(out=out_meta, tmax=200, talktime = 0, runtype = "metapopulation_spatial", sites_sub = grid_sub$sites)
plot_metapop(out_rep, dotot = FALSE, ylim=c(0, 0.4))
plot_metapop(out_rep, dotot = FALSE, ylim=c(0, 0.4), sites = 1)



