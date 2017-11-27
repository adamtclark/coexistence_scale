error
rm(list=ls())
setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")


#load functions
source("run_metapopulation_wrapper.R")

gridout<-makegrid(xlng = 100, ylng = 100)
#population<-populate(gridout, clst = c(0.15, 0.3, 0.8, 3, 15), radlst = Inf)
population<-populate(gridout, clst = c(0.15, 0.3), radlst = 2)
getceq(population)

out<-run_metapopulation(tmax=500, nsteps = 1000, gridout, population, talktime = 10)

plot_metapop(out)
plot_map(out, gridout)
out1<-out



population1<-rewrap_pop(out1, population)
out1<-run_metapopulation(tmax=500, nsteps = 100, gridout, population1, talktime = 0)
plot_metapop(out1)
plot_map(out1, gridout)
