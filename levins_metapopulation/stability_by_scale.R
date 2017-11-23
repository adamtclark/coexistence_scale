error
rm(list=ls())
setwd("~/Dropbox/Projects/032_Coexistence_table/src/levins_metapopulation/")


#load functions
source("run_metapopulation_wrapper.R")

gridout<-makegrid(xlng = 100, ylng = 100)
population<-populate(gridout, nlst = rep(round(0.1*prod(gridout$lng)), 4), clst = c(0.15, 1, 10, 200), mlst = rep(0.1, 4), radlst = Inf)

out<-run_metapopulation(tmax=200, nsteps = 1000, gridout, population, talktime = 0)

plot_metapop(out)



plot(gridout$xpos, gridout$ypos, col=out$full$speciesid, pch=16)




out1<-run_metapopulation(tmax=10, nsteps = 10, gridout, population, talktime = 0)
