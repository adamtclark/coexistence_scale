#effect of column on row
#standard levins
intmat_lev<-rbind(c(1,0,0,0),
                  c(1,1,0,0),
                  c(1,1,1,0),
                  c(1,1,1,1))

#rps
intmat_rps<-rbind(c(1,0,0,1),
                  c(1,1,0,0),
                  c(0,1,1,0),
                  c(0,0,1,1))

#make grid and population
gridout<-makegrid(xlng = 100, ylng = 100)
population<-populate(gridout, nlst = c(0.2, 0.2, 0.2, 0.2), clst = rep(0.3, 4), mlst = rep(0.1, 4), radlst = Inf)


out0<-run_metapopulation(tmax=1000, nsteps = 1000, gridout, population, talktime = 0, runtype="metapopulation", compmat = intmat_lev)
plot_metapop(out0)

out1<-run_metapopulation(tmax=1000, nsteps = 1000, gridout, population, talktime = 0, runtype="rps", compmat = intmat_lev)
plot_metapop(out1)

out2<-run_metapopulation(tmax=1000, nsteps = 1000, gridout, population, talktime = 0, runtype="rps", compmat = intmat_rps)
plot_metapop(out2)


population1<-populate(gridout, nlst = c(0, 0.2, 0.2, 0.2), clst = rep(0.3, 4), mlst = rep(0.1, 4), radlst = Inf)
out3<-run_metapopulation(tmax=1000, nsteps = 1000, gridout, population=population1, talktime = 0, runtype="rps", compmat = intmat_rps)
plot_metapop(out3)


population2<-populate(gridout, nlst = (out3$output[nrow(out3$output),2:5]+c(500, 0, 0, 0))*(1/out$plotdata$ngrid), clst = rep(0.3, 4), mlst = rep(0.1, 4), radlst = Inf)
out4<-run_metapopulation(tmax=1000, nsteps = 1000, gridout, population=population2, talktime = 0, runtype="rps", compmat = intmat_rps)
plot_metapop(out4)



#Check rerun
out<-run_metapopulation(tmax=1000, nsteps = 1000, gridout, population, talktime = 0, runtype="rps", compmat = intmat_rps)
plot_metapop(out)

#rerun simulation
out2<-rerunrun_metapopulation(out=out, tmax=1000, talktime = 0, runtype="rps")
plot_metapop(out2)


