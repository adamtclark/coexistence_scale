#int *psp1dis, int *psp2dis, int *psp1fb, int *psp2fb, int *psp1m, int *psp2m, int *pscenario,
#int *pinitabund, float *pseed, int *pedge, int *pdim, int *ptmax, double outmat[], double outmap0[], double outmap[])  {
#  //sp1dis, sp2dis;		species dispersal: 0=local, 1=global
#  //sp1fb, sp2fb;	 species feedback strength: 0=none, 100=overwhelmingly strong; sign determines sign of feedback
#  //sp1m, sp2m;		 species percent mortality: 0=immortal, <100=perennial, 100=annual
#  //scenario;			 initial soil feedback state: 0=neutral, 1=exotic (spp1), 2=native (spp2)
#  //initabund; 		 initial percent natives
#  //seed; 				   ratio of exotic to native seed production: 1= default (values: 1000, 100, 10, 1, 0.1, 0.01, 0.001)
#  //edge;					 0=absorbing boundaries, 1=one edge all exotics (spp1)
#  //dim;					   dimension of matrix (not counting the edges)
#  //tmax;					 time steps
#  //pr_nocol;       reduction in probability of colonization event - allows for empty cells
#  //stepsize;       step size in by which changes in soil occur (in percent)

#Paramter settings...
sp1dis=1; sp2dis=1          #dispersal
sp1fb=-80; sp2fb=-80        #feedback
sp1m=5; sp2m=5              #stochastic mortality
scenario=0                  #begin with neutral soils
initabund=50                #initial native abundance (sp1)
seed=1.5                    #ratio of seed production (exotic:native)
edge=0                      #absorbing edge conditions
dim=100                     #grid edge size
tmax=1000                   #maximum time
outmat=numeric((tmax+1)*3)  #matrix for storing species abundances
outmap0=numeric((dim+2)^2)  #matrix for storing initial conditions
outmap=outmap0              #matrix for storing end conditions
pr_nocol=0.5                #reduction in probability of colonization event - allows for empty cells
stepsize=1                  #step size in by which changes in soil occur (in percent)

setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

#compile and load C code
system("R CMD SHLIB tmp/sis_R.c")
if(is.loaded("sis_R")) {
  dyn.unload("tmp/sis_R.so")
}

dyn.load("tmp/sis_R.so")

#Run C code
out<-.C("sis_R",
        psp1dis=as.integer(sp1dis), psp2dis=as.integer(sp2dis),
        ps1fb=as.integer(sp1fb), psp2fb=as.integer(sp2fb), 
        psp1m=as.integer(sp1m), psp2m=as.integer(sp2m),
        pscenario=as.integer(scenario), 
        pinitabund=as.integer(initabund), pseed=as.double(seed),
        pedge=as.integer(edge), pdim=as.integer(dim), ptmax=as.integer(tmax),
        ppr_nocol=as.double(pr_nocol), pstepsize=as.integer(stepsize),
        outmat=as.integer(outmat), outmap0=as.integer(outmap0), outmap=as.integer(outmap))

#Store abundance v.s time output in a matrix
m<-matrix(out$outmat, ncol=3)
#Exclude empty cells or invaded cells from the periphery
if(edge==0) {
  m[,1]<-m[,1]-(4*(dim+2)-4)
} else {
  m[,1]<-m[,1]-(3*(dim+2)-2)
  m[,2]<-m[,2]-(dim)
}

#Plot result
matplot(0:(tmax), cbind(m[,-1], rowSums(m[,-1]))/(dim^2), type="l", lty=1, lwd=2, col=c(2,3,1), xlab="time", ylab="p", ylim=c(0,1))
matplot(0:(tmax), cbind(m_sub[,-1], rowSums(m_sub[,-1]))/(length(grid_sub$sites)), type="l", lty=1, lwd=2, col=c(2,3,1), xlab="time", ylab="p", ylim=c(0,1))

abline(h=c(0, 1), v=0, lty=3)



#Sanity check -- 
mean(rowSums(m[,-1])/(dim^2)) # if pr_nocol=0, then this should = 1
mean(rowSums(m)/(dim^2))      # this should always = 1

#Show starting and ending states
matrix(out$outmap0, nrow=(dim+2))

tmp<-matrix(out$outmap, nrow=(dim+2))
plot(row(tmp), col(tmp), col=tmp, pch=16, cex=0.5)
