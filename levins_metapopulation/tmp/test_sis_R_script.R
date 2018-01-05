#int *psp1dis, int *psp2dis, int *psp1fb, int *psp2fb, int *psp1m, int *psp2m, int *pscenario,
#int *pinitabund, float *pseed, int *pedge, int *pdim, int *ptmax, double outmat[], double outmap0[], double outmap[])  {
#  //sp1dis, sp2dis;		species dispersal: 0=local, 1=global
#  //sp1fb, sp2fb;			species feedback strength: 0=none, 100=overwhelmingly strong
#  //sp1m, sp2m;			species percent mortality: 0=immortal, <100=perennial, 100=annual
#  //scenario;				initial soil feedback state: 0=neutral, 1=exotic (spp1), 2=native (spp2)
#  //initabund; 			initial percent natives
#  //seed; 				ratio of exotic to native seed production: 1= default (values: 1000, 100, 10, 1, 0.1, 0.01, 0.001)
#  //edge;					0=absorbing boundaries, 1=one edge all exotics (spp1)
#  //dim;					dimension of matrix (not counting the edges)
#  //tmax;					time steps

sp1dis=1; sp2dis=1    #global dispersal
sp1fb=20; sp2fb=20      #moderate feedback
sp1m=20; sp2m=20        #20% mortality
scenario=2              #neutral feedback
initabund=30            #initial native abundance (sp1)
seed=1                  #ratio of seed production
edge=0                  #edge conditions
dim=100                 #grid edge size
tmax=1000                #maximum time
outmat=numeric((tmax+1)*3)  #matrix for storing species abundances
outmap0=numeric((dim+2)^2)  #matrix for storing initial conditions
outmap=outmap0          #matrix for storing end conditions
pr_nocol=0.8              #reduction in probability of colonization event - allows for empty cells

setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")


system("R CMD SHLIB sis_R.c")
if(is.loaded("sis_R")) {
  dyn.unload("sis_R.so")
}

dyn.load("sis_R.so")


out<-.C("sis_R",
        psp1dis=as.integer(sp1dis), psp2dis=as.integer(sp2dis),
        ps1fb=as.integer(sp1fb), psp2fb=as.integer(sp2fb), 
        psp1m=as.integer(sp1m), psp2m=as.integer(sp2m),
        pscenario=as.integer(scenario), 
        pinitabund=as.integer(initabund), pseed=as.double(seed),
        pedge=as.integer(edge), pdim=as.integer(dim), ptmax=as.integer(tmax),
        ppr_nocol=as.double(pr_nocol),
        outmat=as.integer(outmat), outmap0=as.integer(outmap0), outmap=as.integer(outmap))

m<-matrix(out$outmat, ncol=3)
m[,1]<-m[,1]-(4*(dim+2)-4)
#m

matplot(0:(tmax), m[,-1]/(dim^2), type="l", lty=1, lwd=2, col=1:2, xlab="time", ylab="p", ylim=c(0,1))
abline(h=c(0, 1), v=0, lty=3)




rowSums(m[,-1])/(dim^2)
rowSums(m)/(dim^2)

matrix(out$outmap0, nrow=(dim+2))
matrix(out$outmap, nrow=(dim+2))

