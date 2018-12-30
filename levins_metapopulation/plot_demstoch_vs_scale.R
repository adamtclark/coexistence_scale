#!/usr/bin/env Rscript
#error
rm(list=ls())
setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

#load scripts
source("run_metapopulation_wrapper.R")

#set up for runs
niterations<-10   #CHANGE TO ALTER NUMBER OF ITERATIONS
scalelst<-c(0.005, 0.01, 0.05, 0.1, 0.5, 0.75, 1)

radlst<-Inf

#set up simulations
xfac<-3

tmax<-300 #timeseries length
burnin<-200 #burning for growth rate when rare method

clst_meta = c(0.15, 0.35)*xfac
mlst_meta = rep(0.1, length(clst_meta))*xfac
gridout<-makegrid(xlng = 100, ylng = 100)
population_meta<-populate(gridout, nlst = floor(getceq(clst_meta, mlst_meta)*prod(gridout$lng)),
                          clst = clst_meta, radlst = Inf, mlst = mlst_meta)

varmat1<-matrix(nrow=niterations, ncol=length(scalelst))
varmat2<-matrix(nrow=niterations, ncol=length(scalelst))

for(i in 1:niterations) {
  out<-run_metapopulation(tmax=tmax, nsteps = tmax, gridout, population_meta, talktime = 0, runtype="metapopulation")
  for(j in 1:length(scalelst)) {
    grid_sub<-grid_subset(gridout, size = scalelst[j])
    out2<-rerunrun_metapopulation(out=out, tmax=50, runtype="metapopulation", sites_sub = grid_sub$sites, talktime = 0)
    
    tmp<-apply(out2$output_spatial[-c(1:burnin),-1], 2, function(x) sd(x)/mean(x))
    
    varmat1[i,j]<-tmp[1]
    varmat2[i,j]<-tmp[2]
  }
  print(i/niterations)
}  
  
qtmat1<-t(apply(varmat1, 2, function(x) quantile(x, pnorm(-1:1), na.rm=T)))
qtmat2<-t(apply(varmat2, 2, function(x) quantile(x, pnorm(-1:1), na.rm=T)))

pdf("figures/SUP_FIGURE_demstoch_vs_scale.pdf", colormodel = "cmyk", width=5, height=4)
  par(mar=c(4,4,2,2))
  matplot(scalelst*(100^2), qtmat1, lty=c(2,1,2), lwd=2, type="l", col=1, xlab="Grid Size", ylab="CV(N)", ylim=c(0, 0.25), axes=F, log="x")
  matlines(scalelst*(100^2), qtmat2, lty=c(2,1,2), lwd=2, col=2)
  legend("topright", c("Species 1", "Species 2"), col=c(1,2), lty=1, lwd=2, bty="n")
  
  axis(1, las=2)
  axis(2, las=2)
  box()
dev.off()  



