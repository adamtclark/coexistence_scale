#!/usr/bin/env Rscript
#error
rm(list=ls())

setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

#Load functions
source("run_metapopulation_wrapper.R")
source("util/filled.contour3.R")
source("util/figure_functions.R")
source("util/plot_grid_functions.R")

require(RColorBrewer)

##############################
#Load
##############################
if(length(grep("array_sim_quant_small2.rda", dir("output/tmp")))==0) {
  #Collate data from HPC simulation
  
  flst<-dir("/work/clarka/")
  
  dir_sim<-flst[grep("e014_sim", flst)]
  dir_emp<-flst[grep("e014_emp", flst)]
  
  load(paste("/work/clarka/", dir_sim[1], sep=""))
  load(paste("/work/clarka/", dir_emp[1], sep=""))
  #arrayout; array_sim
  niter<-dim(arrayout)[4]*max(c(length(dir_sim), length(dir_emp)))
  
  array_sim_full_neutral<-array(dim=c(dim(array_sim)[1:3], niter))
  array_sim_full<-array(dim=c(dim(array_sim)[1:3], niter))
  arrayout_full<-array(dim=c(dim(arrayout)[1:3], niter))
  
  n<-1
  for(i in 1:length(dir_sim)) {
    load(paste("/work/clarka/", dir_sim[1], sep=""))
    nnew<-(n+dim(array_sim)[4]-1)
    array_sim_full[,,,n:nnew]<-array_sim
    array_sim_full_neutral[,,,n:nnew]<-array_sim_neut
    n<-nnew+1
  }
  
  n<-1
  for(i in 1:length(dir_emp)) {
    load(paste("/work/clarka/", dir_sim[1], sep=""))
    nnew<-(n+dim(arrayout)[4]-1)
    arrayout_full[,,,n:nnew]<-arrayout
    n<-nnew+1
  }
  
  array_sim_quant_small2_neutral<-apply(array_sim_full_neutral, 1:3, function(x) quantile(x, 0.5, na.rm=T))
  array_sim_quant_small2<-apply(array_sim_full, 1:3, function(x) quantile(x, 0.5, na.rm=T))
  array_quant_small2<-apply(arrayout_full, 1:3, function(x) quantile(x, 0.5, na.rm=T))
  
  save(list = c("array_sim_quant_small2", "array_sim_quant_small2_neutral"), file = "output/tmp/array_sim_quant_small2.rda")
  save(list = "array_quant_small2", file = "output/tmp/array_quant_small2.rda")
}


####################################
# Try to load likelihooods
####################################
tmp<-dir("output/tmp")
havelike<-length(grep("emp_densplot.rda", tmp))>0

if(havelike) {
  load("output/tmp/emp_densplot.rda")
} else {
  print("Note: run 'plot_matchcor_out_array_emp.R' first to generate and plot likelihoods.")
}





##############################
#Plot output
##############################
load("output/tmp/array_sim_quant_small2.rda")
load("output/tmp/array_quant_small2.rda")

load("output/tmp/xscl_small.rda")
load("output/tmp/timebands_small.rda")

scalelst<-c(1,2,4,6,9,12,20,30,49,72,100,196,306,400,812)/c(100^2)

svg("figures/FIGURE_e014_comparison.svg", width=7, height=6)
ofs<-c(0.24, -0.002)

m<-(matrix(nrow=3, 1:9))
m<-cbind(m, 10)
layout(m, widths=c(1,1,1,0.35))
par(mar=c(2,3,1,0), oma=c(2.5,2,2,3))
squse<-c(-8,-4,-2,-1,0,1,2,4,8)#c(-100, seq(-4,0, by=1), seq(1,4,by=1), 100)

if(!havelike) {
  plotout<-plot_cont_emp(arrayout=array_quant_small2, xscalslst=log10(xscl_small), xlst=log10(timebands_small), nlevels=10, logx=TRUE, logy=TRUE, logz=FALSE, nstart=1, sqlst=squse, ofs1=ofs, dolines = FALSE)
  plotout<-plot_cont_emp(arrayout=array_sim_quant_small2, xscalslst=log10(scalelst*100^2), xlst=log10(timebands_small), nlevels=10, logx=TRUE, logy=TRUE, logz=FALSE, nstart=4, sqlst=squse, ofs1=ofs, ylim=range(log10(xscl_small)), dolines = FALSE)
  plotout<-plot_cont_emp(arrayout=array_sim_quant_small2_neutral, xscalslst=log10(scalelst*100^2), xlst=log10(timebands_small), nlevels=10, logx=TRUE, logy=TRUE, logz=FALSE, nstart=7, sqlst=squse, ofs1=ofs, ylim=range(log10(xscl_small)), dolines = FALSE)
} else {
  plotout<-plot_cont_emp(arrayout=array_quant_small2, xscalslst=log10(xscl_small), xlst=log10(timebands_small), nlevels=10, logx=TRUE, logy=TRUE, logz=FALSE, nstart=1, sqlst=squse, ofs1=ofs, dolines = FALSE)
  plotout<-plot_cont_emp(arrayout=array_sim_quant_small2, xscalslst=log10(scalelst*100^2), xlst=log10(timebands_small), nlevels=10, logx=TRUE, logy=TRUE, logz=FALSE, nstart=4, sqlst=squse, ofs1=ofs, ylim=range(log10(xscl_small)), dolines = FALSE, plot_dens_cum=plot_dens_cum)
  plotout<-plot_cont_emp(arrayout=array_sim_quant_small2_neutral, xscalslst=log10(scalelst*100^2), xlst=log10(timebands_small), nlevels=10, logx=TRUE, logy=TRUE, logz=FALSE, nstart=7, sqlst=squse, ofs1=ofs, ylim=range(log10(xscl_small)), dolines = FALSE, plot_dens_cum=plot_dens_cum, revdenscum=TRUE)
}

par(mar=c(2,3,1,0))
sqplot<-c(-8,-4,-2,-1,0,1,2,4,8)#c(-5, seq(-4,0, by=1), seq(1,4,by=1), 5)
filled.legend(z=matrix(sqplot), levels=sqplot, col=plotout$collst2, key.axes = axis(4, at = sqplot, labels = sqplot, las=2))

mtext(text = "Annuals", side = 4, outer = TRUE, line = -4.2, adj = 0.905, cex=1.2)
mtext(text = "C3 Grasses", side = 4, outer = TRUE, line = -4.2, adj = 0.515, cex=1.2)
mtext(text = "C4 Grasses", side = 4, outer = TRUE, line = -4.2, adj = 0.1, cex=1.2)

mtext(text = expression(paste("empirical")), side = 3, outer = TRUE, line = 0, adj = 0.14, cex=1.2)
mtext(text = expression(paste("Levins-OF")), side = 3, outer = TRUE, line = 0, adj = 0.48, cex=1.2)
mtext(text = expression(paste("neutral-OF")), side = 3, outer = TRUE, line = 0, adj = 0.81, cex=1.2)

mtext(text = expression(paste("temporal extent, years")), side = 1, outer = TRUE, line = 1.3, cex=1.2, adj = 0.48)
mtext(text = expression(paste("spatial extent, number of plots")), side = 2, outer = TRUE, line = -0.1, cex=1.2, adj = 0.5)


mtext(text = expression(paste(r[0], t)), side = 3, outer = TRUE, line = -0.5, adj = 1, cex=1.5, xpd=NA)

dev.off()


#Check correlation
if(FALSE) {
  mps<-apply(t(t(scalelst*100^2)), 1, function(x) which.min(abs(x-xscl_small)))
  cor.test(c(array_quant_small2[,mps,1]), c(array_sim_quant_small2[,,1]), use = "complete")
  cor.test(c(array_quant_small2[,mps,2]), c(array_sim_quant_small2[,,2]), use = "complete")
  cor.test(c(array_quant_small2[,mps,3]), c(array_sim_quant_small2[,,3]), use = "complete")
  
  
  summary(lm(c(array_quant_small2[,mps,1])~c(array_sim_quant_small2[,,1])))
  summary(lm(c(array_quant_small2[,mps,2])~c(array_sim_quant_small2[,,2])))
  summary(lm(c(array_quant_small2[,mps,3])~c(array_sim_quant_small2[,,3])))
  
  
  cor.test(c(array_quant_small2[,mps,1]), c(array_sim_quant_small2_neutral[,,1]), use = "complete")
  cor.test(c(array_quant_small2[,mps,2]), c(array_sim_quant_small2_neutral[,,2]), use = "complete")
  cor.test(c(array_quant_small2[,mps,3]), c(array_sim_quant_small2_neutral[,,3]), use = "complete")
  
  summary(lm(c(array_quant_small2[,mps,1])~c(array_sim_quant_small2_neutral[,,1])))
  summary(lm(c(array_quant_small2[,mps,2])~c(array_sim_quant_small2_neutral[,,2])))
  summary(lm(c(array_quant_small2[,mps,3])~c(array_sim_quant_small2_neutral[,,3])))
}

##############################
#Plot temporal trend
##############################
load("output/tmp/state0.rda")
collst<-adjustcolor(c("grey51", brewer.pal(4, "Set1")), alpha.f = 0.7)

#Get average trajectory
abunds_obs<-exp(array_quant_small2[,which.min(abs(xscl_small-100)),])*rep(state0[1:3], each=dim(array_quant_small2)[1])-0.0001
abunds_sim<-exp(array_sim_quant_small2[,which.min(abs(scalelst*100^2-100)),])*rep(state0[1:3], each=dim(array_sim_quant_small2)[1])
abunds_sim_neutral<-exp(array_sim_quant_small2_neutral[,which.min(abs(scalelst*100^2-100)),])*rep(state0[1:3], each=dim(array_sim_quant_small2_neutral)[1])

ofs2<-c(0.06, -0.1)

pdf("figures/SUP_FIGURE_e014_trends.pdf", width=5, height=6)
par(mfrow=c(3,1), mar=c(2,2,1,1), oma=c(2,2,0.5,0))
matplot(timebands_small, abunds_sim, type="l", col=collst[-1], lty=1, lwd=2, xlab="Age", ylab="pcover", xlim=c(1,80), xaxs="i", log="x"); abline(h=0, lty=3)
mtext("Abundance", 2, line=2.5, cex=1.2)
put.fig.letter("a.", offset=ofs2, cex=1.4)

matplot(timebands_small, abunds_sim_neutral, type="l", col=collst[-1], lty=1, lwd=2, xlab="Age", ylab="pcover", xlim=c(1,80), xaxs="i", log="x"); abline(h=0, lty=3)
mtext("Abundance", 2, line=2.5, cex=1.2)
put.fig.letter("b.", offset=ofs2, cex=1.4)

matplot(timebands_small, abunds_obs*100, type="l", col=collst[-1], lty=1, lwd=2, xlab="Age", ylab="pcover", xlim=c(1,80), xaxs="i", log="x"); abline(h=0, lty=3)
mtext("Percent Cover", 2, line=2.5, cex=1.2)
put.fig.letter("c.", offset=ofs2, cex=1.4)

mtext("Field Age", 1, line=2.5, cex=1.2)
dev.off()

