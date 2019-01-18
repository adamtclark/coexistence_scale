#!/usr/bin/env Rscript
#error
rm(list=ls())
#setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

require(data.table)
require(RColorBrewer)

source("util/figure_functions.R")
source("util/filled.contour3.R")

###################################
# Load output data
###################################
if(length(grep("save_processed_data_array_psf.RData", dir("output")))==0) {
  nsp<-2
  
  flst<-dir("output/")
  flst_dyn<-flst[grep("dyn", flst)]
  flst_dyn<-flst_dyn[grep(".rda", flst_dyn, fixed = TRUE)]
  
  flst_dyn<-flst_dyn[grep(".rda", flst_dyn, fixed = TRUE)]
  
  flst_dyn<-flst_dyn[grep("psf", flst_dyn, fixed=T)]
  
  load(paste("output/", flst_dyn[1], sep=""))
  matout_dyn<-data.frame(matout_dyn)
  
  spn1<-paste(rep(c("eig1_pop", "eig2_pop", "r0_pop", "eig1_tot", "eig2_tot", "r0_tot"), each=2), c("sp1", "sp2"), sep=".")
  spn2<-paste(rep(c("eig1_pop", "eig2_pop", "r0_pop", "eig1_tot", "eig2_tot", "r0_tot"), each=4), c("sp1", "sp2", "sp3", "sp4"), sep=".")
  
  clnm_dyn<-c("scale", "tscale", "iter",
              paste(spn1, "psf", sep="."))
  
  #check that dimensions are right:
  colnames(matout_dyn)<-clnm_dyn
  
  #Locations          
  e1pop<-grep("eig1_pop", clnm_dyn)
  e2pop<-grep("eig2_pop", clnm_dyn)
  r0pop<-grep("r0_pop", clnm_dyn)
  e1tot<-grep("eig1_tot", clnm_dyn)
  e2tot<-grep("eig2_tot", clnm_dyn)
  r0tot<-grep("r0_tot", clnm_dyn)
  
  #Modeltype
  modlst<-unique(unlist(matrix(nrow=3, unlist(strsplit(clnm_dyn[e1pop], ".", fixed=T)))[3,]))
  
  #lists
  scalslst<-as.numeric(as.character(gsub(".rda", "", gsub("matout_dyn_psf_", "", flst_dyn, fixed=TRUE), fixed=TRUE)))
  tscalelst<-sort(unique(matout_dyn$tscale))
  iterlst<-sort(unique(matout_dyn$iter))
  
  qtl_lims<-c(0.025, pnorm(-1, 0, 1), 0.5, pnorm(1, 0, 1), 0.975)
  
  ###################################
  # Set up plotting data
  ###################################
  #Two types of plots:
  #1. discrete subset of scales
  #2. countour plot
  #For all cases, calculate:
  #max(eigen); min(r0); mean(beta); mean(cv)
  
  
  matout_eig_pop<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(qtl_lims)))
  matout_r0_pop<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(qtl_lims)))
  
  matout_eig_tot<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(qtl_lims)))
  matout_r0_tot<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(qtl_lims)))
  
  mxfun<-function(x) {
    if(sum(is.na(x))==0) {
      return(max(x, na.rm=T))
    } else {
      return(NA)
    }
  }
  
  mnfun<-function(x) {
    if(sum(is.na(x))==0) {
      return(min(x, na.rm=T))
    } else {
      return(NA)
    }
  }
  
  qtfun<-function(x, ...) {
    #only take quantile function if at least 2/3 of entries are non-na
    if((sum(is.na(x))/length(x))>(1/3)) {
      x[]<-NA
    }
    quantile(x, ...)
  }
  
  for(i in 1:length(scalslst)) {
    print(paste("scale =", scalslst[i]))
    
    load(paste("output/", flst_dyn[i], sep=""))
    matout_dyn<-data.frame(matout_dyn)
    
    print("running dyn...")
    for(j in 1:length(tscalelst)) {
      sbs<-which(matout_dyn$scale==scalslst[i] & matout_dyn$tscale==tscalelst[j])
      for(k in 1:length(modlst)) {
        subscol<-grep(modlst[k], clnm_dyn)
        
        matout_eig_pop[i,j,k,]<-unname(qtfun(apply(matout_dyn[sbs,intersect(e2pop, subscol)], 1, mxfun), qtl_lims, na.rm=TRUE))
        matout_eig_tot[i,j,k,]<-unname(qtfun(apply(matout_dyn[sbs,intersect(e2tot, subscol)], 1, mxfun), qtl_lims, na.rm=TRUE))
        matout_r0_pop[i,j,k,]<-unname(qtfun(apply(matout_dyn[sbs,intersect(r0pop, subscol)], 1, mnfun), qtl_lims, na.rm=TRUE))
        matout_r0_tot[i,j,k,]<-unname(qtfun(apply(matout_dyn[sbs,intersect(r0tot, subscol)], 1, mnfun), qtl_lims, na.rm=TRUE))
      }
      
      if(j/20 == floor(j/20)) {
        print(round(j/length(tscalelst),2))
      }
    }
    
    print("total progress:")
    print(round(i/length(scalslst), 2))
  }
  
  save.image("output/save_processed_data_array_psf.RData")
} else {
  load("output/save_processed_data_array_psf.RData")
}




###################################
# Make Plots
###################################
source("util/plot_grid_functions.R")

###############
# Contour
###############

matout_r0_pop[!is.na(matout_r0_pop) & matout_r0_pop<(-1e6)]<-(-1e3)

#pdf("figures/SUP_FIGURE_sim_continuous_dyn_results_psf.pdf", width=6, height=4, colormodel = "cmyk")
svg("figures/SUP_FIGURE_sim_continuous_dyn_results_psf.svg", width=6, height=4)
sqtmp<-c(-1.5, -1, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 1.5)
sqtmp_plot<-c(-1e4, -1, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 1e4)
logxpos<-c(1,2,5,10,20,50,150,200)

m<-matrix(nrow=1, 1:2)
m<-cbind(m, 3)
layout(m, widths=c(1,1,0.55))

par(mar=c(2,3,1,0), oma=c(2.5,2,2,3))
ofs1<-c(0.21, 0.01)
tmp<-plot_cont(matout_eig_pop, log(scalslst,10), log(tscalelst, 10), nlevels=10, logx=TRUE, logy=TRUE, sqlst = sqtmp_plot, logxps = logxpos, nstart = 1, dolines=FALSE)
tmp<-plot_cont(matout_r0_pop, log(scalslst,10), log(tscalelst, 10), nlevels=10, logx=TRUE, logy=TRUE, sqlst = sqtmp_plot, logxps = logxpos, nstart = 2, revcol = TRUE, dolines=FALSE)

par(mar=c(2,5.5,1,1))
filled.legend(z=matrix(sqtmp), levels=-sqtmp, col=adjustcolor(c(rev(rainbow(sum(sqtmp<0), start=0.55, end=.70)), rev(rainbow(sum(sqtmp>0), start=0, end=0.1))), alpha.f = 0.6), key.axes = axis(4, at = sqtmp, las=2))
axis(2, at=sqtmp, sprintf("%.2f", -sqtmp), las=2)

mtext(text = "PSF", side = 4, outer = TRUE, line = -7.8, adj = 0.515, cex=1.2)

mtext(text = expression(paste(r[e], ", population")), side = 3, outer = TRUE, line = 0, adj = .14, cex=1.2)
mtext(text = expression(paste(r[0], ", population")), side = 3, outer = TRUE, line = 0, adj = 0.64, cex=1.2)

mtext(text = expression(paste("temporal extent, time steps")), side = 1, outer = TRUE, line = 1.3, cex=1.2, adj = 0.4)
mtext(text = expression(paste("spatial extent, fraction of maximum")), side = 2, outer = TRUE, line = -0.1, cex=1.2, adj = 0.46)

mtext(text = expression(paste(r[e])), side = 3, outer = TRUE, line = 0, adj = .9, cex=1.5)
mtext(text = expression(paste(r[0])), side = 3, outer = TRUE, line = 0, adj = 1.035, cex=1.5, xpd=NA)

dev.off()
