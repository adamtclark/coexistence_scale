#!/usr/bin/env Rscript
error
rm(list=ls())
#setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

require(data.table)
require(RColorBrewer)

source("util/figure_functions.R")
source("util/filled.contour3.R")

###################################
# Load output data
###################################
if(FALSE) {
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
  
  ##################################################
  #run bootstrapped identification algorithm.
  ##################################################
  #Array order is:
  #length(scalslst), length(tscalelst), length(modlst), length(iterlst)
  #matout_eig_tot
  #matout_eig_pop
  #matout_r0_pop
  load("output/tmp/xscl_small.rda")
  load("output/tmp/timebands_small.rda")
  scalelst_short<-c(1,2,4,6,9,12,20,30,49,72,100,196,306,400,812)
  scalslst<-xscl_small
  scl_match<-numeric(length(scalelst_short))
  for(i in 1:length(scalelst_short)) {
    scl_match[i]<-which.min(abs(scalelst_short[i]-scalslst))
  }
  
  tscallst_small<-timebands_small
  iterlst<-1:dim(arrayout_full)[4]

  #density estimation function
  compdens_emp<-function(stat, dnsout) {
    if((sum(stat, na.rm=T)!=0) & !is.character(dnsout)) {
      dns_res<-numeric(length(dnsout))
      
      for(i in 1:length(dnsout)) {
        ps<-which.min(abs(stat[i]-dnsout[[i]]$x))
        dns_res[i]<-dnsout[[i]]$y[ps]
      }
      
      return(log10(dns_res))
    } else {
      return(NA)
    }
  }
  
  #use smaller subset of time scales
  densout<-array(dim=c(length(scalelst_short), length(tscallst_small), 2, dim(arrayout_full)[3]))
  
  #extract likelihood estimates
  for(k in 1:length(scalelst_short)) {
    for(l in 1:length(tscallst_small)) {
      for(m in 1:dim(arrayout_full)[3]) {
        stat_emp<-c(arrayout_full[l,scl_match[k],m,])
        stat_sim<-c(array_sim_full[l,k,m,])
        stat_sim_neutral<-c(array_sim_full_neutral[l,k,m,])
        
        if(length(stat_emp)>3) {
          densout[k,l,1,m]<-mean(as.numeric(sign(stat_emp)==sign(stat_sim)))
          densout[k,l,2,m]<-mean(as.numeric(sign(stat_emp)==sign(stat_sim_neutral)))
          
        }
      }
    }
    print(round(k/length(scalelst_short),2))
  }
  
  #get cumulative likelihoods
  arrayify3<-function(x,k,l,spnum,i=1) {
    x<-log(x)
    x[!is.finite(x)]<-NA
    
    if(k==1 & l==1) {
      tmp<-x[k,l,,spnum]
    } else if(k==1 | l==1) {
      tmp<-colSums(x[1:k,1:l,,spnum])
    } else {
      tmp<-apply(x[1:k,1:l,,spnum], 3, function(x) sum(x, na.rm=T))
    }
    tmp<-exp(tmp)
    tmp[i]/sum(tmp, na.rm=T)
  }
  
  splst<-1:dim(densout)[4]
  
  densout_cum<-array(dim=c(length(scalelst_short), length(tscallst_small), length(splst)))
  
  for(k in 1:length(scalelst_short)) {
    for(l in 1:length(tscallst_small)) {
      cumstat_An<-arrayify3(densout,k,l,spnum=1)
      cumstat_C3<-arrayify3(densout,k,l,spnum=2)
      cumstat_C4<-arrayify3(densout,k,l,spnum=3)

      densout_cum[k,l,1]<-cumstat_An
      densout_cum[k,l,2]<-cumstat_C3
      densout_cum[k,l,3]<-cumstat_C4
    }
    print(round(k/length(scalelst_short),2))
  }
  
  #mask na's:
  mask<-which(apply(arrayout_full[,scl_match,,], 1:2, function(x) sum(!is.na(x)))==0)
  
  #clean up
  rm(arrayout_full)
  rm(array_sim_full)
  rm(array_sim_full_neutral)
  
  rm(densout)
 
  save.image("output/save_processed_data_FULL_array_2_emp.RData")
} else {
  load("output/save_processed_data_FULL_array_2_emp.RData")
}


######################
# Plot results
######################
source("util/plot_grid_functions.R")


plot_dens_cum<-densout_cum

plot_dens_cum[!is.finite(plot_dens_cum)]<-NA
for(i in 1:3) {
  tmp<-t(plot_dens_cum[,,i])
  tmp[mask]<-NA
  plot_dens_cum[,,i]<-t(tmp)
}


sqtmp<-c(0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 1.01)
logxpos<-c(1,2,5,10,20,50,150)

ofs1<-c(0.19, -0.002)


svg("figures/FIGURE_match_models_emp.svg", width=4, height=6)
m<-matrix(nrow=3, 1:3)
m<-cbind(m, 4)
layout(m, widths=c(1,0.35))

par(mar=c(2,3,1,0), oma=c(2.5,2,.5,4))

modlst<-c("Annual", "C3", "C4")
scalelst_short_round<-c(1,2,4,5,10,15,20,30,50,75,100,200,300,400,800)
densout_cum_plot<-array(dim=c(dim(plot_dens_cum),3), plot_dens_cum)
plotout<-plot_cont(arrayout=densout_cum_plot, xscalslst=log(scalelst_short,10), xlst=log(tscallst_small, 10), splitcol=0, nlevels=10, sqlst = sqtmp, logx=TRUE, logy=TRUE, logz=FALSE, logxps = logxpos, coltype=4, nstart=1, ciplot=FALSE, cimat=0, revcol=FALSE, dops_subset=FALSE, override_tmpsq=TRUE, ypos = c(1,5,10,25,50,100,200,400,800))


par(mar=c(2,3.5,1,1))
sqtmp2<-sqtmp
sqtmp2[sqtmp2==1.01]<-1
filled.legend(z=(matrix(1:length(sqtmp2))), levels=1:length(sqtmp2), col=adjustcolor(grey.colors(length(sqtmp2)-1), alpha.f = 0.8), key.axes = axis(4, at = 1:length(sqtmp2), labels = sqtmp2, las=2))

mtext(text = "Annuals", side = 4, outer = TRUE, line = -5.5, adj = 0.905, cex=1.2)
mtext(text = "C3 Grasses", side = 4, outer = TRUE, line = -5.5, adj = 0.515, cex=1.2)
mtext(text = "C4 Grasses", side = 4, outer = TRUE, line = -5.5, adj = 0.1, cex=1.2)

mtext(text = expression(paste("temporal span, time steps")), side = 1, outer = TRUE, line = 1.3, cex=1.2, adj = 0.24)
mtext(text = expression(paste("spatial span, number of patches")), side = 2, outer = TRUE, line = -0.1, cex=1.2, adj = 0.5)

mtext(text = expression(paste("Levins-OF Model Identification Success, ", italic(RCL)[italic("m|OF")])), side = 4, outer = TRUE, line = 2.5, cex=1.2)

dev.off()




save(list=c("plot_dens_cum", "scalelst_short", "tscallst_small"), file = "output/tmp/emp_densplot.rda")

