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
if(length(grep("save_processed_data_FULL_array.RData", dir("output")))==0) {
  nsp<-2
  
  flst<-dir("output/")
  flst_dyn<-flst[grep("dyn", flst)]
  flst_cv<-flst[grep("cv", flst)]
  
  flst_dyn<-flst_dyn[grep(".rda", flst_dyn, fixed = TRUE)]
  flst_cv<-flst_cv[grep(".rda", flst_cv, fixed = TRUE)]
  
  load(paste("output/", flst_dyn[1], sep=""))
  load(paste("output/", flst_cv[1], sep=""))
  
  matout_dyn<-data.frame(matout_dyn)
  matout_cv<-data.frame(matout_cv)
  
  spn1<-paste(rep(c("eig1_pop", "eig2_pop", "r0_pop", "eig1_tot", "eig2_tot", "r0_tot", "beta1_e", "beta2_e", "beta_r", "beta_0"), each=2), c("sp1", "sp2"), sep=".")
  spn2<-paste(rep(c("eig1_pop", "eig2_pop", "r0_pop", "eig1_tot", "eig2_tot", "r0_tot", "beta1_e", "beta2_e", "beta_r", "beta_0"), each=4), c("sp1", "sp2", "sp3", "sp4"), sep=".")
  
  spn1<-c(spn1, "beta_0_tot")
  spn2<-c(spn2, "beta_0_tot")
  
  spcv1<-paste(rep(c("CV_pop", "CV_tot"), each=2), c("sp1", "sp2"), sep=".")
  spcv2<-paste(rep(c("CV_pop", "CV_tot"), each=4), c("sp1", "sp2", "sp3", "sp4"), sep=".")
  
  clnm_dyn<-c("scale", "tscale", "iter",
              paste(spn1, "meta", sep="."),
              paste(spn1, "dist", sep="."),
              paste(spn1, "psf", sep="."),
              paste(spn2, "rps", sep="."),
              paste(spn1, "neut", sep="."))
  clnm_cv<-c("lag", "scale", "iter",
             paste(spcv1, "meta", sep="."),
             paste(spcv1, "dist", sep="."),
             paste(spcv1, "psf", sep="."),
             paste(spcv2, "rps", sep="."),
             paste(spcv1, "neut", sep="."))
  
  #check that dimensions are right:
  colnames(matout_dyn)<-clnm_dyn
  colnames(matout_cv)<-clnm_cv
  
  #Locations          
  e1pop<-grep("eig1_pop", clnm_dyn)
  e2pop<-grep("eig2_pop", clnm_dyn)
  r0pop<-grep("r0_pop", clnm_dyn)
  e1tot<-grep("eig1_tot", clnm_dyn)
  e2tot<-grep("eig2_tot", clnm_dyn)
  r0tot<-grep("r0_tot", clnm_dyn)
  bt1e<-grep("beta1_e", clnm_dyn)
  bt2e<-grep("beta2_e", clnm_dyn)
  btr<-grep("beta_r", clnm_dyn)
  bt0<-grep("beta_0", clnm_dyn)
  bt0_tot<-grep("beta_0_tot", clnm_dyn)
  
  cvpop<-grep("CV_pop", clnm_cv)
  cvtot<-grep("CV_tot", clnm_cv)
  
  #Modeltype
  modlst<-unique(unlist(matrix(nrow=3, unlist(strsplit(clnm_dyn[e1pop], ".", fixed=T)))[3,]))
  
  #lists
  scalslst<-as.numeric(as.character(gsub(".rda", "", gsub("matout_dyn_", "", flst_dyn, fixed=TRUE), fixed=TRUE)))
  tscalelst<-sort(unique(matout_dyn$tscale))
  iterlst<-sort(unique(matout_dyn$iter))
  
  laglst<-sort(unique(matout_cv$lag))
  
  qtl_lims<-c(0.025, pnorm(-1, 0, 1), 0.5, pnorm(1, 0, 1), 0.975)
  
  ###################################
  # Set up plotting data
  ###################################
  #Two types of plots:
  #1. discrete subset of scales
  #2. countour plot
  #For all cases, calculate:
  #max(eigen); min(r0); mean(beta); mean(cv)
  
  
  matout_eig_pop<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(iterlst)))
  matout_r0_pop<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(iterlst)))
  
  matout_eig_tot<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(iterlst)))
  matout_r0_tot<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(iterlst)))
  
  #matout_beta_e<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(iterlst)))
  #matout_beta_r<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(iterlst)))
  #matout_beta_0<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(iterlst)))
  #matout_beta_0_tot<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(iterlst)))
  
  matout_invar_pop<-array(dim=c(length(scalslst), length(laglst), length(modlst), length(iterlst)))
  matout_invar_tot<-array(dim=c(length(scalslst), length(laglst), length(modlst), length(iterlst)))
  
  mxfun<-function(x) {
    if(sum(!is.na(x))>0) {
      return(max(x, na.rm=T))
    } else {
      return(NA)
    }
  }
  
  mnfun<-function(x) {
    if(sum(!is.na(x))>0) {
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
    load(paste("output/", flst_cv[i], sep=""))
    matout_dyn<-data.frame(matout_dyn)
    matout_cv<-data.frame(matout_cv)
    

    print("running dyn...")
    for(j in 1:length(tscalelst)) {
      sbs<-which(matout_dyn$scale==scalslst[i] & matout_dyn$tscale==tscalelst[j])
      for(k in 1:length(modlst)) {
        subscol<-grep(modlst[k], clnm_dyn)
        
        matout_eig_pop[i,j,k,]<-apply(matout_dyn[sbs,intersect(e2pop, subscol)], 1, mxfun)
        matout_eig_tot[i,j,k,]<-apply(matout_dyn[sbs,intersect(e2tot, subscol)], 1, mxfun)
        matout_r0_pop[i,j,k,]<-apply(matout_dyn[sbs,intersect(r0pop, subscol)], 1, mnfun)
        matout_r0_tot[i,j,k,]<-apply(matout_dyn[sbs,intersect(r0tot, subscol)], 1, mnfun)
        
        #matout_beta_e[i,j,k,]<-rowMeans(matout_dyn[sbs,intersect(bt2e, subscol)], na.rm=T)
        #matout_beta_r[i,j,k,]<-rowMeans(matout_dyn[sbs,intersect(btr, subscol)], na.rm=T)
        #matout_beta_0[i,j,k,]<-rowMeans(matout_dyn[sbs,intersect(bt0, subscol)], na.rm=T)
        #matout_beta_0_tot[i,j,k,]<-matout_dyn[sbs,intersect(bt0_tot, subscol)]
      }
      
      if(j/20 == floor(j/20)) {
        print(round(j/length(tscalelst),2))
      }
    }
    
    print("running cv...")
    for(j in 1:length(laglst)) {
      sbs<-which(matout_cv$scale==scalslst[i] & matout_cv$lag==laglst[j])
      
      for(k in 1:length(modlst)) {
        subscol<-grep(modlst[k], clnm_cv)
        
        matout_invar_pop[i,j,k,]<-rowMeans(matout_cv[sbs,intersect(cvpop, subscol)], na.rm=T)
        matout_invar_tot[i,j,k,]<-rowMeans(matout_cv[sbs,intersect(cvtot, subscol)], na.rm=T)
      }
      
      if(j/5 == floor(j/5)) {
        print(round(j/length(laglst),2))
      }
    }
    
    print("total progress:")
    print(round(i/length(scalslst), 2))
  }
  
  rm(matout_dyn)
  rm(matout_cv)
  
  save.image("output/save_processed_data_FULL_array.RData")
} else {
  load("output/save_processed_data_FULL_array.RData")
}


##################################################
#NOW run bootstrapped identification algorithm.
##################################################
#Array order is:
#length(scalslst), length(tscalelst), length(modlst), length(iterlst)
#matout_eig_tot
#matout_eig_pop
#matout_r0_pop

if(length(grep("save_processed_data_FULL_array_2.RData", dir("output")))==0) {
  #density estimation function
  compdens<-function(stat, dnsout) {
    if(is.finite(stat) & !is.character(dnsout)) {
      dns_res<-numeric(length(dnsout))
      
      for(i in 1:length(dnsout)) {
        ps<-which.min(abs(stat-dnsout[[i]]$x))
        dns_res[i]<-dnsout[[i]]$y[ps]
      }
      
      return(log10(dns_res))
    } else {
      return(NA)
    }
  }
  
  #use smaller subset of time scales
  tscallst_small<-sort(unique(c(1:26, seq(28, 50, by=2), seq(55, length(tscalelst), by=5), length(tscalelst)-1, length(tscalelst))))
  
  densout<-array(dim=c(length(scalslst), length(tscallst_small), length(modlst), length(iterlst), 4, length(modlst)))
  #last two dimensions:
  #four metrics
  #five models (density for each)
  
  #extract likelihood estimates
  for(k in 1:length(scalslst)) {
    for(l in 1:length(tscallst_small)) {
      dnsout_eig_tot<-apply(matout_eig_tot[k,tscallst_small[l],,], 1, function(x) try(density(x,na.rm=T), silent=T))
      dnsout_eig_pop<-apply(matout_eig_pop[k,tscallst_small[l],,], 1, function(x) try(density(x,na.rm=T), silent=T))
      dnsout_r0_pop<-apply(matout_r0_pop[k,tscallst_small[l],,], 1, function(x) try(density(x,na.rm=T), silent=T))
      
      for(i in 1:length(modlst)) {
        for(j in 1:length(iterlst)) {
          stat_eig_tot<-matout_eig_tot[k,tscallst_small[l],i,j]
          stat_eig_pop<-matout_eig_pop[k,tscallst_small[l],i,j]
          stat_r0_pop<-matout_r0_pop[k,tscallst_small[l],i,j]
        
          densout[k,l,i,j,1,]<-compdens(stat_eig_tot, dnsout_eig_tot)
          densout[k,l,i,j,2,]<-compdens(stat_eig_pop, dnsout_eig_pop)
          densout[k,l,i,j,3,]<-compdens(stat_r0_pop, dnsout_r0_pop)
          densout[k,l,i,j,4,]<-(densout[k,l,i,j,1,])+(densout[k,l,i,j,2,])+(densout[k,l,i,j,3,])
        }
      }
    }
    print(round(k/length(scalslst),2))
  }
  
  
  
  
  
  
  
  #get cumulative likelihoods
  arrayify2<-function(x,i,k,l,modnum) {
    if(k==1 & l==1) {
      tmp<-10^(array(dim=dim(x[k,l,i,,modnum,]), 
                     c(x[k,l,i,,modnum,])))
      
    } else if(k==1 | l==1) {
      tmp<-10^(apply(array(dim=dim(x[1:k,1:l,i,,modnum,]), 
                           c(x[1:k,1:l,i,,modnum,])), 2:3, function(x) sum(x, na.rm=T)))
    } else {
      tmp<-10^apply(array(dim=dim(x[1:k,1:l,i,,modnum,]), 
                          c(x[1:k,1:l,i,,modnum,])), 3:4, function(x) sum(x, na.rm=T))
    }
    tmp[,i]/rowSums(tmp, na.rm=T)
  }
  
  
  
  
  arrayify<-function(x,i,k,l,modnum) {
    if(k==1 & l==1) {
      rowMeans(array(dim=dim(x[k,l,i,,modnum,-i]), 
                           c(x[k,l,i,,modnum,i])-c(x[k,l,i,,modnum,-i])))
    } else if(k==1 | l==1) {
      rowMeans(apply(array(dim=dim(x[1:k,1:l,i,,modnum,-i]), 
                           c(x[1:k,1:l,i,,modnum,i])-c(x[1:k,1:l,i,,modnum,-i])), 2:3, function(x) sum(x, na.rm=T)))
    } else {
      rowMeans(apply(array(dim=dim(x[1:k,1:l,i,,modnum,-i]), 
                           c(x[1:k,1:l,i,,modnum,i])-c(x[1:k,1:l,i,,modnum,-i])), 3:4, function(x) sum(x, na.rm=T)))
    }
  }
  
  
  qlvls<-c(0.025, pnorm(-1, 0, 1), 0.5, pnorm(1, 0, 1), 0.975)
  densout_cum<-array(dim=c(length(scalslst), length(tscallst_small), length(modlst), 4, length(qlvls)))
  
  for(k in 1:length(scalslst)) {
    for(l in 1:length(tscallst_small)) {
      for(i in 1:length(modlst)) {
        cumstat_eig_tot<-arrayify2(densout,i,k,l,modnum=1)
        cumstat_eig_pop<-arrayify2(densout,i,k,l,modnum=2)
        cumstat_r0_pop<-arrayify2(densout,i,k,l,modnum=3)
        cumstat_tot<-arrayify2(densout,i,k,l,modnum=4)
        
        densout_cum[k,l,i,1,]<-quantile(cumstat_eig_tot, qlvls, na.rm=T)
        densout_cum[k,l,i,2,]<-quantile(cumstat_eig_pop, qlvls, na.rm=T)
        densout_cum[k,l,i,3,]<-quantile(cumstat_r0_pop, qlvls, na.rm=T)
        densout_cum[k,l,i,4,]<-quantile(cumstat_tot, qlvls, na.rm=T)
      }
    }
    print(round(k/length(scalslst),2))
  }
  
  
  rm(matout_eig_tot)
  rm(matout_eig_pop)
  rm(matout_r0_pop)
  rm(matout_r0_tot)
  
  #rm(matout_beta_e)
  #rm(matout_beta_r)
  #rm(matout_beta_0)
  #rm(matout_beta_0_tot)
  
  rm(matout_dyn)
  rm(matout_cv)
  rm(matout_invar_pop)
  rm(matout_invar_tot)
  #rm(x)
  
  rm(densout)
  
  save.image("output/save_processed_data_FULL_array_2.RData")
} else {
  load("output/save_processed_data_FULL_array_2.RData")
}


######################
# Plot results
######################
source("util/plot_grid_functions.R")

densout_cum[!is.finite(densout_cum)]<-NA

#sqtmp<-c(0.01, 0.1, 0.5, 1, 5, 10, 50, 100, 500, 1000, 2000)
sqtmp<-c(0, seq(0.2, 0.8, by=0.1), 0.9, 0.99, 1.01)
logxpos<-c(1,2,5,10,20,50,150)
#suppressWarnings(tmp<-log(densout_cum[,,,1,], 10)); tmp[!is.finite(tmp)]<-NA
#arrayout=tmp; xscalslst=log(scalslst,10); xlst=log(tscallst_small, 10); splitcol=0; nlevels=10; sqlst = sqtmp; logx=TRUE; logy=TRUE; logz=TRUE; logxps = logxpos; coltype=3; logxps=0; nstart=1; ciplot=FALSE; cimat=0; revcol=FALSE; dops_subset=FALSE; override_tmpsq=TRUE

ofs1<-c(0.24, -0.002)


svg("figures/FIGURE_match_models_full.svg", width=5.5, height=8)
m<-matrix(nrow=5, 1:10)
m<-cbind(m, 11)
layout(m, widths=c(1,1,0.5))

par(mar=c(2,3,1,0), oma=c(2.5,2,2,4))

for(i in 2:3) {
  tmp<-densout_cum[,,,i,]
  plotout<-plot_cont(arrayout=tmp, xscalslst=log(scalslst,10), xlst=log(tscallst_small, 10), splitcol=0, nlevels=10, sqlst = sqtmp, logx=TRUE, logy=TRUE, logz=FALSE, logxps = logxpos, coltype=4, nstart=1+5*(i-2), ciplot=FALSE, cimat=0, revcol=FALSE, dops_subset=FALSE, override_tmpsq=TRUE)
}

par(mar=c(2,3.5,1,1))
sqtmp2<-sqtmp
sqtmp2[sqtmp2==1.01]<-1
filled.legend(z=(matrix(1:length(sqtmp2))), levels=1:length(sqtmp2), col=adjustcolor(grey.colors(length(sqtmp2)-1), alpha.f = 0.8), key.axes = axis(4, at = 1:length(sqtmp2), labels = sqtmp2, las=2))

mtext(text = "levins", side = 4, outer = TRUE, line = -7.5+1.2, adj = .94, cex=1.2)
mtext(text = "disturbance", side = 4, outer = TRUE, line = -7.5+1.2, adj = 0.74, cex=1.2)
mtext(text = "PSF", side = 4, outer = TRUE, line = -7.5+1.2, adj = 0.515, cex=1.2)
mtext(text = "RPS", side = 4, outer = TRUE, line = -7.5+1.2, adj = 0.305, cex=1.2)
mtext(text = "neutral", side = 4, outer = TRUE, line = -7.5+1.2, adj = .08, cex=1.2)

mtext(text = expression(paste(r[e], ", population")), side = 3, outer = TRUE, line = 0, adj = .14+0.025, cex=1.2)
mtext(text = expression(paste(r[0], ", population")), side = 3, outer = TRUE, line = 0, adj = 0.64+0.065, cex=1.2)

mtext(text = expression(paste("temporal span, time steps")), side = 1, outer = TRUE, line = 1.3, cex=1.2, adj = 0.3+0.06)
mtext(text = expression(paste("spatial span, fraction of maximum")), side = 2, outer = TRUE, line = -0.1, cex=1.2, adj = 0.46)

mtext(text = expression(paste("Model Identification Success, ", italic(RCL)[italic("m|m")])), side = 4, outer = TRUE, line = 2.5, cex=1.2)

dev.off()




svg("figures/SUP_FIGURE_match_models_population.svg", width=3.2, height=8)
m<-matrix(nrow=5, 1:5)
m<-cbind(m, 6)
layout(m, widths=c(1,0.35))

par(mar=c(2,3,1,0), oma=c(2.5,2,2,5))
ofs1<-c(0.255, -0.002)
for(i in 1:1) {
  tmp<-densout_cum[,,,i,]
  plotout<-plot_cont(arrayout=tmp, xscalslst=log(scalslst,10), xlst=log(tscallst_small, 10), splitcol=0, nlevels=10, sqlst = sqtmp, logx=TRUE, logy=TRUE, logz=FALSE, logxps = logxpos, coltype=4, nstart=1+6*(i-1), ciplot=FALSE, cimat=0, revcol=FALSE, dops_subset=FALSE, override_tmpsq=TRUE)
}
par(mar=c(2,3,1,0))
filled.legend(z=(matrix(1:length(sqtmp2))), levels=1:length(sqtmp2), col=adjustcolor(grey.colors(length(sqtmp2)-1), alpha.f = 0.8), key.axes = axis(4, at = 1:length(sqtmp2), labels = sqtmp2, las=2))

mtext(text = "levins", side = 4, outer = TRUE, line = -3.8, adj = .94, cex=1.2)
mtext(text = "disturbance", side = 4, outer = TRUE, line = -3.8, adj = 0.74, cex=1.2)
mtext(text = "PSF", side = 4, outer = TRUE, line = -3.8, adj = 0.515, cex=1.2)
mtext(text = "RPS", side = 4, outer = TRUE, line = -3.8, adj = 0.305, cex=1.2)
mtext(text = "neutral", side = 4, outer = TRUE, line = -3.8, adj = .08, cex=1.2)

mtext(text = expression(paste(r[e], ", community")), side = 3, outer = TRUE, line = 0, adj = .42, cex=1.2)

mtext(text = expression(paste("temporal span, time steps")), side = 1, outer = TRUE, line = 1.3, cex=1.2, adj = 0.1)
mtext(text = expression(paste("spatial span, fraction of maximum")), side = 2, outer = TRUE, line = -0.1, cex=1.2, adj = 0.46)


mtext(text = expression(paste("Model Identification Success, ", italic(RCL)[italic("m|m")])), side = 4, outer = TRUE, line = 3.5, cex=1.2)

dev.off()

 





