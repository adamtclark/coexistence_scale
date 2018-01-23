error
rm(list=ls())
setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

require(data.table)
require(RColorBrewer)

source("~/Dropbox/Rfunctions/figure_functions.R")
source("~/Dropbox/Rfunctions/filled.contour3.R")

###################################
# Load output data
###################################
if(FALSE) {
  nsp<-2
  
  flst<-dir("output/")
  flst_dyn<-flst[grep("dyn", flst)]
  flst_cv<-flst[grep("cv", flst)]
  
  matout_dyn<-data.frame(fread(paste("output/", flst_dyn[1], sep="")))
  matout_cv<-data.frame(fread(paste("output/", flst_cv[1], sep="")))
  
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
  scalslst<-as.numeric(as.character(gsub(".csv", "", gsub("matout_dyn_", "", flst_dyn, fixed=TRUE), fixed=TRUE)))
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
  

  matout_eig_pop<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(qtl_lims)))
  matout_r0_pop<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(qtl_lims)))
  
  matout_eig_tot<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(qtl_lims)))
  matout_r0_tot<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(qtl_lims)))
  
  matout_beta_e<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(qtl_lims)))
  matout_beta_r<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(qtl_lims)))
  matout_beta_0<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(qtl_lims)))
  matout_beta_0_tot<-array(dim=c(length(scalslst), length(tscalelst), length(modlst), length(qtl_lims)))
  
  matout_invar_pop<-array(dim=c(length(scalslst), length(laglst), length(modlst), length(qtl_lims)))
  matout_invar_tot<-array(dim=c(length(scalslst), length(laglst), length(modlst), length(qtl_lims)))
  
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
    matout_dyn<-data.frame(fread(paste("output/", flst_dyn[i], sep=""), verbose=FALSE))
    matout_cv<-data.frame(fread(paste("output/", flst_cv[i], sep=""), verbose=FALSE))
    
    print("running dyn...")
    for(j in 1:length(tscalelst)) {
      sbs<-which(matout_dyn$scale==scalslst[i] & matout_dyn$tscale==tscalelst[j])
      for(k in 1:length(modlst)) {
        subscol<-grep(modlst[k], clnm_dyn)
        
        matout_eig_pop[i,j,k,]<-unname(qtfun(apply(matout_dyn[sbs,intersect(e2pop, subscol)], 1, mxfun), qtl_lims, na.rm=TRUE))
          #unname(qtfun(rowMeans(matout_dyn[sbs,intersect(e2pop, subscol)], na.rm=T), qtl_lims, na.rm=TRUE))
          #unname(qtfun(apply(matout_dyn[sbs,intersect(e2pop, subscol)], 1, mxfun), qtl_lims, na.rm=TRUE))
        matout_eig_tot[i,j,k,]<-unname(qtfun(apply(matout_dyn[sbs,intersect(e2tot, subscol)], 1, mxfun), qtl_lims, na.rm=TRUE))
          #unname(qtfun(rowMeans(matout_dyn[sbs,intersect(e2tot, subscol)], na.rm=T), qtl_lims, na.rm=TRUE))
          #unname(qtfun(apply(matout_dyn[sbs,intersect(e2tot, subscol)], 1, mxfun), qtl_lims, na.rm=TRUE))
        matout_r0_pop[i,j,k,]<-unname(qtfun(apply(matout_dyn[sbs,intersect(r0pop, subscol)], 1, mnfun), qtl_lims, na.rm=TRUE))
          #unname(qtfun(rowMeans(matout_dyn[sbs,intersect(r0pop, subscol)], na.rm=T), qtl_lims, na.rm=TRUE))
          #unname(qtfun(apply(matout_dyn[sbs,intersect(r0pop, subscol)], 1, mnfun), qtl_lims, na.rm=TRUE))
        matout_r0_tot[i,j,k,]<-unname(qtfun(apply(matout_dyn[sbs,intersect(r0tot, subscol)], 1, mnfun), qtl_lims, na.rm=TRUE))
          #unname(qtfun(rowMeans(matout_dyn[sbs,intersect(r0tot, subscol)], na.rm=T), qtl_lims, na.rm=TRUE))
          #unname(qtfun(apply(matout_dyn[sbs,intersect(r0tot, subscol)], 1, mnfun), qtl_lims, na.rm=TRUE))
        
        matout_beta_e[i,j,k,]<-unname(qtfun(rowMeans(matout_dyn[sbs,intersect(bt2e, subscol)], na.rm=T), qtl_lims, na.rm=TRUE))
        matout_beta_r[i,j,k,]<-unname(qtfun(rowMeans(matout_dyn[sbs,intersect(btr, subscol)], na.rm=T), qtl_lims, na.rm=TRUE))
        matout_beta_0[i,j,k,]<-unname(qtfun(rowMeans(matout_dyn[sbs,intersect(bt0, subscol)], na.rm=T), qtl_lims, na.rm=TRUE))
        matout_beta_0_tot[i,j,k,]<-unname(qtfun(matout_dyn[sbs,intersect(bt0_tot, subscol)], qtl_lims, na.rm=TRUE))
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
        
        matout_invar_pop[i,j,k,]<-unname(qtfun(rowMeans(matout_cv[sbs,intersect(cvpop, subscol)], na.rm=T), qtl_lims, na.rm=TRUE))
        matout_invar_tot[i,j,k,]<-unname(qtfun(rowMeans(matout_cv[sbs,intersect(cvtot, subscol)], na.rm=T), qtl_lims, na.rm=TRUE))
      }
      
      if(j/5 == floor(j/5)) {
        print(round(j/length(laglst),2))
      }
    }
    
    print("total progress:")
    print(round(i/length(scalslst), 2))
  }
  
  save.image("output/save_processed_data.RData")
} else {
  load("output/save_processed_data.RData")
}






###################################
# Make Plots
###################################
#c("red2", "darkorange2", "gold1", "forestgreen", "dodgerblue", "dodgerblue4", "orchid4")

#dyn is ordered as: scalslst, tscalelst, modlst, qtl_lims
#cv is orderd as:   scalslst, laglst, modlst, qtl_lims

###############
# Discrete
###############
plot_disc<-function(arrayout, xlst, scalesuse=c(1,2,4,6), cifun=function(x,...) is.finite(x), funcol=3, tmlst=TRUE, smooth=0, nstart=0, ...) {
  for(j in 1:length(modlst)) {
    if(tmlst) {
      rng<-range(c(0, t(arrayout[scalesuse,,j,3])*xlst), na.rm=T)
    } else {
      rng<-range(c(0, arrayout[scalesuse,,j,3]), na.rm=T)
    }
    
    plot(range(xlst[colSums(abs(arrayout[scalesuse,,j,3]), na.rm=T)!=0],na.rm=T), rng, type="n", xlab="", ylab="", xaxs="i", axes=FALSE, ...)
    axis(1); axis(2, las=2); box()
    abline(h=0, lty=3)
    put.fig.letter(paste(letters[nstart], ".", sep=""), "topleft", offset=ofs1, cex=1.6)
    nstart<-nstart+1
    
    n<-1
    for(i in scalesuse) {
      if(tmlst) {
        tmp<-arrayout[i,,j,]*xlst
      } else {
        tmp<-arrayout[i,,j,]
      }
      sigsubs<-cifun(x=tmp[,funcol],i=i,j=j)
      
      scl1<-scl2<-tmp[,3]
      scl1[!sigsubs]<-NA
      scl2[sigsubs]<-NA
      
      if(smooth!=0) {
        if(sum(scl1,na.rm=T)>0) {
          scl1<-predict(loess(scl1~xlst, enp.target = smooth), newdata=data.frame(xlst=xlst))
        }
        if(sum(scl2,na.rm=T)>0) {
          scl2<-predict(loess(scl2~xlst, enp.target = smooth), newdata=data.frame(xlst=xlst))
        }
      }
        
      lines(xlst, scl1, lty=1, lwd=1.5, col=collst[n])
      lines(xlst, scl2, lty=3, lwd=1.5, col=collst[n])
      n<-n+1
    }
  }
}

scalesuse<-c(2,4,5,7)
collst<-adjustcolor(rev(c(rainbow(2, start=0.7, end=0.55), rainbow(2, start=0.1, end=0))), alpha.f = 0.85)
ofs1<-c(0.3, -0.002)
  
pdf("figures/SUP_FIGURE_sim_categorical_dyn_results.pdf", width=6, height=8, colormodel = "cmyk")
par(mfcol=c(length(modlst), 3), mar=c(2,4,1,1), oma=c(2.5,1,1,1.5))

plot_disc(matout_eig_tot, tscalelst, scalesuse, cifun = function(x,...) {x<0}, funcol = 4, smooth=0, nstart=1)
mtext(text = expression(paste("community perturbation response, ", lambda)), side = 2, outer = TRUE, line = -1.5, cex=1.2)

plot_disc(matout_eig_pop, tscalelst, scalesuse, cifun = function(x,...) {x<0}, funcol = 4, smooth=0, nstart=6)
mtext(text = expression(paste("population perturbation response, ", lambda)), side = 2, outer = TRUE, line = -1.8-14.2, cex=1.2)
mtext(text = expression(paste("temporal span, time steps")), side = 1, outer = FALSE, line = 3.5, cex=1.2)

plot_disc(matout_r0_pop, tscalelst, scalesuse, cifun = function(x,...) {x>0}, funcol = 2, smooth=0, nstart=11)
mtext(text = expression(paste("invasion rate when rare, ", italic(r[0]))), side = 2, outer = TRUE, line = -1.8-14.2-14.5, cex=1.2)

mtext(text = "levins", side = 4, outer = TRUE, line = 0, adj = .94, cex=1.2)
mtext(text = "disturbance", side = 4, outer = TRUE, line = 0, adj = 0.74, cex=1.2)
mtext(text = "PSF", side = 4, outer = TRUE, line = 0, adj = 0.515, cex=1.2)
mtext(text = "RPS", side = 4, outer = TRUE, line = 0, adj = 0.305, cex=1.2)
mtext(text = "neutral", side = 4, outer = TRUE, line = 0, adj = .08, cex=1.2)
dev.off()


pdf("figures/SUP_FIGURE_sim_categorical_cv_results.pdf", width=8, height=8, colormodel = "cmyk")
par(mfcol=c(length(modlst), 4), mar=c(2,4,1,1), oma=c(2.5,1,1,1.5))

#ofs1<-c(0.3, -0.01)
#par(mar=c(2,4,1,0))
#plot_disc(matout_invar_tot, laglst, scalesuse, tmlst = FALSE, nstart=1)
plot_disc(matout_beta_0_tot, tscalelst, scalesuse, tmlst = FALSE, nstart=1)
mtext(text = expression(paste("lagged community prediction decay, ", italic("CV"))), side = 2, outer = TRUE, line = -1.1, cex=1.2)

#ofs1<-c(0.23, -0.01)
#par(mar=c(2,3,1,1))
#plot_disc(matout_invar_pop, laglst, scalesuse, tmlst = FALSE, nstart=6)
plot_disc(matout_beta_0, tscalelst, scalesuse, tmlst = FALSE, nstart=6)
mtext(text = expression(paste("lagged population prediction decay, ", italic("CV"))), side = 2, outer = TRUE, line = -15.8, cex=1.2)

#ofs1<-c(0.33, -0.01)
#par(mar=c(2,4.5,1,0))
#plot_disc(matout_beta_0, tscalelst, scalesuse, tmlst = FALSE)
plot_disc(matout_beta_e, tscalelst, scalesuse, cifun=function(x,i,j) {x<matout_beta_0[i,,j,4]}, tmlst = FALSE, nstart=11)
mtext(text = expression(paste("post-perturbation dissimilarity, ", italic("CV"))), side = 2, outer = TRUE, line = -30.5, cex=1.2)

#ofs1<-c(0.23, -0.01)
#par(mar=c(2,3,1,1))
plot_disc(matout_beta_r, tscalelst, scalesuse, cifun=function(x,i,j) {x<matout_beta_0[i,,j,4]}, tmlst = FALSE, nstart=16)
mtext(text = expression(paste("post-invasion dissimilarity, ", italic("CV"))), side = 2, outer = TRUE, line = -45.2, cex=1.2)


mtext(text = expression(paste("temporal lag, time steps")), side = 1, outer = TRUE, line = 1.3, adj = 0.25, cex=1.2)
mtext(text = expression(paste("time steps since event")), side = 1, outer = TRUE, line = 1.3, adj = 0.85, cex=1.2)


mtext(text = "levins", side = 4, outer = TRUE, line = 0, adj = .94, cex=1.2)
mtext(text = "disturbance", side = 4, outer = TRUE, line = 0, adj = 0.74, cex=1.2)
mtext(text = "PSF", side = 4, outer = TRUE, line = 0, adj = 0.515, cex=1.2)
mtext(text = "RPS", side = 4, outer = TRUE, line = 0, adj = 0.305, cex=1.2)
mtext(text = "neutral", side = 4, outer = TRUE, line = 0, adj = .08, cex=1.2)
dev.off()




###############
# Contour
###############
#arrayout=matout_eig_pop; xlst=tscalelst; cifun=function(x,...) is.finite(x); funcol=3; tmlst=TRUE
plot_cont<-function(arrayout, xscalslst, xlst, splitcol=0, nlevels=10, sqlst=0, logx=FALSE, logy=FALSE, logz=FALSE, coltype=1, logxps=0, nstart=1, ciplot=FALSE, cimat=0, ...) {

  if(sum(abs(sqlst), na.rm=T)==0) {
    rng<-range(arrayout[,,,3], na.rm=T)
    rng_rnd<-c(floor(rng[1]*10)/10,
               ceiling(rng[2]*10)/10)
    
    sqlst<-pretty(rng_rnd, nlevels)
  }
  
  if(splitcol==0) {
    dm<-0
    if(sum(sqlst==0)==0) {
      sqlst<-c(sqlst[sqlst<0], 0, sqlst[sqlst>0])
    }
  } else {
    if(splitcol=="mean") {
      dm<-mean(sqlst)
    } else {
      dm<-splitcol
    }
    if(sum(sqlst==dm)==0) {
      sqlst<-c(sqlst[sqlst<dm], dm, sqlst[sqlst>dm])
    }
  }
  
  if(logz) {
    tmpsq<-c(floor(10^sqlst[1]*1000)/1000,
             round(10^sqlst[-c(1, length(sqlst))],3),
             ceiling(10^sqlst[length(sqlst)]*1000)/1000)
    if(sum(tmpsq<=0)>0) {
      tmpsq<-c(1e-6, 0.0001, tmpsq[tmpsq>0])
    }
    tmpsq<-sort(unique(tmpsq))
    sqlst<-log(tmpsq[tmpsq>0],10)
  }
  
  for(j in 1:length(modlst)) {
    if(coltype==1) {
      collst2<-adjustcolor(c(rev(rainbow(sum(sqlst<dm), start=0.55, end=.70)), rev(rainbow(sum(sqlst>dm), start=0, end=0.1))), alpha.f = 0.6)
    } else {
      collst2<-adjustcolor(c(rainbow(sum(sqlst<dm), start=0.15, end=.4), rev(rainbow(sum(sqlst>dm), start=0.75, end=0.85))), alpha.f = 0.6)
    }
    
    if(ciplot) {
      sqlst<-seq(-3, 3, by=1)
      collst2<-c("white", "blue", "lightblue", "grey", "pink", "red", "white")
      
      tmpz<-array(dim=dim(arrayout[,,j,3]))
      
      if(sum(abs(cimat), na.rm=T)==0) {
        tmpz[which(is.finite(arrayout[,,j,3]))]<-0
        
        tmpz[which(arrayout[,,j,4]<0)]<-(-1)
        tmpz[which(arrayout[,,j,5]<0)]<-(-2)
        
        tmpz[which(arrayout[,,j,2]>0)]<-(1)
        tmpz[which(arrayout[,,j,1]>0)]<-(2)
      } else {
        tmpz[which(is.finite(arrayout[,,j,3]))]<-0
        
        tmpz[which(arrayout[,,j,4]<cimat[,,j])]<-(-1)
        tmpz[which(arrayout[,,j,5]<cimat[,,j])]<-(-2)
        
        tmpz[which(arrayout[,,j,2]>cimat[,,j])]<-(1)
        tmpz[which(arrayout[,,j,1]>cimat[,,j])]<-(2)
      }
      
      tmpz<-t(tmpz)
    } else {
      tmpz<-t(arrayout[,,j,3])
    }
    
    tmpps<-colSums(abs(arrayout[,,j,3]), na.rm=T)!=0
    filled.contour3(x = xlst, 
                    y = xscalslst, 
                    z = tmpz, levels = sqlst, col=collst2,axes=F,
                    xlim=range(xlst[tmpps]))
    put.fig.letter(paste(letters[nstart], ".", sep=""), "topleft", offset=ofs1, cex=1.6)
    nstart<-nstart+1
    
    if(logz) {
      contour(x = xlst, 
              y = xscalslst, 
              z = tmpz,
              levels = sqlst,
              labels=10^sqlst,
              add=TRUE,axes=F,
              xlim=range(xlst[tmpps]))
    } else {
      contour(x = xlst, 
              y = xscalslst, 
              z = tmpz,
              levels = sqlst,
              labels=round(sqlst,2),
              add=TRUE,axes=F,
              xlim=range(xlst[tmpps]))
    }
    
    
    if(logx) {
      if(sum(abs(logxps))==0) {
        tmp<-round(10^seq(min(xlst), max(xlst), length=6)/50, 1)*50
        tmp<-tmp[tmp>0]
        tmp<-sort(unique(c(1, tmp)))
      } else {
        tmp<-logxps
      }
      
      axis(1, at=log(tmp,10), labels = tmp, las=2)
    } else {
      axis(1, las=2)
    }
    
    if(logy) {
      tmp<-round(10^seq(min(xscalslst), max(xscalslst), length=6)/0.5, 1)*0.5
      #tmp<-10^xscalslst
      tmp<-tmp[tmp>0]
      tmp<-sort(unique(c(0.01, tmp)))
      
      axis(2, at=log(tmp,10), labels = tmp, las=2)
    } else {
      axis(2, las=2)
    }
    
    box()
    
  }
  
  return(list(sqlst=sqlst, collst2=collst2))
}

pdf("figures/FIGURE_sim_continuous_dyn_results.pdf", width=6.5, height=8, colormodel = "cmyk")
sqtmp<-c(-1.5, -1, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 1.5)
logxpos<-c(1,2,5,10,20,50,150,200)

m<-matrix(nrow=5, 1:15)
m<-cbind(m, 16)
layout(m, widths=c(1,1,1,0.5))

par(mar=c(2,3,1,0), oma=c(2.5,2,2,2.5))
ofs1<-c(0.255, -0.002)
tmp<-plot_cont(matout_eig_tot, log(scalslst,10), log(tscalelst, 10), nlevels=10, logx=TRUE, logy=TRUE, sqlst = sqtmp, logxps = logxpos, nstart = 1)
tmp<-plot_cont(matout_eig_pop, log(scalslst,10), log(tscalelst, 10), nlevels=10, logx=TRUE, logy=TRUE, sqlst = sqtmp, logxps = logxpos, nstart = 6)
tmp<-plot_cont(matout_r0_pop, log(scalslst,10), log(tscalelst, 10), nlevels=10, logx=TRUE, logy=TRUE, sqlst = sqtmp, logxps = logxpos, nstart = 11)

par(mar=c(2,3.5,1,1))
filled.legend(z=matrix(sqtmp), levels=sqtmp, col=adjustcolor(c(rev(rainbow(sum(sqtmp<0), start=0.55, end=.70)), rev(rainbow(sum(sqtmp>0), start=0, end=0.1))), alpha.f = 0.6), key.axes = axis(4, at = sqtmp, las=2))

mtext(text = "levins", side = 4, outer = TRUE, line = -5.5, adj = .94, cex=1.2)
mtext(text = "disturbance", side = 4, outer = TRUE, line = -5.5, adj = 0.74, cex=1.2)
mtext(text = "PSF", side = 4, outer = TRUE, line = -5.5, adj = 0.515, cex=1.2)
mtext(text = "RPS", side = 4, outer = TRUE, line = -5.5, adj = 0.305, cex=1.2)
mtext(text = "neutral", side = 4, outer = TRUE, line = -5.5, adj = .08, cex=1.2)

mtext(text = expression(paste(lambda, ", community")), side = 3, outer = TRUE, line = 0, adj = .095, cex=1.2)
mtext(text = expression(paste(lambda, ", population")), side = 3, outer = TRUE, line = 0, adj = .4505, cex=1.2)
mtext(text = expression(paste(r[0], ", population")), side = 3, outer = TRUE, line = 0, adj = 0.805, cex=1.2)

mtext(text = expression(paste("temporal span, time steps")), side = 1, outer = TRUE, line = 1.3, cex=1.2, adj = 0.45)
mtext(text = expression(paste("spatial span, fraction of maximum")), side = 2, outer = TRUE, line = -0.1, cex=1.2, adj = 0.46)
dev.off()



pdf("figures/FIGURE_sim_continuous_cv_results.pdf", width=8.5, height=8, colormodel = "cmyk")
m<-matrix(nrow=5, 1:20)
m<-cbind(m, 21)
layout(m, widths=c(1,1,1,1,0.5))

logxpos<-c(1,2,5,10,20,50,150,500, 1000)

ofs1<-c(0.25, -0.002)
#sqtmp<-c(-5, -4, -3, -2, log10(0.03), -1, log10(0.3), log10(0.5), 0, log10(1.5), log10(3))
#sqtmp<-c(log10(0.005), log10(0.01), log10(0.02), log10(0.05), log10(0.1), log10(0.2), log10(0.5), log10(1), log10(1.2), log10(1.5), log10(2))
sqtmp<-log10(c(0.005, 0.01, 0.015, 0.025, 0.05, 0.075, 0.15, 0.25, 0.5, 0.75, 1, 2))

par(mar=c(2,3,1,0), oma=c(2.5,2,2,3))
#tmp<-plot_cont(log(matout_invar_tot,10), log(scalslst,10), log(laglst+1,10), nlevels=5, splitcol = log(0.5, 10), logx=TRUE, logy=TRUE, logz=TRUE, coltype = 2, sqlst = sqtmp, nstart = 1, logxps = logxpos)
#tmp<-plot_cont(log(matout_invar_pop,10), log(scalslst,10), log(laglst+1,10), nlevels=5, splitcol = log(0.5, 10), logx=TRUE, logy=TRUE, logz=TRUE, coltype = 2, sqlst = sqtmp, nstart = 6, logxps = logxpos)
#arrayout=log(matout_beta_e,10); xscalslst=log(scalslst, 10); xlst=log(tscalelst, 10); tmlst=TRUE; splitcol=0; nlevels=5; sqlst=0; logx=TRUE; logy=TRUE; logz=TRUE; coltype=2


logxpos<-c(1,2,5,10,20,50,150,200)
tmp<-plot_cont(log(matout_beta_0_tot,10), log(scalslst, 10), log(tscalelst, 10), nlevels=5, splitcol = log(0.5, 10), logx=TRUE, logy=TRUE, logz=TRUE, coltype = 2, sqlst = sqtmp, nstart = 1, logxps = logxpos)
tmp<-plot_cont(log(matout_beta_0,10), log(scalslst, 10), log(tscalelst, 10), nlevels=5, splitcol = log(0.5, 10), logx=TRUE, logy=TRUE, logz=TRUE, coltype = 2, sqlst = sqtmp, nstart = 6, logxps = logxpos)
tmp<-plot_cont(log(matout_beta_e,10), log(scalslst, 10), log(tscalelst, 10), nlevels=5, splitcol = log(0.5, 10), logx=TRUE, logy=TRUE, logz=TRUE, coltype = 2, sqlst = sqtmp, nstart = 11, logxps = logxpos)
tmp<-plot_cont(log(matout_beta_r,10), log(scalslst, 10), log(tscalelst, 10), nlevels=5, splitcol = log(0.5, 10), logx=TRUE, logy=TRUE, logz=TRUE, coltype = 2, sqlst = sqtmp, nstart = 16, logxps = logxpos)


par(mar=c(2,3.5,1,1))
filled.legend(z=matrix(sqtmp), levels=sqtmp, col=adjustcolor(c(rainbow(sum(sqtmp<log(0.5, 10)), start=0.15, end=.4), rev(rainbow(sum(sqtmp>log(0.5, 10)), start=0.75, end=0.85))), alpha.f = 0.6), key.axes = axis(4, at = sqtmp, labels = 10^sqtmp, las=2))

mtext(text = "levins", side = 4, outer = TRUE, line = -5.5, adj = .94, cex=1.2)
mtext(text = "disturbance", side = 4, outer = TRUE, line = -5.5, adj = 0.74, cex=1.2)
mtext(text = "PSF", side = 4, outer = TRUE, line = -5.5, adj = 0.515, cex=1.2)
mtext(text = "RPS", side = 4, outer = TRUE, line = -5.5, adj = 0.305, cex=1.2)
mtext(text = "neutral", side = 4, outer = TRUE, line = -5.5, adj = .08, cex=1.2)

mtext(text = expression(paste(italic(CV), ", community")), side = 3, outer = TRUE, line = 0, adj = .055, cex=1.2)
mtext(text = expression(paste(italic(CV), ", population")), side = 3, outer = TRUE, line = 0, adj = 0.327, cex=1.2)

mtext(text = expression(paste(italic(CV), ", ", italic(lambda))), side = 3, outer = TRUE, line = 0, adj = 0.59, cex=1.2)
mtext(text = expression(paste(italic(CV), ", ", italic(r[0]))), side = 3, outer = TRUE, line = 0, adj = .828, cex=1.2)

mtext(text = expression(paste("temporal lag, time steps")), side = 1, outer = TRUE, line = 1.3, cex=1.2, adj = 0.15)
mtext(text = expression(paste("time steps since event")), side = 1, outer = TRUE, line = 1.3, cex=1.2, adj = 0.75)


mtext(text = expression(paste("spatial span, fraction of maximum")), side = 2, outer = TRUE, line = -0.1, cex=1.2, adj = 0.46)
dev.off()








#CI plots...

pdf("figures/SUP_FIGURE_sim_CI_continuous_dyn_results.pdf", width=6.5, height=8, colormodel = "cmyk")
sqtmp<-c(-1.5, -1, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1, 1.5)
logxpos<-c(1,2,5,10,20,50,150,200)

m<-matrix(nrow=5, 1:15)
m<-cbind(m, 16)
layout(m, widths=c(1,1,1,0.5))

par(mar=c(2,3,1,0), oma=c(2.5,2,2,4.5))
ofs1<-c(0.255, -0.002)
tmp<-plot_cont(matout_eig_tot, log(scalslst,10), log(tscalelst, 10), nlevels=10, logx=TRUE, logy=TRUE, sqlst = sqtmp, logxps = logxpos, nstart = 1, ciplot = TRUE)
tmp<-plot_cont(matout_eig_pop, log(scalslst,10), log(tscalelst, 10), nlevels=10, logx=TRUE, logy=TRUE, sqlst = sqtmp, logxps = logxpos, nstart = 6, ciplot = TRUE)
tmp<-plot_cont(matout_r0_pop, log(scalslst,10), log(tscalelst, 10), nlevels=10, logx=TRUE, logy=TRUE, sqlst = sqtmp, logxps = logxpos, nstart = 11, ciplot = TRUE)

par(mar=c(2,3.5,1,1))
sqtmp<-seq(-2, 2)
filled.legend(z=matrix(sqtmp), levels=sqtmp, col=c("blue", "lightblue", "pink", "red"), key.axes = axis(4, at = seq(-1.5, 1.5), labels = c("97.5% < 0", "+1SD < 0", "-1SD > 0", "2.5% > 0"), las=2))

mtext(text = "levins", side = 4, outer = TRUE, line = -5.5, adj = .94, cex=1.2)
mtext(text = "disturbance", side = 4, outer = TRUE, line = -5.5, adj = 0.74, cex=1.2)
mtext(text = "PSF", side = 4, outer = TRUE, line = -5.5, adj = 0.515, cex=1.2)
mtext(text = "RPS", side = 4, outer = TRUE, line = -5.5, adj = 0.305, cex=1.2)
mtext(text = "neutral", side = 4, outer = TRUE, line = -5.5, adj = .08, cex=1.2)

mtext(text = expression(paste(lambda, ", community")), side = 3, outer = TRUE, line = 0, adj = .095, cex=1.2)
mtext(text = expression(paste(lambda, ", population")), side = 3, outer = TRUE, line = 0, adj = .4505, cex=1.2)
mtext(text = expression(paste(r[0], ", population")), side = 3, outer = TRUE, line = 0, adj = 0.805, cex=1.2)

mtext(text = expression(paste("temporal span, time steps")), side = 1, outer = TRUE, line = 1.3, cex=1.2, adj = 0.45)
mtext(text = expression(paste("spatial span, fraction of maximum")), side = 2, outer = TRUE, line = -0.1, cex=1.2, adj = 0.46)
dev.off()


pdf("figures/SUP_FIGURE_sim_CI_continuous_cv_results.pdf", width=5.5, height=8, colormodel = "cmyk")
m<-matrix(nrow=5, 1:10)
m<-cbind(m, 11)
layout(m, widths=c(1,1,0.5))

logxpos<-c(1,2,5,10,20,50,150,500, 1000)

ofs1<-c(0.24, -0.002)
sqtmp<-c(-5, -4, -3, -2, log10(0.03), -1, log10(0.3), log10(0.5), 0, log10(1.5), log10(3))
par(mar=c(2,3,1,0), oma=c(2.5,2,2,4.5))

logxpos<-c(1,2,5,10,20,50,150,200)
tmp<-plot_cont(log(matout_beta_e,10), log(scalslst, 10), log(tscalelst, 10), nlevels=5, splitcol = log(0.5, 10), logx=TRUE, logy=TRUE, logz=FALSE, coltype = 2, sqlst = sqtmp, nstart = 1, logxps = logxpos, ciplot = TRUE, cimat=log(matout_beta_0,10)[,,,3])
tmp<-plot_cont(log(matout_beta_r,10), log(scalslst, 10), log(tscalelst, 10), nlevels=5, splitcol = log(0.5, 10), logx=TRUE, logy=TRUE, logz=FALSE, coltype = 2, sqlst = sqtmp, nstart = 6, logxps = logxpos, ciplot = TRUE, cimat=log(matout_beta_0,10)[,,,3])
#plot_cont(matout_beta_0, log(tscalelst, 10), nlevels=5, splitcol = FALSE, rng=c(0,1))

par(mar=c(2,3.5,1,1))
sqtmp<-seq(-2, 2)
filled.legend(z=matrix(sqtmp), levels=sqtmp, col=c("blue", "lightblue", "pink", "red"), key.axes = axis(4, at = seq(-1.5, 1.5), labels = c(expression(paste("97.5% < ", beta[0])), expression(paste("+1SD < ", beta[0])), expression(paste("-1SD > ", beta[0])), expression(paste("2.5% > ", beta[0]))), las=2))

mtext(text = "levins", side = 4, outer = TRUE, line = -5.5, adj = .94, cex=1.2)
mtext(text = "disturbance", side = 4, outer = TRUE, line = -5.5, adj = 0.74, cex=1.2)
mtext(text = "PSF", side = 4, outer = TRUE, line = -5.5, adj = 0.515, cex=1.2)
mtext(text = "RPS", side = 4, outer = TRUE, line = -5.5, adj = 0.305, cex=1.2)
mtext(text = "neutral", side = 4, outer = TRUE, line = -5.5, adj = .08, cex=1.2)

mtext(text = expression(paste(italic(CV[beta]), ", ", italic(lambda))), side = 3, outer = TRUE, line = 0, adj = 0.22, cex=1.2)
mtext(text = expression(paste(italic(CV[beta]), ", ", italic(r[0]))), side = 3, outer = TRUE, line = 0, adj = .66, cex=1.2)

mtext(text = expression(paste("time steps since event")), side = 1, outer = TRUE, line = 1.3, cex=1.2, adj = 0.42)


mtext(text = expression(paste("spatial span, fraction of maximum")), side = 2, outer = TRUE, line = -0.1, cex=1.2, adj = 0.46)
dev.off()










###############
# Plot specific model details
###############
source("run_metapopulation_wrapper.R")

collst3<-adjustcolor(c("grey51", brewer.pal(4, "Set1")), alpha.f = 0.7)

##### 1. Show Peturbation vs. stochastic changes at small scales
set.seed(22012018)
runtype = "metapopulation"
gridout<-makegrid(xlng = 100, ylng = 100)
population<-populate(gridout, nlst = 0.1, clst = c(0.15, 0.35)*3, mlst = rep(0.3, 2), radlst = Inf)
grid_sub<-grid_subset(gridout, size = 0.01)
out<-run_metapopulation(tmax=300, nsteps = 300, gridout, population, talktime = 0, sites_sub = grid_sub$sites)
eqret1<-estimate_eqreturn(out, replace_perturb = 0, prtb = 0.2, simtime = 200, perturbsites = grid_sub$sites, sites_sub = grid_sub$sites, doplot = FALSE)
#eqret2<-estimate_eqreturn(out, replace_perturb = 1, simtime = 1000, perturbsites = grid_sub$sites, doplot = FALSE)

tmp<-eqret1$out_lst[[1]]$output_spatial
tmp[,1]<-tmp[,1]+max(out$output_spatial[,1])
tmp<-rbind(out$output_spatial, tmp)


pdf("figures/SUP_FIGURE_dem_stoch_vs_pert.pdf", width=4, height=3, colormodel = "cmyk")
par(mar=c(4,4,1,1), oma=c(0,0,0,0))

sbs<--c(1:250)
plot(tmp[sbs,1]+1-100, tmp[sbs,2], type="l", col=collst3[2], lwd=1.5, xlab="", ylab="", axes=F, xaxs="i")
axis(1); axis(2, las=2); box()
abline(v=300-100, lty=2)
abline(h=out$plotdata$ceq[1]*length(grid_sub$sites), col=collst3[2], lty=3)

mtext(text = expression(paste("simulation time")), side = 1, outer = F, line = 2.5, cex=1.2, adj = 0.42)
mtext(text = expression(paste("number of individuals")), side = 2, outer = F, line = 2.5, cex=1.2, adj = 0.42)

#Add var
lms<-mean(eqret1$out_lst[[1]]$output_spatial[,2])+sd(eqret1$out_lst[[1]]$output_spatial[,2])*c(-1,1)
polygon(c(0, 1000, 1000, 0), c(lms[c(1,1,2,2)]), col=adjustcolor(collst3[2], alpha.f = 0.5))

arrows(tmp[300,1]+1-100, tmp[300,2], tmp[301,1]+1-100, tmp[301,2], length = 0.12, lwd=1.5, angle = 20)
dev.off()


##### 2. Show growth/invasion around disturbance events in dist. model
matout_dyn<-data.frame(fread(paste("output/", flst_dyn[grep("_1.csv", flst_dyn)], sep=""), verbose=FALSE))
sbs<-which(matout_dyn$scale==1)
subscol<-grep(modlst[2], clnm_dyn)

ofs2<-c(0.275, -0.017)
pdf("figures/SUP_FIGURE_disturbance_example.pdf", width=6, height=3, colormodel = "cmyk")
par(mar=c(4,4,1,1), oma=c(0,0,0,0), mfrow=c(1,2))
#Eigen
sp1<-tapply(matout_dyn[sbs,intersect(e2pop, subscol)[1]], matout_dyn[sbs,"tscale"], function(x) median(x, na.rm=T))
sp2<-tapply(matout_dyn[sbs,intersect(e2pop, subscol)[2]], matout_dyn[sbs,"tscale"], function(x) median(x, na.rm=T))

matplot(tscalelst, cbind(sp1, sp2)*tscalelst, type="l", lwd=1.5, col=collst3[3:2], lty=1, xlab="", ylab="", axes=F)
abline(h=0, lty=3)
abline(v=c(0, 50, 100, 150, 200), lty=2)
axis(1); axis(2, las=2); box()

mtext("time since event", 1, line=2.5, cex=1.2)
mtext(expression(paste(italic(lambda))), 2.5, line=2.5, cex=1.2)
put.fig.letter("a.", "topleft", offset=ofs2, cex=1.2)

#r0
sp1<-tapply(matout_dyn[sbs,intersect(r0pop, subscol)[1]], matout_dyn[sbs,"tscale"], function(x) median(x, na.rm=T))
sp2<-tapply(matout_dyn[sbs,intersect(r0pop, subscol)[2]], matout_dyn[sbs,"tscale"], function(x) median(x, na.rm=T))

matplot(tscalelst, cbind(sp1, sp2)*tscalelst, type="l", lwd=1.5, col=collst3[3:2], lty=1, xlab="", ylab="", axes=F)
abline(h=0, lty=3)
abline(v=c(0, 50, 100, 150, 200), lty=2)
axis(1); axis(2, las=2); box()

mtext("time since event", 1, line=2.5, cex=1.2)
mtext(expression(paste(italic(r[0]))), 2.5, line=2.5, cex=1.2)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.2)
dev.off()


##### 3. Show neutral mass effects/die-off at mid-scales
set.seed(23012018)
runtype = "neutral"
gridout<-makegrid(xlng = 100, ylng = 100)
population<-populate(gridout, nlst = 0.1, clst = rep(0.5, 2)*3, mlst = rep(0.7, 2), radlst = Inf)

grid_sub1<-grid_subset(gridout, size = 0.2)
grid_sub2<-grid_subset(gridout, size = 0.01)

set.seed(28012018)
out1<-run_metapopulation(tmax=300, nsteps = 300, runtype = "neutral", gridout, population, talktime = 0, sites_sub = grid_sub1$sites)
set.seed(28012018)
out2<-run_metapopulation(tmax=300, nsteps = 300, runtype = "neutral", gridout, population, talktime = 0, sites_sub = grid_sub2$sites)

r0_neut1<-estimate_rarereturn(out1, simtime=200, burnin=200, runtype="neutral", perturbsites = grid_sub1$sites, sites_sub = grid_sub1$sites, doplot=FALSE)
r0_neut2<-estimate_rarereturn(out2, simtime=200, burnin=200, runtype="neutral", perturbsites = grid_sub2$sites, sites_sub = grid_sub2$sites, doplot=FALSE)

pdf("figures/SUP_FIGURE_neutral_extinction.pdf", width=6, height=3, colormodel = "cmyk")
par(mar=c(4,4,1,1), oma=c(0,0,0,0), mfrow=c(1,2))

tmp1<-r0_neut1$out0_lst[[1]]$output
tmp1[,1]<-tmp1[,1]+max(out1$output[,1])
tmp1<-rbind(out1$output, tmp1)

tmp2<-r0_neut2$out0_lst[[1]]$output
tmp2[,1]<-tmp2[,1]+max(out2$output[,1])
tmp2<-rbind(out2$output, tmp2)

sbs<--c(1:250)
plot(tmp1[sbs,1]+1-100, tmp1[sbs,2]/out1$plotdata$ngrid, type="l", col=collst3[2], lwd=1.5, xlab="", ylab="", axes=F, xaxs="i", ylim=c(0, 0.4))
axis(1); axis(2, las=2); box()
abline(v=300-100, lty=2); abline(h=0, lty=3)
arrows(tmp1[300,1]+1-100, tmp1[300,2]/out1$plotdata$ngrid, tmp1[301,1]+1-100, tmp1[301,2]/out1$plotdata$ngrid, length = 0.12, lwd=1.5, angle = 20)
put.fig.letter("a.", "topleft", offset=ofs2, cex=1.2)

mtext(text = expression(paste("simulation time")), side = 1, outer = F, line = 2.5, cex=1.2, adj = 0.42)
mtext(text = expression(paste("abundance")), side = 2, outer = F, line = 2.5, cex=1.2, adj = 0.42)

plot(tmp2[sbs,1]+1-100, tmp2[sbs,2]/out1$plotdata$ngrid, type="l", col=collst3[2], lwd=1.5, xlab="", ylab="", axes=F, xaxs="i", ylim=c(0, 0.4))
axis(1); axis(2, las=2); box()
abline(v=300-100, lty=2); abline(h=0, lty=3)
arrows(tmp2[300,1]+1-100, tmp2[300,2]/out1$plotdata$ngrid+0.01, tmp2[301,1]+1-100, tmp2[301,2]/out1$plotdata$ngrid-0.02, length = 0.12, lwd=1.5, angle = 20)

mtext(text = expression(paste("simulation time")), side = 1, outer = F, line = 2.5, cex=1.2, adj = 0.42)
mtext(text = expression(paste("abundance")), side = 2, outer = F, line = 2.5, cex=1.2, adj = 0.42)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.2)
dev.off()




##### 4. Show recovery from perturbation at population etc. scale.
set.seed(22012018)
runtype = "neutral"
gridout<-makegrid(xlng = 100, ylng = 100)
population<-populate(gridout, nlst = 0.1, clst = rep(0.5, 2)*3, mlst = rep(0.7, 2), radlst = Inf)
out<-run_metapopulation(tmax=300, nsteps = 300, gridout, population, talktime = 0, runtype = runtype)

eqret1<-estimate_eqreturn(out, runtype = runtype, replace_perturb = 0, prtb = 0.2, simtime = 200, doplot = FALSE)
eqret2<-estimate_eqreturn(out, runtype = runtype, replace_perturb = 1, prtb = 0.2, simtime = 200, doplot = FALSE)

tmp1<-eqret1$out_lst[[1]]$output
tmp1[,1]<-tmp1[,1]+max(out$output[,1])
tmp1<-rbind(out$output, tmp1)
tmp1<-cbind(tmp1[,1], rowSums(tmp1[,-1]), tmp1[,2:3])

tmp2<-eqret2$out_lst[[1]]$output
tmp2[,1]<-tmp2[,1]+max(out$output[,1])
tmp2<-rbind(out$output, tmp2)
tmp2<-cbind(tmp2[,1], rowSums(tmp2[,-1]), tmp2[,2:3])


pdf("figures/SUP_FIGURE_perturb_replace.pdf", width=6, height=3, colormodel = "cmyk")
par(mar=c(4,4,1,1), oma=c(0,0,0,0), mfrow=c(1,2))

sbs<--c(1:250, 350:500)
matplot(tmp1[sbs,1]+1-100, tmp1[sbs,-1]/out$plotdata$ngrid, type="l", col=collst3, lwd=1.5, lty=1, xlab="", ylab="", axes=F, xaxs="i")
axis(1); axis(2, las=2); box()
abline(v=300-100, lty=2)
abline(h=out$plotdata$ceq[1], col=collst3[1], lty=3)
arrows(tmp1[300,1]+1-100, tmp1[300,3]/out$plotdata$ngrid, tmp1[301,1]+1-100, tmp1[301,3]/out$plotdata$ngrid, length = 0.12, lwd=1.5, angle = 20)

mtext(text = expression(paste("simulation time")), side = 1, outer = F, line = 2.5, cex=1.2, adj = 0.42)
mtext(text = expression(paste("abundance")), side = 2, outer = F, line = 2.5, cex=1.2, adj = 0.42)
put.fig.letter("a.", "topleft", offset=ofs2, cex=1.2)


matplot(tmp2[sbs,1]+1-100, tmp2[sbs,-1]/out$plotdata$ngrid, type="l", col=collst3, lwd=1.5, lty=1, xlab="", ylab="", axes=F, xaxs="i")
axis(1); axis(2, las=2); box()
abline(v=300-100, lty=2)
abline(h=out$plotdata$ceq[1], col=collst3[1], lty=3)

arrows(tmp2[300,1]+1-100, tmp2[300,3]/out$plotdata$ngrid, tmp2[301,1]+1-100, tmp2[301,3]/out$plotdata$ngrid, length = 0.12, lwd=1.5, angle = 20)
arrows(tmp2[300,1]+1-100, tmp2[300,4]/out$plotdata$ngrid, tmp2[301,1]+1-100, tmp2[301,4]/out$plotdata$ngrid, length = 0.12, lwd=1.5, angle = 20)

mtext(text = expression(paste("simulation time")), side = 1, outer = F, line = 2.5, cex=1.2, adj = 0.42)
mtext(text = expression(paste("abundance")), side = 2, outer = F, line = 2.5, cex=1.2, adj = 0.42)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.2)

dev.off()









