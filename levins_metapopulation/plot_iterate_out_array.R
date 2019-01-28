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
if(length(grep("save_processed_data_array.RData", dir("output")))==0) {
  nsp<-2
  
  flst<-dir("output/")
  flst_dyn<-flst[grep("dyn", flst)]
  
  flst_dyn<-flst_dyn[grep(".rda", flst_dyn, fixed = TRUE)]
  
  load(paste("output/", flst_dyn[1], sep=""))
  matout_dyn<-data.frame(matout_dyn)
  
  spn1<-paste(rep(c("eig1_pop", "eig2_pop", "r0_pop", "eig1_tot", "eig2_tot", "r0_tot"), each=2), c("sp1", "sp2"), sep=".")
  spn2<-paste(rep(c("eig1_pop", "eig2_pop", "r0_pop", "eig1_tot", "eig2_tot", "r0_tot"), each=4), c("sp1", "sp2", "sp3", "sp4"), sep=".")
  
  clnm_dyn<-c("scale", "tscale", "iter",
          paste(spn1, "meta", sep="."),
          paste(spn1, "dist", sep="."),
          paste(spn1, "psf", sep="."),
          paste(spn2, "rps", sep="."),
          paste(spn1, "neut", sep="."))
  
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
  scalslst<-as.numeric(as.character(gsub(".rda", "", gsub("matout_dyn_", "", flst_dyn, fixed=TRUE), fixed=TRUE)))
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
    matout_dyn<-data.frame(matout_dyn)
    
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
      }
      
      if(j/20 == floor(j/20)) {
        print(round(j/length(tscalelst),2))
      }
    }
    
    print("total progress:")
    print(round(i/length(scalslst), 2))
  }
  
  save.image("output/save_processed_data_array.RData")
} else {
  load("output/save_processed_data_array.RData")
}






###################################
# Make Plots
###################################
source("util/plot_grid_functions.R")

###############
# Contour
###############

#pdf("figures/FIGURE_sim_continuous_dyn_results.pdf", width=5, height=8, colormodel = "cmyk")
svg("figures/FIGURE_sim_continuous_dyn_results.svg", width=5, height=8)
sqtmp<-c(-1.5, -1, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 1.5)
#sqtmp<-c(-1.5,-1,-0.5,-0.25,-0.1,-0.05,0,0.05,0.1,0.25,0.5,1,1.5)
#sqtmp<-c(-1, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 1.5)

logxpos<-c(1,2,5,10,20,50,150,200)

m<-matrix(nrow=5, 1:10)
m<-cbind(m, 11)
layout(m, widths=c(1,1,0.7))

par(mar=c(2,3,1,0), oma=c(2.5,2,2,3))
ofs1<-c(0.255, -0.002)

#Get sd's
matout_eig_pop_SD<-matout_eig_pop
critz<-0.2326174
sdest<-((matout_eig_pop[,,,3]-matout_eig_pop[,,,2])+(matout_eig_pop[,,,4]-matout_eig_pop[,,,3]))/2
matout_eig_pop_SD[,,,3]<-matout_eig_pop_SD[,,,3]/sdest

zcrit<-abs(qnorm(0.05, 0, 1))
sqtmpz<-sort(c(zcrit/sqrt(c(1,2,5,10,30,50)), 0, -zcrit/sqrt(c(1,2,5,10,30,50))))
n_eq<-round((1/(sqtmpz/zcrit))^2*sign(sqtmpz))
sqtmpz<-c(-1000, sqtmpz, 1000)
sqtmpz<-round(sqtmpz, 1)

#plot re
tmp<-plot_cont(arrayout = matout_eig_pop, xscalslst = log(scalslst,10), xlst = log(tscalelst, 10), nlevels=10, logx=TRUE, logy=TRUE, sqlst = sqtmp, logxps = logxpos, nstart = 1, dolines=FALSE, arrayout2 = matout_eig_pop_SD, doz = TRUE, sqz = sqtmpz)

#Get sd's
matout_r0_pop_SD<-matout_r0_pop
critz<-0.2326174
sdest<-((matout_r0_pop_SD[,,,3]-matout_r0_pop_SD[,,,2])+(matout_r0_pop_SD[,,,4]-matout_r0_pop_SD[,,,3]))/2
matout_r0_pop_SD[,,,3]<-matout_r0_pop_SD[,,,3]/sdest

#plot r0
matout_r0_pop2<-matout_r0_pop
matout_r0_pop2[is.na(matout_r0_pop)]<-999 #set extinctions to red
tmp<-plot_cont(arrayout = matout_r0_pop2, xscalslst = log(scalslst,10), xlst = log(tscalelst, 10), nlevels=10, logx=TRUE, logy=TRUE, sqlst = c(sqtmp, 1e6), logxps = logxpos, nstart = 6, revcol = TRUE, dolines=FALSE, doz = TRUE, sqz = sqtmpz, arrayout2 = matout_r0_pop_SD, coltype = 6)

par(mar=c(2,5.5,1,1))
filled.legend(z=matrix(sqtmp), levels=-sqtmp, col=adjustcolor(c(rev(rainbow(sum(sqtmp<0), start=0.55, end=.70)), rev(rainbow(sum(sqtmp>0), start=0, end=0.1))), alpha.f = 0.6), key.axes = axis(4, at = sqtmp, las=2))
axis(2, at=sqtmp, sprintf("%.2f", -sqtmp), las=2)

mtext(text = "levins", side = 4, outer = TRUE, line = -7.5, adj = .94, cex=1.2)
mtext(text = "disturbance", side = 4, outer = TRUE, line = -7.5, adj = 0.74, cex=1.2)
mtext(text = "PSF", side = 4, outer = TRUE, line = -7.5, adj = 0.515, cex=1.2)
mtext(text = "RPS", side = 4, outer = TRUE, line = -7.5, adj = 0.305, cex=1.2)
mtext(text = "neutral", side = 4, outer = TRUE, line = -7.5, adj = .08, cex=1.2)

mtext(text = expression(paste(r[e], ", population")), side = 3, outer = TRUE, line = 0, adj = .14, cex=1.2)
mtext(text = expression(paste(r[0], ", population")), side = 3, outer = TRUE, line = 0, adj = 0.64, cex=1.2)

mtext(text = expression(paste("temporal extent, time steps")), side = 1, outer = TRUE, line = 1.3, cex=1.2, adj = 0.3)
mtext(text = expression(paste("spatial extent, fraction of maximum")), side = 2, outer = TRUE, line = -0.1, cex=1.2, adj = 0.46)

mtext(text = expression(paste(r[e])), side = 3, outer = TRUE, line = 0, adj = .9, cex=1.5)
mtext(text = expression(paste(r[0])), side = 3, outer = TRUE, line = 0, adj = 1.035, cex=1.5, xpd=NA)

dev.off()


#pdf("figures/SUP_FIGURE_sim_continuous_dyn_results_population.pdf", width=3.2, height=8, colormodel = "cmyk")
svg("figures/SUP_FIGURE_sim_continuous_dyn_results_population.svg", width=3.2, height=8)
m<-matrix(nrow=5, 1:5)
m<-cbind(m, 6)
layout(m, widths=c(1,0.35))

par(mar=c(2,3,1,0), oma=c(2.5,2,2,4))
ofs1<-c(0.255, -0.002)
tmp<-plot_cont(matout_eig_tot, log(scalslst,10), log(tscalelst, 10), nlevels=10, logx=TRUE, logy=TRUE, sqlst = sqtmp, logxps = logxpos, nstart = 1)

par(mar=c(2,3,1,0))
filled.legend(z=matrix(sqtmp), levels=-sqtmp, col=adjustcolor(c(rev(rainbow(sum(sqtmp<0), start=0.55, end=.70)), rev(rainbow(sum(sqtmp>0), start=0, end=0.1))), alpha.f = 0.6), key.axes = axis(4, at = sqtmp, -sqtmp, las=2))

mtext(text = "levins", side = 4, outer = TRUE, line = -4.2, adj = .94, cex=1.2)
mtext(text = "disturbance", side = 4, outer = TRUE, line = -4.2, adj = 0.74, cex=1.2)
mtext(text = "PSF", side = 4, outer = TRUE, line = -4.2, adj = 0.515, cex=1.2)
mtext(text = "RPS", side = 4, outer = TRUE, line = -4.2, adj = 0.305, cex=1.2)
mtext(text = "neutral", side = 4, outer = TRUE, line = -4.2, adj = .08, cex=1.2)

mtext(text = expression(paste(r[e], ", community")), side = 3, outer = TRUE, line = 0, adj = .42, cex=1.2)

mtext(text = expression(paste("temporal extent, time steps")), side = 1, outer = TRUE, line = 1.3, cex=1.2, adj = 0.1)
mtext(text = expression(paste("spatial extent, fraction of maximum")), side = 2, outer = TRUE, line = -0.1, cex=1.2, adj = 0.46)


mtext(text = expression(paste(r[e])), side = 3, outer = TRUE, line = 0, adj = .98, cex=1.5)

dev.off()



#CI plots...
#First, get n's needed for significant result
matout_eig_tot_SD<-(abs(matout_eig_tot[,,,3]-matout_eig_tot[,,,4])+abs(matout_eig_tot[,,,3]-matout_eig_tot[,,,2]))/2
matout_eig_pop_SD<-(abs(matout_eig_pop[,,,3]-matout_eig_pop[,,,4])+abs(matout_eig_pop[,,,3]-matout_eig_pop[,,,2]))/2
matout_r0_pop_SD<-(abs(matout_r0_pop[,,,3]-matout_r0_pop[,,,4])+abs(matout_r0_pop[,,,3]-matout_r0_pop[,,,2]))/2

matout_eig_tot_CIplot<-matout_eig_tot
matout_eig_tot_CIplot[,,,3]<-matout_eig_tot_CIplot[,,,3]/matout_eig_tot_SD
matout_eig_pop_CIplot<-matout_eig_pop
matout_eig_pop_CIplot[,,,3]<-matout_eig_pop_CIplot[,,,3]/matout_eig_pop_SD
matout_r0_pop_CIplot<-matout_r0_pop
matout_r0_pop_CIplot[,,,3]<-matout_r0_pop_CIplot[,,,3]/matout_r0_pop_SD

#pdf("figures/SUP_FIGURE_sim_CI_continuous_dyn_results.pdf", width=7, height=8, colormodel = "cmyk")
svg("figures/SUP_FIGURE_sim_CI_continuous_dyn_results.svg", width=6.5, height=8)
zcrit<-abs(qnorm(0.05, 0, 1))
#sqtmp<-c(-5, -2, -1.5, -1, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1, 1.5, 2, 5)
sqtmp<-sort(c(zcrit/sqrt(c(1,2,5,10,30,50)), 0, -zcrit/sqrt(c(1,2,5,10,30,50))))
n_eq<-round((1/(sqtmp/zcrit))^2*sign(sqtmp))
sqtmp<-c(-1000, sqtmp, 1000)
sqtmp<-round(sqtmp, 1)

logxpos<-c(1,2,5,10,20,50,150,200)

m<-matrix(nrow=5, 1:15)
m<-cbind(m, 16)
layout(m, widths=c(1,1,1,0.5))

par(mar=c(2,3,1,2.5), oma=c(2.5,2,2.5,3))
ofs1<-c(0.255, -0.002)
tmp<-plot_cont(matout_eig_tot_CIplot, log(scalslst,10), log(tscalelst, 10), nlevels=length(sqtmp)-1, logx=TRUE, logy=TRUE, sqlst = sqtmp, logxps = logxpos, nstart = 1, dolines = FALSE, coltype = 5)
tmp<-plot_cont(matout_eig_pop_CIplot, log(scalslst,10), log(tscalelst, 10), nlevels=length(sqtmp)-1, logx=TRUE, logy=TRUE, sqlst = sqtmp, logxps = logxpos, nstart = 6, dolines = FALSE, coltype = 5)
tmp<-plot_cont(matout_r0_pop_CIplot, log(scalslst,10), log(tscalelst, 10), nlevels=length(sqtmp)-1, logx=TRUE, logy=TRUE, sqlst = sqtmp, logxps = logxpos, nstart = 11, dolines = FALSE, coltype = 5)

par(mar=c(2,3.5,1,1))
sqtmp_label<-sqtmp
sqtmp_label[1]<-(-2); sqtmp_label[length(sqtmp_label)]<-2
sqtmp_label_text<-sqtmp
sqtmp_label_text[1]<-"< -2"; sqtmp_label_text[length(sqtmp_label_text)]<-"> 2"
n_eq_text<-abs(n_eq)
n_eq_text[is.nan(n_eq_text)]<-"> 50"
n_eq_text<-c("< 1", n_eq_text, "< 1")
filled.legend(z=matrix(sqtmp_label), levels=sqtmp_label, col=adjustcolor(c(rev(rainbow(sum(sqtmp<0)-1, start=0.55, end=.70)), "gray1", "gray1", rev(rainbow(sum(sqtmp>0)-1, start=0, end=0.1))), alpha.f = 0.6), key.axes = axis(4, at = sqtmp_label, sqtmp_label_text, las=2))
axis(2, at = sqtmp_label, n_eq_text, las=2)
text(-0.55, 2.1, expression(paste(italic(n))), xpd=NA, cex=1.5)
text(1.5, 2.1, expression(paste(italic(Z),"-score")), xpd=NA, cex=1.5)

#sqtmp<-seq(-2, 2)
#filled.legend(z=matrix(sqtmp), levels=sqtmp, col=c("blue", "lightblue", "pink", "red"), key.axes = axis(4, at = seq(-1.5, 1.5), labels = c("97.5% < 0", "+1SD < 0", "-1SD > 0", "2.5% > 0"), las=2))

mtext(text = "levins", side = 4, outer = TRUE, line = -8, adj = .94, cex=1.2)
mtext(text = "disturbance", side = 4, outer = TRUE, line = -8, adj = 0.74, cex=1.2)
mtext(text = "PSF", side = 4, outer = TRUE, line = -8, adj = 0.515, cex=1.2)
mtext(text = "RPS", side = 4, outer = TRUE, line = -8, adj = 0.305, cex=1.2)
mtext(text = "neutral", side = 4, outer = TRUE, line = -8, adj = .08, cex=1.2)

mtext(text = expression(paste(r[e], ", community")), side = 3, outer = TRUE, line = 0.2, adj = .065, cex=1.2)
mtext(text = expression(paste(r[e], ", population")), side = 3, outer = TRUE, line = 0.2, adj = .4205, cex=1.2)
mtext(text = expression(paste(r[0], ", population")), side = 3, outer = TRUE, line = 0.2, adj = 0.775, cex=1.2)

mtext(text = expression(paste("temporal extent, time steps")), side = 1, outer = TRUE, line = 1.3, cex=1.2, adj = 0.45)
mtext(text = expression(paste("spatial extent, fraction of maximum")), side = 2, outer = TRUE, line = -0.1, cex=1.2, adj = 0.46)
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
load(paste("output/", flst_dyn[grep("_1.rda", flst_dyn)], sep=""))
matout_dyn<-data.frame(matout_dyn)
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
mtext(expression(paste(italic(r[e]))), 2.5, line=2.5, cex=1.2)
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
tmp<-r0_neut1$out_lst[[1]]$output
tmp[,1]<-tmp[,1]+max(tmp1[,1])
tmp1<-rbind(out1$output, tmp1, tmp)

tmp2<-r0_neut2$out0_lst[[1]]$output
tmp2[,1]<-tmp2[,1]+max(out2$output[,1])
tmp<-r0_neut2$out_lst[[1]]$output
tmp[,1]<-tmp[,1]+max(tmp2[,1])
tmp2<-rbind(out2$output, tmp2, tmp)

sbs<--c(1:250)
plot(tmp1[sbs,1]+1-100, tmp1[sbs,2]/out1$plotdata$ngrid, type="l", col=collst3[2], lwd=1.5, xlab="", ylab="", axes=F, xaxs="i", ylim=c(0, 0.4))
axis(1); axis(2, las=2); box()
abline(v=c(300, 500)-100, lty=2); abline(h=0, lty=3)
arrows(tmp1[300,1]+1-100, tmp1[300,2]/out1$plotdata$ngrid, tmp1[301,1]+1-100, tmp1[301,2]/out1$plotdata$ngrid, length = 0.12, lwd=1.5, angle = 20)
arrows(tmp1[500,1]+1-100, tmp1[500,2]/out1$plotdata$ngrid-0.01, tmp1[501,1]+1-100, tmp1[501,2]/out1$plotdata$ngrid+0.02, length = 0.12, lwd=1.5, angle = 20)

put.fig.letter("a.", "topleft", offset=ofs2, cex=1.2)

mtext(text = expression(paste("simulation time")), side = 1, outer = F, line = 2.5, cex=1.2, adj = 0.42)
mtext(text = expression(paste("abundance")), side = 2, outer = F, line = 2.5, cex=1.2, adj = 0.42)

plot(tmp2[sbs,1]+1-100, tmp2[sbs,2]/out1$plotdata$ngrid, type="l", col=collst3[2], lwd=1.5, xlab="", ylab="", axes=F, xaxs="i", ylim=c(0, 0.4))
axis(1); axis(2, las=2); box()
abline(v=c(300, 500)-100, lty=2); abline(h=0, lty=3)
arrows(tmp2[300,1]+1-100, tmp2[300,2]/out1$plotdata$ngrid+0.01, tmp2[301,1]+1-100, tmp2[301,2]/out1$plotdata$ngrid-0.02, length = 0.12, lwd=1.5, angle = 20)
arrows(tmp2[500,1]+1-100, tmp2[500,2]/out1$plotdata$ngrid-0.01, tmp2[501,1]+1-100, tmp2[501,2]/out1$plotdata$ngrid+0.02, length = 0.12, lwd=1.5, angle = 20)

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









