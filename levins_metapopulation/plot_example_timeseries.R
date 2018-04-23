error
rm(list=ls())
setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

#load functions
require(rEDM)
source("run_metapopulation_wrapper.R")
require(RColorBrewer)
source("~/Dropbox/Rfunctions/figure_functions.R")
source("~/Dropbox/Rfunctions/logit_funs.R")

#plotting offsets
ofs3<-c(0.075, -0.090)

#starting parameters
gridout<-makegrid(xlng = 100, ylng = 100)
xfac<-3
xfac_fast<-7
ptb<-0.33
collst<-adjustcolor(c("grey51", brewer.pal(4, "Set1")), alpha.f = 0.7)
collst2<-adjustcolor(brewer.pal(3, "Dark2"), alpha.f = 0.7)

#simulation lengths
tinit<-300
tsim<-200


############################################################
# Run global models
############################################################
set.seed(180110)

##### Levins metapopulation model
clst_meta = c(0.15, 0.35)*xfac
mlst_meta = rep(0.1, length(clst_meta))*xfac

population_meta<-populate(gridout, nlst = floor(getceq(clst_meta, mlst_meta)*prod(gridout$lng)),
                          clst = clst_meta, radlst = Inf, mlst = mlst_meta)

set.seed(171205)
out_meta<-run_metapopulation(tmax=tinit, gridout = gridout, population = population_meta, talktime = 0)

eig_meta1<-estimate_eqreturn(out_meta, simtime=tsim, runtype="metapopulation", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE)
eig_meta2<-estimate_eqreturn(out_meta, simtime=tsim, runtype="metapopulation", replace_perturb = 0, talktime=0, prtb=ptb, doplot = FALSE)

#r0_meta<-estimate_rarereturn(eig_meta1$out_lst0, simtime=tsim, burnin=tsim, runtype="metapopulation", doplot = FALSE)
r0_meta<-estimate_rarereturn(out_meta, simtime=tsim, burnin=tsim, runtype="metapopulation", doplot = FALSE)

set.seed(171205)
out_meta_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_meta, talktime = 0)
  
getEmeta<-getE(out_meta_long, Elst = 2:10)
E_meta<-getEmeta$Eout
E_meta_tot<-getEmeta$Eout_tot
  
invar_meta<-estimate_invar(out = out_meta_long, E=E_meta, burnin=0, doplot=FALSE, Etot = E_meta_tot, niter = 0)
beta_meta<-beta_estimate(out=out_meta, outlng = out_meta_long, Emat = E_meta, eigout = eig_meta2, r0out = r0_meta, burnin = 100)

##### Hubbell neutal model
set.seed(1712010)
clst_neut<-rep(0.5, 2)*xfac
mlst_neut<-rep(0.1, length(clst_neut))*xfac_fast
tmp<-abs(getceq(clst_neut, mlst_neut))
population_neut<-populate(gridout, nlst = round(rep(unique(tmp[tmp>0])/length(clst_neut), length(clst_neut))*prod(gridout$lng)),
                          clst = clst_neut, radlst = Inf, mlst = mlst_neut)

set.seed(171206)
out_neut<-run_metapopulation(tmax=tinit, gridout = gridout, population = population_neut, talktime = 0, runtype = "neutral")

set.seed(171210)
eig_neut1<-estimate_eqreturn(out_neut, simtime=tsim, runtype="neutral", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE)
eig_neut2<-estimate_eqreturn(out_neut, simtime=tsim, runtype="neutral", replace_perturb = 0, talktime=0, prtb=ptb, doplot = FALSE)

#r0_neut<-estimate_rarereturn(out = eig_neut1$out_lst0, simtime=tsim, burnin=tsim, runtype="neutral", doplot = FALSE)
r0_neut<-estimate_rarereturn(out = out_neut, simtime=tsim, burnin=tsim, runtype="neutral", doplot = FALSE)

set.seed(171206)
out_neut_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_neut, talktime = 0, runtype = "neutral")
getEneut<-getE(out_neut_long, Elst = 2:10)
E_neut<-getEneut$Eout
E_neut_tot<-getEneut$Eout_tot

invar_neut<-estimate_invar(out_neut_long, E=E_neut, burnin=0, doplot=FALSE, Etot = E_neut_tot, niter = 0)

beta_neut<-beta_estimate(out=out_neut, outlng = out_neut_long, Emat = E_neut, eigout = eig_neut2, r0out = r0_neut, burnin = 100)

##### Disturbance model
clst_dist= c(0.145, 0.2)*xfac
mlst_dist = rep(0.1, length(clst_dist))*xfac
distlst<-c(0.95, 0)
prtfrq<-50

population_dist<-populate(gridout, nlst = rep(floor(prod(gridout$lng)/length(clst_dist)*0.8), length(clst_dist)),
                          clst = clst_dist, radlst = Inf, mlst = mlst_dist)
set.seed(171217)
out_dist<-run_metapopulation(tmax=tinit, gridout = gridout, population = population_dist, talktime = 0, runtype = "disturbance", prt = distlst,  prtfrq = prtfrq)

out_dist_0<-rerunrun_metapopulation(out=out_dist, tmax=0, talktime = 0, runtype = "metapopulation", perturb = distlst, replace_perturb = 0)

eig_dist1<-estimate_eqreturn(out_dist_0, simtime=tsim, runtype="disturbance", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE, prt = distlst,  prtfrq = prtfrq)
eig_dist2<-estimate_eqreturn(out_dist_0, simtime=tsim, runtype="disturbance", replace_perturb = 0, talktime=0, prtb=ptb, doplot = FALSE, prt = distlst,  prtfrq = prtfrq)

#out_dist_0<-rerunrun_metapopulation(out=eig_dist1$out_lst0, tmax=0, talktime = 0, runtype = "metapopulation", perturb = distlst, replace_perturb = 0)
out_dist_0<-rerunrun_metapopulation(out=out_dist, tmax=0, talktime = 0, runtype = "metapopulation", perturb = distlst, replace_perturb = 0)
r0_dist<-estimate_rarereturn(out_dist_0, simtime=tsim, burnin=tsim, runtype="disturbance", doplot = FALSE, prt = distlst,  prtfrq = prtfrq)

set.seed(171217)
out_dist_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_dist, talktime = 0, runtype = "disturbance", prt = distlst,  prtfrq = prtfrq)
getEdist<-getE(out_dist_long, Elst = 2:10)
E_dist<-getEdist$Eout
E_dist_tot<-getEdist$Eout_tot

invar_dist<-estimate_invar(out_dist_long, E=E_dist, burnin=100, doplot=FALSE, Etot = E_dist_tot, niter = 0)
out_dist_nodist<-run_metapopulation(tmax=200, gridout = gridout, population = population_dist, talktime = 0, runtype = "metapopulation")
beta_dist<-beta_estimate(out=out_dist, outlng = out_dist_long, Emat = E_dist, eigout = eig_dist2, r0out = r0_dist, burnin = 100)

##### Try psf model
clst_psf= c(1.5, 1)*xfac
mlst_psf = rep(0.1, length(clst_psf))*xfac

population_psf<-populate(gridout, nlst = rep(floor(prod(gridout$lng)/length(clst_psf)*0.8), length(clst_psf)),
                         clst = clst_psf, radlst = Inf, mlst = mlst_psf)
set.seed(180108)
out_psf<-run_metapopulation(tmax=tinit, gridout = gridout, population = population_psf, talktime = 0, runtype = "psf")

eig_psf1<-estimate_eqreturn(out_psf, simtime=tsim, runtype="psf", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE)
eig_psf2<-estimate_eqreturn(out_psf, simtime=tsim, runtype="psf", replace_perturb = 0, talktime=0, prtb=ptb, doplot = FALSE)

#r0_psf<-estimate_rarereturn(eig_psf1$out_lst0, simtime=tsim, burnin=tsim, runtype="psf", doplot = FALSE)
r0_psf<-estimate_rarereturn(out_psf, simtime=tsim, burnin=tsim, runtype="psf", doplot = FALSE)

set.seed(180108)
out_psf_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_psf, talktime = 0, runtype = "psf")
getEpsf<-getE(out_psf_long, Elst = 2:10)
E_psf<-getEpsf$Eout
E_psf_tot<-getEpsf$Eout_tot
invar_psf<-estimate_invar(out_psf_long, E=E_psf, burnin=100, doplot=FALSE, Etot = E_psf_tot, niter = 0)
beta_psf<-beta_estimate(out=out_psf, outlng = out_psf_long, Emat = E_psf, eigout = eig_psf2, r0out = r0_psf, burnin = 10)

##### Try rps model
intmat_rps<-rbind(c(1,0,0,1),
                  c(1,1,0,0),
                  c(0,1,1,0),
                  c(0,0,1,1))

clst_rps<-rep(0.4, 4)*xfac*0.4
mlst_rps<-rep(0.08, length(clst_rps))*xfac_fast*0.4
tmp<-abs(getceq(clst_rps, mlst_rps))
population_rps<-populate(gridout, nlst = round(rep(unique(tmp[tmp>0])/length(clst_rps), length(clst_rps))*prod(gridout$lng)),
                         clst = clst_rps, radlst = Inf, mlst = mlst_rps)
set.seed(180108)
out_rps<-run_metapopulation(tmax=tinit, gridout = gridout, population = population_rps, talktime = 0, runtype = "rps", compmat = intmat_rps)

eig_rps1<-estimate_eqreturn(out_rps, simtime=tsim, runtype="rps", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE)
eig_rps2<-estimate_eqreturn(out_rps, simtime=tsim, runtype="rps", replace_perturb = 0, talktime=0, prtb=ptb, doplot = FALSE)

#r0_rps<-estimate_rarereturn(out = eig_rps1$out_lst0, simtime=tsim, burnin=tsim, runtype="rps", doplot = FALSE)
r0_rps<-estimate_rarereturn(out = out_rps, simtime=tsim, burnin=tsim, runtype="rps", doplot = FALSE)

set.seed(180104)
out_rps_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_rps, talktime = 0, runtype = "rps", compmat = intmat_rps)
getErps<-getE(out_rps_long, Elst = 2:10)
E_rps<-getErps$Eout
E_rps_tot<-getErps$Eout_tot

invar_rps<-estimate_invar(out_rps_long, E=E_rps, burnin=0, doplot=FALSE, Etot = E_rps_tot, niter = 0)

beta_rps<-beta_estimate(out=out_rps, outlng = out_rps_long, Emat = E_rps, eigout = eig_rps2, r0out = r0_rps, burnin = 10)

############################################################
# Plot examples of each model
############################################################
#out<-out_meta; eigout<-eig_meta2; r0out<-r0_meta; collst<-collst; burnin=100; doceq=TRUE; plotpos=1

modplotfun<-function(out, eigout, r0out, collst, burnin=0, doceq=0, plotpos=1, atsq=0, totabund=F, doaxis1=T, figlet=1, ofs1=c(0,0), fcx=0, inbetweenfun=NULL, ...) {
  #original fxn
  abunds<-out$output
  if(burnin>0) {
    abunds<-abunds[-c(1:burnin),]
    abunds[,1]<-abunds[,1]-burnin
  }
  if(totabund) {
    abunds<-cbind(abunds[,1], rowSums(abunds[,-1]), abunds[,-1])
  } else {
    abunds<-cbind(abunds[,1], abunds[,-1])
  }
  
  mxt<-ceiling(max(abunds[,1]))
  mxt_eig<-ceiling(max(eigout$out_lst[[plotpos]]$output[,1]))
  mxt_r0_0<-ceiling(max(r0out$out0_lst[[plotpos]]$output[,1]))
  mxt_r0<-ceiling(max(r0out$out_lst[[plotpos]]$output[,1]))
  
  #eig & r0
  tmp_eig<-eigout$out_lst[[plotpos]]$output
  tmp_eig[,1]<-tmp_eig[,1]+mxt
  if(totabund) {
    tmp_eig<-cbind(tmp_eig[,1], rowSums(tmp_eig[,-1]), tmp_eig[,-1])
  } else {
    tmp_eig<-cbind(tmp_eig[,1], tmp_eig[,-1])
  }
  
  tmp_eig_0<-eigout$out_lst0$output
  tmp_eig_0[,1]<-tmp_eig_0[,1]+mxt
  if(totabund) {
    tmp_eig_0<-cbind(tmp_eig_0[,1], rowSums(tmp_eig_0[,-1]), tmp_eig_0[,-1])
  } else {
    tmp_eig_0<-cbind(tmp_eig_0[,1], tmp_eig_0[,-1])
  }
  
  tmp_r0_0<-r0out$out0_lst[[plotpos]]$output
  tmp_r0_0[,1]<-tmp_r0_0[,1]+mxt#+mxt_eig
  if(totabund) {
    tmp_r0_0<-cbind(tmp_r0_0[,1], rowSums(tmp_r0_0[,-1]), tmp_r0_0[,-1])
  } else {
    tmp_r0_0<-cbind(tmp_r0_0[,1], tmp_r0_0[,-1])
  }
  
  tmp_r0<-r0out$out_lst[[plotpos]]$output
  tmp_r0[,1]<-tmp_r0[,1]+mxt+mxt_r0_0#+mxt_eig
  if(totabund) {
    tmp_r0<-cbind(tmp_r0[,1], rowSums(tmp_r0[,-1]), tmp_r0[,-1])
  } else {
    tmp_r0<-cbind(tmp_r0[,1], tmp_r0[,-1])
  }
  
  #combine
  abunds1<-rbind(abunds, tmp_eig, rep(NA, ncol(abunds)), tmp_eig[nrow(tmp_eig),])
  abunds2<-rbind(abunds, tmp_r0_0, tmp_r0)
  
  pabunds1<-abunds1/out$plotdata$ngrid
  pabunds1[,1]<-abunds1[,1]
  
  pabunds2<-abunds2/out$plotdata$ngrid
  pabunds2[,1]<-abunds2[,1]
  
  ptmp_eig_0<-tmp_eig_0/out$plotdata$ngrid
  ptmp_eig_0[,1]<-tmp_eig_0[,1]
  
  #plot perturb
  suppressWarnings(matplot(pabunds1[,1], pabunds1[,-1], type="l", lty=1, col=collst, lwd=c(1.5, rep(1.5, ncol(pabunds1)-2)), xlab="", ylab="", xaxs="i", axes=F, ylim=range(c(ptmp_eig_0[,-1], pabunds1[,-1]), na.rm=T), ...))
  abline(v=c(mxt), lty=2); abline(h=c(0,1), lty=3, lwd=1)
  if(sum(doceq)==2) {
    if(totabund) {
      abline(h=c(sum(out$plotdata$ceq), out$plotdata$ceq), lty=3, col=collst, lwd=1)
    } else {
      abline(h=c(out$plotdata$ceq), lty=3, col=collst, lwd=1)
    }
  } else if(sum(doceq)==1) {
    if(totabund) {
      abline(h=sum(out$plotdata$ceq), lty=3, col=collst, lwd=1)
    }
  }
  if(doaxis1) {
    if(sum(atsq)==0) {
      axis(1, cex.axis=1.6)
    } else {
      axis(1, at=seq(min(pabunds1[,1], na.rm=T), max(pabunds1[,1]-0.1, na.rm=T), length=3), labels=atsq[1:3], cex.axis=1.6)
    }
  }
  axis(2, las=2, cex.axis=1.6); box()
  
  #add in perturbation null
  matlines(ptmp_eig_0[,1], ptmp_eig_0[,-1], lty=2, col=collst, lwd=c(1, rep(1, ncol(pabunds1)-2)))
  
  #Add in disturbance lines
  ap1<-tmp_eig_0[1,plotpos+1+totabund]/out$plotdata$ngrid
  ap2<-tmp_eig[1,plotpos+1+totabund]/out$plotdata$ngrid
  ap1<-max(c(ap1, ap2+diff(range(pabunds1[,-1],na.rm=T))*0.1))
  arrows(mxt, ap1, mxt, ap2, lend=2, length = 0.06, col=1, lwd=1.5)
  
  put.fig.letter(paste(letters[figlet], ".", sep=""), "topleft", offset=ofs1, cex=fcx)
  
  
  if(!is.null(inbetweenfun)) {
    eval(parse(text=inbetweenfun))
  }
  
  #plot r0
  suppressWarnings(matplot(pabunds2[,1], pabunds2[,-1], type="l", lty=1, col=collst, lwd=c(1.5, rep(1.5, ncol(pabunds2)-2)), xlab="", ylab="", xaxs="i", axes=F, ...))
  
  abline(v=c(mxt, mxt+mxt_r0_0), lty=2); abline(h=c(0,1), lty=3, lwd=1)
  
  
  #Add in disturbance lines
  ap1<-tmp_eig[nrow(tmp_eig),plotpos+1+totabund]/out$plotdata$ngrid
  ap2<-0
  ap1<-max(c(ap1, ap2+diff(range(pabunds2[,-1],na.rm=T))*0.1))
  arrows(mxt, ap1, mxt, ap2, lend=2, length = 0.06, col=1, lwd=1.5)
  
  ap1<-tmp_r0[1,plotpos+1+totabund]/out$plotdata$ngrid
  ap2<-0
  ap1<-max(c(ap1, ap2+diff(range(pabunds2[,-1],na.rm=T))*0.1))
  arrows(mxt+mxt_r0_0, ap2, mxt+mxt_r0_0, ap1, lend=2, length = 0.06, col=1, lwd=1.5)
  
  
  if(sum(doceq)==2) {
    if(totabund) {
      abline(h=c(sum(out$plotdata$ceq), out$plotdata$ceq), lty=3, col=collst, lwd=1)
    } else {
      abline(h=c(out$plotdata$ceq), lty=3, col=collst, lwd=1)
    }
  } else if(sum(doceq)==1) {
    if(totabund) {
      abline(h=sum(out$plotdata$ceq), lty=3, col=collst, lwd=1)
    }
  }
  if(doaxis1) {
    if(sum(atsq)==0) {
      axis(1, cex.axis=1.6)
    } else {
      axis(1, at=seq(min(pabunds2[,1], na.rm=T), max(pabunds2[,1]-0.1, na.rm=T), length=4), labels=atsq[1:4], cex.axis=1.6)
    }
  }
  axis(2, las=2, cex.axis=1.6); box()
  
  put.fig.letter(paste(letters[figlet+1], ".", sep=""), "topleft", offset=ofs1, cex=fcx)
  
  return(list(pabunds1=pabunds1, pabunds2=pabunds2))
}


pdf("figures/FIGURE_model_examples.pdf", width=9, height=10, colormodel = "cmyk")
m<-t(matrix(1:10, nrow=2))
layout(m)
par(mar=c(0.2,3,0,0.2), oma=c(4.2,4,7,2.5))
atsq<-seq(0, 600, by=200)
fcx<-2

ofs1<-c(0.112, -0.04)

tmp<-modplotfun(out=out_meta, eigout=eig_meta2, r0out=r0_meta, collst=collst[-1], burnin=100, doceq=2, atsq=atsq, doaxis1 = F, figlet=1, ofs1=ofs1, fcx=fcx)
tmp<-tmp$pabunds2

#label perturbations
mxt<-max(tmp[,-1], na.rm=T)+diff(range(tmp[,-1], na.rm=T))*0.08
text(-370, mxt, "peturbation", xpd=NA, srt=40, adj = c(0,0), cex=2)
text(200, mxt, "removal", xpd=NA, srt=40, adj = c(0,0), cex=2)
text(400, mxt, "invasion", xpd=NA, srt=40, adj = c(0,0), cex=2)

tmp<-modplotfun(out=out_dist, eigout=eig_dist2, r0out=r0_dist, collst=collst[c(3,2)], burnin=100, doceq=1, plotpos = 2, atsq=atsq, doaxis1 = F, figlet=3, ofs1=ofs1, fcx=fcx, inbetweenfun='abline(v=c(50, 100, 150, 250, 300, 350, 450, 500, 550), lty=3, col="black")')
abline(v=c(50, 100, 150, 250, 300, 350, 450, 500, 550), lty=3, col="black")

tmp<-modplotfun(out=out_psf, eigout=eig_psf2, r0out=r0_psf, collst=collst[-1], burnin=100, doceq=0, atsq=atsq, doaxis1 = F, figlet=5, ofs1=ofs1, fcx=fcx)

tmp<-modplotfun(out=out_rps, eigout=eig_rps2, r0out=r0_rps, collst=collst[-1], burnin=100, doceq=1, atsq=atsq, doaxis1 = F, figlet=7, ofs1=ofs1, fcx=fcx)

ofs1<-c(0.112, -0.055)
tmp<-modplotfun(out=out_neut, eigout=eig_neut2, r0out=r0_neut, collst=collst[-1], burnin=100, doceq=1, atsq=atsq, figlet=9, ofs1=ofs1, fcx=fcx)

mtext("simulation time", 1, line=2.7, cex=1.5, outer = T)
mtext("species or community abundance", 2, line=1.5, cex=1.5, outer = T)

mtext(text = "levins", side = 4, outer = TRUE, line = 0.8, adj = .91+0.025, cex=1.4)
mtext(text = "disturbance", side = 4, outer = TRUE, line = 0.8, adj = 0.715+0.02, cex=1.4)
mtext(text = "PSF", side = 4, outer = TRUE, line = 0.8, adj = 0.49+0.015, cex=1.4)
mtext(text = "RPS", side = 4, outer = TRUE, line = 0.8, adj = 0.278+0.02, cex=1.4)
mtext(text = "neutral", side = 4, outer = TRUE, line = 0.8, adj = .06+0.005, cex=1.4)
dev.off()



############################################################
# Plot examples of eigen and invasion
############################################################
#out<-out_meta; eigout<-eig_meta2; r0out<-r0_meta; collst<-collst; burnin=200; burnine=100; dburnin=100; plotpos=1; doceq=0

statsplotfun<-function(out, eigout, r0out, collst, burnin=0, burnine=0, dburnin=0, plotpos=1, doceq=0, ...) {
  ofs2<-c(0.3, -0.06)
  ofs3<-c(0.3, -0.015)
  fcx<-2
  
  
  #original fxn
  abunds<-out$output
  if(burnin>0) {
    abunds<-abunds[-c(1:burnin),]
    abunds[,1]<-abunds[,1]-burnin
  }
  abunds<-cbind(abunds[,1], rowSums(abunds[,-1]), abunds[,-1])
  
  mxt<-ceiling(max(abunds[,1]))
  mxt_eig<-ceiling(max(eigout$out_lst[[plotpos]]$output[,1]))
  mxt_r0_0<-ceiling(max(r0out$out0_lst[[plotpos]]$output[,1]))
  mxt_r0<-ceiling(max(r0out$out_lst[[plotpos]]$output[,1]))
  
  #eig & r0
  tmp_eig<-eigout$out_lst[[plotpos]]$output
  tmp_eig[,1]<-tmp_eig[,1]+mxt
  tmp_eig<-cbind(tmp_eig[,1], rowSums(tmp_eig[,-1]), tmp_eig[,-1])
  
  tmp_eig_0<-eigout$out_lst0$output
  tmp_eig_0[,1]<-tmp_eig_0[,1]+mxt
  tmp_eig_0<-cbind(tmp_eig_0[,1], rowSums(tmp_eig_0[,-1]), tmp_eig_0[,-1])
  
  tmp_r0_0<-r0out$out0_lst[[plotpos]]$output
  tmp_r0_0[,1]<-tmp_r0_0[,1]+mxt+mxt_eig
  tmp_r0_0<-cbind(tmp_r0_0[,1], rowSums(tmp_r0_0[,-1]), tmp_r0_0[,-1])
  
  tmp_r0<-r0out$out_lst[[plotpos]]$output
  tmp_r0[,1]<-tmp_r0[,1]+mxt+mxt_eig+mxt_r0_0
  tmp_r0<-cbind(tmp_r0[,1], rowSums(tmp_r0[,-1]), tmp_r0[,-1])
  
  #combine
  abunds1<-rbind(abunds, tmp_eig[1:burnine,])
  pabunds1<-abunds1/out$plotdata$ngrid
  pabunds1[,1]<-abunds1[,1]
  
  abunds2<-rbind(tmp_r0_0[-c(1:burnine),], tmp_r0)
  abunds2[,1]<-abunds2[,1]
  pabunds2<-abunds2/out$plotdata$ngrid
  pabunds2[,1]<-abunds2[,1]
  
  ptmp_eig_0<-tmp_eig_0[1:burnine,]/out$plotdata$ngrid
  ptmp_eig_0[,1]<-tmp_eig_0[1:burnine,1]
  
  #Plotting environment
  m<-rbind(c(1,4),
           c(2,4),
           c(3,5))
  layout(m)
  par(mar=c(3,5,2,2), oma=c(2,2,0,0))
  
  ##### plot eigen deviation
  #plot function
  suppressWarnings(matplot(pabunds1[,1]+dburnin, pabunds1[,plotpos+2], type="l", lty=1, col=collst[plotpos+1], lwd=1.5, xlab="", ylab="", xaxs="i", axes=F, xlim=c(floor(min(pabunds1[,1],na.rm=T)), ceiling(max(pabunds1[,1],na.rm=T)))+dburnin, ylim=range(c(pabunds1[,plotpos+2]), ptmp_eig_0[,plotpos+2], na.rm=T), ...))
  put.fig.letter("a.", "topleft", offset=ofs2, cex=fcx)
  
  abline(v=c(mxt)+dburnin, lty=2); abline(h=c(0,1), lty=3, lwd=1)
  abline(h=c(out$plotdata$ceq[plotpos]), lty=3, col=collst[plotpos+1], lwd=1)
  axis(1, cex.axis=1.6)
  axis(2, las=2, cex.axis=1.6); box()
  if(sum(doceq)==2) {
    abline(h=c(sum(out$plotdata$ceq), out$plotdata$ceq[plotpos]), lty=3, col=collst[plotpos+1], lwd=1)
  } else if(sum(doceq)==1) {
    abline(h=sum(out$plotdata$ceq), lty=3, col=collst, lwd=1)
  }
  
  matlines(ptmp_eig_0[,1]+dburnin, ptmp_eig_0[,plotpos+2], lty=2, col=collst[plotpos+1], lwd=1.5)
  
  ap1<-tmp_eig_0[1,plotpos+2]/out$plotdata$ngrid
  ap2<-tmp_eig[1,plotpos+2]/out$plotdata$ngrid
  arrows(mxt+dburnin, ap1, mxt+dburnin, ap2, lend=2, length = 0.06, col=1, lwd=1.5)
  
  mtext("simulation time", 1, line=2.5, cex=1.2, outer = F)
  mtext("abundance", 2, line=4.5, cex=1.2, outer = F)
  
  
  #plot distance
  dst<-abs(pabunds1[-c(1:mxt),plotpos+2]-ptmp_eig_0[,plotpos+2])
  suppressWarnings(matplot(pabunds1[-c(1:mxt),1]+dburnin, dst, type="l", lty=1, col=collst[plotpos+1], lwd=1.5, xlab="", ylab="", xaxs="i", axes=F, xlim=c(floor(mxt), ceiling(max(pabunds1[,1],na.rm=T)))+dburnin, ...))
  put.fig.letter("b.", "topleft", offset=ofs2, cex=fcx)
  
  abline(v=c(mxt)+dburnin, lty=2); abline(h=c(0), lty=3, lwd=1)
  axis(1, cex.axis=1.6)
  axis(2, las=2, cex.axis=1.6); box()
  
  mtext("simulation time", 1, line=2.5, cex=1.2, outer = F)
  mtext("distance", 2, line=4.5, cex=1.2, outer = F)
  
  #plot eigenvalue
  eigest<-log(dst[-1]/dst[-length(dst)])
  sbs<-is.finite(eigest)
  eigest_tot<-cumsum(eigest[sbs])*(1/(1:(sum(sbs))))
  
  suppressWarnings(matplot((pabunds1[-c(1:(mxt+1)),1]+dburnin)[sbs], eigest_tot, type="l", lty=1, col=collst[plotpos+1], lwd=1.5, xlab="", ylab="", xaxs="i", axes=F, xlim=c(floor(mxt), ceiling(max(pabunds1[,1],na.rm=T)))+dburnin, ylim=range(c(eigest_tot, 0), na.rm=T), ...))
  put.fig.letter("c.", "topleft", offset=ofs2, cex=fcx)
  
  abline(v=c(mxt), lty=2); abline(h=c(0), lty=3, lwd=1)
  axis(1, cex.axis=1.6)
  axis(2, las=2, cex.axis=1.6); box()
  
  mtext("simulation time", 1, line=2.5, cex=1.2, outer = F)
  mtext(expression(paste(italic(lambda))), 2, line=3.5, cex=1.2, outer = F)
  
  
  ##### plot invasion
  #plot function
  suppressWarnings(matplot(pabunds2[,1]+dburnin, pabunds2[,-c(1:2)], type="l", lty=1, col=collst[-1], lwd=1.5, xlab="", ylab="", xaxs="i", axes=F, xlim=c(floor(min(pabunds2[,1],na.rm=T)), ceiling(max(pabunds2[,1],na.rm=T)))+dburnin, ...))
  put.fig.letter("d.", "topleft", offset=ofs3, cex=fcx)
  
  abline(v=c(burnine, burnine+mxt_r0_0)+mxt_eig+dburnin, lty=2); abline(h=c(0,1), lty=3, lwd=1)
  if(sum(doceq)==2) {
    abline(h=c(sum(out$plotdata$ceq), out$plotdata$ceq), lty=3, col=collst, lwd=1)
  } else if(sum(doceq)==1) {
    abline(h=sum(out$plotdata$ceq), lty=3, col=collst, lwd=1)
  }
  axis(1, cex.axis=1.6, at=seq(0, 800, by=50), seq(0, 800, by=50)-mxt_eig)
  axis(2, las=2, cex.axis=1.6); box()
  
  ap1<-tmp_eig[nrow(tmp_eig),plotpos+2]/out$plotdata$ngrid
  ap2<-0
  arrows(burnine+mxt_eig+dburnin, ap1, burnine+mxt_eig+dburnin, ap2, lend=2, length = 0.06, col=1, lwd=1.5)
  
  ap1<-tmp_r0[1,plotpos+2]/out$plotdata$ngrid
  ap2<-0
  ap1<-max(c(ap1, ap2+diff(range(pabunds2[,-1],na.rm=T))*0.1))
  arrows(burnine+mxt_eig+mxt_r0_0+dburnin, ap2, burnine+mxt_eig+mxt_r0_0+dburnin, ap1, lend=2, length = 0.06, col=1, lwd=1.5)
  
  mtext("simulation time", 1, line=2.5, cex=1.2, outer = F)
  mtext("abundance", 2, line=3.5, cex=1.2, outer = F)
  
  #plot r0
  grw<-tmp_r0[,plotpos+2]/out$plotdata$ngrid
  grw<-log(grw[-1]/grw[-length(grw)])
  sbs<-is.finite(grw)
  grw_tot<-cumsum(grw[sbs])*(1/(1:sum(sbs)))
  
  suppressWarnings(matplot(((1:length(grw_tot))+burnine+mxt_eig+mxt_r0_0+dburnin)[sbs], grw_tot, type="l", lty=1, col=collst[plotpos+1], lwd=1.5, xlab="", ylab="", xaxs="i", axes=F, ylim=c(0, max(grw_tot)),...))
  put.fig.letter("e.", "topleft", offset=ofs2, cex=fcx)
  
  axis(1, cex.axis=1.6, at=seq(0, 800, by=50), seq(0, 800, by=50)-mxt_eig)
  axis(2, las=2, cex.axis=1.6); box()
  abline(h=0, lty=3)
  
  mtext("simulation time", 1, line=2.5, cex=1.2, outer = F)
  mtext(expression(paste(italic(r[0]))), 2, line=3.5, cex=1.2, outer = F)
  
  return(list(pabunds1=pabunds1, pabunds2=pabunds2))
}



pdf("figures/SUP_FIGURE_eig_r0_examples.pdf", width=5, height=6, colormodel = "cmyk")
tmp<-statsplotfun(out=out_meta, eigout=eig_meta2, r0out=r0_meta, collst=collst, burnin=200, burnine=100, dburnin=100, plotpos=1, doceq = 2)
dev.off()

pdf("figures/SUP_FIGURE_eig_r0_examples_allmodels.pdf", width=5, height=6, colormodel = "cmyk")
tmp<-statsplotfun(out=out_meta, eigout=eig_meta2, r0out=r0_meta, collst=collst, burnin=200, burnine=100, dburnin=100, plotpos=1, doceq = 2)
tmp<-statsplotfun(out=out_dist, eigout=eig_dist2, r0out=r0_dist, collst=collst, burnin=200, burnine=100, dburnin=100, plotpos=1, doceq = 1)  
tmp<-statsplotfun(out=out_psf, eigout=eig_psf2, r0out=r0_psf, collst=collst, burnin=200, burnine=100, dburnin=100, plotpos=1, doceq = 1)  
tmp<-statsplotfun(out=out_rps, eigout=eig_rps2, r0out=r0_rps, collst=collst, burnin=200, burnine=100, dburnin=100, plotpos=1, doceq = 1)  
tmp<-statsplotfun(out=out_neut, eigout=eig_neut2, r0out=r0_neut, collst=collst, burnin=200, burnine=100, dburnin=100, plotpos=1, doceq = 1)  
dev.off()



############################################################
# Plot examples of CV and beta
############################################################
if(FALSE) {
  outlong<-out_meta_long; invarout<-invar_meta; betaout<-beta_meta; collst<-collst; collst2<-collst2; burnin=200; plotpos=1; doceq=2; cleanupbeta<-2
  
  CVplotfun<-function(outlong, invarout, betaout, collst, collst2, burnin=0, plotpos=1, doceq=0, cleanupbeta=FALSE, ...) {
    #Plotting environment
    m<-rbind(c(1,1,1,3,3),
             c(2,2,2,3,3))
    layout(m)
    par(mar=c(4,5,2,1), oma=c(0,0,0,0))
    
    ofs2<-c(0.19, -0.06)
    ofs3<-c(0.29, -0.0225)
    
    ##### Invar
    #plot time series
    invarlst<-outlong$output[-c(1:burnin),plotpos+1]/outlong$plotdata$ngrid
    matplot(outlong$output[-c(1:burnin),1], invarlst, type="l", lty=1, col=collst[-1], lwd=1.5, xlab="", ylab="", xaxs="i",
            ylim=range(invarlst)*c(0.98, 1), axes=FALSE, ...)
    axis(1); axis(2, las=2); box()
    if(sum(doceq)==2) {
      abline(h=c(sum(outlong$plotdata$ceq), outlong$plotdata$ceq), lty=3, col=collst, lwd=1)
    } else if(sum(doceq)==1) {
      abline(h=sum(outlong$plotdata$ceq), lty=3, col=collst, lwd=1)
    }
    mtext("simulation time", 1, line=2.8, cex=1.1)
    mtext("relative abundance", 2, line=3.4, cex=1.1)
    put.fig.letter("a.", "topleft", offset=ofs2, cex=1.6)
    
    tmp<-outlong$output[,2]/outlong$plotdata$ngrid
    tmptm<-outlong$output[,1]
    
    startlst<-(nrow(outlong$output)-burnin)*c(0.2, 0.3, 0.5, 0.7, 0.8)+burnin
    
    y11<-min(tmp[tmptm>startlst[1] & tmptm<startlst[2]]); y12<-max(tmp[tmptm>startlst[1] & tmptm<startlst[2]])
    y21<-min(tmp[tmptm>startlst[4] & tmptm<startlst[5]]); y22<-max(tmp[tmptm>startlst[4] & tmptm<startlst[5]])
    segments(c(startlst[1], startlst[1], startlst[2], startlst[2]), c(y11, y12, y12, y11), c(startlst[1], startlst[2], startlst[2], startlst[1]), c(y12, y12, y11, y11), lwd=2)
    segments(c(startlst[4], startlst[4], startlst[5], startlst[5]), c(y22, y21, y21, y22), c(startlst[4], startlst[5], startlst[5], startlst[4]), c(y21, y21, y22, y22), lwd=2)
    
    arrows(startlst[3], (y11+y21)/2, startlst[4], y21, lwd=2, length = 0.1, lend=2)
    arrows(startlst[3], (y11+y21)/2, startlst[2], y11, lwd=2, length = 0.1, lend=2)
    
    text(startlst[3], (y11+y21)/2, pos=1, labels = "time lag")
    
    text(mean(startlst[1:2]), y11, pos=1, labels = "training set")
    text(mean(startlst[4:5]), y21, pos=1, labels = "testing set")
    
    
    #plot CV
    plot(invarout$pdlag_list[[1]]$laglst, invarout$pdlag_list[[plotpos]]$CVest[,plotpos], xlab="", ylab="", type="l", lty=1, lwd=1.5, col=collst[2], xaxs="i", axes=FALSE, ...); abline(h=0, lty=3)
    axis(1); axis(2, las=2); box()
    
    mtext("time lag", 1, line=2.8, cex=1.1)
    mtext(expression(italic("CV")), 2, line=3.4, cex=1.1)
    put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)
    
    
    ##### Beta
    betalst<-rbind(cbind(rowMeans(betaout$beta_eig), rowMeans(betaout$beta_r0), rowMeans(betaout$beta_0)), rep(NA, 3))
    tlst<-0:(nrow(betaout$beta_eig))
    beta_est<-matrix(nrow=nrow(betalst), ncol=ncol(betalst))
    
    if(cleanupbeta>0) {
      for(i in 1:ncol(beta_est)) {
        beta_est[,i]<-exp(predict(loess(log(betalst[,i])~tlst, enp.target = cleanupbeta), newdata=data.frame(tlst=tlst)))
      }
    } else {
      beta_est<-betalst
    }
    
    hlst<-c(1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001)
    matplot(tlst, beta_est, type="n", lty=1, col=collst2, lwd=1.5, xlab="", ylab="", xaxs="i", log="y", axes=F, ...)
    abline(h=hlst, col=adjustcolor(1, alpha.f = 0.3), lty=1, lwd=1)
    matlines(tlst, beta_est, lty=1, col=collst2, lwd=1.5)
    
    mtext("time since event", 1, line=2.8, cex=1.1)
    mtext(expression(paste("community dissimilarity")), 2, line=3, cex=1.1)
    axis(1); axis(2, las=2); box()
    put.fig.letter("c.", "topleft", offset=ofs3, cex=1.6)
    
    return(invarlst)
  }
  
  if(FALSE) {
    pdf("figures/FIGURE_CV_examples.pdf", width=6, height=4, colormodel = "cmyk")
    tmp<-CVplotfun(outlong=out_meta_long, invarout=invar_meta, betaout=beta_meta, collst=collst, collst2=collst2, burnin=200, plotpos=1, doceq=2, cleanupbeta=20)
    dev.off()
    
    pdf("figures/SUP_FIGURE_CV_examples_allmodels.pdf", width=6, height=4, colormodel = "cmyk")
    tmp<-CVplotfun(outlong=out_meta_long, invarout=invar_meta, betaout=beta_meta, collst=collst, collst2=collst2, burnin=200, plotpos=1, doceq=2, cleanupbeta=20)
    tmp<-CVplotfun(outlong=out_neut_long, invarout=invar_neut, betaout=beta_neut, collst=collst, collst2=collst2, burnin=200, plotpos=1, doceq=0, cleanupbeta=20)
    tmp<-CVplotfun(outlong=out_dist_long, invarout=invar_dist, betaout=beta_dist, collst=collst, collst2=collst2, burnin=200, plotpos=1, doceq=1, cleanupbeta=20)
    tmp<-CVplotfun(outlong=out_psf_long, invarout=invar_psf, betaout=beta_psf, collst=collst, collst2=collst2, burnin=200, plotpos=1, doceq=0, cleanupbeta=20)
    tmp<-CVplotfun(outlong=out_rps_long, invarout=invar_rps, betaout=beta_rps, collst=collst, collst2=collst2, burnin=200, plotpos=1, doceq=0, cleanupbeta=20)
    dev.off()
  }
  
  ############################################################
  # Plot example of map
  ############################################################
  grid_sub<-grid_subset(gridout, size = 0.01)
  grid_sub2<-grid_subset(gridout, size = 0.5)
  
  ofs2<-c(0.3, -0.06)
  ofs4<-c(0.11, -0.045)
  
  pdf("figures/SUP_FIGURE_spatialsubset_map.pdf", width=7, height=7, colormodel = "cmyk")
  par(mar=c(2,2,2,2), oma=c(2,2,0,0), mfcol=c(2,2))
  plot_map(out_meta, gridout = gridout, grid_sub = grid_sub, collst=collst[-1])
  put.fig.letter("a.", "topleft", offset=ofs4, cex=1.6)
  
  segments(grid_sub2$borders[c(1,1,2,2)]+0.5*c(-1,-1,1,1), grid_sub2$borders[c(3,4,4,3)]+0.5*c(-1,1,1,-1), grid_sub2$borders[c(1,2,2,1)]+0.5*c(-1,1,1,-1), grid_sub2$borders[c(4,4,3,3)]+0.5*c(1,1,-1,-1), col="black", lwd=4)
  segments(grid_sub2$borders[c(1,1,2,2)]+0.5*c(-1,-1,1,1), grid_sub2$borders[c(3,4,4,3)]+0.5*c(-1,1,1,-1), grid_sub2$borders[c(1,2,2,1)]+0.5*c(-1,1,1,-1), grid_sub2$borders[c(4,4,3,3)]+0.5*c(1,1,-1,-1), col="white", lwd=1)
  
  shadowtext(51, 59, "1%", cex=1.3)
  shadowtext(51, 90, "50%", cex=1.3)
  
  mtext("x position", 1, line=2.3, cex=1.1)
  mtext("y position", 2, line=2.3, cex=1.1)
  
  
  #Plot spatial lag
  plot_map(out_meta, gridout = gridout, collst=collst[-1])
  put.fig.letter("b.", "topleft", offset=ofs4, cex=1.6)
  
  gseg<-c(20,40,40,60)
  segments(gseg[c(1,1,3,3)], gseg[c(2,4,4,2)], gseg[c(1,3,3,1)], gseg[c(4,4,2,2)], col="black", lwd=4)
  segments(gseg[c(1,1,3,3)], gseg[c(2,4,4,2)], gseg[c(1,3,3,1)], gseg[c(4,4,2,2)], col="white", lwd=1)
  
  gseg<-c(60,40,80,60)
  segments(gseg[c(1,1,3,3)], gseg[c(2,4,4,2)], gseg[c(1,3,3,1)], gseg[c(4,4,2,2)], col="black", lwd=4)
  segments(gseg[c(1,1,3,3)], gseg[c(2,4,4,2)], gseg[c(1,3,3,1)], gseg[c(4,4,2,2)], col="white", lwd=1)
  
  arrows(50, 40, 40, 40, lwd=4, angle = 30, length = 0.15)
  arrows(50, 40, 60, 40, lwd=4, angle = 30, length = 0.15)
  
  arrows(50, 40, 40, 40, lwd=1.5, angle = 30, length = 0.15, col="white")
  arrows(50, 40, 60, 40, lwd=1.5, angle = 30, length = 0.15, col="white")
  
  shadowtext(50, 32, "20 units", cex=1.3)
  
  mtext("x position", 1, line=2.3, cex=1.1)
  mtext("y position", 2, line=2.3, cex=1.1)
  
  #CV lag
  outlong<-out_meta_long; invarout<-invar_meta; betaout<-beta_meta; collst<-collst; collst2<-collst2; burnin=200; plotpos=1; doceq=2; cleanupbeta<-2
  
  invarlst<-outlong$output[-c(1:burnin),plotpos+1]/outlong$plotdata$ngrid
  matplot(outlong$output[-c(1:burnin),1], invarlst, type="l", lty=1, col=collst[-1], lwd=1.5, xlab="", ylab="", xaxs="i",
          ylim=range(invarlst)*c(0.98, 1), axes=FALSE)
  put.fig.letter("c.", "topleft", offset=ofs4, cex=1.6)
  axis(1); axis(2, las=2); box()
  if(sum(doceq)==2) {
    abline(h=c(sum(outlong$plotdata$ceq), outlong$plotdata$ceq), lty=3, col=collst, lwd=1)
  } else if(sum(doceq)==1) {
    abline(h=sum(outlong$plotdata$ceq), lty=3, col=collst, lwd=1)
  }
  mtext("simulation time", 1, line=2.3, cex=1.1)
  mtext("relative abundance", 2, line=2.8, cex=1.1)
  
  tmp<-outlong$output[,2]/outlong$plotdata$ngrid
  tmptm<-outlong$output[,1]
  
  startlst<-(nrow(outlong$output)-burnin)*c(0.2, 0.3, 0.5, 0.7, 0.8)+burnin
  y11<-min(tmp[tmptm>startlst[1] & tmptm<startlst[2]]); y12<-max(tmp[tmptm>startlst[1] & tmptm<startlst[2]])
  y21<-min(tmp[tmptm>startlst[4] & tmptm<startlst[5]]); y22<-max(tmp[tmptm>startlst[4] & tmptm<startlst[5]])
  
  segments(c(startlst[1], startlst[1], startlst[2], startlst[2]), c(y11, y12, y12, y11), c(startlst[1], startlst[2], startlst[2], startlst[1]), c(y12, y12, y11, y11), lwd=2)
  
  arrows((startlst[2]+startlst[1])/2, y11-0.002, startlst[1], y11-0.002, lwd=2, length = 0.1, lend=2)
  arrows((startlst[2]+startlst[1])/2, y11-0.002, startlst[2], y11-0.002, lwd=2, length = 0.1, lend=2)
  text((startlst[2]+startlst[1])/2, y11-0.002, paste("time span =", (startlst[2]-startlst[1])), pos=1)
  
  startlst<-(nrow(outlong$output)-burnin)*c(0.5, 0.9, 0.5, 0.7, 0.8)+burnin
  y11<-min(tmp[tmptm>startlst[1] & tmptm<startlst[2]]); y12<-max(tmp[tmptm>startlst[1] & tmptm<startlst[2]])
  y21<-min(tmp[tmptm>startlst[4] & tmptm<startlst[5]]); y22<-max(tmp[tmptm>startlst[4] & tmptm<startlst[5]])
  
  segments(c(startlst[1], startlst[1], startlst[2], startlst[2]), c(y11, y12, y12, y11), c(startlst[1], startlst[2], startlst[2], startlst[1]), c(y12, y12, y11, y11), lwd=2)
  arrows((startlst[2]+startlst[1])/2, y11-0.002, startlst[1], y11-0.002, lwd=2, length = 0.1, lend=2)
  arrows((startlst[2]+startlst[1])/2, y11-0.002, startlst[2], y11-0.002, lwd=2, length = 0.1, lend=2)
  text((startlst[2]+startlst[1])/2, y11-0.002, paste("time span =", (startlst[2]-startlst[1])), pos=1)
  
  
  
  
  
  
  startlst<-(nrow(outlong$output)-burnin)*c(0.2, 0.3, 0.5, 0.7, 0.8)+burnin
  y11<-min(tmp[tmptm>startlst[1] & tmptm<startlst[2]]); y12<-max(tmp[tmptm>startlst[1] & tmptm<startlst[2]])
  y21<-min(tmp[tmptm>startlst[4] & tmptm<startlst[5]]); y22<-max(tmp[tmptm>startlst[4] & tmptm<startlst[5]])
  
  invarlst<-outlong$output[-c(1:burnin),plotpos+1]/outlong$plotdata$ngrid
  matplot(outlong$output[-c(1:burnin),1], invarlst, type="l", lty=1, col=collst[-1], lwd=1.5, xlab="", ylab="", xaxs="i",
          ylim=range(invarlst)*c(0.98, 1), axes=FALSE)
  put.fig.letter("d.", "topleft", offset=ofs4, cex=1.6)
  axis(1); axis(2, las=2); box()
  if(sum(doceq)==2) {
    abline(h=c(sum(outlong$plotdata$ceq), outlong$plotdata$ceq), lty=3, col=collst, lwd=1)
  } else if(sum(doceq)==1) {
    abline(h=sum(outlong$plotdata$ceq), lty=3, col=collst, lwd=1)
  }
  mtext("simulation time", 1, line=2.3, cex=1.1)
  mtext("relative abundance", 2, line=2.8, cex=1.1)
  
  tmp<-outlong$output[,2]/outlong$plotdata$ngrid
  tmptm<-outlong$output[,1]
  
  startlst<-(nrow(outlong$output)-burnin)*c(0.2, 0.3, 0.5, 0.7, 0.8)+burnin
  
  segments(c(startlst[1], startlst[1], startlst[2], startlst[2]), c(y11, y12, y12, y11), c(startlst[1], startlst[2], startlst[2], startlst[1]), c(y12, y12, y11, y11), lwd=2)
  segments(c(startlst[4], startlst[4], startlst[5], startlst[5]), c(y22, y21, y21, y22), c(startlst[4], startlst[5], startlst[5], startlst[4]), c(y21, y21, y22, y22), lwd=2)
  
  arrows(startlst[3], (y11+y21)/2, startlst[4], y21, lwd=2, length = 0.1, lend=2)
  arrows(startlst[3], (y11+y21)/2, startlst[2], y11, lwd=2, length = 0.1, lend=2)
  
  text(startlst[3], (y11+y21)/2, pos=1, labels = "time lag")
  
  
  dev.off()
}  
