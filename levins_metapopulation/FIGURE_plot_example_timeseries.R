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
ofs1<-c(0.028, -0.17)
ofs2<-c(0.075, -0.09)
ofs3<-c(0.075, -0.090)

#starting parameters
gridout<-makegrid(xlng = 100, ylng = 100)
xfac<-3
xfac_fast<-7
ptb<-0.2
collst<-adjustcolor(c("grey51", brewer.pal(4, "Set1")), alpha.f = 0.7)
collst2<-adjustcolor(c("blue", "red", "black"), alpha.f = 0.7)

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

r0_meta<-estimate_rarereturn(eig_meta2$out_lst0, simtime=tsim, burnin=tsim, runtype="metapopulation", doplot = FALSE)

if(FALSE) {
  set.seed(171205)
  out_meta_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_meta, talktime = 0)
  
  getEmeta<-getE(out_meta_long, Elst = 2:10)
  E_meta<-getEmeta$Eout
  E_meta_tot<-getEmeta$Eout_tot
  
  invar_meta<-estimate_invar(out = out_meta_long, E=E_meta, burnin=0, doplot=FALSE, Etot = E_meta_tot, niter = 0)
  
  beta_meta<-beta_estimate(out=out_meta, outlng = out_meta_long, Emat = E_meta, eigout = eig_meta2, r0out = r0_meta, burnin = 100)
}

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

r0_neut<-estimate_rarereturn(out = eig_neut2$out_lst0, simtime=tsim, burnin=tsim, runtype="neutral", doplot = FALSE)


if(FALSE) {
  set.seed(171206)
  out_neut_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_neut, talktime = 0, runtype = "neutral")
  getEneut<-getE(out_neut_long, Elst = 2:10)
  E_neut<-getEneut$Eout
  E_neut_tot<-getEneut$Eout_tot
  
  invar_neut<-estimate_invar(out_neut_long, E=E_neut, burnin=0, doplot=FALSE, Etot = E_neut_tot, niter = 0)
  
  beta_neut<-beta_estimate(out=out_neut, outlng = out_neut_long, Emat = E_neut, eigout = eig_neut2, r0out = r0_neut, burnin = 100)
}

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

out_dist_0<-rerunrun_metapopulation(out=eig_dist2$out_lst0, tmax=0, talktime = 0, runtype = "metapopulation", perturb = distlst, replace_perturb = 0)
r0_dist<-estimate_rarereturn(out_dist_0, simtime=tsim, burnin=tsim, runtype="disturbance", doplot = FALSE, prt = distlst,  prtfrq = prtfrq)

if(FALSE) {
  set.seed(171217)
  out_dist_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_dist, talktime = 0, runtype = "disturbance", prt = distlst,  prtfrq = prtfrq)
  getEdist<-getE(out_dist_long, Elst = 2:10)
  E_dist<-getEdist$Eout
  E_dist_tot<-getEdist$Eout_tot
  
  invar_dist<-estimate_invar(out_dist_long, E=E_dist, burnin=100, doplot=FALSE, Etot = E_dist_tot, niter = 0)
  
  out_dist_nodist<-run_metapopulation(tmax=200, gridout = gridout, population = population_dist, talktime = 0, runtype = "metapopulation")
  
  beta_dist<-beta_estimate(out=out_dist, outlng = out_dist_long, Emat = E_dist, eigout = eig_dist2, r0out = r0_dist, burnin = 100)
}

##### Try psf model
clst_psf= c(1.5, 1)*xfac
mlst_psf = rep(0.1, length(clst_psf))*xfac

population_psf<-populate(gridout, nlst = rep(floor(prod(gridout$lng)/length(clst_psf)*0.8), length(clst_psf)),
                         clst = clst_psf, radlst = Inf, mlst = mlst_psf)
set.seed(180108)
out_psf<-run_metapopulation(tmax=tinit, gridout = gridout, population = population_psf, talktime = 0, runtype = "psf")

eig_psf1<-estimate_eqreturn(out_psf, simtime=tsim, runtype="psf", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE)
eig_psf2<-estimate_eqreturn(out_psf, simtime=tsim, runtype="psf", replace_perturb = 0, talktime=0, prtb=ptb, doplot = FALSE)

r0_psf<-estimate_rarereturn(eig_psf2$out_lst0, simtime=tsim, burnin=tsim, runtype="psf", doplot = FALSE)

if(FALSE) {
  set.seed(180108)
  out_psf_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_psf, talktime = 0, runtype = "psf")
  getEpsf<-getE(out_psf_long, Elst = 2:10)
  E_psf<-getEpsf$Eout
  E_psf_tot<-getEpsf$Eout_tot
  
  invar_psf<-estimate_invar(out_psf_long, E=E_psf, burnin=100, doplot=FALSE, Etot = E_psf_tot, niter = 0)
  
  beta_psf<-beta_estimate(out=out_psf, outlng = out_psf_long, Emat = E_psf, eigout = eig_psf2, r0out = r0_psf, burnin = 10)
}

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

r0_rps<-estimate_rarereturn(out = eig_rps2$out_lst0, simtime=tsim, burnin=tsim, runtype="rps", doplot = FALSE)

if(FALSE) {
  set.seed(180104)
  out_rps_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_rps, talktime = 0, runtype = "rps", compmat = intmat_rps)
  getErps<-getE(out_rps_long, Elst = 2:10)
  E_rps<-getErps$Eout
  E_rps_tot<-getErps$Eout_tot
  
  invar_rps<-estimate_invar(out_rps_long, E=E_rps, burnin=0, doplot=FALSE, Etot = E_rps_tot, niter = 0)
  
  beta_rps<-beta_estimate(out=out_rps, outlng = out_rps_long, Emat = E_rps, eigout = eig_rps2, r0out = r0_rps, burnin = 10)
}

############################################################
# Plot examples of each model
############################################################
#out<-out_meta; eigout<-eig_meta2; r0out<-r0_meta; collst<-collst; burnin=100; doceq=TRUE; plotpos=1

modplotfun<-function(out, eigout, r0out, collst, burnin=0, doceq=0, plotpos=1, atsq=0, ...) {
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
  abunds<-rbind(abunds, tmp_eig, rep(NA, ncol(abunds)), tmp_eig_0[nrow(tmp_eig_0),], tmp_r0_0, tmp_r0)
  
  pabunds<-abunds/out$plotdata$ngrid
  pabunds[,1]<-abunds[,1]
  
  ptmp_eig_0<-tmp_eig_0/out$plotdata$ngrid
  ptmp_eig_0[,1]<-tmp_eig_0[,1]
  
  #plot
  matplot(pabunds[,1], pabunds[,-1], type="l", lty=1, col=collst, lwd=c(1.5, rep(1.5, ncol(pabunds)-2)), xlab="", ylab="", xaxs="i", axes=F, ...)
  abline(v=c(mxt, mxt+mxt_eig, mxt+mxt_eig+mxt_r0_0), lty=2); abline(h=c(0,1), lty=3, lwd=1)
  if(sum(doceq)==2) {
    abline(h=c(sum(out$plotdata$ceq), out$plotdata$ceq), lty=3, col=collst, lwd=1)
  } else if(sum(doceq)==1) {
    abline(h=sum(out$plotdata$ceq), lty=3, col=collst, lwd=1)
  }
  if(sum(atsq)==0) {
    axis(1)
  } else {
    axis(1, at=c(min(pabunds[,1], na.rm=T), atsq[-c(1, length(atsq))], max(pabunds[,1], na.rm=T)), labels=atsq, xpd=NA)
  }
  axis(2, las=2); box()
  
  #add in perturbation null
  matlines(ptmp_eig_0[,1], ptmp_eig_0[,-1], lty=2, col=collst, lwd=c(1, rep(1, ncol(pabunds)-2)))
  
  #Add in disturbance lines
  ap1<-tmp_eig_0[1,plotpos+2]/out$plotdata$ngrid
  ap2<-tmp_eig[1,plotpos+2]/out$plotdata$ngrid
  ap1<-max(c(ap1, ap2+diff(range(pabunds[,-1],na.rm=T))*0.1))
  arrows(mxt, ap1, mxt, ap2, lend=2, length = 0.06, col=1, lwd=1.5)
  
  ap1<-tmp_eig[nrow(tmp_eig),plotpos+2]/out$plotdata$ngrid
  ap2<-0
  ap1<-max(c(ap1, ap2+diff(range(pabunds[,-1],na.rm=T))*0.1))
  arrows(mxt+mxt_eig, ap1, mxt+mxt_eig, ap2, lend=2, length = 0.06, col=1, lwd=1.5)
  
  ap1<-tmp_r0[1,plotpos+2]/out$plotdata$ngrid
  ap2<-0
  ap1<-max(c(ap1, ap2+diff(range(pabunds[,-1],na.rm=T))*0.1))
  arrows(mxt+mxt_eig+mxt_r0_0, ap2, mxt+mxt_eig+mxt_r0_0, ap1, lend=2, length = 0.06, col=1, lwd=1.5)
  
  return(pabunds)
}


pdf("figures/FIGURE_model_examples.pdf", width=5, height=7, colormodel = "cmyk")
m<-as.matrix(1:5)
layout(m)
par(mar=c(1,1,2.5,1), oma=c(3,3.5,2.8,0))
atsq<-seq(0, 800, by=200)
fcx<-1.4

tmp<-modplotfun(out=out_meta, eigout=eig_meta2, r0out=r0_meta, collst=collst, burnin=100, doceq=2, atsq=atsq)
put.fig.letter("a.", "topleft", offset=ofs1, cex=fcx)

#label perturbations
mxt<-max(tmp[,-1], na.rm=T)+diff(range(tmp[,-1], na.rm=T))*0.08
text(200, mxt, "1. peturbation", xpd=NA, srt=45, adj = c(0,0), cex=1.2)
text(400, mxt, "2. removal", xpd=NA, srt=45, adj = c(0,0), cex=1.2)
text(600, mxt, "3. invasion", xpd=NA, srt=45, adj = c(0,0), cex=1.2)

tmp<-modplotfun(out=out_rps, eigout=eig_rps2, r0out=r0_rps, collst=collst, burnin=100, doceq=1, atsq=atsq)
put.fig.letter("b.", "topleft", offset=ofs1, cex=fcx)

tmp<-modplotfun(out=out_psf, eigout=eig_psf2, r0out=r0_psf, collst=collst, burnin=100, doceq=0, atsq=atsq)
put.fig.letter("c.", "topleft", offset=ofs1, cex=fcx)

tmp<-modplotfun(out=out_dist, eigout=eig_dist2, r0out=r0_dist, collst=collst[c(1,3,2)], burnin=100, doceq=0, plotpos = 2, atsq=atsq)
put.fig.letter("d.", "topleft", offset=ofs1, cex=fcx)

tmp<-modplotfun(out=out_neut, eigout=eig_neut2, r0out=r0_neut, collst=collst, burnin=100, doceq=1, atsq=atsq)
put.fig.letter("e.", "topleft", offset=ofs1, cex=fcx)

mtext("simulation time", 1, line=1.8, cex=1.1, outer = T)
mtext("species or community abundance", 2, line=2, cex=1.1, outer = T)
dev.off()




