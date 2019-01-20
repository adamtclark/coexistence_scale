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
tsim<-50


############################################################
# Run global models
############################################################
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
r0_psf1<-estimate_rarereturn(out_psf, simtime=tsim*10, burnin=20, runtype="psf", doplot = FALSE)
r0_psf2<-estimate_rarereturn(out_psf, simtime=tsim*10, burnin=tsim*4, runtype="psf", doplot = FALSE)

############################################################
# Plot examples of each model
############################################################
#out<-out_meta; eigout<-eig_meta2; r0out<-r0_meta; collst<-collst; burnin=100; doceq=TRUE; plotpos=1

modplotfun<-function(out, eigout, r0out, collst, burnin=0, doceq=0, plotpos=1, atsq=0, totabund=F, doaxis1=T, figlet=1, ofs1=c(0,0), fcx=0, inbetweenfun=NULL, ceq_CUSTOM=NULL, ...) {
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
  #suppressWarnings(matplot(pabunds1[,1], pabunds1[,-1], type="l", lty=1, col=collst, lwd=c(1.5, rep(1.5, ncol(pabunds1)-2)), xlab="", ylab="", xaxs="i", axes=F, ylim=range(c(ptmp_eig_0[,-1], pabunds1[,-1]), na.rm=T), ...))
  #abline(v=c(mxt), lty=2); abline(h=c(0,1), lty=3, lwd=1)
  #if(sum(doceq)==2) {
  #  if(totabund) {
  #    abline(h=c(sum(out$plotdata$ceq), out$plotdata$ceq), lty=3, col=collst, lwd=1)
  #  } else {
  #    abline(h=c(out$plotdata$ceq), lty=3, col=collst, lwd=1)
  #  }
  #} else if(sum(doceq)==1) {
  #  if(totabund) {
  #    abline(h=sum(out$plotdata$ceq), lty=3, col=collst, lwd=1)
  #  }
  #} else if(sum(doceq)==3) {
  #  abline(h=c(ceq_CUSTOM), lty=3, col=collst, lwd=1)
  #}
  #if(doaxis1) {
  #  if(sum(atsq)==0) {
  #    axis(1, cex.axis=1.6)
  #  } else {
  #    axis(1, at=seq(min(pabunds1[,1], na.rm=T), max(pabunds1[,1]-0.1, na.rm=T), length=3), labels=atsq[1:3], cex.axis=1.6)
  #  }
  #}
  #axis(2, las=2, cex.axis=1.6); box()
  #
  ##add in perturbation null
  #matlines(ptmp_eig_0[,1], ptmp_eig_0[,-1], lty=2, col=collst, lwd=c(1, rep(1, ncol(pabunds1)-2)))
  
  #Add in disturbance lines
  #ap1<-tmp_eig_0[1,plotpos+1+totabund]/out$plotdata$ngrid
  #ap2<-tmp_eig[1,plotpos+1+totabund]/out$plotdata$ngrid
  #ap1<-max(c(ap1, ap2+diff(range(pabunds1[,-1],na.rm=T))*0.1))
  #arrows(mxt, ap1, mxt, ap2, lend=2, length = 0.06, col=1, lwd=1.5)
  
  #put.fig.letter(paste(letters[figlet], ".", sep=""), "topleft", offset=ofs1, cex=fcx)
  
  
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
  } else if(sum(doceq)==3) {
    abline(h=c(ceq_CUSTOM), lty=3, col=collst, lwd=1)
  }
  if(doaxis1) {
    if(sum(atsq)==0) {
      axis(1, cex.axis=1.6)
    } else {
      axis(1, at=seq(min(pabunds2[,1], na.rm=T), max(pabunds2[,1]-0.1, na.rm=T), length=4), labels=atsq[1:4], cex.axis=1.6)
    }
  }
  axis(2, las=2, cex.axis=1.6); box()
  
  put.fig.letter(paste(letters[figlet], ".", sep=""), "topleft", offset=ofs1, cex=fcx)
  
  return(list(pabunds1=pabunds1, pabunds2=pabunds2))
}


pdf("figures/SUP_FIGURE_negative_psf_alternatestates.pdf", width=10, height=6, colormodel = "cmyk")
  par(mfrow=c(2,2), mar=c(4,4,2,2), oma=c(2,2,0,0))
  atsq<-seq(0, 600, by=200)
  fcx<-2
  
  ofs1<-c(0.112, -0.04)
  
  tmp<-modplotfun(out=out_psf, eigout=eig_psf2, r0out=r0_psf2, collst=collst[-1], burnin=100, doceq=0, atsq=0, doaxis1 = T, figlet=1, ofs1=ofs1, fcx=fcx, plotpos = 1)
  tmp<-modplotfun(out=out_psf, eigout=eig_psf2, r0out=r0_psf2, collst=collst[-1], burnin=100, doceq=0, atsq=0, doaxis1 = T, figlet=2, ofs1=ofs1, fcx=fcx, plotpos = 2)
  
  tmp<-modplotfun(out=out_psf, eigout=eig_psf2, r0out=r0_psf1, collst=collst[-1], burnin=100, doceq=0, atsq=0, doaxis1 = T, figlet=3, ofs1=ofs1, fcx=fcx, plotpos = 1)
  tmp<-modplotfun(out=out_psf, eigout=eig_psf2, r0out=r0_psf1, collst=collst[-1], burnin=100, doceq=0, atsq=0, doaxis1 = T, figlet=4, ofs1=ofs1, fcx=fcx, plotpos = 2)
  
  mtext("species abundance", 2, line=0, outer=T, cex=2)
  mtext("simulation time", 1, line=0, outer=T, cex=2)
dev.off()

