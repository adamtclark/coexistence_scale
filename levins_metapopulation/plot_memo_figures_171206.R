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
ofs1<-c(0.075, -0.04)
ofs2<-c(0.075, -0.09)
ofs3<-c(0.075, -0.090)

#starting parameters
gridout<-makegrid(xlng = 100, ylng = 100)
xfac<-5
ptb<-0.1
collst<-c("black", brewer.pal(4, "Dark2"))
collst2<-adjustcolor(c("blue", "red", "black"), alpha.f = 0.8)

############################################################
# "Global" run
############################################################

##### Try Tilman metapopulation model
set.seed(171205)
#clst_meta = c(0.15, 0.3, 0.8, 3)*xfac
clst_meta = c(0.15, 0.35)*xfac
mlst_meta = rep(0.1, length(clst_meta))*xfac
#getceq(clst_meta, mlst_meta)

population_meta<-populate(gridout, nlst = floor(getceq(clst_meta, mlst_meta)*prod(gridout$lng)),
                          clst = clst_meta, radlst = Inf, mlst = mlst_meta)
out_meta<-run_metapopulation(tmax=200, gridout = gridout, population = population_meta, talktime = 0)
eig_meta1<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE)
eig_meta2<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation", replace_perturb = 0, talktime=0, prtb=ptb, doplot = FALSE)

r0_meta<-estimate_rarereturn(out_meta, simtime=100, burnin=100, runtype="metapopulation", doplot = FALSE)

out_meta_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_meta, talktime = 0)
getEmeta<-getE(out_meta_long, Elst = 2:10)
E_meta<-getEmeta$Eout

invar_meta<-estimate_invar(out_meta_long, E=E_meta, burnin=0, doplot=FALSE)

beta_meta<-beta_estimate(out=out_meta, outlng = out_meta_long, Emat = E_meta, eigout = eig_meta2, r0out = r0_meta, burnin = 10)

pdf("figures/memo_171206/levis_map.pdf", width=5, height=5, colormodel = "cmyk")
par(mar=c(2,2,2,2), oma=c(2,2,0,0))

plot_map(out_meta_long, gridout = gridout, collst=collst[-1])

mtext("x position", 1, line=2.3, cex=1.1)
mtext("y position", 2, line=2.3, cex=1.1)
dev.off()


##### Try Hubbell neutal model
set.seed(171206)

clst_neut<-rep(0.5, 2)*xfac
mlst_neut<-rep(0.1, length(clst_neut))*xfac
tmp<-abs(getceq(clst_neut, mlst_neut))
population_neut<-populate(gridout, nlst = round(rep(unique(tmp[tmp>0])/length(clst_neut), length(clst_neut))*prod(gridout$lng)),
                          clst = clst_neut, radlst = Inf, mlst = mlst_neut)
out_neut<-run_metapopulation(tmax=200, gridout = gridout, population = population_neut, talktime = 0, runtype = "neutral")

eig_neut1<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE)
eig_neut2<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral", replace_perturb = 0, talktime=0, prtb=ptb, doplot = FALSE)

r0_neut<-estimate_rarereturn(out = out_neut, simtime=100, burnin=100, runtype="neutral", doplot = FALSE)

set.seed(180104)
out_neut_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_neut, talktime = 0, runtype = "neutral")
getEneut<-getE(out_neut_long, Elst = 2:10)
E_neut<-getEneut$Eout
#E_neut[E_neut<4]<-4

invar_neut<-estimate_invar(out_neut_long, E=E_neut, burnin=0, doplot=FALSE)

beta_neut<-beta_estimate(out=out_neut, outlng = out_neut_long, Emat = E_neut, eigout = eig_neut2, r0out = r0_neut, burnin = 10)

##### Try disturbance model
set.seed(171217)

clst_dist= c(0.145, 0.16)*xfac
mlst_dist = rep(0.1, length(clst_dist))*xfac
distlst<-c(0.95, 0)

population_dist<-populate(gridout, nlst = rep(floor(prod(gridout$lng)/length(clst_dist)*0.8), length(clst_dist)),
                          clst = clst_dist, radlst = Inf, mlst = mlst_dist)

out_dist<-run_metapopulation(tmax=200, gridout = gridout, population = population_dist, talktime = 0, runtype = "disturbance", prt = distlst,  prtfrq = 20)
out_dist_0<-rerunrun_metapopulation(out=out_dist, tmax=0, talktime = 0, runtype = "metapopulation", perturb = distlst, replace_perturb = 0)

eig_dist1<-estimate_eqreturn(out_dist_0, simtime=100, runtype="disturbance", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE, prt = distlst,  prtfrq = 20)
eig_dist2<-estimate_eqreturn(out_dist_0, simtime=100, runtype="disturbance", replace_perturb = 0, talktime=0, prtb=ptb, doplot = FALSE, prt = distlst,  prtfrq = 20)

r0_dist<-estimate_rarereturn(out_dist_0, simtime=100, burnin=100, runtype="disturbance", doplot = FALSE, prt = distlst,  prtfrq = 20)

out_dist_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_dist, talktime = 0, runtype = "disturbance", prt = distlst,  prtfrq = 20)
getEdist<-getE(out_dist_long, Elst = 2:10)
E_dist<-getEdist$Eout

invar_dist<-estimate_invar(out_dist_long, E=E_dist, burnin=100, doplot=FALSE)


out_dist_nodist<-run_metapopulation(tmax=200, gridout = gridout, population = population_dist, talktime = 0, runtype = "metapopulation")

beta_dist<-beta_estimate(out=out_dist, outlng = out_dist_long, Emat = E_dist, eigout = eig_dist2, r0out = r0_dist, burnin = 10)

##### Try psf model
set.seed(180108)

clst_psf= c(1.5, 1)*xfac
mlst_psf = rep(0.1, length(clst_psf))*xfac

population_psf<-populate(gridout, nlst = rep(floor(prod(gridout$lng)/length(clst_psf)*0.8), length(clst_psf)),
                          clst = clst_psf, radlst = Inf, mlst = mlst_psf)

out_psf<-run_metapopulation(tmax=200, gridout = gridout, population = population_psf, talktime = 0, runtype = "psf")

eig_psf1<-estimate_eqreturn(out_psf, simtime=100, runtype="psf", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE)
eig_psf2<-estimate_eqreturn(out_psf, simtime=100, runtype="psf", replace_perturb = 0, talktime=0, prtb=ptb, doplot = FALSE)

r0_psf<-estimate_rarereturn(out_psf, simtime=100, burnin=100, runtype="psf", doplot = FALSE)

out_psf_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_psf, talktime = 0, runtype = "psf")
getEpsf<-getE(out_psf_long, Elst = 2:10)
E_psf<-getEpsf$Eout

invar_psf<-estimate_invar(out_psf_long, E=E_psf, burnin=100, doplot=FALSE)

beta_psf<-beta_estimate(out=out_psf, outlng = out_psf_long, Emat = E_psf, eigout = eig_psf2, r0out = r0_psf, burnin = 10)


##### Try rps model
set.seed(171209)

intmat_rps<-rbind(c(1,0,0,1),
                  c(1,1,0,0),
                  c(0,1,1,0),
                  c(0,0,1,1))

clst_rps<-rep(0.4, 4)
mlst_rps<-rep(0.1, length(clst_rps))
tmp<-abs(getceq(clst_rps, mlst_rps))
population_rps<-populate(gridout, nlst = round(rep(unique(tmp[tmp>0])/length(clst_rps), length(clst_rps))*prod(gridout$lng)),
                          clst = clst_rps, radlst = Inf, mlst = mlst_rps)
out_rps<-run_metapopulation(tmax=200, gridout = gridout, population = population_rps, talktime = 0, runtype = "rps", compmat = intmat_rps)

eig_rps1<-estimate_eqreturn(out_rps, simtime=100, runtype="rps", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE)
eig_rps2<-estimate_eqreturn(out_rps, simtime=100, runtype="rps", replace_perturb = 0, talktime=0, prtb=ptb, doplot = FALSE)

r0_rps<-estimate_rarereturn(out = out_rps, simtime=100, burnin=100, runtype="rps", doplot = FALSE)

set.seed(180104)
out_rps_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_rps, talktime = 0, runtype = "rps", compmat = intmat_rps)
getErps<-getE(out_rps_long, Elst = 2:10)
E_rps<-getErps$Eout

invar_rps<-estimate_invar(out_rps_long, E=E_rps, burnin=0, doplot=FALSE)

beta_rps<-beta_estimate(out=out_rps, outlng = out_rps_long, Emat = E_rps, eigout = eig_rps2, r0out = r0_rps, burnin = 10)






#Plot perturbation, Levins
pdf("figures/memo_171206/levins_perturbation.pdf", width=4, height=6, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2), rep(3,2))))
layout(m)

tmp<-eig_meta2$out_lst[[1]]$output
tmp[,1]<-tmp[,1]+200
pmat_meta<-rbind(out_meta$output, tmp)

matplot(pmat_meta[,1], pmat_meta[,-1]/out_meta$plotdata$ngrid, type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=200, lty=2)
ceq<-getceq(clst_meta, mlst_meta)
abline(h=c(sum(ceq), ceq), lty=3, col=collst, lwd=1.5)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
arrows(200, ceq[1]+0.03,
       200, ceq[1],
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_meta2$out_lst[[1]]$output), abs(eig_meta2$out_lst[[1]]$output[,2]-eig_meta2$out_lst0$output[,2])/out_meta$plotdata$ngrid, type="l", ylab="estimated distance", xlab="time", col=collst[2], lwd=2, xaxs="i", ylim=c(0, 0.035), cex.lab=1.5); abline(h=0, lty=3)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

plot(1:nrow(eig_meta2$eigenlst), eig_meta2$eigenlst[,1]*(1:nrow(eig_meta2$eigenlst)), type="l", lwd=2, col=collst[2], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-6, 4), cex.lab=1.5); abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)




#repeat for total
matplot(pmat_meta[,1], rowSums(pmat_meta[,-1])/out_meta$plotdata$ngrid, type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=200, lty=2)
ceq<-getceq(clst_meta, mlst_meta)
abline(h=c(sum(ceq), ceq), lty=3, col=collst, lwd=1.5)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
arrows(200, sum(ceq)+0.03,
       200, sum(ceq),
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_meta2$out_lst[[1]]$output), abs(rowSums(eig_meta2$out_lst[[1]]$output[,-1])-rowSums(eig_meta2$out_lst0$output[,-1]))/out_meta$plotdata$ngrid, type="l", ylab="estimated distance", xlab="time", col=collst[1], lwd=2, xaxs="i", ylim=c(0, 0.035), cex.lab=1.5); abline(h=0, lty=3)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

plot(1:nrow(eig_meta2$eigenlst_tot), eig_meta2$eigenlst_tot[,1]*(1:nrow(eig_meta2$eigenlst_tot)), type="l", lwd=2, col=collst[1], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-6, 0), cex.lab=1.5); abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)
dev.off()


#Plot perturbation, Hubbell
pdf("figures/memo_171206/neutral_perturbation.pdf", width=4, height=6, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2), rep(3,2))))
layout(m)

tmp<-eig_neut1$out_lst[[1]]$output
tmp[,1]<-tmp[,1]+200
pmat_neut1<-rbind(out_neut$output, tmp)

tmp<-eig_neut2$out_lst[[1]]$output
tmp[,1]<-tmp[,1]+200
pmat_neut2<-rbind(out_neut$output, tmp)

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2), rep(3,2))))
layout(m)

#with replacing perturbation
matplot(pmat_neut1[,1], pmat_neut1[,-1]/out_neut$plotdata$ngrid, type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=200, lty=2)
ceq<-getceq(clst_neut, mlst_neut)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
arrows(200, 0.27+0.02,
       200, 0.27,
       length = 0.08, lwd=2, lend=4)


plot(1:nrow(eig_neut1$out_lst[[1]]$output), abs(eig_neut1$out_lst[[1]]$output[,2]-eig_neut1$out_lst0$output[,2])/out_neut$plotdata$ngrid, type="l", ylab="estimated distance", xlab="time", col=collst[2], lwd=2, xaxs="i", ylim=c(0.01, 0.2)); abline(h=0, lty=3)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)

plot(1:nrow(eig_neut1$eigenlst), eig_neut1$eigenlst[,1]*(1:nrow(eig_neut1$eigenlst)), type="l", lwd=2, col=collst[2], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(0, 2)); abline(h=0, lty=3)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)




#repeat for total
matplot(pmat_neut1[,1], rowSums(pmat_neut1[,-1])/out_neut$plotdata$ngrid, type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=200, lty=2)
ceq<-getceq(clst_neut, mlst_neut)
abline(h=ceq, lty=2)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
arrows(200, 0.8+0.02,
       200, 0.8,
       length = 0.08, lwd=2, lend=4)


plot(1:nrow(eig_neut1$out_lst[[1]]$output), abs(rowSums(eig_neut1$out_lst[[1]]$output[,-1])-rowSums(eig_neut1$out_lst0$output[,-1]))/out_neut$plotdata$ngrid, type="l", ylab="estimated distance", xlab="time", col=collst[1], lwd=2, xaxs="i", ylim=c(0, 0.02)); abline(h=0, lty=3)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)

plot(1:nrow(eig_neut1$eigenlst), eig_neut1$eigenlst_tot[,1]*(1:nrow(eig_neut1$eigenlst_tot)), type="l", lwd=2, col=collst[1], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-4, 2)); abline(h=0, lty=3)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)

dev.off()





pdf("figures/memo_171206/neutral_perturbation_noreplace.pdf", width=4, height=6, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2), rep(3,2))))
layout(m)

#without replacing perturbation
matplot(pmat_neut2[,1], pmat_neut2[,-1]/out_neut$plotdata$ngrid, type="l", lty=1, col=collst[c(-1)], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=200, lty=2)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
arrows(200, 0.27+0.02,
       200, 0.27,
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_neut2$out_lst[[1]]$output), abs(eig_neut2$out_lst[[1]]$output[,2]-eig_neut2$out_lst0$output[,2])/out_neut$plotdata$ngrid, type="l", ylab="estimated distance", xlab="time", col=collst[2], lwd=2, xaxs="i", ylim=c(0, 0.15)); abline(h=0, lty=3)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)

plot((1:nrow(eig_neut2$eigenlst))[is.finite(eig_neut2$eigenlst[,1])], eig_neut2$eigenlst[is.finite(eig_neut2$eigenlst[,1]),1]*(1:nrow(eig_neut2$eigenlst))[is.finite(eig_neut2$eigenlst[,1])], type="l", lwd=2, col=collst[2], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-2, 2)); abline(h=0, lty=3)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)

dev.off()


#Hubbell: population-level stability
pdf("figures/memo_171206/neutral_totpop_stability.pdf", width=4, height=3, colormodel = "cmyk")
par(mar=c(2,2,2,2), oma=c(2,2,0,0))

matplot(pmat_neut2[,1], rowSums(pmat_neut2[,-1])/out_neut$plotdata$ngrid, type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=200, lty=2)
abline(h=c(abs(unique(out_neut$plotdata$ceq))), lty=3, col=1, lwd=1.5)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)

dev.off()



#Plot perturbation, dist
pdf("figures/memo_171206/dist_perturbation.pdf", width=4, height=6, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2), rep(3,2))))
layout(m)

tmp<-eig_dist2$out_lst[[1]]$output
tmp[,1]<-tmp[,1]+max(out_dist$output[,1])
pmat_dist<-rbind(out_dist$output, tmp)

matplot(pmat_dist[,1], pmat_dist[,-1]/out_dist$plotdata$ngrid, type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=max(out_dist$output[,1]), lty=2)
ceq<-getceq(clst_dist, mlst_dist)
abline(h=c(sum(ceq), ceq), lty=3, col=collst, lwd=1.5)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
arrows(max(out_dist$output[,1]), 0.22+0.03,
       max(out_dist$output[,1]), 0.22,
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_dist2$out_lst[[1]]$output), abs(eig_dist2$out_lst[[1]]$output[,2]-eig_dist2$out_lst0$output[,2])/out_dist$plotdata$ngrid, type="l", ylab="estimated distance", xlab="time", col=collst[2], lwd=2, xaxs="i", ylim=c(0, 0.08), cex.lab=1.5); abline(h=0, lty=3)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

plot(1:nrow(eig_dist2$eigenlst), eig_dist2$eigenlst[,1]*(1:nrow(eig_dist2$eigenlst)), type="l", lwd=2, col=collst[2], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-3.5, 4), cex.lab=1.5, xlim=c(1, 80)); abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)



#repeat with total
matplot(pmat_dist[,1], rowSums(pmat_dist[,-1])/out_dist$plotdata$ngrid, type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=max(out_dist$output[,1]), lty=2)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
arrows(max(out_dist$output[,1]), 0.4+0.03,
       max(out_dist$output[,1]), 0.4,
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_dist2$out_lst[[1]]$output), abs(rowSums(eig_dist2$out_lst[[1]]$output[,-1])-rowSums(eig_dist2$out_lst0$output[,-1]))/out_dist$plotdata$ngrid, type="l", ylab="estimated distance", xlab="time", col=collst[1], lwd=2, xaxs="i", ylim=c(0, 0.1), cex.lab=1.5); abline(h=0, lty=3)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

plot(1:nrow(eig_dist2$eigenlst), eig_dist2$eigenlst_tot[,1]*(1:nrow(eig_dist2$eigenlst_tot)), type="l", lwd=2, col=collst[1], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-3.5, 4.5), cex.lab=1.5, xlim=c(1, 80)); abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)

dev.off()


#Dist: plot without disturbances
pdf("figures/memo_171206/dist_nodist_stability.pdf", width=4, height=3, colormodel = "cmyk")
par(mar=c(2,2,2,2), oma=c(2,2,0,0))

matplot(out_dist_nodist$output[,1], out_dist_nodist$output[,-1]/out_dist_nodist$plotdata$ngrid, type="l", lty=1, col=collst[2:5], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(h=0, lty=3, col=1, lwd=1.5)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)

dev.off()


#Plot perturbation, psf
pdf("figures/memo_171206/psf_perturbation.pdf", width=4, height=6, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2), rep(3,2))))
layout(m)

tmp<-eig_psf2$out_lst[[1]]$output
tmp[,1]<-tmp[,1]+max(out_psf$output[,1])
pmat_psf<-rbind(out_psf$output, tmp)

matplot(pmat_psf[,1], pmat_psf[,-1]/out_psf$plotdata$ngrid, type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=max(out_psf$output[,1]), lty=2)

mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
arrows(max(out_psf$output[,1]), 0.48+0.03,
       max(out_psf$output[,1]), 0.48,
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_psf2$out_lst[[1]]$output), abs(eig_psf2$out_lst[[1]]$output[,2]-eig_psf2$out_lst0$output[,2])/out_psf$plotdata$ngrid, type="l", ylab="estimated distance", xlab="time", col=collst[2], lwd=2, xaxs="i", ylim=c(0, 0.05), cex.lab=1.5); abline(h=0, lty=3)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

plot(1:nrow(eig_psf2$eigenlst), eig_psf2$eigenlst[,1]*(1:nrow(eig_psf2$eigenlst)), type="l", lwd=2, col=collst[2], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-5, 1), cex.lab=1.5, xlim=c(1, 80)); abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)



#Total
matplot(pmat_psf[,1], rowSums(pmat_psf[,-1])/out_psf$plotdata$ngrid, type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=max(out_psf$output[,1]), lty=2)

mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
arrows(max(out_psf$output[,1]), 0.48+0.03,
       max(out_psf$output[,1]), 0.48,
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_psf2$out_lst[[1]]$output), abs(rowSums(eig_psf2$out_lst[[1]]$output[,-1])-rowSums(eig_psf2$out_lst0$output[,-1]))/out_psf$plotdata$ngrid, type="l", ylab="estimated distance", xlab="time", col=collst[1], lwd=2, xaxs="i", ylim=c(0, 0.05), cex.lab=1.5); abline(h=0, lty=3)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

plot(1:nrow(eig_psf2$eigenlst), eig_psf2$eigenlst_tot[,1]*(1:nrow(eig_psf2$eigenlst_tot)), type="l", lwd=2, col=collst[1], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-8, 0), cex.lab=1.5, xlim=c(1, 80)); abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)
dev.off()




#Plot perturbation, rps
pdf("figures/memo_171206/rps_perturbation.pdf", width=4, height=6, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2), rep(3,2))))
layout(m)

tmp<-eig_rps2$out_lst[[1]]$output
tmp[,1]<-tmp[,1]+max(out_rps$output[,1])
pmat_rps<-rbind(out_rps$output, tmp)

matplot(pmat_rps[,1], pmat_rps[,-1]/out_rps$plotdata$ngrid, type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i", ylim=c(0.1, 0.33))
abline(v=max(out_rps$output[,1]), lty=2)

mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
arrows(max(out_rps$output[,1]), 0.30+0.03,
       max(out_rps$output[,1]), 0.30,
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_rps2$out_lst[[1]]$output), abs(eig_rps2$out_lst[[1]]$output[,2]-eig_rps2$out_lst0$output[,2])/out_rps$plotdata$ngrid, type="l", ylab="estimated distance", xlab="time", col=collst[2], lwd=2, xaxs="i", ylim=c(0, 0.1), cex.lab=1.5); abline(h=0, lty=3)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

plot(1:nrow(eig_rps2$eigenlst), eig_rps2$eigenlst[,1]*(1:nrow(eig_rps2$eigenlst)), type="l", lwd=2, col=collst[2], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-5, 1), cex.lab=1.5, xlim=c(1, 80)); abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)



#Total
matplot(pmat_rps[,1], rowSums(pmat_rps[,-1])/out_rps$plotdata$ngrid, type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=max(out_rps$output[,1]), lty=2)

mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
arrows(max(out_rps$output[,1]), 0.75+0.03,
       max(out_rps$output[,1]), 0.75,
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_rps2$out_lst[[1]]$output), abs(rowSums(eig_rps2$out_lst[[1]]$output[,-1])-rowSums(eig_rps2$out_lst0$output[,-1]))/out_rps$plotdata$ngrid, type="l", ylab="estimated distance", xlab="time", col=collst[1], lwd=2, xaxs="i", ylim=c(0, 0.03), cex.lab=1.5); abline(h=0, lty=3)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

plot(1:nrow(eig_rps2$eigenlst), eig_rps2$eigenlst[,1]*(1:nrow(eig_rps2$eigenlst)), type="l", lwd=2, col=collst[1], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-5.5, 1.5), cex.lab=1.5, xlim=c(1, 80)); abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)
dev.off()














#Plot rare return, Levins
pdf("figures/memo_171206/levis_returnrare.pdf", width=4, height=4.5, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2))))
layout(m)

tmp<-r0_meta$out0_lst[[1]]$output
tmp[,1]<-tmp[,1]+200
pmat_meta<-rbind(out_meta$output, tmp)
tmp<-r0_meta$out_lst[[1]]$output
tmp[,1]<-tmp[,1]+300
pmat_meta<-rbind(pmat_meta, tmp)

matplot(pmat_meta[,1], pmat_meta[,-1]/out_meta$plotdata$ngrid, type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=c(200, 300), lty=2)
ceq<-getceq(clst_meta, mlst_meta)
abline(h=c(sum(ceq), ceq), lty=3, col=collst, lwd=1.5)
abline(h=0, lty=3)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
arrows(200, ceq[1]+0.06,
       200, ceq[1],
       length = 0.08, lwd=2, lend=4)
arrows(300, 0,
       300, 0.06,
       length = 0.08, lwd=2, lend=4)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)

plot(1:nrow(r0_meta$grwrare), r0_meta$grwrare[,1]*(1:nrow(r0_meta$grwrare)), col=collst[2], type="l", lty=1, lwd=1.5, xlab="time span", ylab=expression(italic(paste(r[0], "t"))), xaxs="i")
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(r[0], "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)




#Total
matplot(pmat_meta[,1], rowSums(pmat_meta[,-1])/out_meta$plotdata$ngrid, type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=c(200, 300), lty=2)
ceq<-getceq(clst_meta, mlst_meta)
abline(h=c(sum(ceq)), lty=3, col=collst, lwd=1.5)
abline(h=0, lty=3)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
arrows(200, sum(ceq)+0.06,
       200, sum(ceq),
       length = 0.08, lwd=2, lend=4)
arrows(300, 0.6,
       300, 0.66,
       length = 0.08, lwd=2, lend=4)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)

plot(1:nrow(r0_meta$grwrare), r0_meta$grwrare_tot[,1]*(1:nrow(r0_meta$grwrare_tot)), col=collst[1], type="l", lty=1, lwd=1.5, xlab="time span", ylab=expression(italic(paste(r[0], "t"))), xaxs="i")
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(r[0], "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

dev.off()


#Plot rare return, Hubbell
pdf("figures/memo_171206/neutral_returnrare.pdf", width=4, height=4.5, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2))))
layout(m)

tmp<-r0_neut$out0_lst[[1]]$output
tmp[,1]<-tmp[,1]+200
pmat_neut<-rbind(out_neut$output, tmp)
tmp<-r0_neut$out_lst[[1]]$output
tmp[,1]<-tmp[,1]+300
pmat_neut<-rbind(pmat_neut, tmp)

matplot(pmat_neut[,1], pmat_neut[,-1]/out_neut$plotdata$ngrid, type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=c(200, 300), lty=2)
ceq<-getceq(clst_meta, mlst_meta)
abline(h=c(sum(ceq), ceq), lty=3, col=collst, lwd=1.5)
abline(h=0, lty=3)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
arrows(200, ceq[1]+0.035,
       200, ceq[1],
       length = 0.08, lwd=2, lend=4)
arrows(300, 0,
       300, 0.035,
       length = 0.08, lwd=2, lend=4)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)

plot(1:nrow(r0_neut$grwrare), r0_neut$grwrare[,1]*(1:nrow(r0_neut$grwrare)), col=collst[2], type="l", lty=1, lwd=1.5, xlab="time span", ylab=expression(italic(paste(r[0], "t"))), xaxs="i")
abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(r[0], "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

dev.off()


#Plot rare return, dist
pdf("figures/memo_171206/dist_returnrare.pdf", width=4, height=4.5, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2))))
layout(m)

listpos<-2

mxt<-max(out_dist$output[,1])
mxt2<-max(r0_dist$out0_lst[[listpos]]$output[,1])

tmp<-r0_dist$out0_lst[[listpos]]$output
tmp[,1]<-tmp[,1]+mxt
pmat_dist<-rbind(out_dist$output, tmp)
tmp<-r0_dist$out_lst[[listpos]]$output
tmp[,1]<-tmp[,1]+mxt+mxt2
pmat_dist<-rbind(pmat_dist, tmp)

matplot(pmat_dist[,1], pmat_dist[,-1]/out_dist$plotdata$ngrid, type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=c(mxt, mxt+mxt2), lty=2)
ceq<-getceq(clst_dist, mlst_dist)
abline(h=c(sum(ceq), ceq), lty=3, col=collst, lwd=1.5)
abline(h=0, lty=3)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
arrows(mxt, 0.12+0.06,
       mxt, 0.12,
       length = 0.08, lwd=2, lend=4)
arrows(mxt+mxt2, 0,
       mxt+mxt2, 0.06,
       length = 0.08, lwd=2, lend=4)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)

plot(1:nrow(r0_dist$grwrare), r0_dist$grwrare[,listpos]*(1:nrow(r0_dist$grwrare)), col=collst[listpos+1], type="l", lty=1, lwd=1.5, xlab="time span", ylab=expression(italic(paste(r[0], "t"))), xaxs="i", xlim=c(1, 80))
abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(r[0], "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

dev.off()

#Plot rare return, psf
pdf("figures/memo_171206/psf_returnrare.pdf", width=4, height=4.5, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2))))
layout(m)

tmp<-r0_psf$out0_lst[[1]]$output
tmp[,1]<-tmp[,1]+200
pmat_psf<-rbind(out_psf$output, tmp)
tmp<-r0_psf$out_lst[[1]]$output
tmp[,1]<-tmp[,1]+300
pmat_psf<-rbind(pmat_psf, tmp)

matplot(pmat_psf[,1], pmat_psf[,-1]/out_psf$plotdata$ngrid, type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=c(200, 300), lty=2)
abline(h=0, lty=3)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
arrows(200, 0.5+0.06,
       200, 0.5,
       length = 0.08, lwd=2, lend=4)
arrows(300, 0,
       300, 0.06,
       length = 0.08, lwd=2, lend=4)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)

plot(1:nrow(r0_psf$grwrare), r0_psf$grwrare[,1]*(1:nrow(r0_psf$grwrare)), col=collst[2], type="l", lty=1, lwd=1.5, xlab="time span", ylab=expression(italic(paste(r[0], "t"))), xaxs="i")
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(r[0], "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

dev.off()

#Plot rare return, rps
pdf("figures/memo_171206/rps_returnrare.pdf", width=4, height=4.5, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2))))
layout(m)

tmp<-r0_rps$out0_lst[[1]]$output
tmp[,1]<-tmp[,1]+200
pmat_rps<-rbind(out_rps$output, tmp)
tmp<-r0_rps$out_lst[[1]]$output
tmp[,1]<-tmp[,1]+300
pmat_rps<-rbind(pmat_rps, tmp)

matplot(pmat_rps[,1], pmat_rps[,-1]/out_rps$plotdata$ngrid, type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=c(200, 300), lty=2)
abline(h=0, lty=3)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
arrows(200, 0.3+0.06,
       200, 0.3,
       length = 0.08, lwd=2, lend=4)
arrows(300, 0,
       300, 0.06,
       length = 0.08, lwd=2, lend=4)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)

plot(1:nrow(r0_rps$grwrare), r0_rps$grwrare[,1]*(1:nrow(r0_rps$grwrare)), col=collst[2], type="l", lty=1, lwd=1.5, xlab="time span", ylab=expression(italic(paste(r[0], "t"))), xaxs="i", xlim=c(0, 65))
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(r[0], "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)
abline(h=0, lty=3)

dev.off()









#Plot invar, Levins
pdf("figures/memo_171206/levins_invar.pdf", width=5, height=4, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 2), rep(2,2))))
layout(m)

matplot(out_meta_long$output[,1], out_meta_long$output[,2]/out_meta_long$plotdata$ngrid, type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i", ylim=c(0.3, 0.36))
abline(h=getceq(clst_meta, mlst_meta), col=collst[-1], lty=2, lwd=1.5)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)

segments(c(200, 200, 300, 300), c(0.31, 0.36, 0.36, 0.31), c(200, 300, 300, 200), c(0.36, 0.36, 0.31, 0.31), lwd=2)
segments(c(200, 200, 300, 300)+500, c(0.31, 0.36, 0.36, 0.31), c(200, 300, 300, 200)+500, c(0.36, 0.36, 0.31, 0.31), lwd=2)

arrows(500, 0.31, 690, 0.31, lwd=2, length = 0.1, lend=2)
arrows(500, 0.31, 310, 0.31, lwd=2, length = 0.1, lend=2)

text(500, 0.31, pos=1, labels = "time lag")

text(250, 0.31, pos=1, labels = "training set")
text(750, 0.31, pos=1, labels = "testing set")


plot(invar_meta$pdlag_list[[1]]$laglst, invar_meta$pdlag_list[[1]]$CVest[,1], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[2], xaxs="i", ylim=c(0, 0.05)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)




#Total
matplot(out_meta_long$output[,1], rowSums(out_meta_long$output[,-1])/out_meta_long$plotdata$ngrid, type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(h=sum(getceq(clst_meta, mlst_meta)), col=collst[1], lty=2, lwd=1.5)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)

plot(invar_meta$pdlag_list_tot[[1]]$laglst, invar_meta$pdlag_list_tot[[1]]$CVest[,1], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[1], xaxs="i", ylim=c(0, 0.02)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)

dev.off()




#Plot invar, Hubbell
pdf("figures/memo_171206/neutral_invar.pdf", width=5, height=4, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 2), rep(2,2))))
layout(m)

matplot(out_neut_long$output[,1], out_neut_long$output[,2]/out_neut_long$plotdata$ngrid, type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)

tmp<-out_neut_long$output[,2]/out_neut_long$plotdata$ngrid
tmptm<-out_neut_long$output[,1]

y11<-min(tmp[tmptm>200 & tmptm<300]); y12<-max(tmp[tmptm>200 & tmptm<300])
y21<-min(tmp[tmptm>700 & tmptm<800]); y22<-max(tmp[tmptm>700 & tmptm<800])
segments(c(200, 200, 300, 300), c(y11, y12, y12, y11), c(200, 300, 300, 200), c(y12, y12, y11, y11), lwd=2)
segments(c(200, 200, 300, 300)+500, c(y22, y21, y21, y22), c(200, 300, 300, 200)+500, c(y21, y21, y22, y22), lwd=2)

arrows(500, (y11+y21)/2, 690, y21, lwd=2, length = 0.1, lend=2)
arrows(500, (y11+y21)/2, 310, y11, lwd=2, length = 0.1, lend=2)

text(500, (y11+y21)/2, pos=1, labels = "time lag")

text(250, y11, pos=1, labels = "training set")
text(750, y21, pos=1, labels = "testing set")


plot(invar_neut$pdlag_list[[1]]$laglst, invar_neut$pdlag_list[[1]]$CVest[,1], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[2], xaxs="i", ylim=c(0, 0.3)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)




#total
matplot(out_neut_long$output[,1], rowSums(out_neut_long$output[,-1])/out_neut_long$plotdata$ngrid, type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)
abline(h=getceq(clst_neut, mlst_neut), lty=2)

plot(invar_neut$pdlag_list_tot[[1]]$laglst, invar_neut$pdlag_list_tot[[1]]$CVest[,1], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[1], xaxs="i", ylim=c(0, 0.01)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)

dev.off()



#Plot invar, dist
pdf("figures/memo_171206/dist_invar.pdf", width=5, height=4, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 2), rep(2,2))))
layout(m)

matplot(out_dist_long$output[,1], out_dist_long$output[,3]/out_dist_long$plotdata$ngrid, type="l", lty=1, col=collst[2], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i", ylim=c(-0.05, 0.31))
abline(h=getceq(clst_dist, mlst_dist)[1], col=collst[-1], lty=2, lwd=1.5)
abline(h=0, lty=3)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)
abline(v=100, lty=3)

y11<-0.06; y12<-0.27
y21<-0.09; y22<-0.30
segments(c(200, 200, 300, 300), c(y11, y12, y12, y11), c(200, 300, 300, 200), c(y12, y12, y11, y11), lwd=2)
segments(c(200, 200, 300, 300)+500, c(y21, y22, y22, y21), c(200, 300, 300, 200)+500, c(y22, y22, y21, y21), lwd=2)

arrows(500, (y11+y21)/2, 690, y21, lwd=2, length = 0.1, lend=2)
arrows(500, (y11+y21)/2, 310, y11, lwd=2, length = 0.1, lend=2)

text(500, (y11+y21)/2, pos=1, labels = "time lag")

text(250, y11, pos=1, labels = "training set")
text(750, y21, pos=1, labels = "testing set")


plot(invar_dist$pdlag_list[[2]]$laglst, invar_dist$pdlag_list[[2]]$CVest[,1], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[-1], xaxs="i", ylim=c(0, 0.3)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)



#total
matplot(out_dist_long$output[,1], rowSums(out_dist_long$output[,-1])/out_dist_long$plotdata$ngrid, type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(h=0, lty=3)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)
abline(v=100, lty=3)

plot(invar_dist$pdlag_list_tot[[2]]$laglst, invar_dist$pdlag_list_tot[[2]]$CVest[,1], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[1], xaxs="i", ylim=c(0, 0.3)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)

dev.off()


#Plot invar, psf
pdf("figures/memo_171206/psf_invar.pdf", width=5, height=4, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 2), rep(2,2))))
layout(m)

matplot(out_psf_long$output[,1], out_psf_long$output[,3]/out_psf_long$plotdata$ngrid, type="l", lty=1, col=collst[3], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)
abline(v=100, lty=3)

tmp<-out_psf_long$output[,3]/out_psf_long$plotdata$ngrid
tmptm<-out_psf_long$output[,1]

y11<-min(tmp[tmptm>200 & tmptm<300]); y12<-max(tmp[tmptm>200 & tmptm<300])
y21<-min(tmp[tmptm>700 & tmptm<800]); y22<-max(tmp[tmptm>700 & tmptm<800])
segments(c(200, 200, 300, 300), c(y11, y12, y12, y11), c(200, 300, 300, 200), c(y12, y12, y11, y11), lwd=2)
segments(c(200, 200, 300, 300)+500, c(y22, y21, y21, y22), c(200, 300, 300, 200)+500, c(y21, y21, y22, y22), lwd=2)

arrows(500, (y11+y21)/2, 690, y21, lwd=2, length = 0.1, lend=2)
arrows(500, (y11+y21)/2, 310, y11, lwd=2, length = 0.1, lend=2)

text(500, (y11+y21)/2, pos=1, labels = "time lag")

text(250, y12, pos=3, labels = "training set")
text(750, y22, pos=3, labels = "testing set")


plot(invar_psf$pdlag_list[[1]]$laglst, invar_psf$pdlag_list[[1]]$CVest[,2], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[3], xaxs="i", ylim=c(0, 0.03)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)


#Total
matplot(out_psf_long$output[,1], rowSums(out_psf_long$output[,-1])/out_psf_long$plotdata$ngrid, type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)
abline(v=100, lty=3)

plot(invar_psf$pdlag_list_tot[[1]]$laglst, invar_psf$pdlag_list_tot[[1]]$CVest[,2], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[1], xaxs="i", ylim=c(0, 0.015)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)

dev.off()


#Plot invar, rps
pdf("figures/memo_171206/rps_invar.pdf", width=5, height=4, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 2), rep(2,2))))
layout(m)

matplot(out_rps_long$output[,1], out_rps_long$output[,3]/out_rps_long$plotdata$ngrid, type="l", lty=1, col=collst[3], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i", ylim=c(0.03, 0.4))
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)
abline(v=100, lty=3)

tmp<-out_rps_long$output[,3]/out_rps_long$plotdata$ngrid
tmptm<-out_rps_long$output[,1]

y11<-min(tmp[tmptm>200 & tmptm<300]); y12<-max(tmp[tmptm>200 & tmptm<300])
y21<-min(tmp[tmptm>700 & tmptm<800]); y22<-max(tmp[tmptm>700 & tmptm<800])
segments(c(200, 200, 300, 300), c(y11, y12, y12, y11), c(200, 300, 300, 200), c(y12, y12, y11, y11), lwd=2)
segments(c(200, 200, 300, 300)+500, c(y22, y21, y21, y22), c(200, 300, 300, 200)+500, c(y21, y21, y22, y22), lwd=2)

arrows(500, (y11+y21)/2, 690, y21, lwd=2, length = 0.1, lend=2)
arrows(500, (y11+y21)/2, 310, y11, lwd=2, length = 0.1, lend=2)

text(500, (y11+y21)/2, pos=1, labels = "time lag")

text(250, y12, pos=3, labels = "training set")
text(750, y21, pos=1, labels = "testing set")


plot(invar_rps$pdlag_list[[1]]$laglst, invar_rps$pdlag_list[[1]]$CVest[,2], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[3], xaxs="i", ylim=c(0, 0.4)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)



#Tot
matplot(out_rps_long$output[,1], rowSums(out_rps_long$output[,-1])/out_rps_long$plotdata$ngrid, type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)
abline(v=100, lty=3)

plot(invar_rps$pdlag_list_tot[[1]]$laglst, invar_rps$pdlag_list_tot[[1]]$CVest[,2], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[1], xaxs="i", ylim=c(0, 0.01)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)

dev.off()







#Plot beta, Levins
hlst<-c(1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001)

pdf("figures/memo_171206/beta_levins.pdf", width=5, height=4, colormodel = "cmyk")
par(mar=c(2,2,2,2), oma=c(2,4,0,0))
m<-1
layout(m)
tmp<-beta_meta

matplot(1:nrow(tmp$beta_eig), cbind(rowMeans(tmp$beta_eig), rowMeans(tmp$beta_r0), rowMeans(tmp$beta_0)), type="l", lty=1, col=collst2, lwd=1.5, xlab="", ylab="", xaxs="i", log="y", axes=F)
abline(h=hlst, col=1, lty=3, lwd=1)
mtext("time", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=3.8, cex=1.1)

axis(1); axis(2, las=2); box()
dev.off()

#Plot beta, Neutral
pdf("figures/memo_171206/beta_neut.pdf", width=5, height=4, colormodel = "cmyk")
par(mar=c(2,2,2,2), oma=c(2,4,0,0))
m<-1
layout(m)
tmp<-beta_neut

matplot(1:nrow(tmp$beta_eig), cbind(rowMeans(tmp$beta_eig), rowMeans(tmp$beta_r0), rowMeans(tmp$beta_0)), type="l", lty=1, col=collst2, lwd=1.5, xlab="", ylab="", xaxs="i", log="y", axes=F)
abline(h=hlst, col=1, lty=3, lwd=1)
mtext("time", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=3.8, cex=1.1)

axis(1); axis(2, las=2); box()
dev.off()

#Plot beta, Dist
pdf("figures/memo_171206/beta_dist.pdf", width=5, height=4, colormodel = "cmyk")
par(mar=c(2,2,2,2), oma=c(2,4,0,0))
m<-1
layout(m)
tmp<-beta_dist

matplot(1:nrow(tmp$beta_eig), cbind(rowMeans(tmp$beta_eig), rowMeans(tmp$beta_r0), rowMeans(tmp$beta_0)), type="l", lty=1, col=collst2, lwd=1.5, xlab="", ylab="", xaxs="i", log="y", axes=F)
abline(h=hlst, col=1, lty=3, lwd=1)
mtext("time", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=3.8, cex=1.1)

axis(1); axis(2, las=2); box()
dev.off()

#Plot beta, psf
pdf("figures/memo_171206/beta_psf.pdf", width=5, height=4, colormodel = "cmyk")
par(mar=c(2,2,2,2), oma=c(2,4,0,0))
m<-1
layout(m)
tmp<-beta_psf

matplot(1:nrow(tmp$beta_eig), cbind(rowMeans(tmp$beta_eig), rowMeans(tmp$beta_r0), rowMeans(tmp$beta_0)), type="l", lty=1, col=collst2, lwd=1.5, xlab="", ylab="", xaxs="i", log="y", axes=F)
abline(h=hlst, col=1, lty=3, lwd=1)
mtext("time", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=3.8, cex=1.1)

axis(1); axis(2, las=2); box()
dev.off()

#Plot beta, rps
pdf("figures/memo_171206/beta_rps.pdf", width=5, height=4, colormodel = "cmyk")
par(mar=c(2,2,2,2), oma=c(2,4,0,0))
m<-1
layout(m)
tmp<-beta_rps

matplot(1:nrow(tmp$beta_eig), cbind(rowMeans(tmp$beta_eig), rowMeans(tmp$beta_r0), rowMeans(tmp$beta_0)), type="l", lty=1, col=collst2, lwd=1.5, xlab="", ylab="", xaxs="i", log="y", axes=F)
abline(h=hlst, col=1, lty=3, lwd=1)
mtext("time", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=3.8, cex=1.1)

axis(1); axis(2, las=2); box()
dev.off()












############################################################
# Spatial subset
############################################################
grid_sub<-grid_subset(gridout, size = 0.05)
ptb<-0.2


##### Try Tilman metapopulation model
set.seed(171205)

out_meta<-run_metapopulation(tmax=200, gridout = gridout, population = population_meta, talktime = 0, sites_sub = grid_sub$sites, runtype = "metapopulation")

pdf("figures/memo_171206/levis_spatialsubset_map.pdf", width=5, height=5, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
plot_map(out_meta, gridout = gridout, grid_sub = grid_sub, collst=collst[-1])
shadowtext(50, 50, "5%", cex=2)
mtext("x position", 1, line=2.3, cex=1.1)
mtext("y position", 2, line=2.3, cex=1.1)
dev.off()

eig_meta1<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)
eig_meta2<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation", replace_perturb = 0, talktime=0, prtb=ptb, doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)

r0_meta<-estimate_rarereturn(out_meta, simtime=100, burnin=100, runtype="metapopulation", doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)

out_meta_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_meta, talktime = 0, sites_sub = grid_sub$sites, runtype = "metapopulation")

getEmeta<-getE(out_meta_long, Elst = 2:10, sites_sub = grid_sub$sites)
E_meta<-getEmeta$Eout

invar_meta<-estimate_invar(out_meta_long, E=E_meta, burnin=0, doplot=FALSE, sites_sub = grid_sub$sites)

beta_meta<-beta_estimate(out=out_meta, outlng = out_meta_long, Emat = E_meta, eigout = eig_meta2, r0out = r0_meta, burnin = 10)


##### Try Hubbell neutal model
set.seed(171206)

out_neut<-run_metapopulation(tmax=200, gridout = gridout, population = population_neut, talktime = 0, runtype = "neutral", sites_sub = grid_sub$sites)

eig_neut1<-estimate_eqreturn(out = out_neut, simtime=100, runtype="neutral", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)
eig_neut2<-estimate_eqreturn(out = out_neut, simtime=100, runtype="neutral", replace_perturb = 0, talktime=0, prtb=ptb, doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)

r0_neut<-estimate_rarereturn(out = out_neut, simtime=100, burnin=20, runtype="neutral", doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)

out_neut_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_neut, talktime = 0, runtype = "neutral", sites_sub = grid_sub$sites)

getEneut<-getE(out_neut_long, Elst = 2:10, sites_sub = grid_sub$sites)
E_neut<-getEneut$Eout

invar_neut<-estimate_invar(out_neut_long, E=E_neut, burnin=0, doplot=FALSE, sites_sub = grid_sub$sites)

beta_neut<-beta_estimate(out=out_neut, outlng = out_neut_long, Emat = E_neut, eigout = eig_neut2, r0out = r0_neut, burnin = 10)


##### Try disturbance model
set.seed(171217)

out_dist<-run_metapopulation(tmax=200, gridout = gridout, population = population_dist, talktime = 0, runtype = "disturbance", sites_sub = grid_sub$sites, prt = distlst,  prtfrq = 20)
out_dist_0<-rerunrun_metapopulation(out=out_dist, tmax=0, talktime = 0, runtype = "metapopulation", perturb = distlst, replace_perturb = 0, sites_sub = grid_sub$sites)

eig_dist1<-estimate_eqreturn(out_dist_0, simtime=100, runtype="disturbance", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites, prt = distlst,  prtfrq = 20)
eig_dist2<-estimate_eqreturn(out_dist_0, simtime=100, runtype="disturbance", replace_perturb = 0, talktime=0, prtb=ptb, doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites, prt = distlst,  prtfrq = 20)

r0_dist<-estimate_rarereturn(out_dist_0, simtime=100, burnin=100, runtype="disturbance", doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites, prt = distlst,  prtfrq = 20)

out_dist_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_dist, talktime = 0, runtype = "disturbance", sites_sub = grid_sub$sites, prt = distlst,  prtfrq = 20)

getEdist<-getE(out_dist_long, Elst = 2:10, sites_sub = grid_sub$sites)
E_dist<-getEdist$Eout
#E_dist[E_dist<4]<-4

invar_dist<-estimate_invar(out = out_dist_long, E=E_dist, burnin=100, doplot=FALSE, sites_sub = grid_sub$sites)

beta_dist<-beta_estimate(out=out_dist, outlng = out_dist_long, Emat = E_dist, eigout = eig_dist2, r0out = r0_dist, burnin = 10)



##### Try psf model
set.seed(180108)

out_psf<-run_metapopulation(tmax=200, gridout = gridout, population = population_psf, talktime = 0, runtype = "psf", sites_sub = grid_sub$sites)

eig_psf1<-estimate_eqreturn(out = out_psf, simtime=100, runtype="psf", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)
eig_psf2<-estimate_eqreturn(out = out_psf, simtime=100, runtype="psf", replace_perturb = 0, talktime=0, prtb=ptb, doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)

r0_psf<-estimate_rarereturn(out = out_psf, simtime=100, burnin=100, runtype="psf", doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)

out_psf_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_psf, talktime = 0, runtype = "psf", sites_sub = grid_sub$sites)

getEpsf<-getE(out_psf_long, Elst = 2:10, sites_sub = grid_sub$sites)
E_psf<-getEpsf$Eout

invar_psf<-estimate_invar(out_psf_long, E=E_psf, burnin=100, doplot=FALSE, sites_sub = grid_sub$sites)

beta_psf<-beta_estimate(out=out_psf, outlng = out_psf_long, Emat = E_psf, eigout = eig_psf2, r0out = r0_psf, burnin = 10)



##### Try rps model
set.seed(180108)

out_rps<-run_metapopulation(tmax=200, gridout = gridout, population = population_rps, talktime = 0, runtype = "rps", sites_sub = grid_sub$sites, compmat = intmat_rps)

eig_rps1<-estimate_eqreturn(out = out_rps, simtime=100, runtype="rps", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)
eig_rps2<-estimate_eqreturn(out = out_rps, simtime=100, runtype="rps", replace_perturb = 0, talktime=0, prtb=ptb, doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)

r0_rps<-estimate_rarereturn(out = out_rps, simtime=100, burnin=100, runtype="rps", doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)

out_rps_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_rps, talktime = 0, runtype = "rps", sites_sub = grid_sub$sites, compmat = intmat_rps)

getErps<-getE(out_rps_long, Elst = 2:10, sites_sub = grid_sub$sites)
E_rps<-getErps$Eout

invar_rps<-estimate_invar(out_rps_long, E=E_rps, burnin=100, doplot=FALSE, sites_sub = grid_sub$sites)

beta_rps<-beta_estimate(out=out_rps, outlng = out_rps_long, Emat = E_rps, eigout = eig_rps2, r0out = r0_rps, burnin = 10)















#Plot perturbation, Levins
pdf("figures/memo_171206/levins_perturbation_spsub.pdf", width=4, height=6, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2), rep(3,2))))
layout(m)

tmp<-eig_meta2$out_lst[[1]]$output_spatial
tmp[,1]<-tmp[,1]+200
pmat_meta<-rbind(out_meta$output_spatial, tmp)

matplot(pmat_meta[,1], pmat_meta[,-1]/length(grid_sub$sites), type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=200, lty=2)
ceq<-getceq(clst_meta, mlst_meta)
abline(h=c(sum(ceq), ceq), lty=3, col=collst, lwd=1.5)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
arrows(200, ceq[1]+0.03,
       200, ceq[1],
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_meta2$out_lst[[1]]$output_spatial), abs(eig_meta2$out_lst[[1]]$output_spatial[,2]-eig_meta2$out_lst0$output_spatial[,2])/length(grid_sub$sites), type="l", ylab="estimated distance", xlab="time", col=collst[2], lwd=2, xaxs="i", ylim=c(0, 0.1), cex.lab=1.5); abline(h=0, lty=3)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

plot((1:nrow(eig_meta2$eigenlst))[is.finite(eig_meta2$eigenlst[,1])], (eig_meta2$eigenlst[,1]*(1:nrow(eig_meta2$eigenlst)))[is.finite(eig_meta2$eigenlst[,1])], type="l", lwd=2, col=collst[2], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-4, 1), cex.lab=1.5); abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)




#Total
matplot(pmat_meta[,1], rowSums(pmat_meta[,-1])/length(grid_sub$sites), type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=200, lty=2)
ceq<-getceq(clst_meta, mlst_meta)
abline(h=c(sum(ceq), ceq), lty=3, col=collst, lwd=1.5)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
arrows(200, ceq[1]+0.03,
       200, ceq[1],
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_meta2$out_lst[[1]]$output_spatial), abs(rowSums(eig_meta2$out_lst[[1]]$output_spatial[,-1])-rowSums(eig_meta2$out_lst0$output_spatial[,-1]))/length(grid_sub$sites), type="l", ylab="estimated distance", xlab="time", col=collst[1], lwd=2, xaxs="i", ylim=c(0, 0.08), cex.lab=1.5); abline(h=0, lty=3)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

plot((1:nrow(eig_meta2$eigenlst_tot))[is.finite(eig_meta2$eigenlst_tot[,1])], (eig_meta2$eigenlst_tot[,1]*(1:nrow(eig_meta2$eigenlst_tot)))[is.finite(eig_meta2$eigenlst_tot[,1])], type="l", lwd=2, col=collst[1], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-5, 1), cex.lab=1.5); abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)

dev.off()


#Plot perturbation, disturb
pdf("figures/memo_171206/disturb_perturbation_spsub.pdf", width=4, height=6, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2), rep(3,2))))
layout(m)

tmp<-eig_dist2$out_lst[[1]]$output_spatial
tmp[,1]<-tmp[,1]+max(out_dist$output[,1])
pmat_dist<-rbind(out_dist$output_spatial, tmp)

matplot(pmat_dist[,1], pmat_dist[,-1]/length(grid_sub$sites), type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=max(out_dist$output[,1]), lty=2)
ceq<-getceq(clst_dist, mlst_dist)
abline(h=c(sum(ceq), ceq), lty=3, col=collst, lwd=1.5)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
arrows(max(out_dist$output[,1]), 0.22+0.03,
       max(out_dist$output[,1]), 0.22,
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_dist2$out_lst[[1]]$output_spatial), abs(eig_dist2$out_lst[[1]]$output_spatial[,2]-eig_dist2$out_lst0$output_spatial[,2])/length(grid_sub$sites), type="l", ylab="estimated distance", xlab="time", col=collst[2], lwd=2, xaxs="i", ylim=c(0, 0.15), cex.lab=1.5); abline(h=0, lty=3)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

plot((1:nrow(eig_dist2$eigenlst))[is.finite(eig_dist2$eigenlst[,1])], (eig_dist2$eigenlst[,1]*(1:nrow(eig_dist2$eigenlst)))[is.finite(eig_dist2$eigenlst[,1])], type="l", lwd=2, col=collst[2], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-3, 7), cex.lab=1.5); abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)


#total
matplot(pmat_dist[,1], rowSums(pmat_dist[,-1])/length(grid_sub$sites), type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=max(out_dist$output[,1]), lty=2)
arrows(max(out_dist$output[,1]), 0.4+0.03,
       max(out_dist$output[,1]), 0.4,
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_dist2$out_lst[[1]]$output_spatial), abs(rowSums(eig_dist2$out_lst[[1]]$output_spatial[,-1])-rowSums(eig_dist2$out_lst0$output_spatial[,-1]))/length(grid_sub$sites), type="l", ylab="estimated distance", xlab="time", col=collst[1], lwd=2, xaxs="i", ylim=c(0, 0.15), cex.lab=1.5); abline(h=0, lty=3)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

plot((1:nrow(eig_dist2$eigenlst_tot))[is.finite(eig_dist2$eigenlst_tot[,1])], (eig_dist2$eigenlst_tot[,1]*(1:nrow(eig_dist2$eigenlst_tot)))[is.finite(eig_dist2$eigenlst_tot[,1])], type="l", lwd=2, col=collst[1], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-3, 7), cex.lab=1.5); abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)

dev.off()





#Plot perturbation, Hubbell
pdf("figures/memo_171206/neutral_perturbation_spsub.pdf", width=4, height=6, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2), rep(3,2))))
layout(m)

tmp<-eig_neut1$out_lst[[1]]$output_spatial
tmp[,1]<-tmp[,1]+200
pmat_neut1<-rbind(out_neut$output_spatial, tmp)

tmp<-eig_neut2$out_lst[[1]]$output_spatial
tmp[,1]<-tmp[,1]+200
pmat_neut2<-rbind(out_neut$output_spatial, tmp)

matplot(pmat_neut1[,1], pmat_neut1[,-1]/length(grid_sub$sites), type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=200, lty=2)
ceq<-getceq(clst_neut, mlst_neut)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
arrows(200, 0.28+0.04,
       200, 0.28,
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_neut1$out_lst[[1]]$output_spatial), abs(eig_neut1$out_lst[[1]]$output_spatial[,2]-eig_neut1$out_lst0$output_spatial[,2])/length(grid_sub$sites), type="l", ylab="estimated distance", xlab="time", col=collst[-1], lwd=2, xaxs="i", ylim=c(0, 0.15)); abline(h=0, lty=3)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)

plot((1:nrow(eig_neut1$eigenlst))[is.finite(eig_neut1$eigenlst[,1])], (eig_neut1$eigenlst[,1]*(1:nrow(eig_neut1$eigenlst)))[is.finite(eig_neut1$eigenlst[,1])], type="l", lwd=2, col=collst[-1], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-4, 2)); abline(h=0, lty=3)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)



#Total
matplot(pmat_neut1[,1], rowSums(pmat_neut1[,-1])/length(grid_sub$sites), type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=200, lty=2)
ceq<-getceq(clst_neut, mlst_neut)
abline(h=ceq, lty=2)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
arrows(200, 0.82+0.04,
       200, 0.82,
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_neut1$out_lst[[1]]$output_spatial), abs(rowSums(eig_neut1$out_lst[[1]]$output_spatial[,-1])-rowSums(eig_neut1$out_lst0$output_spatial[,-1]))/length(grid_sub$sites), type="l", ylab="estimated distance", xlab="time", col=collst[1], lwd=2, xaxs="i", ylim=c(0, 0.08)); abline(h=0, lty=3)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)

plot((1:nrow(eig_neut1$eigenlst_tot))[is.finite(eig_neut1$eigenlst_tot[,1])], (eig_neut1$eigenlst_tot[,1]*(1:nrow(eig_neut1$eigenlst_tot)))[is.finite(eig_neut1$eigenlst_tot[,1])], type="l", lwd=2, col=collst[1], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-2, 4)); abline(h=0, lty=3)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)

dev.off()

#Plot perturbation, psf
pdf("figures/memo_171206/psf_perturbation_spsub.pdf", width=4, height=6, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2), rep(3,2))))
layout(m)

tmp<-eig_psf1$out_lst[[1]]$output_spatial
tmp[,1]<-tmp[,1]+200
pmat_psf<-rbind(out_psf$output_spatial, tmp)

matplot(pmat_psf[,1], pmat_psf[,-1]/length(grid_sub$sites), type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=200, lty=2)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
arrows(200, 0.5+0.03,
       200, 0.5,
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_psf2$out_lst[[1]]$output_spatial), abs(eig_psf2$out_lst[[1]]$output_spatial[,2]-eig_psf2$out_lst0$output_spatial[,2])/length(grid_sub$sites), type="l", ylab="estimated distance", xlab="time", col=collst[2], lwd=2, xaxs="i", ylim=c(0, 0.12), cex.lab=1.5); abline(h=0, lty=3)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

plot((1:nrow(eig_psf2$eigenlst))[is.finite(eig_psf2$eigenlst[,1])], (eig_psf2$eigenlst[,1]*(1:nrow(eig_psf2$eigenlst)))[is.finite(eig_psf2$eigenlst[,1])], type="l", lwd=2, col=collst[2], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-4, 1), cex.lab=1.5); abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)


#Total
matplot(pmat_psf[,1], rowSums(pmat_psf[,-1])/length(grid_sub$sites), type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=200, lty=2)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
arrows(200, 0.68+0.03,
       200, 0.68,
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_psf2$out_lst[[1]]$output_spatial), abs(rowSums(eig_psf2$out_lst[[1]]$output_spatial[,-1])-rowSums(eig_psf2$out_lst0$output_spatial[,-1]))/length(grid_sub$sites), type="l", ylab="estimated distance", xlab="time", col=collst[1], lwd=2, xaxs="i", ylim=c(0, 0.12), cex.lab=1.5); abline(h=0, lty=3)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

plot((1:nrow(eig_psf2$eigenlst_tot))[is.finite(eig_psf2$eigenlst_tot[,1])], (eig_psf2$eigenlst_tot[,1]*(1:nrow(eig_psf2$eigenlst_tot)))[is.finite(eig_psf2$eigenlst_tot[,1])], type="l", lwd=2, col=collst[1], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-4, 1), cex.lab=1.5); abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)

dev.off()

#Plot perturbation, rps
pdf("figures/memo_171206/rps_perturbation_spsub.pdf", width=4, height=6, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2), rep(3,2))))
layout(m)

tmp<-eig_rps1$out_lst[[1]]$output_spatial
tmp[,1]<-tmp[,1]+200
pmat_rps<-rbind(out_rps$output_spatial, tmp)

matplot(pmat_rps[,1], pmat_rps[,-1]/length(grid_sub$sites), type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=200, lty=2)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
arrows(200, 0.13+0.03,
       200, 0.13,
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_rps2$out_lst[[1]]$output_spatial), abs(eig_rps2$out_lst[[1]]$output_spatial[,2]-eig_rps2$out_lst0$output_spatial[,2])/length(grid_sub$sites), type="l", ylab="estimated distance", xlab="time", col=collst[2], lwd=2, xaxs="i", ylim=c(0, 0.10), cex.lab=1.5); abline(h=0, lty=3)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

plot((1:nrow(eig_rps2$eigenlst))[is.finite(eig_rps2$eigenlst[,1])], (eig_rps2$eigenlst[,1]*(1:nrow(eig_rps2$eigenlst)))[is.finite(eig_rps2$eigenlst[,1])], type="l", lwd=2, col=collst[2], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-2, 4), cex.lab=1.5); abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)


#Total
matplot(pmat_rps[,1], rowSums(pmat_rps[,-1])/length(grid_sub$sites), type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=200, lty=2)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
arrows(200, 0.76+0.03,
       200, 0.76,
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_rps2$out_lst[[1]]$output_spatial), abs(rowSums(eig_rps2$out_lst[[1]]$output_spatial[,-1])-rowSums(eig_rps2$out_lst0$output_spatial[,-1]))/length(grid_sub$sites), type="l", ylab="estimated distance", xlab="time", col=collst[1], lwd=2, xaxs="i", ylim=c(0, 0.10), cex.lab=1.5); abline(h=0, lty=3)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

plot((1:nrow(eig_rps2$eigenlst_tot))[is.finite(eig_rps2$eigenlst_tot[,1])], (eig_rps2$eigenlst_tot[,1]*(1:nrow(eig_rps2$eigenlst_tot)))[is.finite(eig_rps2$eigenlst_tot[,1])], type="l", lwd=2, col=collst[1], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-3, 1), cex.lab=1.5); abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)

dev.off()












#Plot rare return, Levins
pdf("figures/memo_171206/levis_returnrare_spsub.pdf", width=4, height=4.5, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2))))
layout(m)

tmp<-r0_meta$out0_lst[[1]]$output_spatial
tmp[,1]<-tmp[,1]+200
pmat_meta<-rbind(out_meta$output_spatial, tmp)
tmp<-r0_meta$out_lst[[1]]$output_spatial
tmp[,1]<-tmp[,1]+300
pmat_meta<-rbind(pmat_meta, tmp)

matplot(pmat_meta[,1], pmat_meta[,-1]/length(grid_sub$sites), type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=c(200, 300), lty=2)
ceq<-getceq(clst_meta, mlst_meta)
abline(h=c(sum(ceq), ceq), lty=3, col=collst, lwd=1.5)
abline(h=0, lty=3)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
arrows(200, ceq[1]+0.06,
       200, ceq[1],
       length = 0.08, lwd=2, lend=4)
arrows(300, 0,
       300, 0.06,
       length = 0.08, lwd=2, lend=4)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)

plot(1:nrow(r0_meta$grwrare), r0_meta$grwrare[,1]*(1:nrow(r0_meta$grwrare)), col=collst[2], type="l", lty=1, lwd=1.5, xlab="time span", ylab=expression(italic(paste(r[0], "t"))), xaxs="i")
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(r[0], "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

dev.off()

#Plot rare return, dist
pdf("figures/memo_171206/dist_returnrare_spsub.pdf", width=4, height=4.5, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2))))
layout(m)

mxt<-max(out_dist$output[,1])
mxt2<-max(out_dist$output[,1])+max(r0_dist$out0_lst[[2]]$output[,1])

tmp<-r0_dist$out0_lst[[2]]$output_spatial
tmp[,1]<-tmp[,1]+mxt
pmat_dist<-rbind(out_dist$output_spatial, tmp)
tmp<-r0_dist$out_lst[[2]]$output_spatial
tmp[,1]<-tmp[,1]+mxt2
pmat_dist<-rbind(pmat_dist, tmp)

matplot(pmat_dist[,1], pmat_dist[,-1]/length(grid_sub$sites), type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=c(mxt, mxt2), lty=2)
ceq<-getceq(clst_dist, mlst_dist)
abline(h=c(sum(ceq), ceq), lty=3, col=collst, lwd=1.5)
abline(h=0, lty=3)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
arrows(mxt, 0.11+0.06,
       mxt, 0.11,
       length = 0.08, lwd=2, lend=4)
arrows(mxt2, 0,
       mxt2, 0.06,
       length = 0.08, lwd=2, lend=4)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)

plot(1:nrow(r0_dist$grwrare), r0_dist$grwrare[,2]*(1:nrow(r0_dist$grwrare)), col=collst[3], type="l", lty=1, lwd=1.5, xlab="time span", ylab=expression(italic(paste(r[0], "t"))), xaxs="i", xlim=c(0, 80))
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(r[0], "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)
abline(h=0, lty=3)

dev.off()


#Plot rare return, Hubbell
pdf("figures/memo_171206/neutral_returnrare_spsub.pdf", width=4, height=4.5, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2))))
layout(m)

tmp<-r0_neut$out0_lst[[1]]$output_spatial
tmp[,1]<-tmp[,1]+200
pmat_neut<-rbind(out_neut$output_spatial, tmp)
tmp<-r0_neut$out_lst[[1]]$output_spatial
tmp[,1]<-tmp[,1]+220
pmat_neut<-rbind(pmat_neut, tmp)

matplot(pmat_neut[,1], pmat_neut[,-1]/length(grid_sub$sites), type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=c(200, 220), lty=2)
ceq<-getceq(clst_meta, mlst_meta)
abline(h=c(sum(ceq), ceq), lty=3, col=collst, lwd=1.5)
abline(h=0, lty=3)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
arrows(200, ceq[1]+0.035,
       200, ceq[1],
       length = 0.08, lwd=2, lend=4)
arrows(220, 0,
       220, 0.035,
       length = 0.08, lwd=2, lend=4)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)

plot(1:nrow(r0_neut$grwrare), r0_neut$grwrare[,1]*(1:nrow(r0_neut$grwrare)), col=collst[2], type="l", lty=1, lwd=1.5, xlab="time span", ylab=expression(italic(paste(r[0], "t"))), xaxs="i")
abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(r[0], "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

dev.off()

#Plot rare return, psf
pdf("figures/memo_171206/psf_returnrare_spsub.pdf", width=4, height=4.5, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2))))
layout(m)

tmp<-r0_psf$out0_lst[[1]]$output_spatial
tmp[,1]<-tmp[,1]+200
pmat_psf<-rbind(out_psf$output_spatial, tmp)
tmp<-r0_psf$out_lst[[1]]$output_spatial
tmp[,1]<-tmp[,1]+300
pmat_psf<-rbind(pmat_psf, tmp)

matplot(pmat_psf[,1], pmat_psf[,-1]/length(grid_sub$sites), type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=c(200, 300), lty=2)
abline(h=0, lty=3)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
arrows(200, 0.5+0.06,
       200, 0.5,
       length = 0.08, lwd=2, lend=4)
arrows(300, 0,
       300, 0.06,
       length = 0.08, lwd=2, lend=4)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)

plot(1:nrow(r0_psf$grwrare), r0_psf$grwrare[,1]*(1:nrow(r0_psf$grwrare)), col=collst[2], type="l", lty=1, lwd=1.5, xlab="time span", ylab=expression(italic(paste(r[0], "t"))), xaxs="i")
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(r[0], "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

dev.off()

#Plot rare return, rps
pdf("figures/memo_171206/rps_returnrare_spsub.pdf", width=4, height=4.5, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2))))
layout(m)

tmp<-r0_rps$out0_lst[[1]]$output_spatial
tmp[,1]<-tmp[,1]+200
pmat_rps<-rbind(out_rps$output_spatial, tmp)
tmp<-r0_rps$out_lst[[1]]$output_spatial
tmp[,1]<-tmp[,1]+300
pmat_rps<-rbind(pmat_rps, tmp)

matplot(pmat_rps[,1], pmat_rps[,-1]/length(grid_sub$sites), type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=c(200, 300), lty=2)
abline(h=0, lty=3)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
arrows(200, 0.1+0.06,
       200, 0.1,
       length = 0.08, lwd=2, lend=4)
arrows(300, 0,
       300, 0.06,
       length = 0.08, lwd=2, lend=4)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)

plot(1:nrow(r0_rps$grwrare), r0_rps$grwrare[,1]*(1:nrow(r0_rps$grwrare)), col=collst[2], type="l", lty=1, lwd=1.5, xlab="time span", ylab=expression(italic(paste(r[0], "t"))), xaxs="i")
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(r[0], "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

dev.off()






#Plot invar, Levins
pdf("figures/memo_171206/levins_invar_spsub.pdf", width=5, height=4, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 2), rep(2,2))))
layout(m)

matplot(out_meta_long$output_spatial[,1], out_meta_long$output_spatial[,2]/length(grid_sub$sites), type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i", ylim=c(0.24, 0.41))
abline(h=getceq(clst_meta, mlst_meta)[1], col=collst[-1], lty=2, lwd=1.5)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)

segments(c(200, 200, 300, 300), c(0.27, 0.40, 0.40, 0.27), c(200, 300, 300, 200), c(0.40, 0.40, 0.27, 0.27), lwd=2)
segments(c(200, 200, 300, 300)+500, c(0.27, 0.40, 0.40, 0.27), c(200, 300, 300, 200)+500, c(0.40, 0.40, 0.27, 0.27), lwd=2)

arrows(500, 0.27, 690, 0.27, lwd=2, length = 0.1, lend=2)
arrows(500, 0.27, 310, 0.27, lwd=2, length = 0.1, lend=2)

text(500, 0.27, pos=1, labels = "time lag")

text(250, 0.27, pos=1, labels = "training set")
text(750, 0.27, pos=1, labels = "testing set")


plot(invar_meta$pdlag_list[[1]]$laglst, invar_meta$pdlag_list[[1]]$CVest[,1], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[2], xaxs="i", ylim=c(0, 0.1)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)



#Tot
matplot(out_meta_long$output_spatial[,1], rowSums(out_meta_long$output_spatial[,-1])/length(grid_sub$sites), type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(h=sum(getceq(clst_meta, mlst_meta)), col=collst[1], lty=2, lwd=1.5)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)

plot(invar_meta$pdlag_list_tot[[1]]$laglst, invar_meta$pdlag_list_tot[[1]]$CVest[,1], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[1], xaxs="i", ylim=c(0, 0.1)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)

dev.off()



#Plot invar, dist
pdf("figures/memo_171206/dist_invar_spsub.pdf", width=5, height=4, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 2), rep(2,2))))
layout(m)

matplot(out_dist_long$output_spatial[,1], out_dist_long$output_spatial[,2]/length(grid_sub$sites), type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i", ylim=c(0, 0.35))
abline(h=0, lty=3, lwd=1.5)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)
abline(v=100, lty=3)

segments(c(200, 200, 300, 300), c(0.01, 0.28, 0.28, 0.01), c(200, 300, 300, 200), c(0.28, 0.28, 0.01, 0.01), lwd=2)
segments(c(200, 200, 300, 300)+500, c(0.01, 0.28, 0.28, 0.01), c(200, 300, 300, 200)+500, c(0.28, 0.28, 0.01, 0.01), lwd=2)

arrows(500, 0.28, 690, 0.28, lwd=2, length = 0.1, lend=2)
arrows(500, 0.28, 310, 0.28, lwd=2, length = 0.1, lend=2)

text(500, 0.28, pos=3, labels = "time lag")

text(250, 0.28, pos=3, labels = "training set")
text(750, 0.28, pos=3, labels = "testing set")


plot(invar_dist$pdlag_list[[1]]$laglst, invar_dist$pdlag_list[[1]]$CVest[,1], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[2], xaxs="i", ylim=c(0, 0.8)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)


#Total
matplot(out_dist_long$output_spatial[,1], rowSums(out_dist_long$output_spatial[,-1])/length(grid_sub$sites), type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(h=0, lty=3, lwd=1.5)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)
abline(v=100, lty=3)

plot(invar_dist$pdlag_list_tot[[1]]$laglst, invar_dist$pdlag_list_tot[[1]]$CVest[,1], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[1], xaxs="i", ylim=c(0, 0.8)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)

dev.off()


#Plot invar, Hubbell
pdf("figures/memo_171206/neutral_invar_spsub.pdf", width=5, height=4, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 2), rep(2,2))))
layout(m)

matplot(out_neut_long$output[,1], out_neut_long$output[,2]/out_neut_long$plotdata$ngrid, type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)

tmp<-out_neut_long$output[,2]/out_neut_long$plotdata$ngrid
tmptm<-out_neut_long$output[,1]

y11<-min(tmp[tmptm>200 & tmptm<300]); y12<-max(tmp[tmptm>200 & tmptm<300])
y21<-min(tmp[tmptm>700 & tmptm<800]); y22<-max(tmp[tmptm>700 & tmptm<800])
segments(c(200, 200, 300, 300), c(y11, y12, y12, y11), c(200, 300, 300, 200), c(y12, y12, y11, y11), lwd=2)
segments(c(200, 200, 300, 300)+500, c(y22, y21, y21, y22), c(200, 300, 300, 200)+500, c(y21, y21, y22, y22), lwd=2)

arrows(500, (y11+y21)/2, 690, y21, lwd=2, length = 0.1, lend=2)
arrows(500, (y11+y21)/2, 310, y11, lwd=2, length = 0.1, lend=2)

text(500, (y11+y21)/2, pos=1, labels = "time lag")

text(250, y11, pos=1, labels = "training set")
text(750, y22, pos=3, labels = "testing set")


plot(invar_neut$pdlag_list[[1]]$laglst, invar_neut$pdlag_list[[1]]$CVest[,1], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[2], xaxs="i", ylim=c(0, 0.4)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)


#Total
matplot(out_neut_long$output[,1], rowSums(out_neut_long$output[,-1])/out_neut_long$plotdata$ngrid, type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)
abline(h=getceq(clst_neut, mlst_neut), lty=2)

tmp<-out_neut_long$output[,2]/out_neut_long$plotdata$ngrid
tmptm<-out_neut_long$output[,1]

plot(invar_neut$pdlag_list_tot[[1]]$laglst, invar_neut$pdlag_list_tot[[1]]$CVest[,1], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[1], xaxs="i", ylim=c(0, 0.04)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)

dev.off()


#Plot invar, psf
pdf("figures/memo_171206/psf_invar_spsub.pdf", width=5, height=4, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 2), rep(2,2))))
layout(m)

matplot(out_psf_long$output_spatial[,1], out_psf_long$output_spatial[,2]/length(grid_sub$sites), type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i", ylim=c(0.24, 0.6))
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)
abline(v=100, lty=3)

tmp<-out_psf_long$output_spatial[,2]/length(out_psf_long$sites_sub)
tmptm<-out_psf_long$output_spatial[,1]

y11<-min(tmp[tmptm>200 & tmptm<300]); y12<-max(tmp[tmptm>200 & tmptm<300])
y21<-min(tmp[tmptm>700 & tmptm<800]); y22<-max(tmp[tmptm>700 & tmptm<800])
segments(c(200, 200, 300, 300), c(y11, y12, y12, y11), c(200, 300, 300, 200), c(y12, y12, y11, y11), lwd=2)
segments(c(200, 200, 300, 300)+500, c(y22, y21, y21, y22), c(200, 300, 300, 200)+500, c(y21, y21, y22, y22), lwd=2)

arrows(500, (y11+y21)/2, 690, y21, lwd=2, length = 0.1, lend=2)
arrows(500, (y11+y21)/2, 310, y11, lwd=2, length = 0.1, lend=2)

text(500, (y11+y21)/2, pos=1, labels = "time lag")

text(250, y11, pos=1, labels = "training set")
text(750, y21, pos=1, labels = "testing set")


plot(invar_psf$pdlag_list[[1]]$laglst, invar_psf$pdlag_list[[1]]$CVest[,1], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[2], xaxs="i", ylim=c(0, 0.1)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)



#Total
matplot(out_psf_long$output_spatial[,1], rowSums(out_psf_long$output_spatial[,-1])/length(grid_sub$sites), type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)
abline(v=100, lty=3)

tmp<-out_psf_long$output_spatial[,2]/length(out_psf_long$sites_sub)
tmptm<-out_psf_long$output_spatial[,1]

plot(invar_psf$pdlag_list_tot[[1]]$laglst, invar_psf$pdlag_list_tot[[1]]$CVest[,1], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[1], xaxs="i", ylim=c(0, 0.1)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)
dev.off()




#Plot invar, rps
pdf("figures/memo_171206/rps_invar_spsub.pdf", width=5, height=4, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 2), rep(2,2))))
layout(m)

matplot(out_rps_long$output_spatial[,1], out_rps_long$output_spatial[,2]/length(grid_sub$sites), type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i", ylim=c(0.05, 0.4))
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)
abline(v=100, lty=3)

tmp<-out_rps_long$output_spatial[,2]/length(out_rps_long$sites_sub)
tmptm<-out_rps_long$output_spatial[,1]

y11<-min(tmp[tmptm>200 & tmptm<300]); y12<-max(tmp[tmptm>200 & tmptm<300])
y21<-min(tmp[tmptm>700 & tmptm<800]); y22<-max(tmp[tmptm>700 & tmptm<800])
segments(c(200, 200, 300, 300), c(y11, y12, y12, y11), c(200, 300, 300, 200), c(y12, y12, y11, y11), lwd=2)
segments(c(200, 200, 300, 300)+500, c(y22, y21, y21, y22), c(200, 300, 300, 200)+500, c(y21, y21, y22, y22), lwd=2)

arrows(500, (y11+y21)/2, 690, y21, lwd=2, length = 0.1, lend=2)
arrows(500, (y11+y21)/2, 310, y11, lwd=2, length = 0.1, lend=2)

text(500, (y11+y21)/2, pos=1, labels = "time lag")

text(250, y12, pos=3, labels = "training set")
text(750, y21, pos=1, labels = "testing set")


plot(invar_rps$pdlag_list[[1]]$laglst, invar_rps$pdlag_list[[1]]$CVest[,1], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[2], xaxs="i", ylim=c(0, 0.27)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)



#TOtal
matplot(out_rps_long$output_spatial[,1], rowSums(out_rps_long$output_spatial[,-1])/length(grid_sub$sites), type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
put.fig.letter("a.", "topleft", offset=ofs3, cex=1.6)
abline(v=100, lty=3)

plot(invar_rps$pdlag_list_tot[[1]]$laglst, invar_rps$pdlag_list_tot[[1]]$CVest[,1], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[1], xaxs="i", ylim=c(0, 0.05)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)

dev.off()







#Plot beta, Levins
pdf("figures/memo_171206/beta_levins_subsp.pdf", width=5, height=4, colormodel = "cmyk")
par(mar=c(2,2,2,2), oma=c(2,4,0,0))
m<-1
layout(m)
tmp<-beta_meta

matplot(1:nrow(tmp$beta_eig), cbind(rowMeans(tmp$beta_eig), rowMeans(tmp$beta_r0), rowMeans(tmp$beta_0)), type="l", lty=1, col=collst2, lwd=1.5, xlab="", ylab="", xaxs="i", log="y", axes=F)
abline(h=hlst, col=1, lty=3, lwd=1)
mtext("time", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=3.8, cex=1.1)

axis(1); axis(2, las=2); box()
dev.off()

#Plot beta, Neutral
pdf("figures/memo_171206/beta_neut_subsp.pdf", width=5, height=4, colormodel = "cmyk")
par(mar=c(2,2,2,2), oma=c(2,4,0,0))
m<-1
layout(m)
tmp<-beta_neut

matplot(1:nrow(tmp$beta_eig), cbind(rowMeans(tmp$beta_eig), rowMeans(tmp$beta_r0), rowMeans(tmp$beta_0)), type="l", lty=1, col=collst2, lwd=1.5, xlab="", ylab="", xaxs="i", log="y", axes=F)
abline(h=hlst, col=1, lty=3, lwd=1)
mtext("time", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=3.8, cex=1.1)

axis(1); axis(2, las=2); box()
dev.off()

#Plot beta, Dist
pdf("figures/memo_171206/beta_dist_subsp.pdf", width=5, height=4, colormodel = "cmyk")
par(mar=c(2,2,2,2), oma=c(2,4,0,0))
m<-1
layout(m)
tmp<-beta_dist

matplot(1:nrow(tmp$beta_eig), cbind(rowMeans(tmp$beta_eig), rowMeans(tmp$beta_r0), rowMeans(tmp$beta_0)), type="l", lty=1, col=collst2, lwd=1.5, xlab="", ylab="", xaxs="i", log="y", axes=F)
abline(h=hlst, col=1, lty=3, lwd=1)
mtext("time", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=3.8, cex=1.1)

axis(1); axis(2, las=2); box()
dev.off()

#Plot beta, psf
pdf("figures/memo_171206/beta_psf_subsp.pdf", width=5, height=4, colormodel = "cmyk")
par(mar=c(2,2,2,2), oma=c(2,4,0,0))
m<-1
layout(m)
tmp<-beta_psf

matplot(1:nrow(tmp$beta_eig), cbind(rowMeans(tmp$beta_eig), rowMeans(tmp$beta_r0), rowMeans(tmp$beta_0)), type="l", lty=1, col=collst2, lwd=1.5, xlab="", ylab="", xaxs="i", log="y", axes=F)
abline(h=hlst, col=1, lty=3, lwd=1)
mtext("time", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=3.8, cex=1.1)

axis(1); axis(2, las=2); box()
dev.off()

#Plot beta, rps
pdf("figures/memo_171206/beta_rps_subsp.pdf", width=5, height=4, colormodel = "cmyk")
par(mar=c(2,2,2,2), oma=c(2,4,0,0))
m<-1
layout(m)
tmp<-beta_rps

matplot(1:nrow(tmp$beta_eig), cbind(rowMeans(tmp$beta_eig), rowMeans(tmp$beta_r0), rowMeans(tmp$beta_0)), type="l", lty=1, col=collst2, lwd=1.5, xlab="", ylab="", xaxs="i", log="y", axes=F)
abline(h=hlst, col=1, lty=3, lwd=1)
mtext("time", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=3.8, cex=1.1)

axis(1); axis(2, las=2); box()
dev.off()


