error
rm(list=ls())
setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

#load functions
require(rEDM)
source("run_metapopulation_wrapper.R")
require(RColorBrewer)
source("~/Dropbox/Rfunctions/figure_functions.R")

############################################################
# "Global" run
############################################################
##### Try Tilman metapopulation model
gridout<-makegrid(xlng = 100, ylng = 100)
xfac<-5
ptb<-0.1
collst<-c("black", brewer.pal(4, "Dark2"))

set.seed(171205)
#clst_meta = c(0.15, 0.3, 0.8, 3)*xfac
clst_meta = c(0.15, 0.35, 0.9, 2.8)*xfac
mlst_meta = rep(0.1, length(clst_meta))*xfac
getceq(clst_meta, mlst_meta)

population_meta<-populate(gridout, nlst = floor(getceq(clst_meta, mlst_meta)*prod(gridout$lng)),
                          clst = clst_meta, radlst = Inf, mlst = mlst_meta)
out_meta<-run_metapopulation(tmax=200, gridout = gridout, population = population_meta, talktime = 0)
eig_meta1<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE)

r0_meta<-estimate_rarereturn(out_meta, simtime=100, burnin=100, runtype="metapopulation", doplot = FALSE)

out_meta_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_meta, talktime = 0)
getEmeta<-getE(out_meta_long, Elst = 2:10)
E_meta<-getEmeta$Eout
E_meta[E_meta<4]<-4

invar_meta<-estimate_invar(out_meta_long, E=E_meta, burnin=0, doplot=FALSE)

pdf("figures/memo_171206/levis_spatialsubset_map.pdf", width=5, height=5, colormodel = "cmyk")
par(mar=c(2,2,2,2), oma=c(2,2,0,0))

plot_map(out_meta, gridout = gridout, grid_sub = grid_sub, collst=collst[-1])

shadowtext(50, 50, "5%", cex=2)
mtext("x position", 1, line=2.3, cex=1.1)
mtext("y position", 2, line=2.3, cex=1.1)

dev.off()


##### Try Hubbell neutal model
set.seed(171206)

clst_neut<-rep(0.5, 4)*xfac
mlst_neut<-rep(0.1, length(clst_neut))*xfac
population_neut<-populate(gridout, nlst = round(rep(unique(abs(getceq(clst_neut, mlst_neut)))/length(clst_neut), length(clst_neut))*prod(gridout$lng)),
                          clst = clst_neut, radlst = Inf, mlst = mlst_neut)
out_neut<-run_metapopulation(tmax=200, gridout = gridout, population = population_neut, talktime = 0, runtype = "neutral")

eig_neut1<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE)
eig_neut2<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral", replace_perturb = 0, talktime=0, prtb=ptb, doplot = FALSE)

r0_neut<-estimate_rarereturn(out_neut, simtime=100, burnin=100, runtype="neutral", doplot = FALSE)

out_neut_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_neut, talktime = 0, runtype = "neutral")
getEneut<-getE(out_neut_long, Elst = 2:10)
E_neut<-getEneut$Eout
E_neut[E_neut<4]<-4

invar_neut<-estimate_invar(out_neut_long, E=E_neut, burnin=0, doplot=FALSE)




#Plot perturbation, Levins
pdf("figures/memo_171206/levins_perturbation.pdf", width=4, height=6, colormodel = "cmyk")

ofs1<-c(0.075, -0.04)
ofs2<-c(0.075, -0.09)

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2), rep(3,2))))
layout(m)

tmp<-eig_meta1$out_lst[[1]]$output
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

plot(1:nrow(eig_meta1$out_lst[[1]]$output), abs(eig_meta1$out_lst[[1]]$output[,2]-eig_meta1$out_lst0$output[,2])/out_meta$plotdata$ngrid, type="l", ylab="estimated distance", xlab="time", col=collst[2], lwd=2, xaxs="i", ylim=c(0, 0.035), cex.lab=1.5); abline(h=0, lty=3)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

plot(1:nrow(eig_meta1$eigenlst), eig_meta1$eigenlst[,1]*(1:nrow(eig_meta1$eigenlst)), type="l", lwd=2, col=collst[2], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-6, 0), cex.lab=1.5); abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)

dev.off()


#Plot perturbation, Hubbell
pdf("figures/memo_171206/neutral_perturbation.pdf", width=4, height=6, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2), rep(3,2))))
layout(m)

tmp<-eig_neut1$out_lst[[3]]$output
tmp[,1]<-tmp[,1]+200
pmat_neut1<-rbind(out_neut$output, tmp)

tmp<-eig_neut2$out_lst[[3]]$output
tmp[,1]<-tmp[,1]+200
pmat_neut2<-rbind(out_neut$output, tmp)

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2), rep(3,2))))
layout(m)

#with replacing perturbation
matplot(pmat_neut1[,1], pmat_neut1[,-1]/out_neut$plotdata$ngrid, type="l", lty=1, col=collst[c(4,3,2,5)], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=200, lty=2)
ceq<-getceq(clst_neut, mlst_neut)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
arrows(200, 0.32+0.02,
       200, 0.32,
       length = 0.08, lwd=2, lend=4)


plot(1:nrow(eig_neut1$out_lst[[3]]$output), abs(eig_neut1$out_lst[[3]]$output[,4]-eig_neut1$out_lst0$output[,4])/out_neut$plotdata$ngrid, type="l", ylab="estimated distance", xlab="time", col=collst[2], lwd=2, xaxs="i", ylim=c(0, 0.1)); abline(h=0, lty=3)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)

plot(1:nrow(eig_neut1$eigenlst), eig_neut1$eigenlst[,3]*(1:nrow(eig_neut1$eigenlst)), type="l", lwd=2, col=collst[2], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-4, 2)); abline(h=0, lty=3)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)

dev.off()

pdf("figures/memo_171206/neutral_perturbation_noreplace.pdf", width=4, height=6, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2), rep(3,2))))
layout(m)

#without replacing perturbation
matplot(pmat_neut2[,1], pmat_neut2[,-1]/out_neut$plotdata$ngrid, type="l", lty=1, col=collst[c(4,3,2,5)], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=200, lty=2)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
arrows(200, 0.32+0.02,
       200, 0.32,
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_neut2$out_lst[[3]]$output), abs(eig_neut2$out_lst[[3]]$output[,4]-eig_neut2$out_lst0$output[,4])/out_neut$plotdata$ngrid, type="l", ylab="estimated distance", xlab="time", col=collst[2], lwd=2, xaxs="i", ylim=c(0, 0.1)); abline(h=0, lty=3)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)

plot((1:nrow(eig_neut2$eigenlst))[is.finite(eig_neut2$eigenlst[,3])], eig_neut2$eigenlst[is.finite(eig_neut2$eigenlst[,3]),3]*(1:nrow(eig_neut2$eigenlst))[is.finite(eig_neut2$eigenlst[,3])], type="l", lwd=2, col=collst[2], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-4, 2)); abline(h=0, lty=3)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)

dev.off()


#Hubbell: population-level stability
pdf("figures/memo_171206/neutral_totpop_stability.pdf", width=4, height=3, colormodel = "cmyk")
par(mar=c(2,2,2,2), oma=c(2,2,0,0))

matplot(pmat_neut2[,1], rowSums(pmat_neut2[,-1])/out_neut$plotdata$ngrid, type="l", lty=1, col=collst[1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=200, lty=2)
abline(h=c(mean(abs(ceq))/c(1, 4)), lty=3, col=1, lwd=1.5)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)

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


#Plot invar, Levins
pdf("figures/memo_171206/levins_invar.pdf", width=5, height=4, colormodel = "cmyk")

ofs3<-c(0.075, -0.090)

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

segments(c(200, 200, 300, 300), c(0.22, 0.29, 0.29, 0.22), c(200, 300, 300, 200), c(0.29, 0.29, 0.22, 0.22), lwd=2)
segments(c(200, 200, 300, 300)+500, c(0.054, 0.11, 0.11, 0.054), c(200, 300, 300, 200)+500, c(0.11, 0.11, 0.054, 0.054), lwd=2)

arrows(500, 0.165, 690, 0.11, lwd=2, length = 0.1, lend=2)
arrows(500, 0.165, 310, 0.22, lwd=2, length = 0.1, lend=2)

text(500, 0.175, pos=3, labels = "time lag")

text(250, 0.22, pos=1, labels = "training set")
text(750, 0.11, pos=3, labels = "testing set")


plot(invar_neut$pdlag_list[[1]]$laglst, invar_neut$pdlag_list[[1]]$CVest[,1], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[2], xaxs="i", ylim=c(0, 0.65)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)

dev.off()









############################################################
# Spatial subset
############################################################
##### Try Tilman metapopulation model
grid_sub<-grid_subset(gridout, size = 0.05)
set.seed(171205)
ptb<-0.2

out_meta<-run_metapopulation(tmax=200, gridout = gridout, population = population_meta, talktime = 0, sites_sub = grid_sub$sites, runtype = "metapopulation_spatial")

pdf("figures/memo_171206/levis_spatialsubset_map.pdf", width=5, height=5, colormodel = "cmyk")
par(mar=c(2,2,2,2), oma=c(2,2,0,0))

plot_map(out_meta, gridout = gridout, grid_sub = grid_sub, collst=collst[-1])

shadowtext(50, 50, "5%", cex=2)
mtext("x position", 1, line=2.3, cex=1.1)
mtext("y position", 2, line=2.3, cex=1.1)

dev.off()

eig_meta1<-estimate_eqreturn(out_meta, simtime=100, runtype="metapopulation_spatial", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)
r0_meta<-estimate_rarereturn(out_meta, simtime=100, burnin=100, runtype="metapopulation_spatial", doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)

out_meta_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_meta, talktime = 0, sites_sub = grid_sub$sites, runtype = "metapopulation_spatial")

getEmeta<-getE(out_meta_long, Elst = 2:10, sites_sub = grid_sub$sites)
E_meta<-getEmeta$Eout
E_meta[E_meta<4]<-4

invar_meta<-estimate_invar(out_meta_long, E=E_meta, burnin=0, doplot=FALSE, sites_sub = grid_sub$sites)


##### Try Hubbell neutal model
set.seed(171206)
ptb<-0.2

out_neut<-run_metapopulation(tmax=200, gridout = gridout, population = population_neut, talktime = 0, runtype = "neutral_spatial", sites_sub = grid_sub$sites)

eig_neut1<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral_spatial", replace_perturb = 1, talktime=0, prtb=ptb, doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)
eig_neut2<-estimate_eqreturn(out_neut, simtime=100, runtype="neutral_spatial", replace_perturb = 0, talktime=0, prtb=ptb, doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)

r0_neut<-estimate_rarereturn(out_neut, simtime=100, burnin=100, runtype="neutral_spatial", doplot = FALSE, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites)

out_neut_long<-run_metapopulation(tmax=1000, gridout = gridout, population = population_neut, talktime = 0, runtype = "neutral_spatial", sites_sub = grid_sub$sites)

getEneut<-getE(out_neut_long, Elst = 2:10, sites_sub = grid_sub$sites)
E_neut<-getEneut$Eout
E_neut[E_neut<4]<-4

invar_neut<-estimate_invar(out_neut_long, E=E_neut, burnin=0, doplot=FALSE, sites_sub = grid_sub$sites)






#Plot perturbation, Levins
pdf("figures/memo_171206/levins_perturbation_spsub.pdf", width=4, height=6, colormodel = "cmyk")

ofs1<-c(0.075, -0.04)
ofs2<-c(0.075, -0.09)

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2), rep(3,2))))
layout(m)

tmp<-eig_meta1$out_lst[[1]]$output_spatial
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

plot(1:nrow(eig_meta1$out_lst[[1]]$output_spatial), abs(eig_meta1$out_lst[[1]]$output_spatial[,2]-eig_meta1$out_lst0$output_spatial[,2])/length(grid_sub$sites), type="l", ylab="estimated distance", xlab="time", col=collst[2], lwd=2, xaxs="i", ylim=c(0, 0.1), cex.lab=1.5); abline(h=0, lty=3)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)

plot((1:nrow(eig_meta1$eigenlst))[is.finite(eig_meta1$eigenlst[,1])], (eig_meta1$eigenlst[,1]*(1:nrow(eig_meta1$eigenlst)))[is.finite(eig_meta1$eigenlst[,1])], type="l", lwd=2, col=collst[2], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-4, 2), cex.lab=1.5); abline(h=0, lty=3)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)

dev.off()


#Plot perturbation, Hubbell
pdf("figures/memo_171206/neutral_perturbation_spsub.pdf", width=4, height=6, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2), rep(3,2))))
layout(m)

tmp<-eig_neut1$out_lst[[3]]$output_spatial
tmp[,1]<-tmp[,1]+200
pmat_neut1<-rbind(out_neut$output_spatial, tmp)

tmp<-eig_neut2$out_lst[[3]]$output_spatial
tmp[,1]<-tmp[,1]+200
pmat_neut2<-rbind(out_neut$output_spatial, tmp)

matplot(pmat_neut1[,1], pmat_neut1[,-1]/length(grid_sub$sites), type="l", lty=1, col=collst[c(3,4,2,5)], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
abline(v=200, lty=2)
ceq<-getceq(clst_neut, mlst_neut)
put.fig.letter("a.", "topleft", offset=ofs1, cex=1.6)
mtext("time", 1, line=2.3, cex=1.1)
mtext("relative abundance", 2, line=2.3, cex=1.1)
arrows(200, 0.28+0.04,
       200, 0.28,
       length = 0.08, lwd=2, lend=4)

plot(1:nrow(eig_neut1$out_lst[[3]]$output_spatial), abs(eig_neut1$out_lst[[3]]$output_spatial[,4]-eig_neut1$out_lst0$output_spatial[,4])/length(grid_sub$sites), type="l", ylab="estimated distance", xlab="time", col=collst[2], lwd=2, xaxs="i", ylim=c(0, 0.1)); abline(h=0, lty=3)
put.fig.letter("b.", "topleft", offset=ofs2, cex=1.6)
mtext("time since disturbance", 1, line=2.3, cex=1.1)
mtext("estimated distance", 2, line=2.3, cex=1.1)

plot((1:nrow(eig_neut1$eigenlst))[is.finite(eig_neut1$eigenlst[,3])], (eig_neut1$eigenlst[,3]*(1:nrow(eig_neut1$eigenlst)))[is.finite(eig_neut1$eigenlst[,3])], type="l", lwd=2, col=collst[2], xlab="time span", ylab=expression(italic(paste(lambda, "t"))), xaxs="i", ylim=c(-4, 2)); abline(h=0, lty=3)
put.fig.letter("c.", "topleft", offset=ofs2, cex=1.6)
mtext("time span", 1, line=2.3, cex=1.1)
mtext(expression(italic(paste(lambda, "t"))), 2, line=2.3, cex=1.1)

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


#Plot rare return, Hubbell
pdf("figures/memo_171206/neutral_returnrare_spsub.pdf", width=4, height=4.5, colormodel = "cmyk")

par(mar=c(2,2,2,2), oma=c(2,2,0,0))
m<-t(matrix(nrow=2, c(rep(1, 4), rep(2,2))))
layout(m)

tmp<-r0_neut$out0_lst[[1]]$output_spatial
tmp[,1]<-tmp[,1]+200
pmat_neut<-rbind(out_neut$output_spatial, tmp)
tmp<-r0_neut$out_lst[[1]]$output_spatial
tmp[,1]<-tmp[,1]+300
pmat_neut<-rbind(pmat_neut, tmp)

matplot(pmat_neut[,1], pmat_neut[,-1]/length(grid_sub$sites), type="l", lty=1, col=collst[-1], lwd=1.5, xlab="time", ylab="relative abundance", xaxs="i")
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

segments(c(200, 200, 300, 300), c(0.22, 0.29, 0.29, 0.22), c(200, 300, 300, 200), c(0.29, 0.29, 0.22, 0.22), lwd=2)
segments(c(200, 200, 300, 300)+500, c(0.26, 0.32, 0.32, 0.26), c(200, 300, 300, 200)+500, c(0.32, 0.32, 0.26, 0.26), lwd=2)

arrows(500, 0.24, 690, 0.26, lwd=2, length = 0.1)
arrows(500, 0.24, 310, 0.22, lwd=2, length = 0.1)

text(500, 0.24, pos=1, labels = "time lag")

text(250, 0.22, pos=1, labels = "training set")
text(750, 0.26, pos=1, labels = "testing set")


plot(invar_neut$pdlag_list[[1]]$laglst, invar_neut$pdlag_list[[1]]$CVest[,1], xlab="time lag", ylab=expression(italic(paste("CV"))), type="l", lty=1, lwd=1.5, col=collst[2], xaxs="i", ylim=c(0, 0.4)); abline(h=0, lty=3)
mtext("time lag", 1, line=2.3, cex=1.1)
mtext(expression(italic("CV")), 2, line=2.3, cex=1.1)
put.fig.letter("b.", "topleft", offset=ofs3, cex=1.6)

dev.off()


