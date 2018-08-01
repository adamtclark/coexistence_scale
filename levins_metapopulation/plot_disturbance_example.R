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
# Run disturbance
############################################################
set.seed(1913)

##### Disturbance model
clst_dist= c(0.145, 0.2)*xfac
mlst_dist = rep(0.1, length(clst_dist))*xfac
distlst<-c(0.95, 0)
prtfrq<-50

population_dist<-populate(gridout, nlst = rep(floor(prod(gridout$lng)/length(clst_dist)*0.8), length(clst_dist)),
                          clst = clst_dist, radlst = Inf, mlst = mlst_dist)
set.seed(171217)
out_dist<-run_metapopulation(tmax=2000, gridout = gridout, population = population_dist, talktime = 0, runtype = "disturbance", prt = distlst,  prtfrq = prtfrq)
out_dist_0<-rerunrun_metapopulation(out=out_dist, tmax=300, talktime = 0, runtype = "metapopulation", perturb = distlst, replace_perturb = 0)
out_dist_0$output[,1]<-out_dist_0$output[,1]+max(out_dist$output[,1])

totdat<-rbind(out_dist$output, out_dist_0$output)
totdat<-totdat[totdat[,1]>500,]
totdat[,1]<-totdat[,1]-500
totdat<-totdat[is.finite(rowSums(totdat)),]

pdf("figures/SUP_FIGURE_disturbance_example.pdf", width=9, height=4, colormodel = "cmyk")
  matplot(totdat[,1], totdat[,-1]/prod(gridout$lng), col=collst[-1], type="l", lty=1, lwd=2, xlab="simulation time", ylab="species abundance",
          xlim=c(0, max(totdat[,1])), xaxs="i")
  abline(v=seq(0, 1500, by=50), lty=2)
  abline(h=0, lty=3)
dev.off()

