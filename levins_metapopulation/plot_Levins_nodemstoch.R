error
rm(list=ls())
setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

source("util/figure_functions.R")
source("util/filled.contour3.R")
require(deSolve)


#Functions
LEVmod <- function(t, x, parms, input)  {
  with(as.list(c(parms, x)), {
    p1_col<-c1*(p1O+p1I)

    dp1O <- p1_col*((1-spext)-p1O)-m1*p1O                        #outside of focal patch
    dp1I <- p1_col*(spext-p1I)-m1*p1I                        #inside of focal patch
    
    p2_col<-c2*(p2O+p2I)
    dp2O <- p2_col*((1-spext)-p1O-p2O)-p1_col*p2O-m2*p2O  #outside of focal patch
    dp2I <- p2_col*(spext-p1I-p2I)-p1_col*p2I-m2*p2I  #inside of focal patch
    
    if(hold1) {
      dp1I<-0
    }
    if(hold2) {
      dp2I<-0
    }
    
    res <- c(dp1O, dp1I, dp2O, dp2I)
    list(res)
  })
}


## The parameters 
getceq<-function(clst, mlst=rep(0.1, length(clst))) {
  #calculate analytical equilibrium for Levins model
  #(also works for total community abundance in neutral model)
  
  #get number of species
  nsp<-length(clst)
  #vector for saving results
  ceq<-numeric(nsp)
  
  #follows formula from Tilman 1994 (Ecology 75:2-16)
  for(i in 1:nsp) {
    ceq[i]<-(1-mlst[i]/clst[i])-sum(ceq[0:(i-1)]*(1+clst[0:(i-1)]/clst[i]))
    
    #exclude negative abundances
    if(ceq[i]<0) {
      ceq[i]<-0
    }
  }
  
  return(ceq)
}


## Set up values
times<-seq(0, 150)
spscallst<-10^(seq(log10(0.005), log10(1), length=20))
#spscallst<-c(0.005, 0.01, 0.05, 0.1, 0.5, 0.75, 1)

r0mat1<-r0mat2<-matrix(nrow=length(spscallst), ncol=length(times)-1)
remat1<-remat2<-matrix(nrow=length(spscallst), ncol=length(times)-1)

#ceq
cq<-getceq(clst=c(0.45, 1.05), mlst=c(0.3, 0.3))


## Run simulations
for(i in 1:length(spscallst)) {
  
  #perturbation sp. 1
  parms <- c(c1=0.45, c2=1.05, m1=0.3, m2=0.3, spext=spscallst[i], hold1=FALSE, hold2=FALSE)
  xstart <- c(rep(cq, each=2))*unname(rep(c(1-parms["spext"], parms["spext"]), 2))
  x0<-xstart[2]
  xstart[2]<-xstart[2]*(1-0.2)
  
  names(xstart)<-c("p1O", "p1I", "p2O", "p2I")
  out<-ode(y=xstart, times = times, func = LEVmod, parms = parms)
  
  dx<-abs(x0-out[,3])
  dx[dx==0]<-dx[max(which(dx!=0))]
  remat1[i,]<-cumsum(log(dx[-1]/dx[-length(dx)]))*1/(1:(length(dx)-1))
  
  #perturbation sp. 2
  parms <- c(c1=0.45, c2=1.05, m1=0.3, m2=0.3, spext=spscallst[i], hold1=FALSE, hold2=FALSE)
  xstart <- c(rep(cq, each=2))*unname(rep(c(1-parms["spext"], parms["spext"]), 2))
  x0<-xstart[4]
  xstart[4]<-xstart[4]*(1-0.2)
  
  names(xstart)<-c("p1O", "p1I", "p2O", "p2I")
  out<-ode(y=xstart, times = times, func = LEVmod, parms = parms)
  
  dx<-abs(x0-out[,5])
  dx[dx==0]<-dx[max(which(dx!=0))]
  remat2[i,]<-cumsum(log(dx[-1]/dx[-length(dx)]))*1/(1:(length(dx)-1))
  
  
  
  
  #invasion sp. 1
  parms <- c(c1=0.45, c2=1.05, m1=0.3, m2=0.3, spext=spscallst[i], hold1=TRUE, hold2=FALSE)
  xstart <- c(rep(cq, each=2))*unname(rep(c(1-parms["spext"], parms["spext"]), 2))
  xstart[2]<-0
  
  names(xstart)<-c("p1O", "p1I", "p2O", "p2I")
  out<-ode(y=xstart, times = times, func = LEVmod, parms = parms, method="ode45")
  
  xstart<-out[nrow(out),-1]
  parms["hold1"]<-FALSE
  
  x0<-0
  xstart[2]<-(spscallst[i]-xstart["p2I"])*0.05
  
  out<-ode(y=xstart, times = times, func = LEVmod, parms = parms)
  
  dx<-(out[,3]-x0)
  dx[dx==0]<-dx[max(which(dx!=0))]
  r0mat1[i,]<-cumsum(log(dx[-1]/dx[-length(dx)]))*1/(1:(length(dx)-1))
  
  #invasion sp. 2
  parms <- c(c1=0.45, c2=1.05, m1=0.3, m2=0.3, spext=spscallst[i], hold1=FALSE, hold2=TRUE)
  xstart <- c(rep(cq, each=2))*unname(rep(c(1-parms["spext"], parms["spext"]), 2))
  xstart[4]<-0
  
  names(xstart)<-c("p1O", "p1I", "p2O", "p2I")
  out<-ode(y=xstart, times = times, func = LEVmod, parms = parms)
  
  xstart<-out[nrow(out),-1]
  parms["hold2"]<-FALSE
  
  x0<-0
  xstart[4]<-(spscallst[i]-xstart["p1I"])*0.05
  
  out<-ode(y=xstart, times = times, func = LEVmod, parms = parms, method="ode45")
  
  dx<-(out[,5]-x0)
  dx[dx==0]<-dx[max(which(dx!=0))]
  r0mat2[i,]<-cumsum(log(dx[-1]/dx[-length(dx)]))*1/(1:(length(dx)-1))
}


ofs<-c(0.1, -0.02)
svg("figures/SUP_FIGURE_Levins_nodeom.svg", width=8, height=4)
  m<-c(rep(c(1,2), each=4), 3)
  layout(t(m))
  par(mar=c(2,2,2,2), oma=c(3,3,1,2))
  
  lvl<-seq(-1.6, 0.2, by=0.2)
  #lvl<-c(-1.5, -1, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 1.5)
  collst<-adjustcolor(c(rev(rainbow(sum(lvl<0), start=0.55, end=.70)), rev(rainbow(sum(lvl>0), start=0.099, end=0.1))), alpha.f = 0.6)
  
  filled.contour3(x = log10(times[-1]), y = log10(spscallst), z = t(pmax(remat1, remat2)), levels = lvl, col=collst,
                 xlab="temporal extent, time steps", ylab="spatial extent, fraction of maximum", key.title = title(main = expression(paste(italic(r)[e]))),
                 cex.lab=1.2,
                 plot.axes = { axis(1, at=log10(c(1,2,5,10,20,50,150)), c(1,2,5,10,20,50,150))
                               axis(2, at=log10(c(0.01, 0.05, 0.1, 0.35, 1)), c(0.01, 0.05, 0.1, 0.35, 1)) })
  put.fig.letter("a.", "topleft", offset=ofs, cex=1.8)
  
  lvl2<-sort(-lvl)
  collst2<-adjustcolor(c((rainbow(sum(lvl2<0), start=0.099, end=0.1)),(rainbow(sum(lvl2>0), start=0.55, end=.70))), alpha.f = 0.6)
  filled.contour3(log10(times[-1]), log10(spscallst), t(pmin(r0mat1, r0mat2)), levels = lvl2, col=collst2,
                 xlab="temporal extent, time steps", ylab="spatial extent, fraction of maximum", key.title = title(main = expression(paste(italic(r)[0]))),
                 cex.lab=1.2,
                 plot.axes = { axis(1, at=log10(c(1,2,5,10,20,50,150)), c(1,2,5,10,20,50,150))
                               axis(2, at=log10(c(0.01, 0.05, 0.1, 0.35, 1)), c(0.01, 0.05, 0.1, 0.35, 1)) })
  put.fig.letter("b.", "topleft", offset=ofs, cex=1.8)
  
  filled.legend(z=matrix(lvl), levels=lvl, col=collst, key.axes = axis(4, at = lvl, las=2, labels = sprintf("%.2f", lvl)))
  axis(2, at=lvl, sprintf("%.2f", -lvl), las=2)
  
  mtext(text = expression(paste(r[e])), side = 3, outer = TRUE, line = -1, adj = .9, cex=1.1)
  mtext(text = expression(paste(r[0])), side = 3, outer = TRUE, line = -1, adj = 1.008, cex=1.1, xpd=NA)
  
  mtext(text = expression(paste("temporal extent, time steps")), side = 1, outer = TRUE, line = 1.5, cex=1.2, adj = 0.4)
  mtext(text = expression(paste("spatial extent, fraction of maximum")), side = 2, outer = TRUE, line = 0.3, cex=1.2, adj = 0.5)
dev.off()









