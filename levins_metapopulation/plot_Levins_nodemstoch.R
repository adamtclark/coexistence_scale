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

RPSmod <- function(t, x, parms, input)  {
  with(as.list(c(parms, x)), {
    p1_col<-c1*(p1O+p1I)
    p2_col<-c2*(p2O+p2I)
    p3_col<-c3*(p3O+p3I)
    p4_col<-c4*(p4O+p4I)
    
    dp1O <- p1_col*((1-spext)-p1O-p4O)-p3_col*p1O-p4_col*p1O-m1*p1O  #outside of focal patch
    dp1I <- p1_col*(spext-p1I-p4I)-p3_col*p1I-p4_col*p1I-m1*p1I  #inside of focal patch
    
    dp2O <- p2_col*((1-spext)-p1O-p2O)-p4_col*p2O-p1_col*p2O-m2*p2O  #outside of focal patch
    dp2I <- p2_col*(spext-p1I-p2I)-p4_col*p2I-p1_col*p2I-m2*p2I  #inside of focal patch
    
    dp3O <- p3_col*((1-spext)-p3O-p2O)-p1_col*p3O-p2_col*p3O-m3*p3O  #outside of focal patch
    dp3I <- p3_col*(spext-p3I-p2I)-p1_col*p3I-p2_col*p3I-m3*p3I  #inside of focal patch
    
    dp4O <- p4_col*((1-spext)-p4O-p3O)-p2_col*p4O-p3_col*p4O-m4*p4O  #outside of focal patch
    dp4I <- p4_col*(spext-p4I-p3I)-p2_col*p4I-p3_col*p4I-m4*p4I  #inside of focal patch
    
    if(hold1) {
      dp1I<-0
    }
    if(hold2) {
      dp2I<-0
    }
    if(hold3) {
      dp3I<-0
    }
    if(hold4) {
      dp4I<-0
    }
    
    res <- c(dp1O, dp1I, dp2O, dp2I, dp3O, dp3I, dp4O, dp4I)
    list(res)
  })
}




NUTmod <- function(t, x, parms, input)  {
  with(as.list(c(parms, x)), {
    p1_col<-c1*(p1O+p1I)
    
    dp1O <- p1_col*((1-spext)-p1O-p2O)-m1*p1O                        #outside of focal patch
    dp1I <- p1_col*(spext-p1I-p2I)-m1*p1I                        #inside of focal patch
    
    p2_col<-c2*(p2O+p2I)
    dp2O <- p2_col*((1-spext)-p1O-p2O)-m2*p2O  #outside of focal patch
    dp2I <- p2_col*(spext-p1I-p2I)-m2*p2I  #inside of focal patch
    
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


#############################################
#1. LEVINS:
#############################################
## Set up values
times<-seq(0, 200)
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


#############################################
#2. Disturbance
#############################################
## Set up values
times<-seq(0, 200)
spscallst<-10^(seq(log10(0.005), log10(1), length=20))
#spscallst<-c(0.005, 0.01, 0.05, 0.1, 0.5, 0.75, 1)

r0mat1_D<-r0mat2_D<-matrix(nrow=length(spscallst), ncol=length(times)-1)
remat1_D<-remat2_D<-matrix(nrow=length(spscallst), ncol=length(times)-1)

#ceq
distlst<-c(0.95, 0)
prtfrq<-50
cq<-getceq(clst=c(0.435, 0.600), mlst=c(0.3 - log(1 - distlst[1])/prtfrq, 0.3))


distode<-function(times) {
  nr<-round((length(times)-1)/prtfrq)
  out<-NULL
  for(j in 1:nr) {
    outtmp<-ode(y=xstart, times = 0:prtfrq, func = LEVmod, parms = parms)
    out<-rbind(out, outtmp)
    xstart<-outtmp[nrow(outtmp),-1]
    xstart[c(1,2)]<-xstart[c(1,2)]*(1-distlst[1])
    names(xstart)<-c("p1O", "p1I", "p2O", "p2I")
  }
  out[,1]<-(0:(nrow(out)-1))
  out
}


## Run simulations
for(i in 1:length(spscallst)) {
  
  #perturbation sp. 1
  parms <- c(c1=0.435, c2=0.600, m1=0.3, m2=0.3, spext=spscallst[i], hold1=FALSE, hold2=FALSE)
  xstart <- c(3.019125e-01, 1.517148e-03, 6.145120e-02, 3.088000e-04)
  names(xstart)<-c("p1O", "p1I", "p2O", "p2I")
  xstart[1:2]<-xstart[1:2]*(1-distlst[1])
  out<-distode(times)
  xstart<-out[nrow(out),-1]
  xstart[1:2]<-xstart[1:2]*(1-distlst[1])
  
  x0<-distode(times)
  xstart[2]<-xstart[2]*(1-0.2)
  
  out<-distode(times)
  
  dx<-abs(x0[,3]-out[,3])
  remat1_D[i,]<-(log(dx[-1]/dx[1])/(1:(length(dx)-1)))[1:max(times)]
  
  #perturbation sp. 2
  parms <- c(c1=0.435, c2=0.600, m1=0.3, m2=0.3, spext=spscallst[i], hold1=FALSE, hold2=FALSE)
  xstart <- c(3.019125e-01, 1.517148e-03, 6.145120e-02, 3.088000e-04)
  names(xstart)<-c("p1O", "p1I", "p2O", "p2I")
  xstart[1:2]<-xstart[1:2]*(1-distlst[1])
  out<-distode(times)
  xstart<-out[nrow(out),-1]
  xstart[1:2]<-xstart[1:2]*(1-distlst[1])
  
  x0<-distode(times)
  xstart[4]<-xstart[4]*(1-0.2)
  
  out<-distode(times)
  
  dx<-abs(x0[,5]-out[,5])
  remat2_D[i,]<-(log(dx[-1]/dx[1])/(1:(length(dx)-1)))[1:max(times)]
  
  
  
  
  #invasion sp. 1
  parms <- c(c1=0.435, c2=0.600, m1=0.3, m2=0.3, spext=spscallst[i], hold1=TRUE, hold2=FALSE)
  xstart <- c(3.019125e-01, 1.517148e-03, 6.145120e-02, 3.088000e-04)
  names(xstart)<-c("p1O", "p1I", "p2O", "p2I")
  xstart[1:2]<-xstart[1:2]*(1-distlst[1])
  xstart[2]<-0
  
  names(xstart)<-c("p1O", "p1I", "p2O", "p2I")
  out<-distode(times)
  
  xstart<-out[nrow(out),-1]
  parms["hold1"]<-FALSE
  
  xstart[1:2]<-xstart[1:2]*(1-distlst[1])
  xstart[2]<-(spscallst[i]-xstart["p2I"])*0.05
  x0<-xstart[2]
  
  
  out<-distode(times)
  
  dx<-out[,3]
  dx<-dx[-1]
  r0mat1_D[i,]<-(log(dx/x0)/(1:(length(dx))))[1:max(times)]
  
  
  #invasion sp. 2
  parms <- c(c1=0.435, c2=0.600, m1=0.3, m2=0.3, spext=spscallst[i], hold1=FALSE, hold2=TRUE)
  xstart <- c(3.019125e-01, 1.517148e-03, 6.145120e-02, 3.088000e-04)
  names(xstart)<-c("p1O", "p1I", "p2O", "p2I")
  xstart[1:2]<-xstart[1:2]*(1-distlst[1])
  xstart[4]<-0
  
  names(xstart)<-c("p1O", "p1I", "p2O", "p2I")
  out<-distode(times)
  
  xstart<-out[nrow(out),-1]
  parms["hold2"]<-FALSE
  
  xstart[1:2]<-xstart[1:2]*(1-distlst[1])
  xstart[4]<-(spscallst[i]-xstart["p1I"])*0.05
  x0<-xstart[4]
  
  
  out<-distode(times)
  
  dx<-out[,5]
  dx<-dx[-1]
  r0mat2_D[i,]<-(log(dx/x0)/(1:(length(dx))))[1:max(times)]
  
}






#############################################
#3. RPS:
#############################################
## Set up values
times<-seq(0, 200)
spscallst<-10^(seq(log10(0.005), log10(1), length=20))
#spscallst<-c(0.005, 0.01, 0.05, 0.1, 0.5, 0.75, 1)

r0mat1_RPS<-matrix(nrow=length(spscallst), ncol=length(times)-1)
remat1_RPS<-matrix(nrow=length(spscallst), ncol=length(times)-1)

#ceq
#rnd<-c(1.2, 1.3, 0.8, 0.7)
rnd<-rep(1, 4)
cq<-getceq(clst=c(0.48), mlst=c(0.224))*rnd/(sum(rnd))


## Run simulations
for(i in 1:length(spscallst)) {
  #perturbation sp. 1
  parms <- c(c1=0.48, c2=0.48, c3=0.48, c4=0.48,
             m1=0.224,m2=0.224,m3=0.224,m4=0.224,
             spext=spscallst[i],
             hold1=FALSE, hold2=FALSE, hold3=FALSE, hold4=FALSE)
  xstart <- c(rep(cq, each=2))*unname(rep(c(1-parms["spext"], parms["spext"]), 4))
  names(xstart)<-c("p1O", "p1I", "p2O", "p2I", "p3O", "p3I", "p4O", "p4I")
  
  x0<-ode(y=xstart, times = times, func = RPSmod, parms = parms)
  #matplot(x0[,1], x0[,c(2,4,6,8)+1], type="l")
  
  xstart[2]<-xstart[2]*(1-0.2)
  out<-ode(y=xstart, times = times, func = RPSmod, parms = parms)
  
  dx<-abs(x0[,3]-out[,3])
  dx[dx==0]<-dx[max(which(dx!=0))]
  remat1_RPS[i,]<-cumsum(log(dx[-1]/dx[-length(dx)]))*1/(1:(length(dx)-1))
  
  
  #invasion sp. 1
  xstart <- c(rep(cq, each=2))*unname(rep(c(1-parms["spext"], parms["spext"]), 4))
  xstart[2]<-0
  parms["hold1"]<-TRUE
  
  names(xstart)<-c("p1O", "p1I", "p2O", "p2I", "p3O", "p3I", "p4O", "p4I")
  out<-ode(y=xstart, times = 1:2000, func = RPSmod, parms = parms, method="ode45")
  #matplot(out[,1], out[,-1], type="l")
  
  xstart<-out[nrow(out),-1]
  parms["hold1"]<-FALSE
  
  x0<-0
  xstart[2]<-(spscallst[i]-xstart["p2I"]-xstart["p3I"]-xstart["p4I"])*0.05
  
  out<-ode(y=xstart, times = times, func = RPSmod, parms = parms)
  
  dx<-(out[,3]-x0)
  dx[dx==0]<-dx[max(which(dx!=0))]
  r0mat1_RPS[i,]<-cumsum(log(dx[-1]/dx[-length(dx)]))*1/(1:(length(dx)-1))
}




#############################################
#4. NEUTRAL:
#############################################
## Set up values
times<-seq(0, 200)
spscallst<-10^(seq(log10(0.005), log10(1), length=20))

r0mat1_NUT<-matrix(nrow=length(spscallst), ncol=length(times)-1)
remat1_NUT<-matrix(nrow=length(spscallst), ncol=length(times)-1)

#ceq
cq<-getceq(clst=c(1.5), mlst=c(0.7))


## Run simulations
for(i in 1:length(spscallst)) {
  
  #perturbation sp. 1
  parms <- c(c1=1.5, c2=1.5, m1=0.7, m2=0.7, spext=spscallst[i], hold1=FALSE, hold2=FALSE)
  xstart <- c(rep(cq, each=2))/2*unname(rep(c(1-parms["spext"], parms["spext"]), 2))
  x0<-xstart[2]
  xstart[2]<-xstart[2]*(1-0.2)
  
  names(xstart)<-c("p1O", "p1I", "p2O", "p2I")
  out<-ode(y=xstart, times = times, func = NUTmod, parms = parms)
  #matplot(out[,1], out[,-1], type="l")
  
  dx<-abs(x0-out[,3])
  dx[dx==0]<-dx[max(which(dx!=0))]
  remat1_NUT[i,]<-cumsum(log(dx[-1]/dx[-length(dx)]))*1/(1:(length(dx)-1))
  
  
  #invasion sp. 1
  parms["hold1"]<-TRUE
  xstart <- c(rep(cq, each=2))/2*unname(rep(c(1-parms["spext"], parms["spext"]), 2))
  xstart[2]<-0
  
  names(xstart)<-c("p1O", "p1I", "p2O", "p2I")
  out<-ode(y=xstart, times = times, func = NUTmod, parms = parms, method="ode45")
  
  xstart<-out[nrow(out),-1]
  parms["hold1"]<-FALSE
  
  x0<-0
  xstart[2]<-(spscallst[i]-xstart["p2I"])*0.05
  
  out<-ode(y=xstart, times = times, func = NUTmod, parms = parms)
  
  dx<-(out[,3]-x0)
  dx[dx==0]<-dx[max(which(dx!=0))]
  r0mat1_NUT[i,]<-cumsum(log(dx[-1]/dx[-length(dx)]))*1/(1:(length(dx)-1))
}




#############################################
#PLOT
#############################################
ofs<-c(0.15, -0.05)
svg("figures/SUP_FIGURE_Levins_nodeom.svg", width=6, height=8)
  m<-rbind(c(rep(c(1,2), each=4), 9,9),
           c(rep(c(3,4), each=4), 9,9),
           c(rep(c(5,6), each=4), 9,9),
           c(rep(c(7,8), each=4), 9,9))
  layout((m))
  par(mar=c(2,2,2,2), oma=c(3,3,1,2))
  
  lvl<-seq(-1.6, 1.6, by=0.2)
  collst<-adjustcolor(c(rev(rainbow(sum(lvl<0), start=0.55, end=.70)), rev(rainbow(sum(lvl>0), start=0.0, end=0.1))), alpha.f = 0.6)
  lvl2<-sort(-lvl)
  collst2<-adjustcolor(c((rainbow(sum(lvl2<0), start=0.0, end=0.1)),(rainbow(sum(lvl2>0), start=0.55, end=.70))), alpha.f = 0.6)
  
  
  #Levins
  filled.contour3(x = log10(times[-1]), y = log10(spscallst), z = t(pmax(remat1, remat2)), levels = lvl, col=collst,
                 xlab="temporal extent, time steps", ylab="spatial extent, fraction of maximum", key.title = title(main = expression(paste(italic(r)[e]))),
                 cex.lab=1.2,
                 plot.axes = { axis(1, at=log10(c(1,2,5,10,20,50,150)), c(1,2,5,10,20,50,150), las=2)
                               axis(2, at=log10(c(0.01, 0.05, 0.1, 0.35, 1)), c(0.01, 0.05, 0.1, 0.35, 1), las=2) })
  put.fig.letter("a.", "topleft", offset=ofs, cex=1.8)
  
  filled.contour3(log10(times[-1]), log10(spscallst), t(pmin(r0mat1, r0mat2)), levels = lvl2, col=collst2,
                 xlab="temporal extent, time steps", ylab="spatial extent, fraction of maximum", key.title = title(main = expression(paste(italic(r)[0]))),
                 cex.lab=1.2,
                 plot.axes = { axis(1, at=log10(c(1,2,5,10,20,50,150)), c(1,2,5,10,20,50,150), las=2)
                               axis(2, at=log10(c(0.01, 0.05, 0.1, 0.35, 1)), c(0.01, 0.05, 0.1, 0.35, 1), las=2) })
  put.fig.letter("b.", "topleft", offset=ofs, cex=1.8)
  
  #Disturbance
  filled.contour3(x = log10(times[-1]), y = log10(spscallst), z = t(pmax(remat1_D, remat2_D)), levels = lvl, col=collst,
                  xlab="temporal extent, time steps", ylab="spatial extent, fraction of maximum", key.title = title(main = expression(paste(italic(r)[e]))),
                  cex.lab=1.2,
                  plot.axes = { axis(1, at=log10(c(1,2,5,10,20,50,150)), c(1,2,5,10,20,50,150), las=2)
                    axis(2, at=log10(c(0.01, 0.05, 0.1, 0.35, 1)), c(0.01, 0.05, 0.1, 0.35, 1), las=2) })
  put.fig.letter("c.", "topleft", offset=ofs, cex=1.8)
  
  filled.contour3(log10(times[-1]), log10(spscallst), t(pmin(r0mat1_D, r0mat2_D)), levels = lvl2, col=collst2,
                  xlab="temporal extent, time steps", ylab="spatial extent, fraction of maximum", key.title = title(main = expression(paste(italic(r)[0]))),
                  cex.lab=1.2,
                  plot.axes = { axis(1, at=log10(c(1,2,5,10,20,50,150)), c(1,2,5,10,20,50,150), las=2)
                    axis(2, at=log10(c(0.01, 0.05, 0.1, 0.35, 1)), c(0.01, 0.05, 0.1, 0.35, 1), las=2) })
  put.fig.letter("d.", "topleft", offset=ofs, cex=1.8)
  
  
  
  #RPS
  filled.contour3(x = log10(times[-1]), y = log10(spscallst), z = t((remat1_RPS)), levels = lvl, col=collst,
                  xlab="temporal extent, time steps", ylab="spatial extent, fraction of maximum", key.title = title(main = expression(paste(italic(r)[e]))),
                  cex.lab=1.2,
                  plot.axes = { axis(1, at=log10(c(1,2,5,10,20,50,150)), c(1,2,5,10,20,50,150), las=2)
                    axis(2, at=log10(c(0.01, 0.05, 0.1, 0.35, 1)), c(0.01, 0.05, 0.1, 0.35, 1), las=2) })
  put.fig.letter("e.", "topleft", offset=ofs, cex=1.8)
  
  filled.contour3(log10(times[-1]), log10(spscallst), t((r0mat1_RPS)), levels = lvl2, col=collst2,
                  xlab="temporal extent, time steps", ylab="spatial extent, fraction of maximum", key.title = title(main = expression(paste(italic(r)[0]))),
                  cex.lab=1.2,
                  plot.axes = { axis(1, at=log10(c(1,2,5,10,20,50,150)), c(1,2,5,10,20,50,150), las=2)
                    axis(2, at=log10(c(0.01, 0.05, 0.1, 0.35, 1)), c(0.01, 0.05, 0.1, 0.35, 1), las=2) })
  put.fig.letter("f.", "topleft", offset=ofs, cex=1.8)
  
  
  #Neutral
  filled.contour3(x = log10(times[-1]), y = log10(spscallst), z = t((remat1_NUT)), levels = lvl, col=collst,
                  xlab="temporal extent, time steps", ylab="spatial extent, fraction of maximum", key.title = title(main = expression(paste(italic(r)[e]))),
                  cex.lab=1.2,
                  plot.axes = { axis(1, at=log10(c(1,2,5,10,20,50,150)), c(1,2,5,10,20,50,150), las=2)
                    axis(2, at=log10(c(0.01, 0.05, 0.1, 0.35, 1)), c(0.01, 0.05, 0.1, 0.35, 1), las=2) })
  put.fig.letter("g.", "topleft", offset=ofs, cex=1.8)
  
  filled.contour3(log10(times[-1]), log10(spscallst), t((r0mat1_NUT)), levels = lvl2, col=collst2,
                  xlab="temporal extent, time steps", ylab="spatial extent, fraction of maximum", key.title = title(main = expression(paste(italic(r)[0]))),
                  cex.lab=1.2,
                  plot.axes = { axis(1, at=log10(c(1,2,5,10,20,50,150)), c(1,2,5,10,20,50,150), las=2)
                    axis(2, at=log10(c(0.01, 0.05, 0.1, 0.35, 1)), c(0.01, 0.05, 0.1, 0.35, 1), las=2) })
  put.fig.letter("h.", "topleft", offset=ofs, cex=1.8)
  
  
  
  
  
  
  par(mar=c(2,4,2,2))
  
  filled.legend(z=matrix(lvl), levels=lvl, col=collst, key.axes = axis(4, at = lvl, las=2, labels = sprintf("%.2f", lvl)))
  axis(2, at=lvl, sprintf("%.2f", -lvl), las=2)
  
  mtext(text = expression(paste(italic(r[e]))), side = 3, outer = TRUE, line = -1, adj = .86, cex=1.1)
  mtext(text = expression(paste(italic(r[0]))), side = 3, outer = TRUE, line = -1, adj = 1.008, cex=1.1, xpd=NA)
  
  mtext(text = expression(paste("temporal extent, time steps")), side = 1, outer = TRUE, line = 1.5, cex=1.2, adj = 0.3)
  mtext(text = expression(paste("spatial extent, fraction of maximum")), side = 2, outer = TRUE, line = 0.9, cex=1.2, adj = 0.5)
  
  mtext(text = "levis", side = 2, line = 4, cex=1.2, adj = 0.945)
  mtext(text = "disturbance", side = 2, line = 4, cex=1.2, adj = 0.648)
  mtext(text = "RPS", side = 2, line = 4, cex=1.2, adj = 0.352)
  mtext(text = "neutral", side = 2, line = 4, cex=1.2, adj = 0.055)
dev.off()
  
  
