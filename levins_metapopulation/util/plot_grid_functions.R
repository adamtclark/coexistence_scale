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


#arrayout=tmp; xscalslst=log(scalslst,10); xlst=log(tscalelst, 10); splitcol=0; nlevels=10; sqlst = sqtmp; logx=TRUE; logy=TRUE; logz=FALSE; logxps = logxpos; coltype=1; logxps=0; nstart=1; ciplot=FALSE; cimat=0; revcol=FALSE; dops_subset=FALSE; override_tmpsq=TRUE
plot_cont<-function(arrayout, xscalslst, xlst, splitcol=0, nlevels=10, sqlst=0, logx=FALSE, logy=FALSE, logz=FALSE, coltype=1, logxps=0, nstart=1, ciplot=FALSE, cimat=0, revcol=FALSE, dops_subset=TRUE, override_tmpsq=FALSE, ...) {
  
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
    if(override_tmpsq) {
      tmpsq<-sqlst
    } else {
      tmpsq<-c(floor(10^sqlst[1]*1000)/1000,
               round(10^sqlst[-c(1, length(sqlst))],3),
               ceiling(10^sqlst[length(sqlst)]*1000)/1000)
      if(sum(tmpsq<=0)>0) {
        tmpsq<-c(1e-6, 0.0001, tmpsq[tmpsq>0])
      }
      tmpsq<-sort(unique(tmpsq))
    }
    sqlst<-log(tmpsq[tmpsq>0],10)
  }
  
  for(j in 1:length(modlst)) {
    if(coltype==1) {
      collst2<-adjustcolor(c(rev(rainbow(sum(sqlst<dm), start=0.55, end=.70)), rev(rainbow(sum(sqlst>dm), start=0, end=0.1))), alpha.f = 0.6)
    } else if(coltype==2) {
      collst2<-adjustcolor(c(rainbow(sum(sqlst<dm), start=0.15, end=.4), rev(rainbow(sum(sqlst>dm), start=0.75, end=0.85))), alpha.f = 0.6)
    } else if(coltype==3) {
      #collst2<-adjustcolor(c(rainbow(sum(sqlst<0.5), start=0.15, end=.4), rev(rainbow(sum(sqlst>0.5), start=0.75, end=0.85))), alpha.f = 0.6)
      collst2<-adjustcolor(c((rainbow(sum(sqlst<0.5), start=0, end=0.1)), (rainbow(sum(sqlst>0.5), start=0.55, end=.70))), alpha.f = 0.6)
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
    if(revcol) {
      cl2<-rev(collst2)
    } else {
      cl2<-(collst2)
    }
    
    #reduce number of scales
    #ps_subset<-sort(unique(round(exp(seq(log(1), log(length(xlst)), length=100)),0)))
    if(dops_subset) {
      ps_subset<-sort(unique(c(1:26, seq(28, 50, by=2), seq(55, length(xlst), by=5), length(xlst)-1, length(xlst))))
    } else {
      ps_subset<-1:length(xlst)
    }
    
    filled.contour3(x = xlst[ps_subset], 
                    y = xscalslst, 
                    z = tmpz[ps_subset,], levels = sqlst, col=cl2,axes=F,
                    xlim=range(xlst[tmpps]))
    put.fig.letter(paste(letters[nstart], ".", sep=""), "topleft", offset=ofs1, cex=1.6)
    nstart<-nstart+1
    
    if(logz) {
      contour(x = xlst[ps_subset], 
              y = xscalslst, 
              z = tmpz[ps_subset,],
              levels = sqlst,
              labels=10^sqlst,
              add=TRUE,axes=F,
              xlim=range(xlst[tmpps]))
    } else {
      contour(x = xlst[ps_subset], 
              y = xscalslst, 
              z = tmpz[ps_subset,],
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




#for empirical data
#arrayout=array_quant; xscalslst=log10(xscl); xlst=log10(timebands); nlevels=10; logx=TRUE; logy=TRUE; logz=FALSE; nstart=1; ofs1=c(0, 0); sqlst=squse
plot_cont_emp<-function(arrayout, xscalslst, xlst, nlevels=10, logx=FALSE, logy=FALSE, logz=FALSE, coltype=1, nstart=1, sqlst=0, ofs1=c(0, 0), ...) {
  rng<-range(arrayout[,,], na.rm=T)
  rng_rnd<-c(floor(rng[1]*10)/10,
               ceiling(rng[2]*10)/10)
    
  if(sum(abs(sqlst))==0) {
    sqlst<-pretty(rng_rnd, nlevels)
  }
  
  if(sum(sqlst==0)==0) {
    sqlst<-c(sqlst[sqlst<0], 0, sqlst[sqlst>0])
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
  
  for(j in 1:dim(arrayout)[3]) {
    dm<-0
    collst2<-adjustcolor(c((rainbow(sum(sqlst<0), start=0, end=0.1)), (rainbow(sum(sqlst>=0), start=0.55, end=.70))), alpha.f = 0.6)
    
    tmpps<-colSums(abs(arrayout[,,j]), na.rm=T)!=0

    filled.contour3(x = xlst, 
                    y = xscalslst, 
                    z = arrayout[,,j], levels = sqlst, col=collst2,axes=F,
                    xlim=range(xlst[tmpps], na.rm=T), ...)
    put.fig.letter(paste(letters[nstart], ".", sep=""), "topleft", offset=ofs1, cex=1.6)
    nstart<-nstart+1
    
    if(logz) {
      contour(x = xlst, 
              y = xscalslst, 
              z = arrayout[,,j],
              levels = sqlst,
              labels="", method="edge",#10^sqlst,
              add=TRUE,axes=F,
              xlim=range(xlst[tmpps],na.rm=T))
    } else {
      contour(x = xlst, 
              y = xscalslst, 
              z = arrayout[,,j],
              levels = sqlst,
              labels="", method="edge",#round(sqlst,2),
              add=TRUE,axes=F,
              xlim=range(xlst[tmpps],na.rm=T))
    }
    
    
    if(logx) {
      tmp<-c(1,2,5,10,20,50,80)
      
      axis(1, at=log(tmp,10), labels = tmp, las=2)
    } else {
      axis(1, las=2)
    }
    
    if(logy) {
      tmp<-c(1, 5, 10, 25, 50, 100, 200, 400)
      
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


