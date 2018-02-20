#from https://waterprogramming.wordpress.com/2015/12/
put.fig.letter <- function(label, location="topleft", x=NULL, y=NULL, offset=c(0, 0), ...) {
  if(length(label) > 1) {
    warning("length(label) > 1, using label[1]")
  }
  if(is.null(x) | is.null(y)) {
    coords <- switch(location,
                     topleft = c(0.015,0.98),
                     topcenter = c(0.5525,0.98),
                     topright = c(0.985, 0.98),
                     bottomleft = c(0.015, 0.02), 
                     bottomcenter = c(0.5525, 0.02), 
                     bottomright = c(0.985, 0.02),
                     c(0.015, 0.98) )
  } else {
    coords <- c(x,y)
  }
  this.x <- grconvertX(coords[1] + offset[1], from="nfc", to="user")
  this.y <- grconvertY(coords[2] + offset[2], from="nfc", to="user")
  text(labels=label[1], x=this.x, y=this.y, xpd=NA, ...)
}



matbarplot<-function(heightdat, sedat, namelst, collst, logged=FALSE, yup=NULL, pty=NULL) {
  #heightdat is bar heights
  #sedat is standard errors
  #namelst is names for bars
  #collst is colors for bars
  #logged tells whether heights and SE's are log-transformed (will plot in linear space)
  #yup is a scaling factor for the y-axis (in fold-change relative to mean +2SE)
  
  
  if(is.null(yup)) {
    yup<-1.05
  }
  
  yxlm<-range(c(heightdat-2*sedat, heightdat+2*sedat))
  
  if(logged) {
    if(is.null(pty)) {
      pty<-10^pretty(seq((yxlm[1]), (yxlm[2]), length=4), n = 4)
      pty<-round(pty/10)*10
      pty<-pty[pty>0]
    }
    
    hqtl<-rbind(10^(heightdat+sedat*2),
                10^(heightdat+sedat),
                10^(heightdat),
                10^(heightdat-sedat),
                10^(heightdat-sedat*2))
    
    plot(c(0, 1), c(1, 1), xlab="", ylab="", axes=F, type="n", ylim=c(min(pty), max(pty)*yup), xaxs="i", log="y")
  } else {
    pty<-pretty((seq(yxlm[1], yxlm[2], length=4)), n = 4)
    
    hqtl<-rbind((heightdat+sedat*2),
                (heightdat+sedat),
                (heightdat),
                (heightdat-sedat),
                (heightdat-sedat*2))
    
    plot(c(0, 1), c(1, 1), xlab="", ylab="", axes=F, type="n", ylim=c(min(pty), max(pty)*yup), xaxs="i")
  }
  axis(2, pty, pty, las=2)
  
  xps<-seq(0, 1, length=length(heightdat)*3+2)
  dx<-min(diff(xps))
  
  xps<-xps[-c(1, 2, length(xps)-1, length(xps))]
  xps<-xps[seq(1, length(xps), by=3)]
  
  
  axis(1, at=xps, labels=FALSE)
  if(!logged) {
    text(x=xps, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
       labels=namelst, srt=45, xpd=TRUE, cex=1.2)
  } else {
    text(x=xps, y=exp(log(pty[1])-diff(log(range(pty)))*0.16),
         labels=namelst, srt=45, xpd=TRUE, cex=1.2)
  }
  
  
  for(i in 1:length(heightdat)) {
    polygon(c(xps[i]-dx, xps[i]-dx, xps[i]+dx, xps[i]+dx),
            c(min(pty)*0.01, hqtl[3,i], hqtl[3,i], min(pty)*0.01), col=collst[i])
  }
  
  segments(xps, hqtl[4,], xps, hqtl[2,], lwd=3, lend=2)
  segments(xps, hqtl[5,], xps, hqtl[1,], lwd=1, lend=2)
  
  
  segments(xps-dx/9, hqtl[5,], xps+dx/9, hqtl[5,], lwd=1, lend=2)
  segments(xps-dx/9, hqtl[1,], xps+dx/9, hqtl[1,], lwd=1, lend=2)
  
  
  box()
  return(list(xps=xps, dx=dx, hqtl=hqtl, ylm=c(min(pty), max(pty))))
}


matbarplot_group<-function(heightdat, sedat, namelst, collst, yup=NULL, dopoly=TRUE, dopvalchecks=TRUE, dodiffchecks=TRUE, niter=1000, pty=NULL, yxlm=NULL, xsrt=45, xadj=0.85, xcex=1.2, posttr=NULL, labside=NULL, overwritechar=NULL, doyax=TRUE,  doxax=TRUE) {
  #heightdat is bar heights
  #sedat is standard errors
  #namelst is names for bars
  #collst is colors for bars
  #logged tells whether heights and SE's are log-transformed (will plot in linear space)
  #yup is a scaling factor for the y-axis (in fold-change relative to mean +2SE)
  
  
  if(is.null(yup)) {
    yup<-1
  }
  
  if(!is.null(yxlm[1])==0) {
    yxlm<-range(c(heightdat-2*sedat*yup, heightdat+2*sedat*yup), na.rm=T)
  }
  
  if(!is.null(pty[1])==0) {
    pty<-pretty((seq(yxlm[1], yxlm[2], length=4)), n = 4)
  }
  
  hqtl<-rbind((heightdat+sedat*2),
              (heightdat+sedat),
              (heightdat),
              (heightdat-sedat),
              (heightdat-sedat*2))
  if(!is.null(posttr)) {
    hqtl<-posttr(hqtl)
  }
  
  plot(c(0, 1), c(1, 1), xlab="", ylab="", axes=F, type="n", ylim=c(min(yxlm), max(yxlm)*yup), xaxs="i", yaxs="i")
  abline(h=0, lty=3)
  if(doyax) {
    axis(2, pty, pty, las=2)
  }
  
  xps<-seq(0, 1, length=nrow(heightdat)*3+2)
  dx<-min(diff(xps))
  sdx<-dx/ncol(heightdat)
  
  xps<-xps[-c(1, 2, length(xps)-1, length(xps))]
  xps<-xps[seq(1, length(xps), by=3)]
  
  if(doxax) {
    axis(1, at=xps, labels=FALSE)
    text(x=xps, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
        labels=namelst, srt=xsrt, adj=xadj, xpd=NA, cex=xcex)
  }
  
  dpsq<-seq(-dx, dx, length=ncol(heightdat))
  if(dopoly) {
    for(j in 1:ncol(heightdat)) {
      for(i in 1:nrow(heightdat)) {
        ps3<-1:nrow(heightdat)+nrow(heightdat)*(3-1)
        
        dpos<-dpsq[j]
        polygon(c(xps[i]-sdx+dpos, xps[i]-sdx+dpos, xps[i]+sdx+dpos, xps[i]+sdx+dpos),
                c(0, hqtl[ps3,j][i], hqtl[ps3,j][i], 0), col=collst[i,j])
      }
    }
  }
  tmp<-xps+mean(diff(xps))/2
  abline(v=tmp[-length(tmp)], lty=1, col="grey")
  
  for(j in 1:ncol(heightdat)) {
    ps1<-1:nrow(heightdat)+nrow(heightdat)*(1-1)
    ps2<-1:nrow(heightdat)+nrow(heightdat)*(2-1)
    ps3<-1:nrow(heightdat)+nrow(heightdat)*(3-1)
    ps4<-1:nrow(heightdat)+nrow(heightdat)*(4-1)
    ps5<-1:nrow(heightdat)+nrow(heightdat)*(5-1)
    
    dpos<-dpsq[j]
    
    segments(xps+dpos, hqtl[ps4,j], xps+dpos, hqtl[ps2,j], lwd=3, lend=2)
    segments(xps+dpos, hqtl[ps5,j], xps+dpos, hqtl[ps1,j], lwd=1, lend=2)
    
    
    segments(xps-sdx/9+dpos, hqtl[ps5,j], xps+sdx/9+dpos, hqtl[ps5,j], lwd=1, lend=2)
    segments(xps-sdx/9+dpos, hqtl[ps1,j], xps+sdx/9+dpos, hqtl[ps1,j], lwd=1, lend=2)
  }
  
  if(dopvalchecks|dodiffchecks) {
    simoutdat<-array(dim=c(nrow(heightdat), ncol(heightdat), niter))
    for(i in 1:nrow(heightdat)) {
      for(j in 1:ncol(heightdat)) {
        simoutdat[i,j,]<-suppressWarnings(try(rnorm(niter, heightdat[i,j], sedat[i,j]), silent=TRUE))
      }
    }
    
  }
  
  if(dopvalchecks) {
    pvsym<-c("", "*", "**", "***")
    pvmtch<-c(1, 0.05, 0.01, 0.001)
    dy<-abs(diff(range(pty)))*0.025
    
    for(j in 1:ncol(heightdat)) {
      dpos<-dpsq[j]
      for(i in 1:nrow(heightdat)) {
        if(sum(is.finite(simoutdat[i,j,]), na.rm=T)>0) {
          pobs<-sum((simoutdat[i,j,]<0)/niter, na.rm=T)
          if(heightdat[i,j]<0) {
            pobs<-1-pobs
          } 
          ps<-pvsym[sum(pobs<pvmtch)]
          
          if(heightdat[i,j]>0) {
            ps1<-1:nrow(heightdat)+nrow(heightdat)*(1-1)
            if(is.null(labside) || labside=="up") {
              text(xps[i]+dpos, hqtl[ps1,j][i]+dy, ps, cex=0.75)
            } else {
              text(xps[i]+dpos, hqtl[ps5,j][i]-dy, ps, cex=0.75)
            }
          } else {
            ps5<-1:nrow(heightdat)+nrow(heightdat)*(5-1)
            if(is.null(labside) || labside=="down") {
              text(xps[i]+dpos, hqtl[ps5,j][i]-dy, ps, cex=0.75)
            } else {
              text(xps[i]+dpos, hqtl[ps1,j][i]+dy, ps, cex=0.75)
            }
          }
        }
      }
    }
  }
  
  if(dodiffchecks) {
    if(!dopvalchecks) {
      dy<-abs(diff(range(pty)))*0.025*1.2
    } else {
      dy<-abs(diff(range(pty)))*0.025*2.5
    }
    
    for(i in 1:nrow(heightdat)) {
      isfh<-which(is.finite(heightdat[i,]))
      compmat<-matrix(ncol=ncol(heightdat), nrow=ncol(heightdat))
      
      for(j in 1:(ncol(heightdat)-1)) {
        for(k in (j+1):ncol(heightdat)) {
          pest<-sum(simoutdat[i,j,]>simoutdat[i,k,], na.rm=T)/niter
          pest<-min(c(pest, 1-pest))
          
          compmat[j,k]<-pest
        }
      }
      
      lp<-1
      glabs<-character(ncol(heightdat))
      if(sum(compmat[isfh,isfh]<=0.05, na.rm=T)==sum(is.finite(compmat[isfh,isfh]))) {
        glabs[isfh]<-letters[1:length(isfh)]
      } else if(sum(compmat[isfh,isfh]>0.05, na.rm=T)==sum(is.finite(compmat[isfh,isfh]))) {
        glabs[isfh]<-letters[1]
      } else {
        for(j in 1:ncol(heightdat)) {
          if(j<ncol(heightdat)) {
            tps<-c(rep(FALSE, j-1), TRUE, compmat[j,(j+1):ncol(heightdat)]>0.05)
          } else {
            tps<-c(rep(FALSE, j-1), TRUE)
          }
          
          if(sum(tps)>1) {
            glabs[tps]<-paste(glabs[tps], letters[lp], sep="")
            lp<-lp+1
          } else {
            if(glabs[tps]=="") {
              glabs[tps]<-letters[lp]
              lp<-lp+1
            }
          }
        }
      }
      
      if(!is.null(overwritechar[1])) {
        glabs<-overwritechar[i,]
      }
      
      for(j in isfh) {
        dpos<-dpsq[j]
        if(heightdat[i,j]>0) {
          ps1<-1:nrow(heightdat)+nrow(heightdat)*(1-1)
          
          if(is.null(labside) || labside=="up") {
            text(xps[i]+dpos, hqtl[ps1,j][i]+dy, glabs[j], cex=0.9)
          } else {
            text(xps[i]+dpos, hqtl[ps5,j][i]-dy, glabs[j], cex=0.9)
          }
        } else {
          ps5<-1:nrow(heightdat)+nrow(heightdat)*(5-1)
          if(is.null(labside) || labside=="down") {
            text(xps[i]+dpos, hqtl[ps5,j][i]-dy, glabs[j], cex=0.9)
          } else {
            text(xps[i]+dpos, hqtl[ps1,j][i]+dy, glabs[j], cex=0.9)
          }
        } 
      }
    }
  }
  
  box()
  return(list(xps=xps, dx=dx, hqtl=hqtl, ylm=c(min(pty), max(pty))))
}



shadowtext <- function(x, y=NULL, labels, col='white', bg='black', 
                       theta= seq(0, 2*pi, length.out=50), r=0.1, ... ) {
  
  xy <- xy.coords(x,y)
  xo <- r*strwidth('A')
  yo <- r*strheight('A')
  
  # draw background text with small shift in x and y in background colour
  for (i in theta) {
    text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
  }
  # draw actual text in exact xy position in foreground colour
  text(xy$x, xy$y, labels, col=col, ... )
}
