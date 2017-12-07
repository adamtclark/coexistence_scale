error
rm(list=ls())
setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

require(data.table)

matout_tot<-data.frame(fread("output/matout_tot.csv"))
nsp<-(ncol(matout_tot)-3)/10

colnames(matout_tot)<-c("scale", "tscale", "iter",
                        paste("eig_meta1", 1:nsp, sep="."),
                        paste("eig_meta2", 1:nsp, sep="."),
                        paste("r0_meta", 1:nsp, sep="."),
                        paste("eig_neut1", 1:nsp, sep="."),
                        paste("eig_neut2", 1:nsp, sep="."),
                        paste("r0_neut", 1:nsp, sep="."),
                        paste(paste(c("lag", "CV"), rep(c("meta", "neutral"), each=2), sep="."), rep(1:nsp, each=4), sep="."))

#lists
scalslst<-sort(unique(matout_tot$scale))
tscalelst<-sort(unique(matout_tot$tscale))
iterlst<-sort(unique(matout_tot$iter))
laglst<-sort(unique(unname(unlist(matout_tot[,c(colnames(matout_tot)[grep("lag", colnames(matout_tot))])]))))

spmat<-rep(1:6, each=length(iterlst)*nsp)
spmat_lag<-rep(1:2, each=length(iterlst)*nsp)

matout_tot$lagtot<-c(laglst, rep(NA, length(tscalelst)-length(laglst)))

qtl_lims<-c(0.025, pnorm(-1, 0, 1), 0.5, pnorm(1, 0, 1), 0.975)

#by species
matout_eigr0<-array(dim=c(length(scalslst), length(tscalelst), length(qtl_lims), nsp*6))
matout_invar<-array(dim=c(length(scalslst), length(laglst), length(qtl_lims), nsp*2))
matout_eigr0_noscl<-matout_eigr0

#means
matout_eigr0_mean<-array(dim=c(length(scalslst), length(tscalelst), length(qtl_lims),6))
matout_invar_mean<-array(dim=c(length(scalslst), length(laglst), length(qtl_lims), 2))
matout_eigr0_mean_noscl<-matout_eigr0_mean

eigr0pos<-(1:(nsp*6))+3
invarpos<-grep("CV.", colnames(matout_tot), fixed = T)

names_eigr0pos<-colnames(matout_tot)[eigr0pos]
names_invarpos<-colnames(matout_tot)[invarpos]

#format data
for(i in 1:length(scalslst)) {
  for(j in 1:length(tscalelst)) {
    sbs<-which(matout_tot$scale==scalslst[i] & matout_tot$tscale==tscalelst[j])
    
    #mean by species
    matout_eigr0_noscl[i,j,,]<-apply(matout_tot[sbs,eigr0pos], 2, function(x) quantile(x, qtl_lims, na.rm=T))
    matout_eigr0[i,j,,]<-apply(matout_tot[sbs,eigr0pos], 2, function(x) quantile(x*tscalelst[j], qtl_lims, na.rm=T))
    
    #mean across species
    matout_eigr0_mean_noscl[i,j,,]<-matrix(nrow=5, unlist(tapply(unlist(matout_tot[sbs,eigr0pos]), spmat, function(x) quantile(x, qtl_lims, na.rm=T))))
    matout_eigr0_mean[i,j,,]<-matrix(nrow=5, unlist(tapply(unlist(matout_tot[sbs,eigr0pos]), spmat, function(x) quantile(x*tscalelst[j], qtl_lims, na.rm=T))))
  }
  
  for(j in 1:length(laglst)) {
    sbs<-which(matout_tot$scale==scalslst[i] & matout_tot$lagtot==laglst[j])
    
    #mean by species
    tmp<-apply(matout_tot[sbs,invarpos], 2, function(x) quantile(x, qtl_lims, na.rm=T))
    matout_invar[i,j,,]<-tmp
    
    #mean across species
    matout_invar_mean[i,j,,]<-matrix(nrow=5, unlist(tapply(unlist(matout_tot[sbs,invarpos]), spmat_lag, function(x) quantile(x, qtl_lims, na.rm=T))))
  }
  
  print(round(i/length(scalslst), 2))
}


###################################
# Make Plots
###################################
sclsuse<-c(1:7)
collst<-c("red2", "darkorange2", "gold1", "forestgreen", "dodgerblue", "dodgerblue4", "orchid4")

###############
# By Mean
###############

#lambda t
for(adpos in c(1,2,4,5)) {
  plot(range(tscalelst), range(c(0, matout_eigr0_mean[sclsuse,,c(3),adpos]), na.rm=T), type="n", xlab="time scale", ylab=expression(italic(paste(lambda, "t"))), xaxs="i")
  abline(h=0, lty=3)
  
  n<-1
  for(scl in sclsuse) {
    lines(tscalelst, matout_eigr0_mean[scl,,c(3),adpos], col=collst[n], lty=1, lwd=lwduse)
    #polygon(c(tscalelst, rev(tscalelst)), c(matout_eigr0_mean[scl,,c(2),adpos], rev(matout_eigr0_mean[scl,,c(4),adpos])),
    #          col=adjustcolor(collst[n], alpha.f = 0.2), border = NA)
    n<-n+1
  }
}

#r0 t
for(adpos in c(3,6)) {
  plot(range(tscalelst), range(c(0, matout_eigr0_mean[sclsuse,,c(3),adpos]), na.rm=T), type="n", xlab="time scale", ylab=expression(italic(paste(r[0], "t"))), xaxs="i")
  abline(h=0, lty=3)
  
  n<-1
  for(scl in sclsuse) {
    lines(tscalelst, matout_eigr0_mean[scl,,c(3),adpos], lty=c(1), col=collst[n], lwd=lwduse)
    #polygon(c(tscalelst, rev(tscalelst)), c(matout_eigr0_mean[scl,,c(2),adpos], rev(matout_eigr0_mean[scl,,c(4),adpos])),
    #        col=adjustcolor(collst[n], alpha.f = 0.2), border = NA)
    n<-n+1
  }
}

#invar
for(adpos in c(1,2)) {
  plot(range(laglst), range(c(0, matout_invar_mean[sclsuse,,c(3),adpos]), na.rm=T), type="n", xlab="time lag", ylab=expression(italic(paste("CV"))), xaxs="i")
  abline(h=0, lty=3)
  
  n<-1
  for(scl in sclsuse) {
    lines(laglst, matout_invar_mean[scl,,c(3),adpos], col=collst[n], lty=1, lwd=lwduse)
    #polygon(c(laglst, rev(tscalelst)), c(matout_invar_mean[scl,,c(2),adpos], rev(matout_invar_mean[scl,,c(4),adpos])),
    #          col=adjustcolor(collst[n], alpha.f = 0.2), border = NA)
    n<-n+1
  }
}







###############
# By Species
###############
pdf("figures/stability_by_scale_iter.pdf", width=6, height=4, colormodel = "cmyk")
par(mfrow=c(2,nsp), mar=c(2,2,1,1), oma=c(3,4,2,0))

tmp<-c(intersect(grep("eig_", names_eigr0pos, fixed=T), grep("1.", names_eigr0pos, fixed=T)),
       intersect(grep("eig_", names_eigr0pos, fixed=T), grep("2.", names_eigr0pos, fixed=T)))

lwduse<-1.5

#lambda t
m<-1
for(adpos in tmp) {
  plot(range(tscalelst), range(c(0, matout_eigr0[sclsuse,,c(3),adpos]), na.rm=T), type="n", xlab="time scale", ylab=expression(italic(paste(lambda, "t"))), xaxs="i",
       ylim=c(-2.5, 0.5))
  abline(h=0, lty=3)
  
  if(m<=nsp) {
    mtext(paste("sp.", m), side=3, outer=F, line=1)
  }
  
  if(m==1) {
    mtext("Levins Model", 2, line=2.5)
  } else if(m==(nsp+1)) {
    mtext("Neutral Model", 2, line=2.5)
  }
  
  n<-1
  for(scl in sclsuse) {
    matlines(tscalelst, matout_eigr0[scl,,c(3),adpos], col=collst[n], lty=1, lwd=lwduse)
    n<-n+1
  }
  box()
  
  m<-m+1
  if(m>(2*nsp)) {
    mtext(expression(italic(paste(lambda, "t"))), side=2, outer=T, line=1.5, cex=1.8)
    mtext("time span", side=1, outer=T, line=1.5, cex=1.5)
    m<-1
    
    legend("bottomright", legend = scalslst[sclsuse], col = collst, lty=1, lwd=lwduse, bty="n", title = "spatial extent", ncol=2, cex=0.7)
  }
}


#r0 t
m<-1
for(adpos in grep("r0_", names_eigr0pos, fixed=T)) {
  plot(range(tscalelst), range(c(0, matout_eigr0[sclsuse,,c(3),adpos]), na.rm=T), type="n", xlab="time scale", ylab=expression(italic(paste(r[0], "t"))), xaxs="i",
       ylim=c(-0.5,4.5))
  abline(h=0, lty=3)
  
  if(m<=nsp) {
    mtext(paste("sp.", m), side=3, outer=F, line=1)
  }
  
  if(m==1) {
    mtext("Levins Model", 2, line=2.5)
  } else if(m==(nsp+1)) {
    mtext("Neutral Model", 2, line=2.5)
  }
  
  n<-1
  for(scl in sclsuse) {
    matlines(tscalelst, matout_eigr0[scl,,c(3),adpos], lty=c(1), col=collst[n], lwd=lwduse)
    n<-n+1
  }
  box()
  
  m<-m+1
  if(m>(2*nsp)) {
    mtext(expression(italic(paste(r[0], "t"))), side=2, outer=T, line=1.5, cex=1.8)
    mtext("time span", side=1, outer=T, line=1.5, cex=1.5)
    
    legend("topright", legend = scalslst[sclsuse], col = collst, lty=1, lwd=lwduse, bty="n", title = "spatial extent", ncol=2, cex=0.7)
  }
}


posmat<-matrix(ncol=2, data=c(1,3,5,7,2,4,6,8))
m<-1
for(adpos in 1:length(posmat)) {
  plot(range(laglst[1:16]), range(c(0, matout_invar[sclsuse,1:16,c(3), posmat[adpos]]), na.rm=T), type="n", xlab="time lag", ylab=expression(italic(paste("CV"))), xaxs="i",
       ylim=c(0, 0.35))
  abline(h=0, lty=3)
  
  if(m<=nsp) {
    mtext(paste("sp.", m), side=3, outer=F, line=1)
  }
  
  if(m==1) {
    mtext("Levins Model", 2, line=2.5)
  } else if(m==(nsp+1)) {
    mtext("Neutral Model", 2, line=2.5)
  }
  
  n<-1
  for(scl in sclsuse) {
    matlines(laglst[1:16], matout_invar[scl,1:16,c(3),posmat[adpos]], col=collst[n], lty=1, lwd=lwduse)
    n<-n+1
  }
  box()
  
  if(m==1) {
    legend("topright", legend = scalslst[sclsuse], col = collst, lty=1, lwd=lwduse, bty="n", title = "spatial extent", ncol=2, cex=0.7)
  }
  
  m<-m+1
  if(m>(2*nsp)) {
    mtext(expression(italic("CV")), side=2, outer=T, line=1.5, cex=1.5)
    mtext("time lag", side=1, outer=T, line=1.5, cex=1.5)
  }
}
dev.off()
