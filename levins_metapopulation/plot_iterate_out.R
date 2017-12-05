error
rm(list=ls())
setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")


matout_tot<-read.csv("output/matout_tot.csv")
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

matout_tot$lagtot<-c(laglst, rep(NA, length(tscalelst)-length(laglst)))

qtl_lims<-c(0.025, pnorm(-1, 0, 1), 0.5, pnorm(1, 0, 1), 0.975)

matout_eigr0<-array(dim=c(length(scalslst), length(tscalelst), length(qtl_lims), nsp*6))
matout_invar<-array(dim=c(length(scalslst), length(laglst), length(qtl_lims), nsp*2))

eigr0pos<-(1:(nsp*6))+3
invarpos<-grep("CV.", colnames(matout_tot), fixed = T)

names_eigr0pos<-colnames(matout_tot)[eigr0pos]
names_invarpos<-colnames(matout_tot)[invarpos]

#format data
for(i in 1:length(scalslst)) {
  for(j in 1:length(tscalelst)) {
    sbs<-which(matout_tot$scale==scalslst[i] & matout_tot$tscale==tscalelst[j])
    
    tmp<-apply(matout_tot[sbs,eigr0pos], 2, function(x) quantile(x, qtl_lims, na.rm=T))
    matout_eigr0[i,j,,]<-tmp
  }
  
  for(j in 1:length(laglst)) {
    sbs<-which(matout_tot$scale==scalslst[i] & matout_tot$lagtot==laglst[j])
    
    tmp<-apply(matout_tot[sbs,invarpos], 2, function(x) quantile(x, qtl_lims, na.rm=T))
    matout_invar[i,j,,]<-tmp
  }
}


#make plots

matplot(tscalelst, matout_eigr0[7,,c(1,3,5),1], type="l", lty=c(2,1,2), col=1); abline(h=0, lty=3)
matplot(tscalelst, matout_eigr0[1,,c(1,3,5),1], type="l", lty=c(2,1,2), col=1); abline(h=0, lty=3)

matplot(tscalelst, matout_eigr0[7,,c(1,3,5),13], type="l", lty=c(2,1,2), col=1); abline(h=0, lty=3)
matplot(tscalelst, matout_eigr0[1,,c(1,3,5),13], type="l", lty=c(2,1,2), col=1); abline(h=0, lty=3)

matplot(tscalelst, matout_eigr0[7,,c(1,3,5),17], type="l", lty=c(2,1,2), col=1); abline(h=0, lty=3)
matplot(tscalelst, matout_eigr0[1,,c(1,3,5),17], type="l", lty=c(2,1,2), col=1); abline(h=0, lty=3)







matplot(laglst, matout_invar[7,,c(1,3,5),1], type="l", lty=c(2,1,2), col=1); abline(h=0, lty=3)
matplot(laglst, matout_invar[7,,c(1,3,5),2], type="l", lty=c(2,1,2), col=1); abline(h=0, lty=3)

matplot(laglst, matout_invar[1,,c(1,3,5),1], type="l", lty=c(2,1,2), col=1); abline(h=0, lty=3)
matplot(laglst, matout_invar[1,,c(1,3,5),2], type="l", lty=c(2,1,2), col=1); abline(h=0, lty=3)




