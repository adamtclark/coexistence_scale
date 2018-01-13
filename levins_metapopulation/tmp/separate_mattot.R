for(i in 1:length(scalslst)) {
  sbs<-which(matout_tot$scale==scalslst[i])
  
  #separate
  matout_dyn<-matout_tot[sbs,1:123]
  matout_cv<-matout_tot[sbs,c(1,3,124:147)]
  matout_cv<-cbind(lag=c(lglst, rep(NA, max(matout_tot[,"tscale"])-length(lglst))), matout_cv)
  matout_cv<-matout_cv[!is.na(matout_cv[,"lag"]),]
  
  #save outputs to csv
  write.csv(matout_dyn, paste("output/matout_dyn_", scalslst[i], ".csv", sep=""), row.names=F)
  write.csv(matout_cv, paste("output/matout_cv_", scalslst[i], ".csv", sep=""), row.names=F)
  
  print(i)
}
