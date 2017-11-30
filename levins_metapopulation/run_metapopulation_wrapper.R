#!/usr/bin/env Rscript
#error
#rm(list=ls())
#setwd("~/Dropbox/ActiveWork/Projects/019_Dave habitat loss/src/")

####################
# Functions
####################
makegrid<-function(xlng=50, ylng=50) {
  xpos<-rep(1:xlng, ylng)
  ypos<-rep(1:ylng, each=xlng)
  
  posmat<-matrix(nrow=xlng, ncol=ylng, data=1:length(xpos))
  
  return(list(xpos=xpos, ypos=ypos, lng=c(xlng, ylng), posmat=posmat))
}


populate<-function(gridout, nlst=0.1, clst=c(3,20), mlst=rep(0.1, length(clst)), radlst=3) {
  if(sum(nlst)<1) {
    nlst<-rep(round(prod(gridout$lng)*nlst), length(clst))
    while(sum(nlst)/prod(gridout$lng)>1) {
      nlst<-nlst/(prod(gridout$lng)*1.2)
      nlst<-round(nlst)
    }
  }
  
  spid<-factor(rep(1:length(nlst), nlst)) #species identities
  deathdate<-NULL #scheduled date of death
  #radlst #dispersal distance for each species
  
  #species locations
  pos_sp<-sample(1:length(gridout$xpos), sum(nlst))
  patch_occupied<-numeric(length(gridout$xpos))
  patch_occupied[pos_sp]<-1
  
  return(list(spid=spid, nlst=nlst, clst=clst, mlst=mlst, radlst=radlst, pos_sp=pos_sp, patch_occupied=patch_occupied))
}

#tmax=20; nsteps=20; talktime=1
run_metapopulation<-function(tmax, nsteps=tmax, gridout, population, talktime=1, runtype="metapopulation") {
  
  if(runtype=="metapopulation") {
    system("R CMD SHLIB run_metapopulation.c", ignore.stdout = TRUE)
    if(!is.loaded("run_metapopulation")) {
      dyn.load("run_metapopulation.so")
    } else {
      dyn.unload("run_metapopulation.so")
      dyn.load("run_metapopulation.so")
    }
  } else if(runtype=="neutral") {
    system("R CMD SHLIB run_neutral_metapopulation.c", ignore.stdout = TRUE)
    if(!is.loaded("run_neutral_metapopulation")) {
      dyn.load("run_neutral_metapopulation.so")
    } else {
      dyn.unload("run_neutral_metapopulation.so")
      dyn.load("run_neutral_metapopulation.so")
    }
  } else {
    return("error: run type must be 'metapopulation' or 'neutral'")
  }
  
  gridsize<-prod(gridout$lng); nsp<-length(population$clst); xylim<-gridout$lng; destroyed<-rep(0, gridsize)
  c_sptraits<-population$clst; m_sptraits<-population$mlst; abundances<-population$nlst
  output<-numeric((nsteps+1)*(nsp+1))
  
  #make_colsites
  if(is.infinite(population$radlst)) {
    colsites<-cbind(gridout$xpos, gridout$ypos)-1
    colsites<-colsites[!(colsites[,1]==0 & colsites[,2]==0),]
  } else {
    dists<-seq(-population$radlst, population$radlst)
    colsites<-expand.grid(dists, dists)
    colsites<-colsites[!(colsites[,1]==0 & colsites[2]==0),]
    colsites<-colsites[sqrt(colsites[,1]^2+colsites[,2]^2)<=population$radlst,]
    colsites<-unname(unlist(colsites))
  }
  ncolsites<-length(colsites)/2
  
  #species info
  speciesid<-rep(nsp, gridsize) #empty
  speciesid[population$pos_sp]<-as.numeric(population$spid)-1
  
  #make events
  eventtimes_c<-numeric(gridsize)
  eventtimes_m<-numeric(gridsize)
  
  n<-1
  for(i in 1:length(population$spid)) {
    x<-runif(1)
    eventtimes_c[population$pos_sp[i]]<-log(-x+1)/(-c(population$clst[population$spid[i]]))
    
    x<-runif(1)
    eventtimes_m[population$pos_sp[i]]<-log(-x+1)/(-c(population$mlst[population$spid[i]]))
  }
  
  eventtimes_c[!is.finite(eventtimes_c)]<-tmax+1
  eventtimes_m[!is.finite(eventtimes_m)]<-tmax+1
  
  #print(speciesid)
  #print(round(eventtimes_c, 3))
  #print(round(eventtimes_m, 3))
  
  if(gridsize>1e6) {
    print("error: must compile code with larger buffer size!")
  } else {
    
    if(runtype=="metapopulation") {
      runname<-"run_metapopulation"
    } else if(runtype=="neutral") {
      runname<-"run_neutral_metapopulation"
    }
    cout<-.C(runname,
             ptmax= as.double(tmax), pgridsize=as.integer(gridsize), pnsp=as.integer(nsp), xylim=as.integer(xylim), destroyed=as.integer(destroyed), #grid
             c_sptraits=as.double(c_sptraits), m_sptraits=as.double(m_sptraits), abundances=as.integer(abundances), colsites=as.integer(colsites), pncolsites=as.integer(ncolsites), #traits
             eventtimes_c=as.double(eventtimes_c), eventtimes_m=as.double(eventtimes_m), #events
             speciesid=as.integer(speciesid), #species
             output=as.double(output), pnsteps=as.integer(nsteps),
             ptalktime=as.integer(talktime))
    
    if(runtype=="metapopulation") {
      dyn.unload("run_metapopulation.so")
    } else if(runtype=="neutral") {
      dyn.load("run_neutral_metapopulation.so")
    }
    
    out<-matrix(cout$output, nrow=nsteps+1)
    nsp<-ncol(out)-1
    ngrid<-prod(gridout$lng)
    
    ceq<-numeric(nsp)
    for(i in 1:nsp) {
      ceq[i]<-(1-population$mlst[i]/population$clst[i])-sum(ceq[0:(i-1)]*(1+population$clst[0:(i-1)]/population$clst[i]))
    }
    
    plotdata<-list(ceq=ceq, ngrid=ngrid)
    
    out<-out[out[,1]!=0 | (1:nrow(out))==1,]
    
    return(list(output=out, full=cout, plotdata=plotdata))
  }
}


rerunrun_metapopulation<-function(out, tmax, nsteps=tmax, talktime=1, runtype="metapopulation", perturb=rep(0, out$full$pnsp), perturbsites=1:out$plotdata$ngrid) {
  if(!all(perturb<=1 & perturb>=0)) {
    return("error!: perturb elements must be between 0 and 1")
  }
  
  if(runtype=="metapopulation") {
    if(!is.loaded("run_metapopulation")) {
      dyn.load("run_metapopulation.so")
    } else {
      dyn.unload("run_metapopulation.so")
      dyn.load("run_metapopulation.so")
    }
    
    runname<-"run_metapopulation"
  } else if(runtype=="neutral") {
    if(!is.loaded("run_neutral_metapopulation")) {
      dyn.load("run_neutral_metapopulation.so")
    } else {
      dyn.unload("run_neutral_metapopulation.so")
      dyn.load("run_neutral_metapopulation.so")
    }
    
    runname<-"run_neutral_metapopulation"
  }
  
  eventtimes_c<-out$full$eventtimes_c-out$full$ptmax
  eventtimes_c[eventtimes_c<0]<-0
  eventtimes_m<-out$full$eventtimes_m-out$full$ptmax
  eventtimes_m[eventtimes_m<0]<-0
  
  output<-numeric((nsteps+1)*(out$full$pnsp+1))
  
  abundances<-out$full$abundances
  speciesid<-out$full$speciesid
  
  #conduct perturbation
  if(sum(abs(perturb))>0) {
    pert_abund<-ceiling(table(speciesid[perturbsites])[1:out$full$pnsp]*perturb)
    
    for(i in 1:length(pert_abund)) {
      if(pert_abund[i]>0) {
        #select individuals to remove
        remove<-sample(which(speciesid[perturbsites]==(i-1)), pert_abund[i], rep=FALSE)
        speciesid[perturbsites][remove]<-out$full$pnsp
        
        #remove c and m events from killed off individuals
        eventtimes_c[perturbsites][remove]<-0
        eventtimes_m[perturbsites][remove]<-0
        
        abundances[i]<-abundances[i]-pert_abund[i]
      }
    }
  }
  
  
  cout<-.C(runname,
           ptmax= as.double(tmax), pgridsize=as.integer(out$full$pgridsize), pnsp=as.integer(out$full$pnsp), xylim=as.integer(out$full$xylim), destroyed=as.integer(out$full$destroyed), #grid
           c_sptraits=as.double(out$full$c_sptraits), m_sptraits=as.double(out$full$m_sptraits), abundances=as.integer(abundances), colsites=as.integer(out$full$colsites), pncolsites=as.integer(out$full$pncolsites), #traits
           eventtimes_c=as.double(eventtimes_c), eventtimes_m=as.double(eventtimes_m), #events
           speciesid=as.integer(speciesid), #species
           output=as.double(output), pnsteps=as.integer(nsteps),
           ptalktime=as.integer(talktime))
  
  if(runtype=="metapopulation") {
    dyn.unload("run_metapopulation.so")
  } else if(runtype=="neutral") {
    dyn.load("run_neutral_metapopulation.so")
  }
  
  outnew<-matrix(cout$output, nrow=nsteps+1)
  nsp<-ncol(outnew)-1
  ngrid<-prod(gridout$lng)
  
  ceq<-numeric(nsp)
  for(i in 1:nsp) {
    ceq[i]<-(1-out$full$m_sptraits[i]/out$full$c_sptraits[i])-sum(ceq[0:(i-1)]*(1+out$full$c_sptraits[0:(i-1)]/out$full$c_sptraits[i]))
  }
  
  plotdata<-list(ceq=ceq, ngrid=ngrid)
  
  outnew<-outnew[outnew[,1]!=0 | (1:nrow(outnew))==1,]
  
  return(list(output=outnew, full=cout, plotdata=plotdata))
}








getceq<-function(clst, mlst=rep(0.1, length(clst))) {
  nsp<-length(clst)
  ceq<-numeric(nsp)
  for(i in 1:nsp) {
    ceq[i]<-(1-mlst[i]/clst[i])-sum(ceq[0:(i-1)]*(1+clst[0:(i-1)]/clst[i]))
  }
  ceq
}

plot_metapop<-function(output) {
  out<-output$out
  ceq<-output$plotdata$ceq
  ngrid<-output$plotdata$ngrid
  
  sbs<-which(out[,1]>0)
  if(sbs[1]!=1) {
    sbs<-c(1, sbs)
  }
  #out[,-1][out[,-1]==0]<-NA
  matplot(out[sbs,1], out[sbs,-1]/ngrid, type="l", xlab="time", ylab="p",
          lty=1, col=1:ncol(out), lwd=2, ylim=c(0,1))
  abline(h=c(0,1), lty=3)
  abline(h=ceq,
         lty=2, col=1:ncol(out), lwd=2)
}

plot_map<-function(out, gridout) {
  tmp<-out$full$speciesid; tmp[tmp==out$full$pnsp]<-NA
  plot(gridout$xpos, gridout$ypos, col=tmp+1, pch=16)
}


rewrap_pop<-function(out, population) {
  snp<-length(population$nlst)
  
  spid<-NULL
  pos_sp<-NULL
  
  for(i in 1:snp) {
    spid<-c(spid, rep(i, out$full$abundances[i]))
    pos_sp<-c(pos_sp, which(out$full$speciesid==(i-1)))
  }
  spid<-factor(spid, levels=1:snp)
  
  patch_occupied<-numeric(length(gridout$xpos))
  patch_occupied[pos_sp]<-1
  
  population<-list(spid=spid,
                   nlst=unname(table(spid)),
                   clst=population$clst,
                   mlst=population$mlst,
                   radlst=population$radlst,
                   pos_sp=pos_sp,
                   patch_occupied=patch_occupied)
  
  population
}


#rEDM functions
predict_vs_L<-function(outcol, E=1, burnin=0, Luse=floor((seq((30), (length(outcol)-burnin), length=20))), niter=0,  doplot=TRUE) {
  
  #calculates predictive power vs. library length
  #niter=0 calculates for all possible time windows
  predmat<-NULL
  
  for(i in 1:length(Luse)) {
    #find windows to use for prediction
    window_starts<-seq((burnin+1), (length(outcol)-(Luse[i]-1)), by=Luse[i])
    if(niter!=0 & niter<length(window_starts)) {
      window_use<-sample(window_starts, niter, replace = F)
    } else {
      window_use<-window_starts
    }
    
    predmat_tmp<-matrix(ncol=5, nrow=length(window_use))
    predmat_tmp[,1]<-Luse[i]
    
    #calculate predictive power
    for(j in 1:length(window_use)) {
      tmp<-suppressWarnings(simplex(outcol, E=E, lib=c(window_use[j], (window_use[j]+Luse[i]-1)), pred=c(1+burnin, length(outcol)), stats_only = FALSE))
      
      predmat_tmp[j,2]<-tmp$rho
      predmat_tmp[j,3]<-tmp$rmse#/mean(tmp[[1]]$model_output$pred, na.rm=T)
      
      predmat_tmp[j,4]<-mean(tmp$model_output[[1]]$pred, na.rm=T)
      predmat_tmp[j,5]<-tmp$num_pred
    }
    
    predmat<-rbind(predmat, predmat_tmp)
  }
  
  colnames(predmat)<-c("Luse", "rho", "rmse", "mean", "n")
  
  CVest<-t(matrix(nrow=5, unlist(tapply(predmat[,"rmse"]/predmat[,"mean"], predmat[,"Luse"], function(x) quantile(x[is.finite(x)], c(0.025, pnorm(-1, 0, 1), 0.5, pnorm(1, 0, 1), 0.975), na.rm=T)))))
  
  #minimum library length that gives us 50% of the best predictive power
  Lmin<-Luse[min(which((CVest[,5]-min(CVest[,5]))/(CVest[,5])<0.5))]
  
  if(doplot) {
    plot(Luse, CVest[,3], type="n", xlab="library length", ylab="CV", ylim=c(0, max(c(1, CVest[,3]))), xlim=c(min(c(0, Luse)), max(Luse)), xaxs="i")
    polygon(c(Luse, rev(Luse)), c(CVest[,1], rev(CVest[,5])), col=adjustcolor("black", alpha.f = 0.2), border=NA)
    polygon(c(Luse, rev(Luse)), c(CVest[,2], rev(CVest[,4])), col=adjustcolor("black", alpha.f = 0.2), border=NA)
    lines(Luse, CVest[,3], lwd=2)
    
    abline(v=c(min(Luse), Lmin), h=c(0, 1), lty=3)
  }
  
  return(list(Lmin=Lmin, CVest=CVest, predmat=predmat))
}

test_predict_tlag<-function(outcol, Luse, E=1, burnin=0, laglst=c(floor((seq((0), ((length(outcol)-burnin-Luse)), length=20)))), niter=0, doplot=TRUE) {
  predmat<-NULL
  
  #Cycle through lag lengths
  for(i in 1:length(laglst)) {
    #Find all possible comparisons
    p1full<-seq(burnin+1, length(outcol)-(Luse-1)-laglst[i], by=Luse)
    
    if(niter!=0 & niter<length(p1full)) {
      p1<-sample(p1full, niter, replace = F)
    } else {
      p1<-p1full
    }
    
    p2<-p1+laglst[i]
    
    predmat_tmp<-matrix(ncol=5, nrow=length(p1))
    predmat_tmp[,1]<-laglst[i]
    
    for(j in 1:length(p1)) {
      tmplib<-c(p1[j], (p1[j]+Luse-1))
      tmppred<-c(p2[j], (p2[j]+Luse-1))
      
      tmp<-suppressWarnings(simplex(outcol,E=E, lib=tmplib, pred=tmppred, stats_only = FALSE))
      
      predmat_tmp[j,2]<-tmp$rho
      predmat_tmp[j,3]<-tmp$rmse#/mean(tmp[[1]]$model_output$pred, na.rm=T)
      
      predmat_tmp[j,4]<-mean(tmp$model_output[[1]]$pred, na.rm=T)
      predmat_tmp[j,5]<-tmp$num_pred
    }
    predmat<-rbind(predmat, predmat_tmp)
  }
  
  colnames(predmat)<-c("Luse", "rho", "rmse", "mean", "n")
  
  CVest<-t(matrix(nrow=5, unlist(tapply(predmat[,"rmse"]/predmat[,"mean"], predmat[,"Luse"], function(x) quantile(x[is.finite(x)], c(0.025, pnorm(-1, 0, 1), 0.5, pnorm(1, 0, 1), 0.975), na.rm=T)))))
  
  if(doplot) {
    plot(laglst, CVest[,3], type="n", xlab="time lag", ylab="CV", ylim=c(0, max(c(1, CVest[,3]))), xlim=c(min(c(0, laglst)), max(laglst)), xaxs="i")
    polygon(c(laglst, rev(laglst)), c(CVest[,1], rev(CVest[,5])), col=adjustcolor("black", alpha.f = 0.2), border=NA)
    polygon(c(laglst, rev(laglst)), c(CVest[,2], rev(CVest[,4])), col=adjustcolor("black", alpha.f = 0.2), border=NA)
    lines(laglst, CVest[,3], lwd=2)
    
    abline(v=c(min(laglst)), h=c(0, 1), lty=3)
  }
  
  return(list(CVest=CVest, predmat=predmat))
}


if(FALSE) {
  gridout<-makegrid(xlng = 100, ylng = 100)
  population<-populate(gridout, nlst = rep(round(0.1*prod(gridout$lng)), 4), clst = c(0.15, 1, 10, 100), mlst = rep(0.1, 4), radlst = Inf)
  
  out<-run_metapopulation(tmax=200, nsteps = 1000, gridout, population, talktime = 0)
  
  plot_metapop(out)
}
