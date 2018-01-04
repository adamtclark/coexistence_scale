#!/usr/bin/env Rscript
#setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")

#TODO:
#4. Add Stan's PSF model
#5. Example of rock-paper-scissors



########################################
# Metapopulation simulation functions
########################################
makegrid<-function(xlng=50, ylng=50) {
  #creates spatial grid of size xlng by ylng
  #for use in subsequent functions
  
  xpos<-rep(1:xlng, ylng)
  ypos<-rep(1:ylng, each=xlng)
  
  posmat<-matrix(nrow=xlng, ncol=ylng, data=1:length(xpos))
  
  #xpos and ypos are vectors of x and y coordinates for each grid cell
  #lng shows length of x and y components of grid
  #posmat is a matrix of size xlng by ylng, each elements containing its numeric position
  #(i.e. 1st element is 1, last element is xlng*ylng)
  return(list(xpos=xpos, ypos=ypos, lng=c(xlng, ylng), posmat=posmat))
}

grid_subset<-function(gridout, size=0.1) {
  #creates a subset of an existing grid
  #gridout is a grid object
  #where "size" is the fraction of original area requested
  #note, subset is always in the center of the original grid,
  #and achieved size may not exactly match requested size.
  
  lng<-gridout$lng
  
  #produce a square with integer side lengths
  #that is near requested area
  #note, sz_tmp is the length of a side of that square
  sz_tmp<-round(sqrt(prod(lng)*size))
  
  #if size is zero, but requested size was greather than zero,
  #then make size equal to 1
  if(sz_tmp==0 & size!=0) {sz_tmp<-1}
  #if size is greater than that of smallest side of the original grid,
  #reduce request to maximum possible size
  if(sz_tmp>min(lng)) {sz_tmp<-min(lng)}
  
  #find boundaries of region
  x1<-ceiling(lng[1]/2)-ceiling(sz_tmp/2)+1
  x2<-ceiling(lng[1]/2)+ceiling(sz_tmp/2-1)+1
  
  y1<-ceiling(lng[2]/2)-ceiling(sz_tmp/2)+1
  y2<-ceiling(lng[2]/2)+ceiling(sz_tmp/2-1)+1
  
  #find sites within region
  sites<-c(gridout$posmat)[gridout$xpos%in%c(x1:x2) & gridout$ypos%in%c(y1:y2)]
  
  #frac_tot is achieved size, as fraction of total grid
  return(list(sites=sites, frac_tot=(sz_tmp^2)/prod(lng), borders=c(x1=x1, x2=x2, y1=y1, y2=y2)))
}

populate<-function(gridout, nlst=0.1, clst=c(3,20), mlst=rep(0.1, length(clst)), radlst=3, nsp=length(clst)) {
  #populate a grid with a community
  #gridout is a grid object
  #nlst is the desired starting population for species,
  #(in fraction of total grid size, or in number of grid cells)
  #note, can be single number, or one per species.
  #clst is colonization rates for each species
  #mlst are mortality rates for each species
  #radlst is dispersal radius (in grid cells) - set to "Inf" for global dispersal
  #nsp is optional argument describing number of species
  
  #expand if only one number given for >1 species
  if(length(nlst)<nsp) {
    nlst<-rep(nlst[1], nsp)
  }
  
  #expand nlst if given as fraction
  if(any(nlst<1)) {
    nlst<-round(prod(gridout$lng)*nlst)
  } else {
    #otherwise, force to integer
    nlst<-round(nlst)
  }
  
  #if sizes are too large, reduce until they fit in the grid
  while((sum(nlst)/prod(gridout$lng))>1) {
    nlst<-nlst*0.9
    nlst<-round(nlst)
  }
  
  #expand clst and mlst if only one element and n>1
  if(length(clst)!=nsp) {
    clst<-rep(clst[1], nsp)
  }
  if(length(mlst)!=nsp) {
    mlst<-rep(mlst[1], nsp)
  }
  
  #create a vector with species ID's
  spid<-factor(rep(1:length(nlst), nlst)) #species identities
  deathdate<-NULL #scheduled date of death
  
  #select random starting species locations
  pos_sp<-sample(1:length(gridout$xpos), sum(nlst))
  patch_occupied<-numeric(length(gridout$xpos))
  patch_occupied[pos_sp]<-1
  
  #spid is vector of species identities for each individual
  #nlst is abundances for each species
  #pos_sp is location for each individual
  #patch_occupied identifies each grid cell as occupied or not occupied
  return(list(spid=spid, nlst=nlst, clst=clst, mlst=mlst, radlst=radlst, pos_sp=pos_sp, patch_occupied=patch_occupied))
}

loadrun<-function(runtype) {
  trigger<-NULL
  #loads compiled files for runs
  #NOTE to run in Winows, need to replace ".so" with ".dll" here and is all subsequent files
  
  if(runtype%in%c("metapopulation", "disturbance")) {
    system("R CMD SHLIB run_metapopulation_spatialsub.c", ignore.stdout = TRUE)
    if(!is.loaded("run_metapopulation_spatialsub")) {
      dyn.load("run_metapopulation_spatialsub.so")
    } else {
      dyn.unload("run_metapopulation_spatialsub.so")
      dyn.load("run_metapopulation_spatialsub.so")
    }
  } else if(runtype%in%c("neutral")) {
    system("R CMD SHLIB run_neutral_metapopulation_spatialsub.c", ignore.stdout = TRUE)
    if(!is.loaded("run_neutral_metapopulation_spatialsub")) {
      dyn.load("run_neutral_metapopulation_spatialsub.so")
    } else {
      dyn.unload("run_neutral_metapopulation_spatialsub.so")
      dyn.load("run_neutral_metapopulation_spatialsub.so")
    }
  } else {
    trigger<-"error: run type must be 'metapopulation', 'neutral', or disturbance'"
  }
  return(trigger)
}

run_metapopulation<-function(tmax, nsteps=tmax, gridout, population, talktime=1, runtype="metapopulation", sites_sub=0, prt=0, prtfrq=0) {
  #tmax is maximum time for simulation
  #nsteps is number of time steps for which to save output
  #gridout is a grid object, population is a population object
  #talktime = 1 produces verbose output describing simulation progress, = 0 is silent
  #runtype must be one of the following:
    #1. 'metapopulation' - Levins model
    #2. 'neutral' - NZS Hubbell model
    #3. 'disturbance' - Levins with periodic disturbances
  #sites_sub is a vector of sites desired for spatial subset - defaults to zero when not tracked
  
  #prt and prtfrq are only applied when using the 'disturbance' runtype
  #prt is vector with one element per species, specifying fraction of population lost during perturbation events
  #prtfrq specifies perturbation frequency (in time steps)
  
  #load compiled files - exit if runtype is not valid
  runtp<-loadrun(runtype)
  if(!is.null(runtp)) {
    return(runtp)
  }
  
  #extract information from grid and population
  gridsize<-prod(gridout$lng); nsp<-length(population$clst); xylim<-gridout$lng; destroyed<-rep(0, gridsize)
  c_sptraits<-population$clst; m_sptraits<-population$mlst; abundances<-population$nlst
  
  #vector for output
  output<-numeric((nsteps+1)*(nsp+1))
  
  #make colsites vector
  #shows potential sites available for colonization based on rad
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
  
  #populate sites with individuals
  speciesid<-rep(nsp, gridsize)
  speciesid[population$pos_sp]<-as.numeric(population$spid)-1
  
  #specify times for next birth and death events, using random exponential distribution
  eventtimes_c<-numeric(gridsize)
  eventtimes_m<-numeric(gridsize)
  
  n<-1
  for(i in 1:length(population$spid)) {
    x<-runif(1)
    eventtimes_c[population$pos_sp[i]]<-log(-x+1)/(-c(population$clst[population$spid[i]]))
    
    x<-runif(1)
    eventtimes_m[population$pos_sp[i]]<-log(-x+1)/(-c(population$mlst[population$spid[i]]))
  }
  
  #truncate infinite time spans
  eventtimes_c[!is.finite(eventtimes_c)]<-tmax+1
  eventtimes_m[!is.finite(eventtimes_m)]<-tmax+1
  
  if(gridsize>1e6) {
    #error if memory allocation in compiled file needs to be increased (variable "indiv" in .c code)
    #we use this error, rather than a malloc command, because of a bug in malloc for very large vectors
    #in some older compilers
    print("error: must compile code with larger buffer size, or decrease grid size!")
  } else {
    
    #script names for runs
    if(runtype%in%c("metapopulation", "disturbance")) {
      runname<-"run_metapopulation_spatialsub"
    } else if(runtype%in%c("neutral")) {
      runname<-"run_neutral_metapopulation_spatialsub"
    }
    
    #create extra output data for disturbance
    if(runtype%in%c("disturbance")) {
      tmax_sub<-prtfrq
      nprt_cyc<-round(tmax/prtfrq)
      nsteps_sub<-tmax_sub
      
      output<-numeric((nsteps_sub+1)*(nsp+1))
    }
    
    #note, for now, conditional below includes all potential run types
    #rock/paper/scissors etc. would fall in an "else" loop
    if(runtype%in%c("metapopulation", "neutral", "disturbance")) {
      
      #get lengths, outputs, and abundances for subsets - note, this may be of length zero
      pnsites_sub<-length(sites_sub[sites_sub!=0])
      output_sub<-output
      if(sum(sites_sub)>0) {
        c_sites_sub<-sites_sub-1
        abundances_sub<-unname(table(speciesid[sites_sub])[1:nsp])
      } else {
        c_sites_sub<-0
        abundances_sub<-rep(0, nsp)
      }
      
      if(runtype=="disturbance") {
        #matrices for storing total output
        output_tot<-NULL
        output_spatial_tot<-NULL
        
        #list for storing total output
        cout_lst<-NULL
        
        #complete first run before first disturbance event
        #(see c script for description of variables)
        cout<-.C(runname,
                 ptmax= as.double(tmax_sub), pgridsize=as.integer(gridsize), pnsp=as.integer(nsp), xylim=as.integer(xylim), destroyed=as.integer(destroyed), spdestroy=as.integer(rep(0, nsp)), #grid
                 c_sptraits=as.double(c_sptraits), m_sptraits=as.double(m_sptraits), abundances=as.integer(abundances), colsites=as.integer(colsites), pncolsites=as.integer(ncolsites), #traits
                 eventtimes_c=as.double(eventtimes_c), eventtimes_m=as.double(eventtimes_m), #events
                 speciesid=as.integer(speciesid), #species
                 output=as.double(output), pnsteps=as.integer(nsteps_sub),
                 ptalktime=as.integer(talktime),
                 abundances_sub=as.integer(abundances_sub), sites_sub=as.integer(c_sites_sub), pnsites_sub=as.integer(pnsites_sub), output_sub=as.double(output_sub))
        cout_lst[[1]]<-cout
        
        #conduct subsequent disturbances
        if(nprt_cyc>1) {
          #data for rerun script
          nsp<-length(population$clst)
          ngrid<-prod(gridout$lng)
          ceq<-numeric(nsp)
          plotdata<-list(ceq=ceq, ngrid=ngrid)
          
          for(i in 1:(nprt_cyc-1)) {
            #extract output from previous run
            out<-matrix(cout$output, nrow=nsteps_sub+1)
            out<-out[out[,1]!=0 | (1:nrow(out))==1,]
            
            #extract subset output from previous run
            out_spatial<-matrix(cout$output_sub, nrow=nsteps_sub+1)
            out_spatial<-out_spatial[out_spatial[,1]!=0 | (1:nrow(out_spatial))==1,]
            
            #create temporary list for re-run
            out_tmp<-list(output=out, full=cout, plotdata=plotdata, output_spatial=0, sites_sub=sites_sub)
            
            #add total accumulated time to time vector
            out[,1]<-out[,1]+tmax_sub*(i-1)
            out_spatial[,1]<-out_spatial[,1]+tmax_sub*(i-1)
            #store output into full lists
            output_tot<-rbind(output_tot, out)
            output_spatial_tot<-rbind(output_spatial_tot, out_spatial)
            
            #rerun to next disturbance event
            #(see script for description of variables)
            out_rerun<-rerunrun_metapopulation(out_tmp, tmax_sub, nsteps=tmax_sub, talktime=0, runtype="metapopulation", perturb=prt, sites_sub = sites_sub)
            cout<-out_rerun$full
            
            cout_lst[[1+i]]<-cout
          }
          #save outputs to regular names
          out<-output_tot
          out_spatial<-output_spatial_tot
          cout<-cout_lst
          
        } else {
          #if only one disturbance cycle occurrs, then extract information as usual
          out_spatial<-matrix(cout$output_sub, nrow=nsteps+1)
          out_spatial<-out_spatial[out_spatial[,1]!=0 | (1:nrow(out_spatial))==1,]
        }
        
      } else {
        #run c script
        #(see c script for description of variables)
        cout<-.C(runname,
                 ptmax= as.double(tmax), pgridsize=as.integer(gridsize), pnsp=as.integer(nsp), xylim=as.integer(xylim), destroyed=as.integer(destroyed), spdestroy=as.integer(rep(0, nsp)), #grid
                 c_sptraits=as.double(c_sptraits), m_sptraits=as.double(m_sptraits), abundances=as.integer(abundances), colsites=as.integer(colsites), pncolsites=as.integer(ncolsites), #traits
                 eventtimes_c=as.double(eventtimes_c), eventtimes_m=as.double(eventtimes_m), #events
                 speciesid=as.integer(speciesid), #species
                 output=as.double(output), pnsteps=as.integer(nsteps),
                 ptalktime=as.integer(talktime),
                 abundances_sub=as.integer(abundances_sub), sites_sub=as.integer(c_sites_sub), pnsites_sub=as.integer(pnsites_sub), output_sub=as.double(output_sub))
        
        #save output
        out_spatial<-matrix(cout$output_sub, nrow=nsteps+1)
        out_spatial<-out_spatial[out_spatial[,1]!=0 | (1:nrow(out_spatial))==1,]
      }
    }
    
    #unload c code after run this is necesary because of an aparent bug,
    #that seems to prevent R from succesfully closing down c scripts
    #with large memory allocations
    if(runtype=="metapopulation"  || (runtype=="disturbance" && nprt_cyc==1)) {
      dyn.unload("run_metapopulation_spatialsub.so")
    } else if(runtype=="neutral") {
      dyn.unload("run_neutral_metapopulation_spatialsub.so")
    }
    
    #update ouput and plotting data
    if((!(runtype %in% c("disturbance"))) || (nprt_cyc==1)) {
      out<-matrix(cout$output, nrow=nsteps+1)
      nsp<-ncol(out)-1
      ngrid<-prod(gridout$lng)
    }
    
    #calculate analytical expectation of equilibrium (for Levins model)
    ceq<-getceq(population$clst, population$mlst)
    
    #save data for use in plotting
    plotdata<-list(ceq=ceq, ngrid=ngrid)
    
    #remove empty output data
    out<-out[out[,1]!=0 | (1:nrow(out))==1,]
    
    #out is matrix of times and abundances
    #full is total data from c run
    #output_spatial is times and abundances for spatial subset
    #sites_sub are the locations in the spatial subset
    return(list(output=out, full=cout, plotdata=plotdata, output_spatial=out_spatial, sites_sub=sites_sub))
  }
}




rerunrun_metapopulation<-function(out, tmax, nsteps=tmax, talktime=1, runtype="metapopulation", perturb=rep(0, length(out$plotdata$ceq)), perturbsites=1:out$plotdata$ngrid, addn=0, addsites=perturbsites, replace_perturb=0, sites_sub=0, prt=0, prtfrq=0, habitatdestructsites=0, habitatdestruct_species=rep(0, length(out$plotdata$ceq))) {
  #reruns the simulation models in "run_metapopulation", starting at the last event in the previous simulation.
  #useful for iterating the model - e.g. for simulating disturbance events.
  
  #out is the output from 'run_metapopulation'
  #tmax, nsteps, talktime, runtype, sites_sub, prt, and prtfrq are as described in 'run_metapopulation'
  #perturb is a vector listing size of negative perturbation to be enacted on each species (defaults to zero)
  #perturbsites is a vector of locations at which negative perturbations take place
  #addn is a vector listing size of positive perturbation to be enacted on each species (defaults to zero)
  #addsites is a vector of locations at which positive perturbations take place
  #replace_perturb=1 replaces removed individuals with those from a randomly chosen (other species); =0 does not replace individuals
  
  #habitatdestructsites and habitatdestruct_species are internal variables used in the "invasion when rare" routine, and prevent
  #individual species from colonizing a particular region
  
  ################ Extract data
  #pull out last full observation, if there are multiple in the list
  if(!length(unique(names(out$full)))>1) {
    out$full<-out$full[[length(out$full)]]
  }
  
  #make sure perturbations are of correct magnitude
  if(!all(perturb<=1 & perturb>=0) | !all(addn<=1 & addn>=0)) {
    return("error!: perturb and addn elements must be between 0 and 1")
  }
  
  #load c code
  runtp<-loadrun(runtype)
  if(!is.null(runtp)) {
    return(runtp)
  }
  
  #names of c functions for runs
  if(runtype%in%c("metapopulation", "disturbance")) {
    runname<-"run_metapopulation_spatialsub"
  } else if(runtype%in%c("neutral")) {
    runname<-"run_neutral_metapopulation_spatialsub"
  }
  
  #extract event times from previous simulation (these have not yet occurred)
  eventtimes_c<-out$full$eventtimes_c-out$full$ptmax
  eventtimes_c[eventtimes_c<0]<-0
  eventtimes_m<-out$full$eventtimes_m-out$full$ptmax
  eventtimes_m[eventtimes_m<0]<-0
  
  #extract output from previous simulation
  output<-numeric((nsteps+1)*(out$full$pnsp+1))
  
  #extract abundances
  abundances<-out$full$abundances
  #extract species identities
  speciesid<-out$full$speciesid
  
  ################ Conduct perturbations
  #conduct negative perturbation
  if(sum(abs(perturb))>0) {
    #extract number of individuals of each species to perturb
    pert_abund<-ceiling(table(factor(speciesid[perturbsites], levels=0:out$full$pnsp))[1:out$full$pnsp]*perturb)
    
    for(i in 1:length(pert_abund)) {
      if(pert_abund[i]>0) {
        #select individuals to remove
        pstmp<-which(speciesid[perturbsites]==(i-1))
        remove<-sample(pstmp, min(c(pert_abund[i], length(pstmp))), rep=FALSE)
        speciesid[perturbsites][remove]<-out$full$pnsp
        
        #remove c and m events from killed off individuals
        eventtimes_c[perturbsites][remove]<-0
        eventtimes_m[perturbsites][remove]<-0
        
        #update abundances
        abundances[i]<-abundances[i]-pert_abund[i]
      }
    }
  }
  
  #conduct positive perturbation, or replace removed individuals
  if(sum(abs(addn))>0 | replace_perturb!=0) {
    #extract number of individuals of each species to perturb
    add_abund<-ceiling((pmax(0, length(addsites)-sum(table(factor(speciesid[addsites][speciesid[addsites]!=out$full$pnsp], levels=0:out$full$pnsp))[1:out$full$pnsp]))*addn))
    #if selected number is zero, but desired size is >0, add at least one individual
    add_abund[addn>0 & add_abund==0]<-1
    
    #replace removed individuals, if replace_perturb==1
    if(replace_perturb!=0) {
      tmp<-numeric(length(pert_abund))
      for(i in 1:length(pert_abund)) {
        if(pert_abund[i]>0) {
          if(length(pert_abund)>2) {
            tmp[sample((1:length(pert_abund))[-i],1)]<-pert_abund[i]
          } else {
            tmp[-i]<-pert_abund[i]
          }
        }
      }
      
      #total number of individuals for each species to add
      add_abund<-add_abund+tmp
    }
    
    #add new individuals
    for(i in 1:length(add_abund)) {
      if(add_abund[i]>0) {
        #select locations to add individuals
        add<-sample(which(speciesid[addsites]==out$full$pnsp), add_abund[i], rep=FALSE)
        speciesid[addsites][add]<-(i-1)
        
        #add c and m events for added individuals
        for(j in 1:add_abund[i]) {
          x<-runif(1)
          eventtimes_c[addsites][add[j]]<-log(-x+1)/(-c(out$full$c_sptraits[i]))
          
          x<-runif(1)
          eventtimes_m[addsites][add[j]]<-log(-x+1)/(-c(out$full$m_sptraits[i]))
        }
        
        eventtimes_c[!is.finite(eventtimes_c)]<-tmax+1
        eventtimes_m[!is.finite(eventtimes_m)]<-tmax+1
        
        #update abundances
        abundances[i]<-abundances[i]+add_abund[i]
      }
    }
  }
  
  ################ Run simulations
  #create extra output data for disturbance
  if(runtype%in%c("disturbance")) {
    tmax_sub<-prtfrq
    nprt_cyc<-round(tmax/prtfrq)
    nsteps_sub<-tmax_sub
    
    output<-numeric((nsteps_sub+1)*(out$full$pnsp+1))
  }
  
  #note, for now, conditional below includes all potential run types
  #rock/paper/scissors etc. would fall in an "else" loop
  if(runtype%in%c("metapopulation", "neutral", "disturbance")) {
    
    #get lengths, outputs, and abundances for subsets - note, this may be of length zero
    pnsites_sub<-length(sites_sub[sites_sub!=0])
    output_sub<-output
    if(sum(sites_sub)>0) {
      c_sites_sub<-sites_sub-1
      abundances_sub<-unname(table(factor(speciesid[sites_sub], levels=0:out$full$pnsp))[1:out$full$pnsp])
    } else {
      c_sites_sub<-0
      abundances_sub<-rep(0, out$full$pnsp)
    }
    
    #set up "destoryed" habitat for individual species
    #this routine is automatically controlled by the invasion when rare function
    tmpdestory<-out$full$destroyed
    if(sum(habitatdestructsites)==0) {
      tmpspdestroy<-rep(0, out$full$pnsp)
    } else {
      tmpspdestroy<-habitatdestruct_species
      tmpdestory[habitatdestructsites]<-1
    }
    
    if(runtype=="disturbance") {
      #matrices for storing total output
      output_tot<-NULL
      output_spatial_tot<-NULL
      
      #list for storing total output
      cout_lst<-NULL
      
      #complete first run before first disturbance event
      #(see c script for description of variables)
      cout<-.C(runname,
               ptmax= as.double(tmax_sub), pgridsize=as.integer(out$full$pgridsize), pnsp=as.integer(out$full$pnsp), xylim=as.integer(out$full$xylim), destroyed=as.integer(tmpdestory), spdestroy=as.integer(tmpspdestroy), #grid
               c_sptraits=as.double(out$full$c_sptraits), m_sptraits=as.double(out$full$m_sptraits), abundances=as.integer(abundances), colsites=as.integer(out$full$colsites), pncolsites=as.integer(out$full$pncolsites), #traits
               eventtimes_c=as.double(eventtimes_c), eventtimes_m=as.double(eventtimes_m), #events
               speciesid=as.integer(speciesid), #species
               output=as.double(output), pnsteps=as.integer(nsteps_sub),
               ptalktime=as.integer(talktime),
               abundances_sub=as.integer(abundances_sub), sites_sub=as.integer(c_sites_sub), pnsites_sub=as.integer(pnsites_sub), output_sub=as.double(output_sub))
      #replace destoryed sites
      if(sum(habitatdestructsites)!=0) {
        cout$destroyed[]<-0
      }
      
      cout_lst[[1]]<-cout
      
      #conduct subsequent disturbances
      if(nprt_cyc>1) {
        #data for rerun script
        nsp<-length(out$full$c_sptraits)
        ngrid<-prod(gridout$lng)
        ceq<-numeric(nsp)
        plotdata<-list(ceq=ceq, ngrid=ngrid)
        
        for(i in 1:(nprt_cyc-1)) {
          #extract output from previous run
          outtmp<-matrix(cout$output, nrow=nsteps_sub+1)
          outtmp<-outtmp[outtmp[,1]!=0 | (1:nrow(outtmp))==1,]
          
          #extract subset output from previous run
          out_spatial<-matrix(cout$output_sub, nrow=nsteps_sub+1)
          out_spatial<-out_spatial[out_spatial[,1]!=0 | (1:nrow(out_spatial))==1,]
          
          #create temporary list for re-run
          out_tmp<-list(output=outtmp, full=cout, plotdata=plotdata, output_spatial=0, sites_sub=sites_sub)
          
          #add total accumulated time to time vector
          outtmp[,1]<-outtmp[,1]+tmax_sub*(i-1)
          output_tot<-rbind(output_tot, outtmp)
          
          #store output into full lists
          out_spatial[,1]<-out_spatial[,1]+tmax_sub*(i-1)
          output_spatial_tot<-rbind(output_spatial_tot, out_spatial)
          
          #rerun to next disturbance event
          #(see script for description of variables)
          out_rerun<-rerunrun_metapopulation(out_tmp, tmax_sub, nsteps=tmax_sub, talktime=0, runtype="metapopulation", perturb=prt, sites_sub = sites_sub, habitatdestructsites=habitatdestructsites, habitatdestruct_species=habitatdestruct_species)
          cout<-out_rerun$full
          
          cout_lst[[1+i]]<-cout
        }
        #save outputs to regular names
        outnew<-output_tot
        out_spatial<-output_spatial_tot
        cout<-cout_lst
        
      } else {
        #if only one disturbance cycle occurrs, then extract information as usual
        out_spatial<-matrix(cout$output_sub, nrow=nsteps+1)
        out_spatial<-out_spatial[out_spatial[,1]!=0 | (1:nrow(out_spatial))==1,]
      }
      
    } else {
      #run c script
      #(see c script for description of variables)
      cout<-.C(runname,
               ptmax= as.double(tmax), pgridsize=as.integer(out$full$pgridsize), pnsp=as.integer(out$full$pnsp), xylim=as.integer(out$full$xylim), destroyed=as.integer(tmpdestory), spdestroy=as.integer(tmpspdestroy), #grid
               c_sptraits=as.double(out$full$c_sptraits), m_sptraits=as.double(out$full$m_sptraits), abundances=as.integer(abundances), colsites=as.integer(out$full$colsites), pncolsites=as.integer(out$full$pncolsites), #traits
               eventtimes_c=as.double(eventtimes_c), eventtimes_m=as.double(eventtimes_m), #events
               speciesid=as.integer(speciesid), #species
               output=as.double(output), pnsteps=as.integer(nsteps),
               ptalktime=as.integer(talktime),
               abundances_sub=as.integer(abundances_sub), sites_sub=as.integer(c_sites_sub), pnsites_sub=as.integer(pnsites_sub), output_sub=as.double(output_sub))
      #replace destoryed sites
      if(sum(habitatdestructsites)!=0) {
        cout$destroyed[]<-0
      }
      
      #save output
      out_spatial<-matrix(cout$output_sub, nrow=nsteps+1)
      out_spatial<-out_spatial[out_spatial[,1]!=0 | (1:nrow(out_spatial))==1,]
    }
  }
  
  ################ Clean up
  #unload c code after run this is necesary because of an aparent bug,
  #that seems to prevent R from succesfully closing down c scripts
  #with large memory allocations
  if(runtype=="metapopulation"  || (runtype=="disturbance" && nprt_cyc==1)) {
    dyn.unload("run_metapopulation_spatialsub.so")
  } else if(runtype%in%c("neutral")) {
    dyn.unload("run_neutral_metapopulation_spatialsub.so")
  }
  
  #update ouput and plotting data
  if((!(runtype %in% c("disturbance"))) || (nprt_cyc==1)) {
    outnew<-matrix(cout$output, nrow=nsteps+1)
    nsp<-ncol(outnew)-1
    ngrid<-prod(gridout$lng)
  }
  
  #calculate analytical expectation of equilibrium (for Levins model)
  ceq<-getceq(out$full$c_sptraits, out$full$m_sptraits)
  
  #save data for use in plotting
  plotdata<-list(ceq=ceq, ngrid=ngrid)
  
  #remove empty output data
  outnew<-outnew[outnew[,1]!=0 | (1:nrow(outnew))==1,]
  
  #output is identical to that from 'run_metapopulation'
  
  return(list(output=outnew, full=cout, plotdata=plotdata, output_spatial=out_spatial, sites_sub=sites_sub))
}


########################################
# Plotting functions
########################################

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

plot_metapop<-function(output, sites=0, dotot=TRUE, ylim=c(0,1)) {
  #plot time trend for results from 'run_metapopulation' function
  #sites!=0 plots spatial subset
  #dotot=TRUE means that total community abundance is also plotted
  #ylim allows setting y limits for plot
  
  #extract analytical expectation
  ceq<-output$plotdata$ceq
  
  #extract data
  if(sum(sites)==0) {
    out<-output$output
    ngrid<-output$plotdata$ngrid
  } else {
    out<-output$output_spatial
    ngrid<-length(output$sites_sub)
  }
  
  #exclude any empty time steps (for which time=0)
  sbs<-which(out[,1]>0)
  if(sbs[1]!=1) {
    #make sure to include first time step
    sbs<-c(1, sbs)
  }
  
  #plot
  if(dotot) {
    matplot(out[sbs,1], cbind(rowSums(out[sbs,-1]/ngrid), out[sbs,-1]/ngrid), type="l", xlab="time", ylab="p",
            col=1:(ncol(out)), lty=c(1, rep(1, ncol(out)-1)), lwd=2, ylim=ylim, xaxs="i", xlim=c(0, ceiling(max(out[sbs,1]))))
  } else {
    matplot(out[sbs,1], out[sbs,-1]/ngrid, type="l", xlab="time", ylab="p",
            col=2:(ncol(out)), lty=c(rep(1, ncol(out)-1)), lwd=2, ylim=ylim, xaxs="i", xlim=c(0, ceiling(max(out[sbs,1]))))
  }
  
  abline(h=c(0,1), lty=3)
  
  #add equilibria
  if(length(unique(abs(ceq)))>1) {
    abline(h=ceq,
         lty=2, col=2:ncol(out), lwd=2)
    if(dotot) {
      abline(h=sum(ceq), lwd=2, lty=2)
    }
  } else {
    abline(h=unique(abs(ceq)),
           lty=2, col=1, lwd=2)
  }
}

plot_map<-function(out, gridout, grid_sub=NULL, collst=2:(length(out$plotdata$ceq)+1)) {
  #plots 2D map of locations of individuals at last time step in 'run_metapopulation' simulation
  #out is output from simulatio
  #gridout is grid used to run simulation
  #grid_sub is output from grid subset function, if plotting of subset is desired
  #collst allows control of plotting colors
  
  #extract data
  tmp<-out$full$speciesid; tmp[tmp==out$full$pnsp]<-NA
  
  #plot locations of individuals
  plot(gridout$xpos, gridout$ypos, col=collst[tmp+1], pch=16, xlab="x position", ylab="y position", cex=0.8, xaxs="i", yaxs="i")
  
  #add segments for grid subset
  if(!is.null(grid_sub)) {
    segments(grid_sub$borders[c(1,1,2,2)]+0.5*c(-1,-1,1,1), grid_sub$borders[c(3,4,4,3)]+0.5*c(-1,1,1,-1), grid_sub$borders[c(1,2,2,1)]+0.5*c(-1,1,1,-1), grid_sub$borders[c(4,4,3,3)]+0.5*c(1,1,-1,-1), col="black", lwd=4)
    segments(grid_sub$borders[c(1,1,2,2)]+0.5*c(-1,-1,1,1), grid_sub$borders[c(3,4,4,3)]+0.5*c(-1,1,1,-1), grid_sub$borders[c(1,2,2,1)]+0.5*c(-1,1,1,-1), grid_sub$borders[c(4,4,3,3)]+0.5*c(1,1,-1,-1), col="white", lwd=1)
  }
}

rewrap_pop<-function(out, population) {
  #re-wraps output from 'run_metapopulation' back into
  #a 'population' structure.
  #Note, however, that precise timings of next events are lost when using this function,
  #and that it should therefore not be used to iterate (as otherwise all new events reset).
  #For iterated runs, 'rerunrun_metapopulation' should be used.
  
  #gets number of species
  snp<-length(population$nlst)
  
  #lists for storing abundances and positions
  spid<-NULL
  pos_sp<-NULL
  
  #extract abundances and positions
  for(i in 1:snp) {
    spid<-c(spid, rep(i, out$full$abundances[i]))
    pos_sp<-c(pos_sp, which(out$full$speciesid==(i-1)))
  }
  spid<-factor(spid, levels=1:snp)
  
  #determine which spots are occupied
  patch_occupied<-numeric(length(gridout$xpos))
  patch_occupied[pos_sp]<-1
  
  #save in population structure format
  population<-list(spid=spid,
                   nlst=unname(table(spid)),
                   clst=population$clst,
                   mlst=population$mlst,
                   radlst=population$radlst,
                   pos_sp=pos_sp,
                   patch_occupied=patch_occupied)
  
  return(population)
}


########################################
# Equilibrium analysis functions
########################################

estimate_eqreturn<-function(out, simtime=100, runtype="metapopulation", perturbsites=1:out$plotdata$ngrid, doplot=TRUE, prtb=0.1, replace_perturb=0, useeq=0, talktime=0, sites_sub=0, prt=0,  prtfrq=0) {
  #estimates rate of return to equilibrium after a small (negative) perturbation for each species in a simulated system
  
  #out is output from a run of the 'run_metapopulation' function
  #simtime is length of time to simulate after perturbation (peturbation occurrs at t=0)
  #perturbsites is vector with locations, from which individuals will be randomly selected to perturb
  #prtb indicates the stength of the perturbation (in fraction of existing population size)
  #replace_perturb indicates whether or not perturbed individuals are replaced with individuals from another species
  #useeq=0 indicates that simulated trajectory should be used to estimate divergence - otherwise, rate of return to analytical equilibrium is used
    #NOTE: if useeq!=0, then it must be a vector with one element per species, specifying the 
    #abundance of each species at equilibrium, in terms of fraction of total sites occupied
  #doplot is logical, indicating whether or not results should be plotted
  #runtype, talktime, sites_sub, prt, and prtfrq, (and also perturbsites, prtb, and replace_perturb) are as described in 'run_metapopulation'
  
  #list for storing results
  out_lst<-NULL
  
  #run simulation for each species, with perturbation
  for(i in 1:length(out$plotdata$ceq)) {
    #add perturbation for each species, sequentially
    pt<-rep(0, length(out$plotdata$ceq))
    pt[i]<-prtb
    
    #rerun model
    out_lst[[i]]<-rerunrun_metapopulation(out=out, tmax=simtime, runtype = runtype, perturb=pt, perturbsites=perturbsites, replace_perturb=replace_perturb, talktime=talktime, sites_sub = sites_sub, prt=prt, prtfrq=prtfrq)
    
    #store spatial subset results, if applicable
    if(sum(sites_sub)!=0) {
      out_lst[[i]]$tmp<-out_lst[[i]]$output
      out_lst[[i]]$output<-out_lst[[i]]$output_spatial
      out_lst[[i]]$plotdata$ngrid<-length(out_lst[[i]]$sites_sub)
    }
  }
  
  #re-run each simulation without the disturbance (to calculate divergence)
  if(useeq==0) {
    pt<-rep(0, length(out$plotdata$ceq))
    out_lst0<-rerunrun_metapopulation(out=out, tmax=simtime, runtype = runtype, perturb=pt, perturbsites=perturbsites, talktime=talktime, sites_sub = sites_sub, prt=prt, prtfrq=prtfrq)
    
    if(sum(sites_sub)!=0) {
      out_lst0$tmp<-out_lst0$output
      out_lst0$output<-out_lst0$output_spatial
      out_lst0$plotdata$ngrid<-length(out_lst0$sites_sub)
    }
  }
  
  #vectors for storing pseudo-eigenvalues (i.e. rate of return to equilibrium trajectory for each species)
  eigenlst<-matrix(ncol=length(out_lst), nrow=(simtime-1), data=NA)
  eigenlst_sd<-matrix(ncol=length(out_lst), nrow=(simtime-1), data=NA)
  
  #calculate distance between trajectories
  for(sppos in 1:length(out_lst)) {
    if(sum(useeq)==0) {
      #distance between simulations
      pred_diff<-abs(out_lst0$output[,sppos+1]-out_lst[[sppos]]$output[,sppos+1])/(out$plotdata$ngrid)
    } else {
      #distance between perturbed simulation and analytical expectation
      pred_diff<-abs(useeq[sppos]-(out_lst[[sppos]]$output[,sppos+1])/(out$plotdata$ngrid))
    }
    
    #transform distances into growth rate for each time step
    logdiftmp<-log(pred_diff[-1]/pred_diff[-length(pred_diff)])
    sbs<-is.finite(logdiftmp)
    
    #store mean growth rate for each time scale (i.e. 0 to tmax)
    eigenlst[1:(length(logdiftmp)),sppos][sbs]<-cumsum(logdiftmp[sbs])/(1:sum(sbs))
    eigenlst_sd[1:(length(logdiftmp)),sppos][sbs]<-sqrt(cumsum(logdiftmp[sbs]^2)/(1:sum(sbs))-(eigenlst[1:(length(logdiftmp)),sppos][sbs])^2)
  }
  
  #plot results
  if(doplot) {
    #multiply by time scale to get total growth over time interval
    eigenlst_tmp<-eigenlst*(1:nrow(eigenlst))
    eigenlst_sd_tmp<-eigenlst_sd*(1:nrow(eigenlst_sd))
    
    #plot mean +/- one sd
    matplot(c(1,nrow(eigenlst_tmp)), range(c(0, eigenlst_tmp), na.rm=T), col=1:ncol(eigenlst_tmp)+1, lty=1, xlab="time span", ylab=expression(paste(lambda, "t")), type="n"); abline(h=0, lty=3)
    for(i in 1:ncol(eigenlst_tmp)) {
      sbs<-which(!is.na(eigenlst_tmp[,i]+eigenlst_sd_tmp[,i]))
      if(sum(sbs)>0) {
        polygon(c(1:nrow(eigenlst_tmp[sbs,]), rev(1:nrow(eigenlst_tmp[sbs,]))),
              c(eigenlst_tmp[sbs,i]+eigenlst_sd_tmp[sbs,i], rev(eigenlst_tmp[sbs,i]-eigenlst_sd_tmp[sbs,i])), col=adjustcolor(i+1, alpha.f = 0.1), border=NA)
      }
    }
    
    matlines(1:nrow(eigenlst_tmp), eigenlst_tmp, col=1:ncol(eigenlst_tmp)+1, lty=1, lwd=2)
  }
  
  #if spatial subset was used, unpack results
  if(sum(sites_sub)!=0) {
    #undo shifting of positions
    for(i in 1:length(out_lst)) {
      out_lst[[i]]$output<-out_lst[[i]]$tmp
      out_lst[[i]]$plotdata$ngrid<-prod(out$plotdata$ngrid)
      out_lst[[i]]$tmp<-NULL
    }
    
    if(exists("out_lst0")) {
      out_lst0$output<-out_lst0$tmp
      out_lst0$plotdata$ngrid<-prod(out$plotdata$ngrid)
      out_lst0$tmp<-NULL
    }
  }
  
  #eigenlist is matrix of return rates
  #eigenlst_sd shows standard deviation for elements in eigenlist
  #out_lst includes all simulated peturbed trajectories
  #out_lst0 includes all simulated non-perturbed trajectories (if applicable)
  return(list(eigenlst=eigenlst, eigenlst_sd=eigenlst_sd, out_lst=out_lst, out_lst0=out_lst0))
}


estimate_rarereturn<-function(out, simtime=100, burnin=100, runtype="metapopulation", perturbsites=1:out$plotdata$ngrid, doplot=TRUE, talktime=0, add_amount=0.05, sites_sub=0, prt=0,  prtfrq=0) {
  #estimates growth rate when rare for each spieces in a simulation model
  
  #out is result from 'run_metapopulation'
  #simtime is time to simulate after re-introduction
  #burnin is time to simulate after population is reduced to zero (to allow rest of the community to reach equilibrium)
  #perturbsites is vector with locations, from which individuals will be randomly selected to perturb
  #doplot is logical, indicating whether or not results should be plotted
  #add_amount is fraction of total sites that should be populated by species during re-introduction
  #runtype, talktime, sites_sub, prt,and prtfrq are as described in 'run_metapopulation'
  
  #lists for output
  out_lst<-NULL
  out0_lst<-NULL
  
  #for each species drop to zero abundance and simulate for burnin steps
  for(i in 1:length(out$plotdata$ceq)) {
    #add perturbation to each species, sequentially
    pt<-rep(0, length(out$plotdata$ceq))
    pt[i]<-1
    
    #rerun with species excluded
    tmp<-rerunrun_metapopulation(out=out, tmax=burnin, runtype = runtype, perturb=pt, perturbsites=perturbsites, talktime=0, sites_sub=sites_sub, prt=prt, prtfrq=prtfrq, habitatdestructsites = perturbsites, habitatdestruct_species = pt)
    out0_lst[[i]]<-tmp
    
    #now, add species back in
    at<-rep(0, length(out$plotdata$ceq))
    at[i]<-add_amount
    
    #if simulation is a disturbance simulation, update briefly to standard metapopulation
    #this will prevent a disturbance from wiping out added individuals right after introduction
    if(runtype %in% c("disturbance")) {
      tmptyp<-"metapopulation"
      tmp<-rerunrun_metapopulation(out=tmp, tmax=0, runtype = tmptyp, perturb=prt, perturbsites=1:tmp$plotdata$ngrid, talktime=0, sites_sub=sites_sub)
    }
    
    #run with species added back in
    out_lst[[i]]<-rerunrun_metapopulation(out=tmp, tmax=simtime, runtype = runtype, addn=at,  perturb=pt, perturbsites=perturbsites, addsites=perturbsites, talktime=0, sites_sub=sites_sub, prt=prt, prtfrq=prtfrq)
  }
  
  #organize spatial subset data
  if(sum(sites_sub)>0) {
    for(i in 1:length(out$plotdata$ceq)) {
      out_lst[[i]]$tmp<-out_lst[[i]]$output
      out_lst[[i]]$output<-out_lst[[i]]$output_spatial
      out_lst[[i]]$plotdata$ngrid<-length(out_lst[[i]]$sites_sub)
    }
  }
  
  #matrices for storing growth rate data
  grwrare<-matrix(ncol=length(out_lst), nrow=(simtime-1))
  grwrare_sd<-matrix(ncol=length(out_lst), nrow=(simtime-1))
  
  #calculate growth rate for each species
  for(sppos in 1:length(out_lst)) {
    pred_grw<-out_lst[[sppos]]$output[,sppos+1]/out$plotdata$ngrid
    
    #mean growth rate over time interval
    logdiftmp<-log(pred_grw[-1]/pred_grw[-length(pred_grw)])
    grwrare[,sppos][1:length(logdiftmp)]<-cumsum(logdiftmp)/(1:(length(pred_grw)-1))
    
    grwrare_sd[,sppos][1:length(logdiftmp)]<-sqrt(cumsum(logdiftmp^2)/(1:(length(pred_grw)-1))-(grwrare[,sppos]^2)[1:length(logdiftmp)])
  }
  
  if(doplot) {
    #plot results, if desired
    
    #covert into total growth observed over time period
    grwrare_tmp<-grwrare*(1:nrow(grwrare))
    grwrare_sd_tmp<-grwrare_sd*(1:nrow(grwrare_sd))
    
    matplot(c(1, nrow(grwrare_tmp)), range(c(grwrare_tmp[is.finite(grwrare_tmp)], 0, na.rm=T)), col=1:ncol(grwrare_tmp)+1, lty=1, xlab="time span", ylab=expression(paste("r"[0], "t")), type="n"); abline(h=0, lty=3)
    for(i in 1:ncol(grwrare_tmp)) {
      sbs<-is.finite(grwrare_tmp[,i]+grwrare_sd_tmp[,i])
      polygon(c(1:nrow(grwrare_tmp[sbs,]), rev(1:nrow(grwrare_tmp[sbs,]))),
              c(grwrare_tmp[sbs,i]+grwrare_sd_tmp[sbs,i], rev(grwrare_tmp[sbs,i]-grwrare_sd_tmp[sbs,i])), col=adjustcolor(i+1, alpha.f = 0.1), border=NA)
    }
    
    matlines(1:nrow(grwrare_tmp), grwrare_tmp, col=(1:ncol(grwrare_tmp))+1, lty=1, lwd=2)
  }
  
  #if spatial subset, rearrange output
  if(sum(sites_sub)!=0) {
    #undo shifting of positions
    for(i in 1:length(out_lst)) {
      out_lst[[i]]$output<-out_lst[[i]]$tmp
      out_lst[[i]]$plotdata$ngrid<-prod(out$plotdata$ngrid)
      out_lst[[i]]$tmp<-NULL
    }
  }
  
  #grwrare is matrix of growth rates
  #grwrare_sd is matrix with standard deviations for each element in grwrare
  #out_lst is list with all simulations of growth after reintroduction
  #out0_lst is list with all simulations after population is excluded
  return(list(grwrare=grwrare, grwrare_sd=grwrare_sd, out_lst=out_lst, out0_lst=out0_lst))
}


estimate_invar<-function(out, E=0, burnin=0, Luse=0, laglst=0, niter=0, doplot=TRUE, sites_sub=0) {
  #Calculates coefficient of variation comparing simulation dynamics to predictions based on subsets of past observations
  #Basic idea is to use a historical library of observations to predict observations n timesteps into the future
  #Theoretically, if the historical library contains enough points to define an equilibrium state, and if the system is at equilibrium,
    #then predictive power should be high regardless of time lag between training and testing set.
  
  #out is result from 'run_metapopulation' function
  #E is embedding dimension to use for simplex algorithm. If zero, then "best" E is automatically selected
  #burnin is number of initial time steps to exclude
  #Luse is list of library lengths to test. Defaults to zero, which attempts to choose sensible limits.
    #These are used to try to find a window size that provides good predictions of system dynamics
  #laglst is list of temporal lags to test. Default is zero, which automatically selects lag sizes based on training and total library sizes
  #niter is number of different training windows to test - when set to zero, all possible windows are tested
  #doplot is binary variable that determines whether or not results are plotted
  #sites_sub includes subset of sites to use. Defaults to zero, which includes all sites.
  
  if(sum(E)==0) {
    tmp<-getE(out)
    E<-tmp$Eout
  }
  
  if(sum(Luse)==0) {
    Luse<-floor((seq((30), min(c(ceiling(nrow(out$output)/5), nrow(out$output)-burnin)), length=10)))
  }
  
  #exclude windows that are too long to test any reasonable lag lengths
  Luse<-Luse[Luse<((nrow(out$output)-burnin-max(E))/2)]
  
  #if only one E is given, expand to include one entry per species
  if(length(E)==1) {
    E<-rep(E, length(out$plotdata$ceq))
  }
  
  #extract data if spatial subset is being analyzed
  if(sum(sites_sub)>0) {
    out$tmp_data<-out$output
    out$tmp_lng<-out$plotdata$ngrid
    
    out$output<-out$output_spatial
    out$plotdata$ngrid<-length(out$sites_sub)
  }
  
  #lists for storing results
  pdL_list<-NULL
  pdlag_list<-NULL
  
  for(i in 1:length(out$plotdata$ceq)) {
    #for each species, first identify best window size to use
    pdL_list[[i]]<-predict_vs_L(out$output[,i+1], E=E[i], burnin=burnin, Luse=Luse, niter=niter, doplot=FALSE)
    Lusetmp<-pdL_list[[i]]$Lmin
    
    #select appropriate lag list, or test that user provided bounds are reasonable
    if(sum(laglst)==0) {
      laglst_use<-c(floor((seq((0), ((nrow(out$output)-burnin-Lusetmp)), length=20))))
    } else {
      laglst_use<-laglst
      laglst_use<-laglst_use[laglst_use<=((nrow(out$output)-burnin-Lusetmp))]
    }
    
    #use selected window size and lags to test predictive power of training and testing sets
    pdlag_list[[i]]<-test_predict_tlag(out$output[,i+1], Luse=Lusetmp, E=E[i], burnin=burnin, laglst=laglst_use, niter=niter, doplot=FALSE)
  }
  
  if(doplot) {
    #if plotting is desired, plot results
    
    #find range to use for plotting axes
    mx<-0
    tl<-0
    for(i in 1:length(pdlag_list)) {
      tmpmx<-max(range(pdlag_list[[i]]$CVest[,2:4], na.rm=T))
      
      if(is.finite(tmpmx)) {
        mx<-max(c(mx, tmpmx, na.rm=T))
      }
      tl<-max(c(tl, pdlag_list[[i]]$laglst), na.rm=T)
    }
    
    #plot mean CV vs. lag
    plot(c(0, tl), c(0, max(c(mx), na.rm=T)), xlab="time lag", ylab="CV", type="n", xaxs="i")
    abline(h=0, lty=3)
    
    #add CI's based on standard deviation of among-window variability
    #note that for very long lags, only a few windows exist, and standard deviations are therefore not
    #indicative of true variability
    for(i in 1:length(pdlag_list)) {
      sbs<-which(is.finite(rowSums(pdlag_list[[i]]$CVest)))
      polygon(c(pdlag_list[[i]]$laglst[sbs], rev(pdlag_list[[i]]$laglst[sbs])),
              c(pdlag_list[[i]]$CVest[sbs,2], rev(pdlag_list[[i]]$CVest[sbs,4])),
              col=adjustcolor(i+1, 0.2), border=NA)
      lines(pdlag_list[[i]]$laglst, pdlag_list[[i]]$CVest[,3], col=i+1, lwd=2, lty=1)
    }
  }

  #pdL_list includees results from window size selection
  #pdlag_list incudes results for CV for each lag
  return(list(pdL_list=pdL_list, pdlag_list=pdlag_list))
}




########################################
# Invariance helper functions
########################################

getE<-function(out, Elst=2:10, doplot=FALSE, sites_sub=0) {
  #finds optimal embedding dimension based on predictive power (rho)
  
  #extract data if spatial subset is being analyzed
  if(sum(sites_sub)>0) {
    out$output<-out$output_spatial
  }
  
  #run simplex algorithm for each species
  Eout<-numeric(length(out$plotdata$ceq))
  simplout<-NULL
  for(i in 1:length(Eout)) {
    simplout[[i]]<-suppressWarnings(simplex(out$output[,i+1], E=Elst))
    
    #selects "best" E, which is the smallest E for which the difference between rho and the maximum
    #rho observed is smaller than 10% of the total range observed across all rho values.
    #this is an ad-hoc method meant to help find the "saturation" point
    Eout[i]<-Elst[min(which((max(simplout[[i]]$rho)-simplout[[i]]$rho)/diff(range(simplout[[i]]$rho))<0.1))]
  }
  
  if(doplot) {
    #if plotting is desired, plot actual relationshisp for rho vs. E for each species
    rng<-c(NA, NA)
    for(i in 1:length(simplout)) {
      tmp<-range(simplout[[i]]$rho, na.rm=T)
      rng[1]<-min(c(tmp[1], rng[1]), na.rm=T)
      rng[2]<-max(c(tmp[2], rng[2]), na.rm=T)
    }
    
    plot(range(Elst), rng, xlab="E", ylab=expression(paste(rho)), type="n")
    for(i in 1:length(simplout)) {
      lines(simplout[[i]]$E, simplout[[i]]$rho, lty=1, lwd=2, col=i+1)
      abline(v=Eout[i], lty=3, col=i+1)
    }
  }
  
  #Eout is a vector with the best E for each species
  #simplout includes all rho for all E and species
  return(list(Eout=Eout, simplout=simplout))
}


predict_vs_L<-function(outcol, E=1, burnin=0, Luse=0, niter=0, doplot=TRUE) {
  #calculates predictive power vs. library length across window sizes.
  #is meant to identify an optimal window size
  
  #outcol is a vector of observations, for which predictions are to be made
  #E, burnin, Luse, niter, and doplot are as described in 'estimate_invar'
  
  #if Luse=0, select a reasonable range of window sizes
  if(sum(Luse)==0) {
    Luse<-floor((seq((30), min(c(ceiling(length(outcol)/5), length(outcol)-burnin)), length=10)))
  }
  
  #variable for saving results
  predmat<-NULL
  
  for(i in 1:length(Luse)) {
    #for each window size, find all possible windows to use for prediction
    window_starts<-seq((burnin+1), (length(outcol)-(Luse[i]-1)), by=Luse[i])
    
    #select a subset of window starts to use, or use all if niter==0
    if(niter!=0 & niter<length(window_starts)) {
      window_use<-sample(window_starts, niter, replace = F)
    } else {
      window_use<-window_starts
    }
    
    #matrix for storing prediction results
    predmat_tmp<-matrix(ncol=5, nrow=length(window_use))
    predmat_tmp[,1]<-Luse[i]
    
    #calculate predictive power using each window
    for(j in 1:length(window_use)) {
      #for each window, apply the simplex algorithm to predict the whole time series using leave one out cross validation
      tmp<-suppressWarnings(simplex(outcol, E=E, lib=c(window_use[j], (window_use[j]+Luse[i]-1)), pred=c(1+burnin, length(outcol)), stats_only = FALSE))
      
      #store summary statistics
      predmat_tmp[j,2]<-tmp$rho
      predmat_tmp[j,3]<-tmp$rmse
      
      predmat_tmp[j,4]<-mean(tmp$model_output[[1]]$pred, na.rm=T)
      predmat_tmp[j,5]<-tmp$num_pred
    }
    
    #concatenate results
    predmat<-rbind(predmat, predmat_tmp)
  }
  
  #column names for output
  colnames(predmat)<-c("Luse", "rho", "rmse", "mean", "n")
  
  #calculate mean CV, and confidence intervals, for each window size
  #note, CV is simply root mean square error/mean observation
  CVest<-t(matrix(nrow=5, unlist(tapply(predmat[,"rmse"]/predmat[,"mean"], predmat[,"Luse"], function(x) quantile(x[is.finite(x)], c(0.025, pnorm(-1, 0, 1), 0.5, pnorm(1, 0, 1), 0.975), na.rm=T)))))
  
  #identify minimum library length that gives us 50% of the best predictive power
  Lmin<-Luse[min(which((CVest[,5]-min(CVest[,5],na.rm=T))/(CVest[,5])<0.5),na.rm=T)]
  
  if(doplot) {
    #if plotting is desired, plot results (widow size vs. CV and CI's)
    
    plot(Luse, CVest[,3], type="n", xlab="library length", ylab="CV", ylim=c(0, max(c(1, CVest[,3]))), xlim=c(min(c(0, Luse)), max(Luse)), xaxs="i")
    polygon(c(Luse, rev(Luse)), c(CVest[,1], rev(CVest[,5])), col=adjustcolor("black", alpha.f = 0.2), border=NA)
    polygon(c(Luse, rev(Luse)), c(CVest[,2], rev(CVest[,4])), col=adjustcolor("black", alpha.f = 0.2), border=NA)
    lines(Luse, CVest[,3], lwd=2)
    
    abline(v=c(min(Luse), Lmin), h=c(0, 1), lty=3)
  }
  
  #Lmin is minimum library length that gives 50% of best predictive power
  #CVest is estimated CV (0.025% quantile, mean-sd, mean, mean+sd, and 0.975% quantile) for each window size
  #Luse is list of window lengths tested
  return(list(Lmin=Lmin, CVest=CVest, predmat=predmat, Luse=Luse))
}


test_predict_tlag<-function(outcol, Luse, E=1, burnin=0, laglst=0, niter=0, doplot=TRUE) {
  #calculates predictive power vs. time distance between training and testing sets
  #given a particular time window size, Luse.

  #outcol is a vector of observations, for which predictions are to be made
  #Luse should be calculated using the 'predict_vs_L' function
  #E, burnin, laglst, niter, and doplot are as described in 'estimate_invar'
  
  #variable for storing results
  predmat<-NULL
  
  #if laglst is set to zero, choose a sensible set of lags to test
  if(sum(laglst)==0) {
    laglst<-c(floor((seq((0), ((length(outcol)-burnin-Luse)), length=20))))
  }
  
  #Cycle through lag lengths
  for(i in 1:length(laglst)) {
    #Find all possible comparisons between windows
    p1full<-seq(burnin+1, length(outcol)-(Luse-1)-laglst[i], by=Luse)
    
    #either choose subset of comparisons, or use all if niter=0
    if(niter!=0 & niter<length(p1full)) {
      p1<-sample(p1full, niter, replace = F)
    } else {
      p1<-p1full
    }
    
    #position of test set
    p2<-p1+laglst[i]
    
    #matrix for storing results
    predmat_tmp<-matrix(ncol=5, nrow=length(p1))
    predmat_tmp[,1]<-laglst[i]
    
    for(j in 1:length(p1)) {
      #for each training and test set, run simplex algorithm to predict p2 using data from p1
      tmplib<-c(p1[j], (p1[j]+Luse-1))
      tmppred<-c(p2[j], (p2[j]+Luse-1))
      
      tmp<-suppressWarnings(simplex(outcol,E=E, lib=tmplib, pred=tmppred, stats_only = FALSE))
      
      #extract summary statistics
      predmat_tmp[j,2]<-tmp$rho
      predmat_tmp[j,3]<-tmp$rmse
      
      predmat_tmp[j,4]<-mean(tmp$model_output[[1]]$pred, na.rm=T)
      predmat_tmp[j,5]<-tmp$num_pred
    }
    
    #concatenate results
    predmat<-rbind(predmat, predmat_tmp)
  }
  
  #names of result columns
  colnames(predmat)<-c("Luse", "rho", "rmse", "mean", "n")
  
  #calculate mean CV, and confidence intervals, for each window size
  #note, CV is simply root mean square error/mean observation
  CVest<-t(matrix(nrow=5, unlist(tapply(predmat[,"rmse"]/predmat[,"mean"], predmat[,"Luse"], function(x) quantile(x[is.finite(x)], c(0.025, pnorm(-1, 0, 1), 0.5, pnorm(1, 0, 1), 0.975), na.rm=T)))))
  
  if(doplot) {
    #if plotting is desired, plot results (widow size vs. CV and CI's)
    
    plot(laglst, CVest[,3], type="n", xlab="time lag", ylab="CV", ylim=c(0, max(c(1, CVest[,3]))), xlim=c(min(c(0, laglst)), max(laglst)), xaxs="i")
    polygon(c(laglst, rev(laglst)), c(CVest[,1], rev(CVest[,5])), col=adjustcolor("black", alpha.f = 0.2), border=NA)
    polygon(c(laglst, rev(laglst)), c(CVest[,2], rev(CVest[,4])), col=adjustcolor("black", alpha.f = 0.2), border=NA)
    lines(laglst, CVest[,3], lwd=2)
    
    abline(v=c(min(laglst)), h=c(0, 1), lty=3)
  }
  
  #CVest is estimated CV (0.025% quantile, mean-sd, mean, mean+sd, and 0.975% quantile) for each window size
  #predmat is full outputs for all lags and windows tested
  #laglst is vector of lag lengths tested
  #Luse is window length used
  return(list(CVest=CVest, predmat=predmat, laglst=laglst, Luse=Luse))
}


########################################
# Parallel functions
########################################

runpar<-function(...) {
  #Function to automate stability tests across models
  #Is embarrassingly parallel, and can be used to iterate runs
  
  require(rEDM)
  
  #run simulations
  if(FALSE) {
  #if(grid_sub$frac_tot==1) {   #global sclae
    #run levins
    out_meta<-run_metapopulation(tmax=tmax, gridout = gridout, population = population_meta, talktime = 0)
    
    getEmeta<-getE(out_meta, Elst = 2:10)
    E_meta<-getEmeta$Eout
    E_meta[E_meta<minE]<-minE
    
    eig_meta1<-estimate_eqreturn(out_meta, simtime=simtime, runtype="metapopulation", replace_perturb = 1, talktime=0, prtb=ptb, doplot=FALSE)
    eig_meta2<-estimate_eqreturn(out_meta, simtime=simtime, runtype="metapopulation", replace_perturb = 0, talktime=0, prtb=ptb, doplot=FALSE)
    r0_meta<-estimate_rarereturn(out_meta, simtime=simtime, burnin=burnin, runtype="metapopulation", doplot=FALSE)
    invar_meta<-estimate_invar(out_meta, E=E_meta, burnin=invarburn, doplot=FALSE, laglst = lglst)
    
    #run dist
    out_dist<-run_metapopulation(tmax=tmax, gridout = gridout, population = population_dist, talktime = 0, runtype = "disturbance", prt = distlst, prtfrq = distfrq)
    
    getEdist<-getE(out_dist, Elst = 2:10)
    E_dist<-getEdist$Eout
    E_dist[E_dist<minE]<-minE
    
    eig_dist1<-estimate_eqreturn(out_dist, simtime=simtime, runtype="disturbance", replace_perturb = 1, talktime=0, prtb=ptb, doplot=FALSE, prt = distlst, prtfrq = distfrq)
    eig_dist2<-estimate_eqreturn(out_dist, simtime=simtime, runtype="disturbance", replace_perturb = 0, talktime=0, prtb=ptb, doplot=FALSE, prt = distlst, prtfrq = distfrq)
    r0_dist<-estimate_rarereturn(out_dist, simtime=simtime, burnin=burnin, runtype="disturbance", doplot=FALSE, prt = distlst, prtfrq = distfrq)
    invar_dist<-estimate_invar(out_dist, E=E_dist, burnin=invarburn, doplot=FALSE, laglst = lglst)
    
    #run netural
    out_neut<-run_metapopulation(tmax=tmax, gridout = gridout, population = population_neut, talktime = 0, runtype = "neutral")
    
    getEneut<-getE(out_neut, Elst = 2:10)
    E_neut<-getEneut$Eout
    E_neut[E_neut<minE]<-minE
    
    eig_neut1<-estimate_eqreturn(out_neut, simtime=simtime, runtype="neutral", replace_perturb = 1, talktime=0, prtb=ptb, doplot=FALSE)
    eig_neut2<-estimate_eqreturn(out_neut, simtime=simtime, runtype="neutral", replace_perturb = 0, talktime=0, prtb=ptb, doplot=FALSE)
    r0_neut<-estimate_rarereturn(out_neut, simtime=simtime, burnin=burnin, runtype="neutral", doplot=FALSE)
    invar_neut<-estimate_invar(out_neut, E=E_neut, burnin=invarburn, doplot=FALSE, laglst = lglst)
    
  } else { #spatial subset
    #run levins
    out_meta<-run_metapopulation(tmax=tmax, gridout = gridout, population = population_meta, talktime = 0, runtype = "metapopulation_spatial", sites_sub = grid_sub$sites)
    
    getEmeta<-getE(out_meta, Elst = 2:10, sites_sub = grid_sub$sites)
    E_meta<-getEmeta$Eout
    E_meta[E_meta<minE]<-minE
    
    eig_meta1<-estimate_eqreturn(out_meta, simtime=simtime, runtype="metapopulation_spatial", replace_perturb = 1, talktime=0, prtb=ptb, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites, doplot=FALSE)
    eig_meta2<-estimate_eqreturn(out_meta, simtime=simtime, runtype="metapopulation_spatial", replace_perturb = 0, talktime=0, prtb=ptb, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites, doplot=FALSE)
    r0_meta<-estimate_rarereturn(out_meta, simtime=simtime, burnin=burnin, runtype="metapopulation_spatial", perturbsites = grid_sub$sites, sites_sub = grid_sub$sites, doplot=FALSE)
    invar_meta<-estimate_invar(out_meta, E=E_meta, burnin=invarburn, doplot=FALSE, sites_sub = grid_sub$sites, laglst = lglst)
    
    #run disturbance
    out_dist<-run_metapopulation(tmax=tmax, gridout = gridout, population = population_dist, talktime = 0, runtype = "disturbance_spatial", sites_sub = grid_sub$sites, prt = distlst, prtfrq = distfrq)
    
    getEdist<-getE(out_dist, Elst = 2:10, sites_sub = grid_sub$sites)
    E_dist<-getEdist$Eout
    E_dist[E_dist<minE]<-minE
    
    eig_dist1<-estimate_eqreturn(out_dist, simtime=simtime, runtype="disturbance_spatial", replace_perturb = 1, talktime=0, prtb=ptb, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites, doplot=FALSE, prt = distlst, prtfrq = distfrq)
    eig_dist2<-estimate_eqreturn(out_dist, simtime=simtime, runtype="disturbance_spatial", replace_perturb = 0, talktime=0, prtb=ptb, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites, doplot=FALSE, prt = distlst, prtfrq = distfrq)
    r0_dist<-estimate_rarereturn(out_dist, simtime=simtime, burnin=burnin, runtype="disturbance_spatial", perturbsites = grid_sub$sites, sites_sub = grid_sub$sites, doplot=FALSE, prt = distlst, prtfrq = distfrq)
    invar_dist<-estimate_invar(out_dist, E=E_dist, burnin=invarburn, doplot=FALSE, sites_sub = grid_sub$sites, laglst = lglst)
    
    
    #run neutral
    out_neut<-run_metapopulation(tmax=tmax, gridout = gridout, population = population_neut, talktime = 0, runtype = "neutral_spatial", sites_sub = grid_sub$sites)
    
    getEneut<-getE(out_neut, Elst = 2:10)
    E_neut<-getEneut$Eout
    E_neut[E_neut<minE]<-minE
    
    eig_neut1<-estimate_eqreturn(out_neut, simtime=simtime, runtype="neutral_spatial", replace_perturb = 1, talktime=0, prtb=ptb, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites, doplot=FALSE)
    eig_neut2<-estimate_eqreturn(out_neut, simtime=simtime, runtype="neutral_spatial", replace_perturb = 0, talktime=0, prtb=ptb, sites_sub = grid_sub$sites, perturbsites = grid_sub$sites, doplot=FALSE)
    r0_neut<-estimate_rarereturn(out_neut, simtime=simtime, burnin=burnin, runtype="neutral_spatial", perturbsites = grid_sub$sites, sites_sub = grid_sub$sites, doplot=FALSE)
    invar_neut<-estimate_invar(out_neut, E=E_neut, burnin=invarburn, doplot=FALSE, sites_sub = grid_sub$sites, laglst = lglst)
  }
  
  #make output matrix
  nsp<-length(clst_meta)
  neq<-nrow(eig_meta1$eigenlst)
  nlg<-length(invar_meta$pdlag_list[[1]]$laglst)
  mlng<-max(c(neq, nlg))
  matout<-matrix(nrow=mlng, ncol=nsp*15, data=NA)
  
  #collect and output data
  matout[1:neq,(1:nsp)+nsp*(0)]<-eig_meta1$eigenlst
  matout[1:neq,(1:nsp)+nsp*(1)]<-eig_meta2$eigenlst
  matout[1:neq,(1:nsp)+nsp*(2)]<-r0_meta$grwrare
  
  matout[1:neq,(1:nsp)+nsp*(3)]<-eig_dist1$eigenlst
  matout[1:neq,(1:nsp)+nsp*(4)]<-eig_dist2$eigenlst
  matout[1:neq,(1:nsp)+nsp*(5)]<-r0_dist$grwrare
  
  matout[1:neq,(1:nsp)+nsp*(6)]<-eig_neut1$eigenlst
  matout[1:neq,(1:nsp)+nsp*(7)]<-eig_neut2$eigenlst
  matout[1:neq,(1:nsp)+nsp*(8)]<-r0_neut$grwrare
  
  for(ii in 1:length(invar_meta$pdlag_list)) {
    matout[1:length(invar_meta$pdlag_list[[ii]]$laglst),nsp*9+(ii)+5*(ii-1)]<-invar_meta$pdlag_list[[ii]]$laglst
    matout[1:length(invar_meta$pdlag_list[[ii]]$laglst),nsp*9+(ii)+5*(ii-1)+1]<-invar_meta$pdlag_list[[ii]]$CVest[,3]
    
    matout[1:length(invar_dist$pdlag_list[[ii]]$laglst),nsp*9+(ii)+5*(ii-1)+2]<-invar_dist$pdlag_list[[ii]]$laglst
    matout[1:length(invar_dist$pdlag_list[[ii]]$laglst),nsp*9+(ii)+5*(ii-1)+3]<-invar_dist$pdlag_list[[ii]]$CVest[,3]
    
    matout[1:length(invar_neut$pdlag_list[[ii]]$laglst),nsp*9+(ii)+5*(ii-1)+4]<-invar_neut$pdlag_list[[ii]]$laglst
    matout[1:length(invar_neut$pdlag_list[[ii]]$laglst),nsp*9+(ii)+5*(ii-1)+5]<-invar_neut$pdlag_list[[ii]]$CVest[,3]
  }
  
  return(t(matout))
}


########################################
# Example of running functions
########################################

if(FALSE) {
  #make grid and population
  gridout<-makegrid(xlng = 100, ylng = 100)
  population<-populate(gridout, nlst = 0.1, clst = c(0.12, 0.8), mlst = rep(0.1, 4), radlst = Inf)
  
  #run initial simulation
  out<-run_metapopulation(tmax=200, nsteps = 1000, gridout, population, talktime = 0)
  plot_metapop(out)
  
  #rerun simulation
  out2<-rerunrun_metapopulation(out=out, tmax=50, talktime = 0)
  plot_metapop(out2)
  
  #calculate stability metrics
  eqret<-estimate_eqreturn(out)
  rareret<-estimate_rarereturn(out)
  invar<-estimate_invar(out)
}

