#load 
e014dat<-read.csv("~/Dropbox/Projects/old/019_Dave habitat loss/src/data/e014_abund_data.csv")

#Get c values
m<-0.1
clst<-numeric(nrow(e014dat))
plst<-numeric(nrow(e014dat))

for(i in 1:nrow(e014dat)) {
  preq<-e014dat$cumSum[i]
  clst[i]<-m/(1-preq)
  
  plst[i]<-preq-sum(plst)
}

gridout<-makegrid(xlng = 1000, ylng = 1000)
n<-200
population<-populate(gridout, nlst = round(plst[1:n]*prod(gridout$lng)), clst = clst[1:n], mlst = rep(m, n), radlst = 100)

x1<-Sys.time()
out<-run_metapopulation(tmax=100, nsteps = 100, gridout, population, talktime = 1)
x2<-Sys.time()
x2-x1


sum(out$output[nrow(out$output),-1]>0)
