error
rm(list=ls())
setwd("~/Dropbox/Projects/032_Coexistence_mechanisms/src/levins_metapopulation/")
require(data.table)
require(RColorBrewer)
source("~/Dropbox/Rfunctions/figure_functions.R")



if(FALSE) {
d001<-read.table("data/e001_Plant aboveground biomass data.txt", header=T, sep="\t")
d002<-read.table("data/e002_Plant aboveground biomass data.txt", header=T, sep="\t")
d014<-data.frame(fread("data/e014_long_covars_161017.csv"))
d120<-read.table("data/e120_Plant aboveground biomass data.txt", header=T, sep="\t", fill = TRUE)
dnut<-data.frame(fread("../../../035_Yann_Species loss/src/data/comb-by-plot-clim-soil-diversity-04-Dec-2017.csv"))



#e001/e002
d002$chr<-as.character(d002$Sampling.date.mm.dd.yyyy.)
d002$year<-as.numeric(matrix(nrow=3, unlist(strsplit(d002$chr, "/", fixed=T)))[3,])
d002<-d002[d002$Subplot=="Whole",]

d001<-d001[d001$Field!="D",]

i001<-unique(data.frame(year=c(d001$Year, d002$year), exp=c(d001$Exp, rep(2, nrow(d002))),
                        field=c(d001$Field, d002$Field), plot=c(d001$Plot, d002$plot.number)))
i001$index<-paste(i001$exp, i001$field, i001$plot)



cumsum(rev(table(colSums(table(i001$year, i001$index)))))
cumsum(rev(table(tapply(i001$year, i001$index, function(x) diff(range(x))))))

#e014
i014<-unique(data.frame(year=d014$year, field=d014$field, transect=d014$transect, plot=d014$plot))
i014$index<-paste(i014$field, i014$transect, i014$plot)

cumsum(rev(table(colSums(table(i014$year, i014$index)))))
cumsum(rev(table(tapply(i014$year, i014$index, function(x) diff(range(x))))))

#e120
i120<-unique(data.frame(year=d120$Year, plot=d120$Plot, strip=d120$Strip))
i120$index<-paste(i120$plot, i120$strip)

cumsum(rev(table(colSums(table(i120$year, i120$index)))))
cumsum(rev(table(tapply(i120$year, i120$plot, function(x) diff(range(x))))))

#xcoord<-c(70, 323, 64, 316)
#ycoord<-c(604, 608, 866, 872)
#width: 252-253
#height: 262-264

#Nutnet
inut<-unique(data.frame(year=dnut$year, site=dnut$site_code, block=dnut$block, plot=dnut$plot))
inut$index<-paste(inut$site, inut$block, inut$plot)

cumsum(rev(table(colSums(table(inut$year, inut$index)))))
cumsum(rev(table(tapply(inut$year, inut$index, function(x) diff(range(x))))))




#Jena:

#"Small"
#2002-2008: 82*(3.5^2)
#"Large"
#2002-2010: 96*(20*20)
#2010-present: 96*(104.75)

#0.2-0.5 m 2-4

}


#Make plots:
plotdat<-read.table("data/realworld_extents.csv", header=T, sep=";")
collst<-(c("lightblue", "darkblue", "blue", "red", "darkgreen"))
#collst<-c((brewer.pal(4, "Blues"))[-1], "red", "darkgreen")
anglst<-seq(0, 360, by=30)

explst<-sort(unique(plotdat$experiment))
typelst<-sort(unique(plotdat$type))

plotdat$summed_sample_size<-plotdat$spatial_plots*plotdat$sample_size_m2
plotdat$summed_plot_size<-plotdat$spatial_plots*plotdat$plot_size


sbs1<-(plotdat$type==typelst[1])
sbs2<-(plotdat$type==typelst[2])

pdf("figures/FIGURE_realworld_scales.pdf", width=8, height=4, colormodel = "cmyk")
par(mfrow=c(1,2), mar=c(2,2,1,2), oma=c(2,3,0.5,0))
fcx<-1.4
ofs1<-c(0.12, -0.005)


plot(log10(plotdat[sbs1,]$temporal_years), log10(plotdat[sbs1,]$summed_sample_size), xlab="", ylab="", type="n", axes=F, xlim=log10(c(1, 32)), xaxs="i")
mtext(expression(paste("temporal span, years")), side = 1, line = 1, outer=TRUE, cex=1.2)
mtext(expression(paste("spatial span, m"^2)), side = 2, line = 1, outer=TRUE, cex=1.2)
xaxl<-c(1,2,4,8,16,32)
yaxl<-c(1,10,25,50,100,250,500,1000,3000)
axis(1, at=log10(xaxl), xaxl)
axis(2, at=log10(yaxl), yaxl, las=2)
box()

expord<-rev(order(tapply(plotdat[sbs1,]$summed_sample_size, plotdat[sbs1,]$experiment, mean)))
for(i in expord) {
  spsize<-log10(plotdat[sbs1,]$summed_sample_size[plotdat[sbs1,]$experiment==explst[i]])
  tmsize<-log10(plotdat[sbs1,]$temporal_years[plotdat[sbs1,]$experiment==explst[i]])
  
  polygon(c(tmsize, -1, -1, max(tmsize), max(tmsize)), c(spsize, max(spsize), -1, -1, min(spsize)), col=collst[i], density = 10, angle = anglst[i], lwd=1.5)
}
put.fig.letter("a.", "topleft", offset=ofs1, cex=fcx)




plot(log10(plotdat[sbs2,]$temporal_years), log10(plotdat[sbs2,]$summed_plot_size), xlab="", ylab="", type="n", axes=F, xlim=log10(c(1, 32)), xaxs="i")
xaxl<-c(1,2,4,8,16,32)
yaxl<-c(1,10,25,50,100,250,500,1000,3000,5000,10000,20000,50000)
axis(1, at=log10(xaxl), xaxl)
axis(2, at=log10(yaxl), yaxl, las=2)
box()

expord<-rev(order(tapply(plotdat[sbs2,]$summed_sample_size, plotdat[sbs2,]$experiment, mean)))
for(i in expord) {
  spsize<-log10(plotdat[sbs2,]$summed_plot_size[plotdat[sbs2,]$experiment==explst[i]])
  tmsize<-log10(plotdat[sbs2,]$temporal_years[plotdat[sbs2,]$experiment==explst[i]])
  
  polygon(c(tmsize, -1, -1, max(tmsize), max(tmsize)), c(spsize, max(spsize), -1, -1, min(spsize)), col=collst[i], density = 10, angle = anglst[i], lwd=1.5)
}
put.fig.letter("b.", "topleft", offset=ofs1, cex=fcx)


dev.off()

collst
explst
  