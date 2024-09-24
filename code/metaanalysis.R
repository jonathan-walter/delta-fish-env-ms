library(trend)
library(MARSS)
library(viridis)

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

sppmods <- readRDS("../outputs/sppmods_fallAge0_regionwide.RDS")


#drop selected species and species with poor fits
sppmods <- sppmods[!names(sppmods) %in% c("Ameiurus catus"
                                          ,"Cymatogaster aggregata"
                                          ,"Hysterocarpus traskii"
                                          ,"Ictalurus punctatus"
                                          ,"Menidia audens"
                                          ,"Microgadus proximus"
                                          ,"Parophrys vetulus"
                                          ,"Platichthys stellatus"
                                          ,"Pimephales promelas"
                                          ,"Pogonichthys macrolepidotus"
                                          ,"Porichthys notatus"
                                          ,"Sebastes auriculatus"
                                          ,"Syngnathus leptorhynchus"
)]

nullMods <- readRDS("../outputs/nullmod_fallAge0_regionwide.RDS")
nullMods <- nullMods[names(nullMods) %in% names(sppmods)]

#nullBICs <- read.csv("../data/nullBICs.csv")
#nullBICs <- nullBICs[10:34,]

sppDetails <- read.csv("../data/SpeciesDetails.csv")
sppDetails <- sppDetails[1:36,]
sppDetails$Salinity <- factor(sppDetails$Salinity, 
                                 levels=c("Marine","Marine/Brackish","Marine/Brackish/Freshwater",
                                          "Marine/Brackish/Freshwater (anadromous)","Brackish","Brackish/Freshwater",
                                          "Brackish/Freshwater (anadromous)","Freshwater"))
sppDetails$Water.Column <- factor(sppDetails$Water.Column,
                                  levels=c("Demersal","Benthopelagic","Pelagic"))
sppDetails <- sppDetails[order(sppDetails$Salinity, sppDetails$Water.Column, sppDetails$Family),]
sppDetails <- sppDetails[sppDetails$Scientific.name %in% names(sppmods),]

years <- 1980:2020

sppmats <- readRDS("../data/sppmats_fallAge0_regionwide.RDS")
sppmats <- sppmats[names(sppmats) %in% names(nullMods)]

sppmods <- sppmods[names(sppmods) %in% names(nullMods)]


## Extract AICc from null models ------------------------------------------------------------------

nullAICc <- rep(NA, length(nullMods))

for(ii in 1:length(nullMods)){
  
  nullAICc[ii] <- nullMods[[ii]]$AICc
  
}

nullAICc <- data.frame(taxon=names(nullMods), AICc=nullAICc)



## Look at trends in model states -----------------------------------------------------------------

states <- list()
topmods <- list()
states.se <- list()

for(ii in 1:length(sppmods)){
  
  yy <- years
  
  if(nullAICc$AICc[ii] < min(sppmods[[ii]]$AICc)){
    states.ii <- nullMods[[ii]]$states
    states[[paste0(names(sppmods)[ii])]] <- states.ii
    states.se[[paste0(names(sppmods)[ii])]] <- nullMods[[ii]]$states.se
    lc.ii <- exp(nullMods[[ii]]$states - nullMods[[ii]]$states.se)
    uc.ii <- exp(nullMods[[ii]]$states + nullMods[[ii]]$states.se)
    
    topmods[[paste0(names(sppmods)[ii])]] <- nullMods[[ii]]
    
  }
  else{
    states.ii <- sppmods[[ii]]$topmod$states
    topmods[[paste0(names(sppmods)[ii])]] <- sppmods[[ii]]$topmod
    states[[paste0(names(sppmods)[ii])]] <- states.ii
    states.se[[paste0(names(sppmods)[ii])]] <- sppmods[[ii]]$topmod$states.se
    lc.ii <- exp(sppmods[[ii]]$topmod$states - sppmods[[ii]]$topmod$states.se)
    uc.ii <- exp(sppmods[[ii]]$topmod$states + sppmods[[ii]]$topmod$states.se)
  }
  
  if(!is.null(sppmods[[ii]]$dropyr)){
    yy <- years[-sppmods[[ii]]$dropyr]
  }

  plot(yy, exp(states.ii), main=names(sppmods)[ii], type="l")
  lines(yy, lc.ii, lty=2)
  lines(yy, uc.ii, lty=2)

}

## test for trends in states 
state_trends <- list()  
for(ii in 1:length(states)){
  states.ii <- exp(states[[ii]])
  state_trends[[paste0(names(sppmods)[ii])]] <- sens.slope(states.ii)
}


### Plot states for manuscript --------------------------------------------------------------------

regpal <- data.frame(region=c("Central Bay","Confluence","Napa River","North","San Pablo","South","South Bay","Suisun Bay and Marsh"),
                     color=c("#8a3ffc","#ff7eb6", "#33b1ff","#007d79","#fa4d56","#6fdc8c","#d2a106","#08bdba"))

alpha <- c("a)","b)","c)","d)","e)","f)","g)","h)","i)","j)","k)","l)","m)","n)","o)","p)","q)","r)","s)","t)","u)","v)","w)","x)")


pdf("../outputs/fig_sppStates_rw.pdf", width=6.5, height=8)

par(mfrow=c(7,3), mar=c(0.5,2.1,1.5,0.5), oma=c(3.1,1.5,0,0.5))

for(ii in 1:length(states)){
  
  tmp1 <- states[[paste0(sppDetails$Scientific.name[ii])]]
  tmp2 <- exp(sppmats[[paste0(sppDetails$Scientific.name[ii])]])
  
  yy <- years
  if(!is.null(sppmods[[paste0(sppDetails$Scientific.name[ii])]]$dropyr)){
    yy <- years[-sppmods[[paste0(sppDetails$Scientific.name[ii])]]$dropyr]
    if(nrow(tmp2)==1){
      tmp2 <- tmp2[-sppmods[[paste0(sppDetails$Scientific.name[ii])]]$dropyr]
    }
    else{
      tmp2 <- tmp2[,-sppmods[[paste0(sppDetails$Scientific.name[ii])]]$dropyr]
    }
  }
  
  plot(rep(yy, each=nrow(tmp2)), c(tmp2), pch=".", ylab="", xlab="",
       col="darkgrey", xaxt="n", xlim=range(years), log="y")
  # plot(yy, exp(tmp1), type="l", lwd=1.5, ylab="CPUE", xlab="Year",
  #      col="black", xaxt="n", xlim=range(years),
  #      ylim=range(c(exp(tmp1+tmp2)),exp(tmp1-tmp2)))
  axis(1, at=c(1980,1990,2000,2010,2020), labels=FALSE)
  mtext(sppDetails$Common.name[ii], cex=3/4, line=0.1)
  mtext(alpha[ii], at=1980, cex=3/4, line=0.1)
  
  tmp3 <- state_trends[[paste0(sppDetails$Scientific.name[ii])]]
  if(tmp3$p.value >= 0.05){
    lines(yy, exp(tmp1), lwd=1.5, col="black")
  }
  if(tmp3$p.value < 0.05 & tmp3$estimates < 0){
    lines(yy, exp(tmp1), lwd=1.5, col="red")
  }
  if(tmp3$p.value < 0.05 & tmp3$estimates > 0){
    lines(yy, exp(tmp1), lwd=1.5, col="blue")
  }
  if(ii %in% c(19,20,21)){
    axis(1, at=c(1980,1990,2000,2010,2020))
  }
  
}
mtext("Year",1,outer=T, cex=3/4, line=1.9)
mtext("Age-0 abundance (scaled CPE)",2,outer=T, cex=3/4, line=0.2)

dev.off()


## Make variable importances for the all predictors 

AICcMat <- NULL

for(ii in 1:length(sppmods)){
  AICcMat <- rbind(AICcMat, sppmods[[ii]]$AICc)
}

AICcMat <- cbind(nullAICc$AICc, AICcMat)
colnames(AICcMat) <- c("None","qTot","X2","qOut","deltaT","secchi","chl","meso","micro","upwelling","sst")

rownames(AICcMat) <- names(sppmods)

akaike.weights <- function(x){
  
  dAIC <- x-min(x)
  wAIC <- exp(-0.5*dAIC)/sum(exp(-0.5*dAIC))
  return(wAIC)
}

wAICcMat <- t(apply(AICcMat, MARGIN=1, FUN=akaike.weights))


wAICcMat <- wAICcMat[match(sppDetails$Scientific.name, rownames(wAICcMat)),]


effectDir <- c("-","-","-","*","*","-","*","*","+","+","*","+","-","-","-","-","+","-","+","+","+")

pdf("../outputs/fig_wAICc.pdf", width=4.25, height=6.5)
par(mar=c(6.5,9.8,4,0.5), mgp=c(2,0.7,0))
image(x=1:11, y=1:nrow(wAICcMat), t(wAICcMat), col=viridis(50), xaxt="n", yaxt="n", xlab="", ylab="", zlim=c(0,1))
axis(1, at=1:11, labels=c("None", "Delta inflow","X2","Delta outflow", "Delta temp", "Secchi","Chlorophyll-a","Mesozoop.","Microzoop.",
                         "Upwelling","SST"), las=2)
axis(2, at=1:nrow(wAICcMat), labels=sppDetails$Common.name, las=2)
for(ii in 1:nrow(wAICcMat)){
  text(y=ii, x=which.max(wAICcMat[ii,]),effectDir[ii], col="black")
}
#scalebar
par(fig=c(0,1,0.91,1), mar=c(1,10,1.1,1), new=TRUE)
image(t(matrix(1:50, byrow=TRUE, nrow=1)), col=viridis(50), yaxt="n", xaxt="n")
axis(1, at=c(0,0.5,1))
mtext("AICc weight", line=0.1, cex=0.9)

dev.off()



## How strongly are model predicted states correlated with raw data from different sampling programs?

#make matrix to hold outputs

progcombs <- c("Bay Study-Midwater trawl-South Bay",
               "Bay Study-Otter trawl-South Bay",
               "Bay Study-Midwater trawl-Central Bay",
               "Bay Study-Otter trawl-Central Bay",
               "Bay Study-Midwater trawl-San Pablo",
               "Bay Study-Otter trawl-San Pablo",
               "Bay Study-Midwater trawl-Suisun Bay and Marsh",
               "Bay Study-Otter trawl-Suisun Bay and Marsh",
               "Bay Study-Midwater trawl-Confluence",
               "Bay Study-Otter trawl-Confluence",
               "DJFMP-Beach seine-Central Bay",
               "DJFMP-Midwater trawl-Central Bay",
               "DJFMP-Beach seine-San Pablo",
               "DJFMP-Midwater trawl-San Pablo",
               "DJFMP-Beach seine-Suisun Bay and Marsh",
               "DJFMP-Midwater trawl-Suisun Bay and Marsh",
               "DJFMP-Beach seine-Confluence",
               "DJFMP-Midwater trawl-Confluence",
               "DJFMP-Beach seine-North",
               "DJFMP-Midwater trawl-North",
               "DJFMP-Beach seine-South",
               "DJFMP-Midwater trawl-South",
               "FMWT-Midwater trawl-Central Bay",
               "FMWT-Midwater trawl-San Pablo",
               "FMWT-Midwater trawl-Napa River",
               "FMWT-Midwater trawl-Suisun Bay and Marsh",
               "FMWT-Midwater trawl-Confluence",
               "FMWT-Midwater trawl-North",
               "FMWT-Midwater trawl-South",
               "Suisun-Beach seine-Suisun Bay and Marsh",
               "Suisun-Otter trawl-Suisun Bay and Marsh")


spp <- sppDetails$Scientific.name

programCorr <- matrix(NA, nrow=length(spp), ncol=length(progcombs))
rownames(programCorr) <- spp
colnames(programCorr) <- progcombs

for(ii in 1:length(spp)){
  
  raw.spp <- sppmats[[paste0(spp[ii])]]
  state.spp <- states[[paste0(spp[ii])]]
  dropyr <- sppmods[[paste0(spp[ii])]]$dropyr
  
  if(!is.null(dropyr)){
    if(length(raw.spp)<43){
      raw.spp <- raw.spp[-dropyr]
    }
    else{
      raw.spp <- raw.spp[,-dropyr]
    }
    
  }
  
  if(length(raw.spp)<43){
    programCorr[ii, which(progcombs==rownames(raw.spp))] <- cor(x=c(raw.spp), y=c(state.spp), use="pairwise.complete.obs")
  }
  else{
    for(jj in 1:nrow(raw.spp)){
      programCorr[ii, which(progcombs==rownames(raw.spp)[jj])] <- cor(x=c(raw.spp[jj,]), y=c(state.spp), use="pairwise.complete.obs")
    }
  }
  
}

programCorr2 <- programCorr[,colSums(is.na(programCorr))<nrow(programCorr)]


pal <- colorRampPalette(c("red","grey85","blue"))


progcombs.nice <- c("Bay Study-MWT-South Bay",
               "Bay Study-OT-South Bay",
               "Bay Study-MWT-Central Bay",
               "Bay Study-OT-Central Bay",
               "Bay Study-MWT-San Pablo",
               "Bay Study-OT-San Pablo",
               "Bay Study-MWT-Suisun",
               "Bay Study-OT-Suisun",
               "Bay Study-MWT-Confluence",
               "Bay Study-OT-Confluence",
               "DJFMP-BS-Central Bay",
               "DJFMP-MWT-Central Bay",
               "DJFMP-BS-San Pablo",
               "DJFMP-MWT-San Pablo",
               "DJFMP-BS-Suisun",
               "DJFMP-MWT-Suisun",
               "DJFMP-BS-Confluence",
               "DJFMP-MWT-Confluence",
               "DJFMP-BS-North Delta",
               "DJFMP-MWT-North Delta",
               "DJFMP-BS-South Delta",
               "DJFMP-MWT-South Delta",
               "FMWT-MWT-Central Bay",
               "FMWT-MWT-San Pablo",
               "FMWT-MWT-Napa River",
               "FMWT-MWT-Suisun",
               "FMWT-MWT-Confluence",
               "FMWT-MWT-North Delta",
               "FMWT-MWT-South Delta",
               "Suisun-BS-Suisun",
               "Suisun-OT-Suisun")

progcombs.nice <- progcombs.nice[colSums(is.na(programCorr))<nrow(programCorr)]

pdf("../outputs/fig_program_correlations.pdf", height=7.75, width=6.5)
par(mar=c(12.5,10,4,1))
image(y=1:length(spp), x=1:ncol(programCorr2), t(programCorr2), xaxt="n", yaxt="n", xlab="", ylab="", 
      zlim=c(-1,1), col=pal(50))
axis(1, at=1:ncol(programCorr2), labels=progcombs.nice, las=2)
axis(2, at=1:length(spp), labels=sppDetails$Common.name, las=2)
abline(v=c(10.5,17.5,23.5))
#scalebar
par(fig=c(0,1,0.925,1), mar=c(1,10,1.1,1), new=TRUE)
image(t(matrix(1:50, byrow=TRUE, nrow=1)), col=pal(50), yaxt="n", xaxt="n")
axis(1, at=c(0,0.5,1), labels=c(-1,0,1))
mtext("Pearson correlation", line=0.1, cex=0.9)

dev.off()



### Look at estimated sampling variances of fitted models

sampvarMat <- matrix(NA, nrow=length(spp), ncol=length(progcombs))
rownames(sampvarMat) <- spp
colnames(sampvarMat) <- progcombs

for(ii in 1:length(spp)){
  
  tmp <- topmods[[paste0(spp[ii])]]$par$R
 
  for(jj in 1:length(tmp)){
    tmp2 <- rownames(tmp)[jj]
    tmp2 <- strsplit(tmp2, split=",")[[1]][1]
    tmp2 <- substr(tmp2, 2, 500)
    sampvarMat[ii, colnames(sampvarMat)==tmp2] <- tmp[jj,1]
  }
  
}

sampvarMat2 <- sampvarMat[,colSums(is.na(sampvarMat))<nrow(sampvarMat)]


pdf("../outputs/fig_sampling_variance.pdf", height=7.75, width=6.5)
par(mar=c(12.5,10,4,1))
image(y=1:length(spp), x=1:ncol(sampvarMat2), t(sampvarMat2), xaxt="n", yaxt="n", xlab="", ylab="", col=viridis(50))
axis(1, at=1:ncol(programCorr2), labels=progcombs.nice, las=2)
axis(2, at=1:length(spp), labels=sppDetails$Common.name, las=2)
abline(v=c(10.5,17.5,23.5))
#scalebar
par(fig=c(0,1,0.925,1), mar=c(1,10,1.1,1), new=TRUE)
image(t(matrix(1:50, byrow=TRUE, nrow=1)), col=viridis(50), yaxt="n", xaxt="n")
axis(1, at=c(0,0.5,1), labels=c(0.17,0.58,1.33))
mtext("Observation variance", line=0.1, cex=0.9)
dev.off()

cor(c(programCorr2),c(sampvarMat2), use="pairwise.complete.obs") #These are strongly negatively correlated!

sampvarLong <- data.frame(species = rep(rownames(sampvarMat2), times=ncol(sampvarMat2)),
                          sampgrp = rep(colnames(sampvarMat2), each=nrow(sampvarMat2)),
                          sampvar = c(sampvarMat2)
                          )

sampvarLong$program <- NA
sampvarLong$method <- NA
sampvarLong$region <- NA

for(ii in 1:nrow(sampvarLong)){
  tmp <- strsplit(sampvarLong$sampgrp[ii], split="-")[[1]]
  sampvarLong$program[ii] <- tmp[1]
  sampvarLong$method[ii] <- tmp[2]
  sampvarLong$region[ii] <- tmp[3]
}

boxplot(sampvar ~ program, data=sampvarLong, xlab="")
boxplot(sampvar ~ method, data=sampvarLong, las=2, xlab="")
boxplot(sampvar ~ region, data=sampvarLong, las=2, xlab="")

pdf("../outputs/fig_sampvar_programMethod.pdf")
par(mar=c(12,4,1,1))
boxplot(sampvar ~ program*method, data=sampvarLong, las=2, drop=TRUE, xlab="", ylab="Observation variance")
dev.off()

pdf("../outputs/fig_sampvar_programMethodRegion.pdf", width=11)
par(mar=c(20,4,1,1))
boxplot(sampvar ~ program*region*method, data=sampvarLong, las=2, drop=TRUE, xlab="", ylab="Observation variance")
dev.off()


# combined figure
names1 <- c("DJFMP-BS","Bay Study-MWT","DJFMP-MWT","FMWT-MWT","Bay Study-OT","Suisun-OT")
names2 <- c("DJFMP-Central Bay-BS","DJFMP-Confluence-BS","DJFMP-North Delta-BS","DJFMP-San Pablo Bay-BS",
            "DJFMP-South Delta-BS","Bay Study-Central Bay-MWT","Bay Study-Confluence-MWT",
            "DJFMP-Confluence-MWT","FMWT-Confluence-MWT","FMWT-Napa River-MWT","DJFMP-North Delta-MWT",
            "FMWT-North Delta-MWT","Bay Study-San Pablo Bay-MWT","FMWT-San Pablo Bay-MWT",
            "FMWT-South Delta-MWT","Bay Study-South Bay-MWT","Bay Study-Suisun Bay & Marsh-MWT",
            "FMWT-Suisun Bay & Marsh-MWT","Bay Study-Central Bay-OT","Bay Study-Confluence-OT",
            "Bay Study-San Pablo Bay-OT","Bay Study-South Bay-OT","Bay Study-Suisun Bay & Marsh-OT",
            "Suisun-Suisun Bay & Marsh-OT")


pdf("../outputs/fig_sampvar_2panel.pdf", width=5, height=8)
layout(matrix(1:2,nrow=2,byrow=TRUE),heights=c(length(names1),length(names2)))
par(mar=c(0.1,15.5,0.5,1), mgp=c(2,0.7,0))
boxplot(sampvar ~ program*method, data=sampvarLong, las=2, drop=TRUE, 
        ylab="", xlab="Observation variance", horizontal = TRUE,
        names=names1, xaxt="n")
text(par("usr")[1]+0.05*diff(par("usr")[1:2]), par("usr")[4]-0.06*diff(par("usr")[3:4]),"a)", cex=0.9)
par(mar=c(3.5,15.5,0.1,1))
boxplot(sampvar ~ program*region*method, data=sampvarLong, las=2, 
        drop=TRUE, ylab="", xlab="Observation variance", horizontal=TRUE, names=names2)
text(par("usr")[1]+0.05*diff(par("usr")[1:2]), par("usr")[4]-0.015*diff(par("usr")[3:4]),"b)", cex=0.9)
dev.off()

