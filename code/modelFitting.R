## Modeling population dynamics for Age-0 fish surveys

rm(list=ls())

library(tidyverse)
library(lubridate)
library(mgcv)
library(MARSS)
library(fields)
library(viridis)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("fn_modsel_regionwide.R") #helper function

##Initial processing of fish data -----------------------------------------------------------------
dat_raw <- read.csv("../data/fish_cleaned_fallSurveys.csv")

dat_raw$Date <- as.Date(dat_raw$Date)
dat_raw$Month <- factor(month(dat_raw$Date))
dat_raw$Year <- factor(year(dat_raw$Date))

taxa <- unique(dat_raw$Taxon)
years <- 1980:2020
tt <- length(years)

dat_raw <- dat_raw[!is.na(dat_raw$region),]
dat_raw <- dat_raw[!is.na(dat_raw$Method),]

##Organize environmental covariate data -----------------------------------------------------------

env_comb <- read.csv("../data/env_drivers_combined.csv")

covars <- list(matrix(env_comb$qTot, nrow=1, byrow=TRUE)
               ,matrix(env_comb$X2, nrow=1, byrow=TRUE)
               ,matrix(env_comb$qOut, nrow=1, byrow=TRUE)
               ,matrix(env_comb$deltaT, nrow=1, byrow=TRUE)
               ,matrix(env_comb$secchi, nrow=1, byrow=TRUE)
               ,matrix(env_comb$chl, nrow=1, byrow=TRUE)
               ,matrix(env_comb$mesozoop, nrow=1, byrow=TRUE)
               ,matrix(env_comb$microzoop, nrow=1, byrow=TRUE)
               ,matrix(env_comb$upwelling, nrow=1, byrow=TRUE)
               ,matrix(env_comb$sst, nrow=1, byrow=TRUE)
)


covarmat <- cbind(qTot, X2, qOut, deltaT, secchi, chl, mesozoop, microzoop, upwelling, sst)
cor(covarmat)
## combos with high colinearity (abs(cor) >= 0.6) to exclude are:
#qTot x X2, secchi x microzoop, microzoop x upwelling


myscale <- function(x){
  x2 <- x - mean(x, na.rm=TRUE)
  return(x2/sd(x2, na.rm=TRUE))
}

## Loop over species and run model

dmin = 0.5 #must be present (non-zero, non-na) in 75% of timesteps

sppmods <- list()
sppmats <- list()

for(ss in 1:36){

  print(ss)
  dropyr <- NULL

  dat_spp <- dat_raw[dat_raw$Taxon == taxa[ss],]
  if(nrow(dat_spp)==0){next}
  dat_spp <- aggregate(CPUE ~  Source + Method + region + Year, data=dat_spp, FUN=mean)
  dat_spp$Year <- as.numeric(as.character(dat_spp$Year))
  dat_spp <- dat_spp[!is.na(dat_spp$CPUE),]
  dat_spp <- dat_spp[!is.infinite(dat_spp$CPUE),]
  dat_spp$SourceMethodRegion <- paste(dat_spp$Source, dat_spp$Method, dat_spp$region, sep="-")

  munit <- unique(dat_spp$SourceMethodRegion)
  nn <- length(munit)

  sppmat <- matrix(NA, nrow=nn, ncol=tt)

  for(ii in 1:nn){
    for(jj in 1:tt){
      if(any(dat_spp$SourceMethodRegion==munit[ii]
             & dat_spp$Year==years[jj]))
      {
        sppmat[ii,jj] <- dat_spp$CPUE[dat_spp$SourceMethodRegion==munit[ii]
                                      & dat_spp$Year==years[jj]]
      }
    }
  }

  a <- min(sppmat[sppmat > 0], na.rm=T)
  sppmat <- log(sppmat + a)

  nasum <- apply(sppmat, 1, function(x){sum(is.na(x))})
  zsum <- apply(sppmat, 1, function(x){sum(x==log(a), na.rm=T)})
  pmiss <- (nasum + zsum)/tt

  if(!any(pmiss <= (1-dmin))){next}

  sppmat <- sppmat[pmiss <= (1-dmin),]

  if(is.null(dim(sppmat))){
    sppmat <- matrix(sppmat, nrow=1, byrow=TRUE)
  }

  sppmat <- t(apply(sppmat, 1, myscale))
  rownames(sppmat) <- munit[pmiss <= (1-dmin)]
  colnames(sppmat) <- years

  ydata <- sppmat
  allNA <- apply(ydata, MARGIN=2, FUN=function(x){all(is.na(x))})
  if(allNA[1]){
    NArun <- rle(allNA)
    dropyr <- 1:NArun$lengths[1]
    if(nrow(ydata)==1){
      ydata <- ydata[-dropyr]
    }
    else{
      ydata <- ydata[,-dropyr]
    }
  }

  ## make model
  regions <- NULL
  for(jj in 1:nrow(sppmat)){
    regions <- c(regions, strsplit(rownames(sppmat)[jj], "-")[[1]][3])
  }
  Z=factor(rep(1, nrow(sppmat)))
  sppmats[[paste(taxa[ss])]] <- sppmat
  sppmods[[paste(taxa[ss])]] <- modsel_regionwide(ydata, covars, dropyr)
  sppmods[[paste(taxa[ss])]]$dropyr <- dropyr

}

saveRDS(sppmats, file="../outputs/sppmats_fallAge0_regionwide.RDS")
saveRDS(sppmods, file="../outputs/sppmods_fallAge0_regionwide.RDS")



### Make no covariate model

mods_noCovar <- list()

for(ss in 1:length(taxa)){

  print(ss)
  dropyr <- NULL

  dat_spp <- dat_raw[dat_raw$Taxon == taxa[ss],]
  if(nrow(dat_spp)==0){next}
  dat_spp <- aggregate(CPUE ~  Source + Method + region + Year, data=dat_spp, FUN=mean)
  dat_spp$Year <- as.numeric(as.character(dat_spp$Year))
  dat_spp <- dat_spp[!is.na(dat_spp$CPUE),]
  dat_spp <- dat_spp[!is.infinite(dat_spp$CPUE),]
  dat_spp$SourceMethodRegion <- paste(dat_spp$Source, dat_spp$Method, dat_spp$region, sep="-")

  munit <- unique(dat_spp$SourceMethodRegion)
  nn <- length(munit)

  sppmat <- matrix(NA, nrow=nn, ncol=tt)

  for(ii in 1:nn){
    for(jj in 1:tt){
      if(any(dat_spp$SourceMethodRegion==munit[ii]
             & dat_spp$Year==years[jj]))
      {
        sppmat[ii,jj] <- dat_spp$CPUE[dat_spp$SourceMethodRegion==munit[ii]
                                      & dat_spp$Year==years[jj]]
      }
    }
  }

  a <- min(sppmat[sppmat > 0], na.rm=T)
  sppmat <- log(sppmat + a)

  nasum <- apply(sppmat, 1, function(x){sum(is.na(x))})
  zsum <- apply(sppmat, 1, function(x){sum(x==log(a), na.rm=T)})
  pmiss <- (nasum + zsum)/tt

  if(!any(pmiss <= (1-dmin))){next}

  sppmat <- sppmat[pmiss <= (1-dmin),]

  if(is.null(dim(sppmat))){
    sppmat <- matrix(sppmat, nrow=1, byrow=TRUE)
  }

  sppmat <- t(apply(sppmat, 1, myscale))
  rownames(sppmat) <- munit[pmiss <= (1-dmin)]
  colnames(sppmat) <- years
  
  ydata <- sppmat
  allNA <- apply(ydata, MARGIN=2, FUN=function(x){all(is.na(x))})
  if(allNA[1]){
    NArun <- rle(allNA)
    dropyr <- 1:NArun$lengths[1]
    if(nrow(ydata)==1){
      ydata <- ydata[-dropyr]
    }
    else{
      ydata <- ydata[,-dropyr]
    }
    sppyears <- years[-dropyr]
  }

  ## make model

  Z=factor(rep(1, nrow(sppmat)))
  sppmats[[paste(taxa[ss])]] <- sppmat
  mods_noCovar[[paste(taxa[ss])]] <- MARSS(ydata, model=list(U="unequal", B="diagonal and unequal", Z=Z, R="diagonal and unequal", tinitx=1),
                                           control=list(maxit=50000, safe=TRUE), silent=TRUE)
  mods_noCovar[[paste(taxa[ss])]]$dropyr <- dropyr

}


saveRDS(mods_noCovar, "../outputs/nullmods_fallAge0_regionwide.RDS")
