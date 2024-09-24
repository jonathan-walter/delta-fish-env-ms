## function for doing covariate selection using MARSS

## The model structure is going to be the same for all runs and is hard-coded.
## The inputs are 
##  1) a response matrix (designed to be fish species), rows are region-sampling program combos and cols are time (months)
##        - it's assumed this has already been cleaned, and that the row names give the name of the sampling program and region
##  2) a list of covariate matrices
##        - it's assumed that the row names of the matrices give the covariate names
##  3) a control object for MARSS models for passing fitting/opimitization parameters


modsel_regionwide <- function(resp, covar, dropyr, control=list(maxit=50000, safe=TRUE)){
  
  allMod <- list()
  AICc <- rep(NA, length(covar))
  
  if(is.null(dropyr)){
    for(ii in 1:length(covar)){
      
      #mod <- MARSS(resp, model=list(U="unequal", Z=Z, Q="unconstrained", C="unconstrained", R="diagonal and unequal", c=covar[[ii]]),
      #             control=control, silent=TRUE)
      mod <- try(MARSS(resp, model=list(U="unequal", B="diagonal and unequal", Z=Z, C="unconstrained", R="diagonal and unequal", c=covar[[ii]]),
                       control=control, silent=TRUE))
      if(grepl("try-error", class(mod))){
        allMod[[ii]] <- NULL
        AICc[ii] <- NA
      }
      else{
        allMod[[ii]] <- mod
        AICc[ii] <- mod$AICc
      }
    }
  }
  else{
    for(ii in 1:length(covar)){
      
      #mod <- MARSS(resp, model=list(U="unequal", Z=Z, Q="unconstrained", C="unconstrained", R="diagonal and unequal", c=covar[[ii]]),
      #             control=control, silent=TRUE)
      mod <- try(MARSS(resp, model=list(U="unequal", B="diagonal and unequal", Z=Z, C="unconstrained", R="diagonal and unequal", 
                                        c=matrix(covar[[ii]][-dropyr], nrow=1, byrow=TRUE)),
                       control=control, silent=TRUE))
      if(grepl("try-error", class(mod))){
        allMod[[ii]] <- NULL
        AICc[ii] <- NA
      }
      else{
        allMod[[ii]] <- mod
        AICc[ii] <- mod$AICc
      }
    }
  }

  
  return(list(AICc=AICc, allMod=allMod, topmod = allMod[[which.min(AICc)]]))
  
  
}
