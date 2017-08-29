# PDobay for AG Bourquin, KISPI
# pamela.dobay@kispi.uzh.ch

#@...
# Routines that extract information for a drug from a list
# of drug response profiles for the purpose of comparison
#@...


getMax <- function(exp.res.ave, dl){
  # given an exp.res.ave object from the averageRes() function
  # and a list of drugs, extract the fits
  lapply(exp.res.ave, function(x){
    lapply(x$max, function(y){
      extractMax(y)
    }) -> fits
    if (length(which(is.na(dl)))>0){
      fits[-which(is.na(dl))] -> fits
      dl[-which(is.na(dl))] -> dl
    }
    names(fits) <- dl
    return(fits)
  }) -> fits
  return(fits)
}


extractFit <- function(summary, drug){
  lapply(summary, function(x){
    if(drug %in% names(x)){
      return(x[[which(names(x) %in% drug)[1]]]$x)
    }
  }) -> x.fits
  
  lapply(summary, function(x){
    if(drug %in% names(x)){
      return(x[[which(names(x) %in% drug)[1]]]$y)
    }
  }) -> y.fits
  
  # create ggplot-compatible data frame
  res.df <- matrix(0, nrow=0, ncol=4)
  colnames(res.df) <- c("x", "y", "patient", "drug")
  for (i in 1:length(x.fits)){
    if (!is.null(x.fits[[i]])){
      cbind(x.fits[[i]], y.fits[[i]], 
            rep(names(x.fits)[i], 100), 
            rep(drug, 100)) -> t
      colnames(t) <- colnames(res.df)
      rbind(res.df, t) -> res.df
    }
  }
  
  if (nrow(res.df) > 0){
    as.data.frame(res.df) -> res.df
    as.numeric(as.character(res.df$x)) -> res.df$x
    as.numeric(as.character(res.df$y)) -> res.df$y
  }
  return(res.df)
}

extractParam <- function(exp.res.ave, param){
  # Extract a given a param (IC50, EC50, Emax, AUC) from an exp.res.ave object
  t(sapply(exp.res.ave, function(x) x[[which(names(x) %in% param)]])) -> res
  return(res)
}