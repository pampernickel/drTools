library(RCurl)
source_https('https://raw.githubusercontent.com/pampernickel/chemblr/pampernickel-patch-1/package/R/query.R')

.matchDrugs <- function(res.df, assembled, drug.list.all){
  # designed to match drugs created from a data frame generated with the drTools averageRes()
  # call, can be expanded to deal with strings; in new version, also requires
  # drug.list.all to correct for drugs with a misspelled name
  res <- common.drugs <- unmatched.drugs <- NA
  if (is.matrix(res.df) || is.data.frame(res.df)){
    res.df[which(res.df >= 4.5)] <- 4.5
    res.df[which(res.df <= min(assembled, na.rm=T))] <- min(assembled, na.rm=T)
    toupper(colnames(res.df)) -> colnames(res.df)
    toupper(colnames(assembled)) -> colnames(assembled)
    
    # cases where nothing was fitted
    apply(res.df, 2, function(x) length(which(is.na(x)))) -> nas
    as.numeric(which(nas == nrow(res.df))) -> ind
    
    if (length(ind) > 0){
      # remove drugs where no fit was performed, then check
      # again if the matrix has been converted to a named numeric
      res.df[,-ind] -> res.df
      nrow(res.df) -> e
      as.matrix(res.df) -> conv
      if (is.null(dim(res.df)) & ncol(conv) < nrow(conv)){
        # case of named numeric vector
        t(conv) -> res.df
      }
    }
    gsub("\\.", "-", colnames(res.df)) -> colnames(res.df)
    sapply(colnames(res.df), function(x) 
      ifelse(length(grep(x, drug.list.all$all.names))==1, 
             drug.list.all$final.name[grep(x, drug.list.all$all.names)],
             x)) -> colnames(res.df)
    intersect(colnames(res.df), colnames(assembled)) -> common.drugs
    
    if (length(grep("_",colnames(res.df))) > 0 && length(common.drugs) == 0){
      # split & sub
      gsub("\\.", "-", as.character(sapply(strsplit(colnames(res.df), "_"), 
                                           function(x) x[1]))) -> colnames(res.df)
      intersect(colnames(res.df), colnames(assembled)) -> common.drugs
    }
    
    as.character(colnames(res.df)[which(colnames(res.df) %ni% common.drugs)]) -> unmatched
    findClosestMatch(unmatched, drug.list.all) -> corrected
    
    # then check if the corrected form of the drug name is among the column names of assembled;
    # if they are, then remove from corrected$new.drug
    if (length(which(corrected$match %ni% colnames(assembled))) > 0){
      unmatched <- corrected$match[which(corrected$match %ni% colnames(assembled))]
    }
    
    corrected$match[which(corrected$match %in% colnames(assembled))] -> corrected.fin
  
    if (length(corrected.fin) > 0){
      for (i in 1:length(corrected.fin)){
        as.character(corrected.fin[i]) -> colnames(res.df)[which(colnames(res.df) %in% corrected.drug[i])]
      }
    }
    common.drugs <- c(common.drugs, corrected.fin[which(corrected.fin %ni% NA)])
  }
  
  list(common.drugs, unmatched, res.df) -> res
  names(res) <- c("common", "unmatched", "corrected")
  return(res)
}

.mapResponse <- function(res.df, assembled, drug.list.all, pat.name=NULL){
  # function under dataFormattingFuncs.r
  # Function that maps patient data onto currently available population data
  # This function is NOT meant to add new drug data, which will be performed
  # in the collectFit() function (under dataCollectionsFuncs.r)
  nrow(res.df) -> e
  as.matrix(res.df) -> conv
  if (is.null(dim(res.df)) & ncol(conv) < nrow(conv)){
    # case of named numeric vector
    t(conv) -> res.df
  }
  
  .matchDrugs(res.df, assembled, drug.list.all) -> resm
  resm$common -> common.drugs
  resm$corrected -> res.df
  res.df[,which(colnames(res.df) %in% common.drugs)] -> res.df.sub
  assembled[,which(colnames(assembled) %in% common.drugs)] -> assembled.sub
  
  c.ind <- NA
  if (is.matrix(res.df.sub) || is.data.frame(res.df.sub)){
    sapply(colnames(assembled.sub), function(x) 
      which(toupper(colnames(res.df.sub)) %in% x)) -> c.ind
  } else {
    sapply(colnames(assembled.sub), function(x) 
      which(toupper(names(res.df.sub)) %in% x)[1]) -> c.ind
  }
  
  if (is.matrix(res.df.sub)){
    nrow(res.df.sub) -> ml
  } else {
    max(sapply(c.ind, function(x) length(x))) -> ml
  }
  
  respMat <- matrix(NA, nrow=ml, ncol=length(c.ind))
  for (i in 1:length(c.ind)){
    if (is.matrix(res.df.sub)){
      for (j in 1:ml){
        respMat[j,i] <- res.df.sub[j,c.ind[[i]]]
      }
    } else {
      respMat[1:length(c.ind[[i]]),i] <- res.df.sub[c.ind[[i]]]
    }
  }
  colnames(respMat) <- colnames(assembled.sub)
  
  if (is.matrix(res.df.sub)){
    rownames(respMat) <- rownames(res.df.sub)
  }
  
  rbind(assembled.sub, respMat) -> assembled.sub
  which(rownames(assembled.sub) %ni% 
          rownames(assembled)) -> r.ind
  
  if (!is.matrix(res.df.sub) && is.null(pat.name)){
    rownames(assembled.sub)[which(rownames(assembled.sub) %ni% rownames(assembled))] <- 
      paste("POI",r.ind, sep=",")
  } else if (length(which(rownames(assembled.sub) %ni% rownames(assembled))) == length(pat.name)) {
    rownames(assembled.sub)[which(rownames(assembled.sub) %ni% rownames(assembled))] <- 
      pat.name
  } else {
    stop("Number of patient names specified does not match available data.")
  }
  
  cbind(rownames(assembled.sub), assembled.sub) -> assembled.sub
  colnames(assembled.sub)[1] <- "id"
  return(assembled.sub)
}

condenseToDF <- function(l){
  # create a data frame from a set of lists with unequal length; this
  # function is specifically designed for lists of patient parameters, e.g.
  # IC50, AUC, Emax lists; l may contain either data frames
  # or lists; in the case that it is a list of data frames
  # l has to be preprocessed such that the drug names have already been converted
  # to the standard drug names, i.e. can be compared across different entries
  # in l
  if (is.list(l[[1]]) || is.numeric(l[[1]])){
    unique(unlist(sapply(l, function(x) 
      toupper(names(x))))) -> dl
    res <- matrix(NA, nrow=length(l), ncol=length(dl))
    colnames(res) <- dl
    rownames(res) <- names(l)
    for (i in 1:length(l)){
      res[i, match(names(l[[i]]), colnames(res))] <- l[[i]]
    }
  } else if (is.matrix(l[[1]]) || is.matrix[l[[1]]]){
    unique(unlist(sapply(l, function(x) 
      toupper(colnames(x))))) -> dl
    res <- matrix(NA, nrow=sum(sapply(l, function(x) nrow(x))), ncol=length(dl))
    colnames(res) <- dl
    rownames(res) <- as.character(unlist(sapply(l, function(x) rownames(x))))
    ctr <- 1
    for (i in 1:length(l)){
      for (j in 1:nrow(l[[i]])){
        res[ctr, match(colnames(l[[i]]), colnames(res))] <- l[[i]][j,]
        ctr+1 -> ctr
      }
    }
  }
  return(res)
}
