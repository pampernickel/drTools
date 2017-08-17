.matchDrugs <- function(res.df){
  # designed to match drugs created from a data frame generated with the drTools averageRes()
  # call, can be expanded to deal with strings
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
    
    colnames(res.df)[which(colnames(res.df) %ni% common.drugs)] -> unmatched
  }
  
  list(common.drugs, unmatched) -> res
  names(res) <- c("common", "unmatched")
  return(res)
}

.mapResponse <- function(res.df, assembled, drug.list.all=NULL, pat.name=NULL){
  # conversion
  nrow(res.df) -> e
  as.matrix(res.df) -> conv
  if (is.null(dim(res.df)) & ncol(conv) < nrow(conv)){
    # case of named numeric vector
    t(conv) -> res.df
  }
  
  .matchDrugs(res.df) -> resm
  resm$common -> common.drugs
  res.df[,which(colnames(res.df) %in% common.drugs)] -> res.df.sub
  assembled[,which(colnames(assembled) %in% common.drugs)] -> assembled.sub
  
  c.ind <- NA
  if (is.matrix(res.df.sub)){
    sapply(colnames(assembled.sub), function(x) 
      which(toupper(colnames(res.df.sub)) %in% x)) -> c.ind
  } else {
    sapply(colnames(assembled.sub), function(x) 
      which(toupper(names(res.df.sub)) %in% x)) -> c.ind
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