# PDobay for AG Bourquin, KISPI
# pamela.dobay@kispi.uzh.ch

#@...
# Routine that collects and saves drug response information from all runs
# using patient material (i.e. with LK200x tag); would currently go for the
# solution of storing in an r data file (per patient) and in the 20161012.heatmap.all.fits.rda
# matrix, which contains the response data of the blood paper cohort
#@...

getNewestFile <- function(){
  list.files("./r.data.files", pattern=".rda", full.names=T) -> files
  if (length(files) == 1 && length(grep("heatmap", files)) == 1){
    load(files[grep("heatmap", files)], envir = parent.frame()) #assembled
  } else {
    # check if there are files that follow the format Day MM DD Time YY
    paths <- dir("./r.data.files", pattern=".rda", full.names=TRUE)
    paths[which(tail(file.info(paths)$ctime) %in% max(tail(file.info(paths)$ctime)))] -> f
    load(f, envir = parent.frame()) # assembled
  }
  
  # do the same for the fit summary files
  list.files("./r.data.files/rawFits", pattern=".rda", full.names=T) -> files
  if (length(files) == 1){
    load(files, envir = parent.frame()) #summary
  } else {
    # check if there are files that follow the format Day MM DD Time YY
    paths <- dir("./r.data.files/rawFits", pattern=".rda", full.names=TRUE)
    paths[which(tail(file.info(paths)$ctime) %in% max(tail(file.info(paths)$ctime)))] -> f
    load(f, envir = parent.frame()) #summary
  }
}

addNewDrugs <- function(res.df, assembled, drug.list.all, patient.name){
  # Need R interface to ChEMBL API; as this is a dev package, need to add
  # auto-install options for users. Currently, the ChEMBL API does not have
  # full functionality, so the solution is to linked to a local forked version
  # of chemblr
  if ("devtools" %in% rownames(installed.packages()) && !is.loaded("devtools")){
    library(devtools)
    if ("chemblr" %in% rownames(installed.packages()) && !is.loaded("chemblr")){
      library(chemblr)  
    }
  } else if ("devtools" %ni% rownames(installed.packages())){
    # if devtools is not installed, it's guaranteed the dev package is not installed
    print("Installing required packages...")
    install.packages("devtools")
    if ("chemblr" %ni% rownames(installed.packages())){
      install_github("rajarshi/chemblr/package")
    }
  }
  
  .matchDrugs(res.df, assembled, drug.list.all) -> res
  newDrugs <- matrix(NA, nrow=nrow(assembled), ncol=length(res$unmatched))
  for (i in 1:length(res$unmatched)){
    newDrugs[which(rownames(assembled) %in% patient.name),i] <- 
      res$corrected[1,which(colnames(res$corrected) %in% res$unmatched[i])]
  }
  colnames(newDrugs) <- res$unmatched
  cbind(assembled, newDrugs) -> assembledNew
  return(assembledNew)
}

collectFit <- function(ares, drug.list.all, patient.name=""){
  # automatically check the available files in the ./r.data.files subdirectory; change
  # routine to create a subdirectory with derived, i.e. data frame files created using
  # the condenseToDF() function and also save fit files (rawFits).
  # This function takes the average of fits, and saves the whole ares file together
  # with the summary files; also appends the newest fit to the heatmap file
  # `assembled'
  print("Getting newest fit files...")
  getNewestFile()
  t(sapply(ares, function(x) x$logIC50)) -> res.df
  
  # explicitly add patient name and check if has the LK code
  if (length(grep("LK", patient.name, ignore.case = T)) == 0){
    stop("Only primary patient data can be added to the ALL record. 
         The name that you specified is inconsistent with patient information.")
  }
  
  if (patient.name %ni% rownames(assembled)){
    .mapResponse(res.df, assembled, drug.list.all, patient.name) -> assembled.sub
    
    # merge assembled, assembled.sub
    assembled.sub[which(rownames(assembled.sub) %ni% rownames(assembled)),] -> pat
    match(colnames(pat), colnames(assembled)) -> ind.match
    ind.match[which(ind.match %ni% NA)] -> ind.match
    rep(NA, ncol(assembled)) -> fin
    fin[ind.match[which(ind.match %ni% NA)]] <- as.numeric(as.character(pat[1,2:ncol(pat)]))
    rbind(assembled, fin) -> assembled
    rownames(assembled)[nrow(assembled)] <- patient.name
    
    # then append any information from new drugs/unmatched drugs
    addNewDrugs(res.df, assembled, drug.list.all, patient.name) -> assembled
    save(assembled, file=paste("./r.data.files/", date(), ".rda", sep=""))
    
    # additionally, keep a copy of the fit in a summary file
    # of the form summary$`PATIENT`$DRUG; ensure that drug names match
    # the standard name when available; deals with cases where only one patient is added
    if (length(ares) > 1){
      stop("Add patient fits one by one.")
    } else {
      lapply(ares, function(x) 
        lapply(x$max, function(y) extractMax(y))) -> pat.fits
      gsub("\\.", "-", toupper(names(pat.fits[[1]]))) -> names(pat.fits[[1]])
      findClosestMatch(toupper(names(pat.fits[[1]])), drug.list.all) -> drugMatches
      drugMatches$match -> names(pat.fits[[1]])
      
      # add to summary
      if (names(pat.fits) %ni% names(summary)){
        c(summary, pat.fits) -> summary
      }
    }
  } else if (patient.name %in% rownames(assembled)){
    warning("Patient record already exists in this file. No record added.")
  }
}