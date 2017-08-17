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
    load(files[grep("heatmap", files)]) #assembled
  } else {
    # check if there are files that follow the format Day MM DD Time YY
    paths <- dir("./r.data.files", pattern=".rda", full.names=TRUE)
    paths[which(tail(file.info(paths)$ctime) %in% max(tail(file.info(paths)$ctime)))] -> f
    load(f) # assembled
  }
}

addFit <- function(res.df, drug.list.all, patient.name=""){
  # automatically check the available files in the ./r.data.files subdirectory
  getNewestFile()
  
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
    match(colnames(pat), colnames(assembled)) -> ind.match
    ind.match[which(ind.match %in% NA)] -> ind.match
    
    save(assembled, file=paste("./r.data.files/", date(), ".rda", sep=""))
  } else if (patient.name %ni% rownames(assembled)){
    warning("Patient record already exists in this file. No record added.")
  }
}