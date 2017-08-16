# PDobay for AG Bourquin, KISPI
# pamela.dobay@kispi.uzh.ch

#@...
# Routine that collects and saves drug response information from all runs
# using patient material (i.e. with LK200x tag); would currently go for the
# solution of storing in an r data file (per patient) and in the 20161012.heatmap.all.fits.rda
# matrix, which contains the response data of the blood paper cohort
#@...

addFit <- function(assembled, res.df, patient.name=""){
  # explicitly add patient name and check if has the LK code
  if (length(grep("LK", patient.name, ignore.case = T) == 0)){
    stop("Only primary patient data can be added to the ALL record. 
         The name that you specified is inconsistent with patient information.")
  }
  
  .mapResponse(res.df, assembled, drug.list.all, patient.name) -> assembled.sub
  if (getData ==T && patient.name %ni% rownames(assembled)){
    # merge assembled, assembled.sub
    assembled.sub[which(rownames(assembled.sub) %ni% rownames(assembled)),] -> pat
    match(colnames(pat), colnames(assembled)) -> ind.match
    ind.match[which(ind.match %ni% NA)] -> ind.match
    rep(NA, ncol(assembled)) -> fin
    fin[ind.match[which(ind.match %ni% NA)]] <- as.numeric(as.character(pat[1,2:ncol(pat)]))
    rbind(assembled, fin) -> assembledNew
    rownames(assembled) <- patient.name
    save(assembledNew, file=paste("./r.data.files/", date(), ".rda", sep=""))
  } if (getData ==T && pat.name %ni% rownames(assembled)){
    warning("Patient record already exists in this file. No record added.")
  }
}