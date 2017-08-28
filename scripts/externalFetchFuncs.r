# PDobay for AG Bourquin, KISPI
# pamela.dobay@kispi.uzh.ch

#@...
# Routines that extract information for a drug from resources
# linked to clinical trials
#@...

# Use of ClinicalTrials.gov data is subject to these Terms and Conditions. 
# If you link to this site, please provide proper attribution to NLM and 
# ClinicalTrials.gov.

library(rclinicaltrials) #devtools

getTrials <- function(query, otherTerms="", maxK=10){
  # wrapper for .getTrials
  .getTrials(query, otherTerms, maxK) -> tids
  formatNCTIDs(tids) -> res
  return(res)
}

.getTrials <- function(query, otherTerms="", maxK=10){
  # given a query (a drug or a vector of drugs), return information on clinical trials, where available; 
  # in case the `otherTerms' parameter is not blank, add this option to the queries
  sapply(query, function(x){
    # note that clinicaltrials_search on a case where no results are available results in an
    # error; to handle this, the query will be restricted to the drug name and will be filtered
    # based on the other search terms at a later point
    y <- clinicaltrials_search(query = x)
    tid <- NA
    if (nrow(y) > 0){
      #y <- clinicaltrials_download(query = paste(x, paste(otherTerms, collapse="+"), sep="+"), 
      #                          count = maxK, include_results = TRUE)  
      y <- clinicaltrials_download(query = x, 
              count = maxK, include_results = TRUE)
      if (nrow(y[[1]][[1]]) > 0){
        # focus on hematological malignancies, currently hardcoded, but will be changed
        # when a cloned version of the package is used instead
        unique(c(grep("leukemia", y[[1]][[1]]$primary_condition, ignore.case = T),
          grep("lymphoma", y[[1]][[1]]$primary_condition, ignore.case = T))) -> ind
        y[[1]][[1]][ind,] -> fin
        if (length(ind) > 0){
          tid <- fin$nct_id
        } else {
          # return all results anyway
          tid <- y[[1]][[1]]$nct_id
        }
      }
    }
    return(tid)
  }) -> tids
  return(tids)
}

formatNCTIDs <- function(tids){
  # given some ids from clinicaltrials.gov
  sapply(tids, function(x)
    paste(as.character(sapply(x, function(y)
      paste("<a href=", paste("https://clinicaltrials.gov/ct2/show/", y, sep=""), ">", y, "<\a>", sep="")
    )), collapse=",") -> ids
  ) -> ids
  return(ids)
}