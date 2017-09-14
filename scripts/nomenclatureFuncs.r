# set of functions for making names uniform
# standard Format: data frame with samples on rows and drugs on columns

require(stringr)
require(stringdist)

'%ni%' <- Negate('%in%')

format <- function(exp.res, attribute=c("logIC50", "logEC10", "logEC50", "logEC90",
                                        "AUC", "max"), drug.list){
    if (attribute %ni% c("logIC50", "logEC10", "logEC50", "logEC90",
                         "AUC", "max")){
      stop("Unknown attribute.")
    } else {
      which(names(exp.res$res[[1]]) %in% attribute) -> ind
      t(as.data.frame(exp.res$res[[1]][[ind]])) -> df
      print(dim(df))
      colnames(df) <- drug.list
      rownames(df) <- names(exp.res$res)
      return(df)
    }
}

formatDrugNames <- function(res.df, drug.list.all){
  # convert names to standard nomenclature
  sapply(colnames(res.df), function(x) 
    ifelse(length(grep(x, drug.list.all$all.names, ignore.case=T)) > 0, 
           as.character(drug.list.all$final.name[grep(x, drug.list.all$all.names, ignore.case=T)]), 
           NA)) -> drug.names.fin
  if (length(which(is.na(drug.names.fin))) > 0){
    warning("Several drug names not matched.")
  }
  
  drug.names.fin[which(is.na(drug.names.fin))] <- toupper(colnames(res.df)[which(is.na(drug.names.fin))])
  colnames(res.df) <- drug.names.fin
  return(res.df)
}

findClosestMatch <- function(drug.name, drug.list.all){
  # input: character vector of drug names
  # output: standard name of the drug, as indicated in the "final.name" column of drug.list.all
  strsplit(drug.list.all$all.names, ";") -> all.names
  no.match.fin <- NA
  sapply(drug.name, function(x){
    match <- NA
    toupper(x) -> x
    if (toupper(x) %in% toupper(drug.list.all$final.name)){
      x -> match
    } else if (length(grep(x, drug.list.all$all.names, ignore.case=T)) > 0) {
      if (length(grep(x, drug.list.all$all.names, ignore.case=T)) == 1){
        grep(x, drug.list.all$all.names, ignore.case=T) -> ind
        drug.list.all$final.name[ind] -> match
      } else {
        # get closer match
        #str_match(string, pattern)
      }
    } else {
      # get closest match, with a minimum distance requirement; again for
      # addressing typos
      sapply(drug.list.all$final.name, function(y) adist(y, x, counts=F)) -> d
      # check top candidates, based on minimum distance
      drug.list.all$final.name[which(d %in% min(d))] -> pmatch
      if (length(pmatch) > 0 && min(d) <= 2){
        unlist(strsplit(x, "")) -> all.char
        sapply(pmatch, function(z) 
          which(length(which(all.char %in% unlist(strsplit(z, "")))) == length(all.char))) -> check
        unlist(check) -> check
        if (length(check) > 0){
          pmatch[which(check == 1)] -> match
        }
      }
    }
  }) -> match
  
  # then tag cases where the drug might be a new drug; these will be added as a new drug in the critical files
  # as well as the drug.list flat file
  which(sapply(match, function(x) length(x)) == 0) -> no.match
  if (length(no.match) > 0){
    names(no.match) -> no.match.fin
  }
  
  as.character(match) -> match
  match[no.match] <- drug.name[no.match]
  res <- list(match, no.match.fin)
  names(res) <- c("match", "new.drug")
  return(res)
}