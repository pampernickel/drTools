# set of functions for making names uniform
# standard Format: data frame with samples on rows and drugs on columns

require(stringr)

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
  sapply(drug.name, function(x){
    match <- NA
    if (toupper(x) %in% drug.list.all$final.name){
      toupper(x) -> match
    } else if (length(grep(x, drug.list.all$all.names, ignore.case=T)) > 0) {
      if (length(grep(x, drug.list.all$all.names, ignore.case=T)) == 1){
        grep(x, drug.list.all$all.names, ignore.case=T) -> ind
        drug.list.all$final.name[ind] -> match
      } else {
        # get closer match
      }
    } else {
      # get closest match, but give a warning that it's second guessing
    }
  }) -> match
  return(match)
}