library(xlsx)

# try to retrieve any annotations from mic_backup first, before moving on to local
# copy

loadDrugList <- function(){
  drug.list.all <- NA
  if (dir.exists("/Volumes/data/drTools/drugList")){
    files <- list.files("/Volumes/data/drTools/drugList", 
                        pattern=".xlsx", full.names = T)
    if (length(files) > 1){
      # get largest file
      files[which(file.size(files) %in% max(file.size(files)))] -> files
      warning("There is more than one excel file in the directory. 
           The largest file in the directory was loaded into the workspace.")
    }
    read.xlsx(files, sheetIndex = 1, header=T, stringsAsFactors = F) -> drug.list.all
  } else {
    # try from local
    drug.list.all <- read.csv("./annotations/drug.list.final_format.csv", stringsAsFactors = F)
  }
  return(drug.list.all)
}
# add files to backup r data files to mic_backup or any appropriate backup

getNewestFile <- function(){
  list.files("./r.data.files", pattern=".rda", full.names=T) -> files
  if (length(files) == 1 && length(grep("heatmap", files)) == 1){
    load(files[grep("heatmap", files)], envir = parent.frame()) #assembled
  } else {
    # check if there are files that follow the format Day MM DD Time YY
    # modify to extract most recent based on file creation timestamp
    paths <- dir("./r.data.files", pattern=".rda", full.names=TRUE)
    paths[which(file.info(paths)$mtime %in% max(file.info(paths)$mtime))] -> f
    load(f, envir = parent.frame()) # assembled
  }

  # do the same for the fit summary files
  list.files("./r.data.files/rawFits", pattern=".rda", full.names=T) -> files
  if (length(files) == 1){
    load(files, envir = parent.frame()) #summary
  } else {
    # check if there are files that follow the format Day MM DD Time YY
    paths <- dir("./r.data.files/rawFits", pattern=".rda", full.names=TRUE)
    paths[which(file.info(paths)$mtime %in% max(file.info(paths)$mtime))] -> f
    load(f, envir = parent.frame()) #summary
  }
}