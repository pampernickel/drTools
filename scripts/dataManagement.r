library(xlsx)

# try to retrieve any annotations from mic_backup first, before moving on to local
# copy

loadDrugList <- function(){
  drug.list.all <- NA
  if (dir.exists("/Volumes/data/drTools/drugList")){
    files <- list.files("/Volumes/data/drTools/drugList", 
                        pattern=".xlsx", full.names = T)
    if (length(files) > 1){
      stop("There is more than one excel file in the directory. 
           Please remove whatever file has no business here.")
    }
    read.xlsx(files, sheetIndex = 1, header=T, stringsAsFactors = F) -> drug.list.all
  } else {
    # try from local
    drug.list.all <- read.csv("./annotations/drug.list.final_format.csv", stringsAsFactors = F)
  }
  return(drug.list.all)
}
# add files to backup r data files to mic_backup or any appropriate backup
