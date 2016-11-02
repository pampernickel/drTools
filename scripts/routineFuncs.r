`%ni%` <- Negate(`%in%`)
library(RCurl)

start <- function(){
  # start saving all results from a session in a single folder
  #file.path("./results", format(Sys.time(), "%F %H-%M")) ->> rd
  file.path("./results", Sys.time()) ->> rd
  dir.create(rd)
}

toFile <- function(ic50){
  write.csv(ic50, paste(rd, "/raw_vals.csv", sep=""))
}

end <- function(){
  rm(rd)
}

loadAllDependencies <- function(){
  print("Loading scripts from pampernickel/drTools...")
  getURL('https://raw.githubusercontent.com/pampernickel/drTools/master/scripts/analysis.Funcs.R', ssl.verifypeer = F) -> script
  eval(parse(text=script))
  getURL('https://raw.githubusercontent.com/pampernickel/drTools/master/scripts/fitting.functions.R', ssl.verifypeer = F) -> script
  eval(parse(text=script))
  getURL('https://raw.githubusercontent.com/pampernickel/drTools/master/scripts/plateReader.R', ssl.verifypeer = F) -> script
  eval(parse(text=script))
  getURL('https://raw.githubusercontent.com/pampernickel/drTools/master/scripts/visFuncs.R', ssl.verifypeer = F) -> script
  eval(parse(text=script))
  getURL('https://raw.githubusercontent.com/pampernickel/drTools/master/scripts/nomenclatureFuncs.r', ssl.verifypeer = F) -> script
  eval(parse(text=script))
  
  library(ggplot2)
  library(reshape2)
  print("Done.")
}