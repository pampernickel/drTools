`%ni%` <- Negate(`%in%`)
library(RCurl)

start <- function(){
  # start saving all results from a session in a single folder
  #file.path("./results", Sys.time()) ->> rd
  
  # for the report mode, use a fixed name file name
  dir.create("./temp")
}

toFile <- function(ic50){
  write.csv(ic50, paste(rd, "/raw_vals.csv", sep=""))
}


# https://www.r-bloggers.com/source_https-sourcing-an-r-script-from-github-over-https/
source_https <- function(url, ...) {
  # load package
  # parse and evaluate each .R script
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}

loadAllDependencies <- function(){
  
  print("Loading scripts from pampernickel/drTools...")
  source_https('https://raw.githubusercontent.com/pampernickel/drTools/master/scripts/plateReader.R')
  print("Loading I/O tools...")
  source_https('https://raw.githubusercontent.com/pampernickel/drTools/master/scripts/fitting.functions.R')
  print("Loading fitting tools...")
  source_https('https://raw.githubusercontent.com/pampernickel/drTools/master/scripts/analysis.Funcs.R')
  print("Loading analysis tools...")
  source_https('https://raw.githubusercontent.com/pampernickel/drTools/master/scripts/visFuncs.R')
  print("Loading visualization tools...")
  source_https('https://raw.githubusercontent.com/pampernickel/drTools/master/scripts/nomenclatureFuncs.r')
  print("Loading nomenclature tools...")
  
  library(ggplot2)
  library(reshape2)
  print("Done.")
}