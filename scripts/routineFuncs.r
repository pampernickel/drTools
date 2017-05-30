library(RCurl)

`%ni%` <- Negate(`%in%`)

getSlopes <- function(y.dat){
  # get slopes for all pairs of points in a fit
  rbind(y.dat[1:(length(y.dat)-1)], y.dat[2:length(y.dat)]) -> pts
  slopes <- apply(pts, 2, 
                  function(x) lm(c(x[1],x[2])~c(1,2))$coefficients[2])
  return(slopes)
}

isDescending <- function(x){
  res <- F
  if (all(diff(x) >= 0) == F){
    res <- T
  }
  return(res)
}

start <- function(){
  # start saving all results from a session in a single folder
  #file.path("./results", Sys.time()) ->> rd
  
  # for the report mode, use a fixed name file name
  #dir.create("./temp")
}

toFile <- function(ic50, filename=""){
  if (filename %in% ""){
    write.csv(ic50, paste(rd, "/raw_vals.csv", sep=""))
  } else {
    write.csv(ic50, paste(rd, "/", filename, ".csv", sep=""))
  }
}

# https://www.r-bloggers.com/source_https-sourcing-an-r-script-from-github-over-https/
source_https <- function(url, ...) {
  # load package
  # parse and evaluate each .R script
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}

loadComboDependencies <- function(){
  source_https("https://raw.githubusercontent.com/pampernickel/drTools/master/scripts/comboReader.r")
  print("Loading scripts from xtmgah/DDCV...")
  source_https("https://raw.githubusercontent.com/pampernickel/DDCV/master/DDCV_function/IC50.R")
  source_https("https://raw.githubusercontent.com/pampernickel/DDCV/master/DDCV_function/isobologram.R")
  source_https("https://raw.githubusercontent.com/pampernickel/DDCV/master/DDCV_function/shapeA.R")
  source_https("https://raw.githubusercontent.com/pampernickel/DDCV/master/DDCV_function/cIndex.R")
  source_https("https://raw.githubusercontent.com/pampernickel/DDCV/master/DDCV_function/cIndex2.R")
  source_https("https://raw.githubusercontent.com/pampernickel/DDCV/master/DDCV_function/dContour.R")
  print("Done.")
}

loadAllDependencies <- function(){
  print("Loading scripts from pampernickel/drTools...")
  source_https('https://raw.githubusercontent.com/pampernickel/drTools/master/scripts/plateReader.R')
  print("Loading I/O tools...")
  source_https('https://raw.githubusercontent.com/pampernickel/drTools/master/scripts/fitting.functions.R')
  source_https('https://raw.githubusercontent.com/pampernickel/drTools/master/scripts/outputFuncs.r')
  print("Loading fitting tools...")
  source_https('https://raw.githubusercontent.com/pampernickel/drTools/master/scripts/analysis.Funcs.R')
  print("Loading analysis tools...")
  source_https('https://raw.githubusercontent.com/pampernickel/drTools/master/scripts/visFuncs.R')
  print("Loading visualization tools...")
  source_https('https://raw.githubusercontent.com/pampernickel/drTools/master/scripts/nomenclatureFuncs.r')
  print("Loading nomenclature tools...")
  loadComboDependencies()
  library(ggplot2)
  library(reshape2)
  print("Done.")
}