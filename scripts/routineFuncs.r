`%ni%` <- Negate(`%in%`)

start <- function(){
  # start saving all results from a session in a single folder
  file.path("./results", format(Sys.time(), "%F %H-%M")) ->> rd
  dir.create(rd)
}

toFile <- function(ic50){
  write.csv(ic50, paste(rd, "/raw_vals.csv", sep=""))
}

end <- function(){
  rm(rd)
}