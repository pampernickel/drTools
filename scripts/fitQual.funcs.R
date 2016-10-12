# --- functions to check the quality of R data fits vs PRISM data
compareIC50 <- function(prismIC50, exp.res.ave, drug.list, match.type){
  
  # --- create 2-col matrix with prism and drc values
  prism.c <- c()
  drc.c <- c()
  patients <- c()
  drugs <- c()
  if (match.type == "greedy"){
    for (i in 1:nrow(prismIC50)){
      # --- first, check if patient is on exp.res.ave
      gsub("X", "", rownames(prismIC50)) -> names
      names[i] -> curr.patient
      # --- find IC50 values per drug in prismIC50 and in exp.res.ave
      grep(curr.patient, names(exp.res.ave)) -> matches
      if (length(matches) == 1){
        exp.res.ave[[matches]][[1]] -> IC50.calc
        # go through drug.list, which has shorter names,
        # and find equivalent entries in prismIC50
        for (j in 1:length(drug.list)){
          if (length(grep(drug.list[j], colnames(prismIC50))) == 1){
            # --- just do for exact matches first
            prism.c <- c(prismIC50[i,grep(drug.list[j], colnames(prismIC50))], prism.c)
            drc.c <- c(IC50.calc[j], drc.c)
            patients <- c(as.character(curr.patient), patients)
            drugs <- c(as.character(drug.list[j]), drugs)
          }
        }
      }
    }
  } else if (match.type == "exact"){
    for (i in 1:nrow(prismIC50)){
      #gsub("X", "", rownames(prismIC50)) -> names
      prismIC50$X[i] -> curr.patient
      grep(prismIC50$X[i], names(exp.res.ave)) -> matches
      if (length(matches) == 1){
        exp.res.ave[[matches]][[1]] -> IC50.calc
        # go through drug.list, which has shorter names,
        # and find equivalent entries in prismIC50
        for (j in 1:length(drug.list)){
          if (length(grep(drug.list[j], colnames(prismIC50))) == 1){
            # --- just do for exact matches first
            prism.c <- c(as.numeric(as.character(prismIC50[i,grep(drug.list[j], colnames(prismIC50))])), prism.c)
            drc.c <- c(IC50.calc[j], drc.c)
            patients <- c(as.character(curr.patient), patients)
            drugs <- c(as.character(drug.list[j]), drugs)
          }
        }
      }
    }
  }
  
  abs(prism.c - drc.c) -> diff
  cbind(prism.c, drc.c, diff, patients, drugs) -> fit.comparison
  # sort according to diff
  fit.comparison[order(-diff),] -> fit.comparison
  #View(cbind(prism.c, drc.c, patients, drugs))
  png("Correl.plot.png", w=460, h=460)
  plot(drc.c~prism.c, xlab="Ext IC50", ylab="DRC IC50")
  dev.off()
  return(fit.comparison)
}