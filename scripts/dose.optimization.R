ALL.60 <- read.csv("~/Documents/SIB/ALL.project/ALL.pipeline/ALL.60.csv", stringsAsFactors = TRUE)
ALL.60 -> IC50

for (i in 1:ncol(IC50)){
  as.numeric(IC50[,i]) -> IC50[,i]
}

ALL.60$Name -> rownames(IC50)
IC50[,-1] -> IC50 
optimizeDoses(IC50, 4, 8) -> doses.4

