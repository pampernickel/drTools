# optimize dose response by performing consensus clustering for the maximum
# no. of k (currently dependent on the number of plates/unique dose combinations)
# that can be handled

# group samples according to the similarity of their response profiles,
library(pracma)

# then get the spread of the data
# the number of equidistant dose points can then be determined based on
# the data spread of responders.

# also contains the visualization function
optimizeDoses <- function(IC50, k, n){
  # IC50: contains matrix of IC50s associated per sample per drug
  # assumed to contain only a list of drugs for which the response
  # patterns are interesting for a certain subset of patients
  # do not contain drugs that are not reactive in 90% of the cohort
  
  # k: no of plates
  # n: no of doses
  # first, find groups based on IC50
  IC50 <- t(IC50)
  breaks <- as.data.frame(matrix(nrow=ncol(IC50), ncol=n))
  for (i in 1:ncol(IC50)){
    hist(IC50[,i], plot =  FALSE, breaks=seq(min(IC50[,i]), max(IC50[,i]), length=8)) -> temp
    breaks[i,] <- temp$breaks
  }
  rownames(breaks) <- colnames(IC50)
  getClusters(breaks, k) -> breaks
  
  # for each category, get range of response, equidistant points for response
  unique(breaks$fit.cluster) -> groups
  doses <- matrix(nrow=length(groups),ncol=n)
  for (i in 1:length(unique(breaks$fit.cluster))){
      as.numeric(as.character(unlist(c(breaks[which(breaks$fit.cluster %in% groups[i]),c(1:8)])))) -> range
      #hist(range, plot=TRUE, breaks=seq(min(range), max(range), length=8)) -> temp
      (max(range)-min(range))/(n-1) -> increment
      doses[i,] <- round(seq(from=min(range), to=max(range), by=increment),1)
  }
  
  as.data.frame(doses) -> doses
  for (i in 1:n){
    colnames(doses)[i] <- paste("dose",i,sep=".")
  }
  
  for (i in 1:k){
    rownames(doses)[i] <- paste("group",groups[i],sep=".")
  }
  cbind(as.character(rownames(breaks)), as.character(breaks$fit.cluster)) -> drug.groupings
  doseSets <- list(drug.groupings, doses)
  return(doseSets)
}

getClusters <- function(mydata,k){
  wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
  for (i in 2:k) wss[i] <- sum(kmeans(mydata,centers=i)$withinss)
  #plot(1:6, wss, type="b", xlab="Number of Clusters",
  #ylab="Within groups sum of squares")
  fit <- kmeans(mydata, k)
  aggregate(mydata,by=list(fit$cluster),FUN=mean)
  mydata <- data.frame(mydata, fit$cluster)
  mydata[order(mydata$fit.cluster),] -> mydata
  return(mydata)
}

getICL <- function(matrix, k, lim, dir){
  d = sweep(matrix,1, apply(matrix,1,median,na.rm=T))
  results = ConsensusClusterPlus(d,maxK=k,reps=lim,pItem=0.8,pFeature=1, title=dir, clusterAlg="hc", distance="pearson",seed=1262118388.71279,plot="pdf", writeTable =TRUE)
  icl = calcICL(results,title=paste(dir,"icl",sep="."),plot="pdf")
  list(results, icl) -> res
  return(res)
}

getOptK <- function(resICL, lim){
  cdf <- c()
  for (i in 2:lim){
    # --- get CDFs, consensus index based on the histogram of the
    # -- consensus matrices for i in 2 to 10
    hist(resICL[[1]][[i]]$consensusMatrix, 100) -> temp
    auc <- trapz(cumsum(temp$density), c(1:length(cumsum(temp$density))))
    c(cdf, auc) -> cdf
  }
  
  cdf/max(cdf) -> cdf
  abs(diff(cdf)) -> diff.cdf
  abs(diff(diff.cdf)) -> diff.cdf
  opt.k <- min(which(diff.cdf < 0.05))+1
  return(opt.k)
}

assignClust <- function(matrix,resICL,significance,lim){
  cbind(lapply(split(resICL$itemConsensus[which(resICL$itemConsensus$k==lim),4], resICL$itemConsensus[which(resICL$itemConsensus$k==lim),3]), which.max), lapply(split(resICL$itemConsensus[which(resICL$itemConsensus$k==lim),4], resICL$itemConsensus[which(resICL$itemConsensus$k==lim),3]), max))->clustTable
  
  unlist(clustTable[which(unlist(clustTable[,2])>0.0),1])->all.clust
  rownames(clustTable)[which(unlist(clustTable[,2])>0.0)] -> all.names
  unlist(clustTable[-which(unlist(clustTable[,2])<significance),1])->sub.members
  rownames(clustTable)[-which(unlist(clustTable[,2])<significance)] -> sub.names
  
  ind.core <- which(colnames(matrix) %in% sub.names)
  groups <- c(1:length(colnames(matrix)))
  groups[1:length(colnames(matrix))] <- lim+1
  
  for (i in 1:lim){
    groups[which(colnames(matrix) %in% all.names[which(all.names %in% sub.names[which(sub.members %in% i)])])] <- i
  }
  
  return(groups)
}
