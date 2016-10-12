
signif.test <- function(x, drug, test=c("cor", "t.test"), patients){
  test.res <- data.frame(matrix(NA, ncol = 5, nrow = 0))
  combos <- combn(5,2)
  for (i in 1:length(x)){
    x[[i]][,grep(drug,colnames(x[[i]]))] -> t
    # --- perform pairwise test of controls
    for (j in 1:ncol(combos)){
      grep(paste(drug,combos[1,j],sep="."), colnames(t)) -> ind.1
      grep(paste(drug,combos[2,j],sep="."), colnames(t)) -> ind.2
      as.numeric(as.matrix(t[,ind.1])) -> a
      as.numeric(as.matrix(t[,ind.2])) -> b
      
      if (test == "cor"){
        if (length(which(is.na(a))) || length(which(is.na(b)))){
          print("Warning: NAs detected. Correlation test being skipped for entries...")
        } else {
          cor.test(a, b, method="s") -> cor
          res <- cbind(as.character(patients[i]), as.character(combos[1,j]), as.character(combos[2,j]), as.character(cor$estimate), as.character(cor$p.value))
          rbind(test.res, res) -> test.res
        }
      } else {
        t.test(a, b) -> t.test.res
        res <- cbind(as.character(patients[i]), as.character(combos[1,j]), as.character(combos[2,j]), as.character(t.test.res$statistic), as.character(t.test.res$p.value))
        rbind(test.res, res) -> test.res
      }
    }
  }
  
  if (test == "cor"){
    colnames(test.res) <- c("Patient", "Plate x", "Plate y", "rho", "p.val")
  } else {
    colnames(test.res) <- c("Patient", "Plate x", "Plate y", "t-stat", "p.val")
  }
  return(test.res)
}

plot.cor <- function(test.res, drug, patients){
  for (i in 1:nrow(test.res)){
    j <- which(patients %in% test.res$Patient[i])
    x[[j]][,grep(drug,colnames(x[[j]]))] -> t
    grep(paste(drug,test.res[i,2],sep="."), colnames(t)) -> ind.1
    grep(paste(drug,test.res[i,3],sep="."), colnames(t)) -> ind.2
    as.numeric(as.matrix(t[,ind.1])) -> a
    as.numeric(as.matrix(t[,ind.2])) -> b
    plot.res <- list(a,b)
    return(plot.res)
  }
}

plot.selection <- function(test.res, list, drug, patients){
  for (i in 1:nrow(test.res)){
    j <- which(patients %in% test.res$Patient[i])
    list[[j]][,grep(drug,colnames(list[[j]]))] -> t
    grep(paste(drug,test.res[i,2],sep="."), colnames(t)) -> ind.1
    grep(paste(drug,test.res[i,3],sep="."), colnames(t)) -> ind.2
    as.numeric(as.matrix(t[,ind.1])) -> a
    as.numeric(as.matrix(t[,ind.2])) -> b
    par(mfrow=c(1,2))
    print(a)
    print(b)
    x.ax <- c(10000,1000,100,10,1, 10000,1000,100,10,1)
    plot(a~x.ax, log="x")
    plot(b~x.ax, log="x")
    #return(plot.res)
  }
}

checkReplicates <- function(list, drugs, patients, rho.min){
  res <- list()
  for (i in 1:length(list)){
    # per patient, check the correlation between the dose responses per drug
    print(paste("Patient", patients[i], sep=" "))
    names <- c()
    if (is.data.frame(list[[i]])){
      checkPerDrug(df, drugs) -> res[[i]]
      #res[[i]] <- names
    } else if (is.list(list[[i]])){
      # --- perform extra layer
      names <- list()
      for (j in 1:length(list[[i]])){
        list[[i]][[j]] -> df
        checkPerDrug(df, drugs) -> names[[j]]
      }
      names -> res[[i]]
    }
  }
  return(res)
}

# --- extra function to handle 5-pt and 8-pt data organization of experiments
checkPerDrug <- function(df, drugs){
  # --- restrict drugs to those that match between the data frame
  # --- and the list of drugs
  substr(colnames(df), 1, nchar(colnames(df))-2) -> colnames.sub
  unique(drugs[which(drugs %in% colnames.sub)]) -> drugs
  names <- c()
  for (j in 1:length(drugs)){
    df[,grep(drugs[j],colnames(df))] -> t
    as.numeric(gsub(",", "", t[,1])) -> t[,1]
    as.numeric(gsub(",", "", t[,2])) -> t[,2]
    if (length(which(is.na(t))) == 0){
      cor.test(as.numeric(t[,1]), as.numeric(t[,2]),method="p") -> cor.res
      cor.res$p.value -> p.val
      cor.res$estimate -> rho
      if (p.val > 0.05 & abs(rho) < rho.min){
        print(drugs[j])
        names <- c(as.character(drugs[j]), names)
      }
    }
  }
  return(names)
}

# --- visualize a result
visres <- function(x, res, patients){
  dir.create(paste(getwd(),"/quality.control", sep=""))
  getwd() -> old.dir
  setwd(paste(getwd(),"/quality.control", sep=""))
  for (i in 1:length(res)){
    res[[i]] -> drugs.i
    print(paste("Length of drugs for patient ", patients[i], ":", length(drugs.i)))
    
    logd <- c(4,3,2,1,0)
    x1 <- logd[rev(order(logd))]
    rows <- ceiling(length(drugs.i)/5)
    if (rows > 0){
      png(paste("Patient",patients[i],".png",sep=""), w=3000, h=3000)
      par(mfrow=c(rows,5))
      
      for (j in 1:length(drugs.i)){
        x[[i]][,grep(drugs.i[j],colnames(x[[i]]))] -> t
        as.numeric(gsub(",", "", t[,1])) -> t[,1]
        as.numeric(gsub(",", "", t[,2])) -> t[,2]
        plot(x1, as.numeric(t[,1])/max(as.numeric(t[,1])), col="red", type="l", main=drugs.i[j], xlab="log conc", ylab="% live, norm to 1", xlim=c(0,4), ylim=c(0,1), cex.main=5, cex.axis=4, cex.lab=4, lwd=3)
        lines(x1, as.numeric(t[,2])/max(as.numeric(t[,2])), col="blue", lwd=3, cex.main=5, cex.axis=4, cex.lab=4, ylim=c(0,1))
      }
      dev.off()
    }
  }
  setwd(old.dir)
}

# --- visualize data scatter
scatter <- function(x, drug){
  t.all <- c()
  patient.list <- c()
  plate.list <- c()
  for (i in 1:length(x)){
    # --- within each patient, get the distribution of values for a drug
    patients[i] -> curr.patient
    x[[i]][,grep(drug,colnames(x[[i]]))] -> t
    as.numeric(gsub(",", "", t[,1])) -> t[,1]
    as.numeric(gsub(",", "", t[,2])) -> t[,2]
    c(t[,1], t[,2]) -> temp
    t.all <- c(temp, t.all)
    c(colnames(t), plate.list) -> plate.list
    c(rep(curr.patient, length(temp)), patient.list) -> patient.list
  }
  
  for (i in 1:length(plate.list)){
    plate.list[i] <- substr(plate.list[i],6,6)
  }
  
  rbind(t.all,patient.list,plate.list) -> df
  t(df) -> df
  as.data.frame(df) -> df
  return(df)
}

# --- visualize the plot for the dose response of a particular drug
# --- for all of the patients ("all") or for a list of patients:
# --- vector of patient names
showResponse <- function(drug.name, experiments, patients){
  if (patients %in% "all"){
    drug.response <- matrix(nrow=0, ncol=3)
    for (i in 1:length(experiments$normalized)){
      curr.exp <- experiments$normalized[[i]]
      grep(drug.name, colnames(curr.exp)) -> ind
      curr.exp[,ind] -> dr
      apply(dr,1,mean)->mean.dr
      round(mean.dr,2) -> mean.dr
      rep(names(experiments$normalized)[i], length(mean.dr)) -> patient
      cbind(curr.exp[,1], mean.dr, patient) -> temp
      rbind(temp, drug.response) -> drug.response
      as.data.frame(drug.response) -> drug.response
    }
    
    log10(as.numeric(as.character(drug.response$V1))) -> drug.response$V1
    as.numeric(as.character(drug.response$mean.dr)) -> drug.response$mean.dr
    ggplot(drug.response, aes(x=V1, y=mean.dr, group=patient))+geom_line(aes(colour=patient))+geom_point(aes(colour=patient))+ggtitle(drug.name)
  } else {
    
  }
}
# ---
# http://www.stat.sc.edu/~hitchcock/chapter8_R_examples.txt
cancor2<-function(x,y,dec=4){
  #Canonical Correlation Analysis to mimic SAS PROC CANCOR output.
  #Basic formulas can be found in Chapter 10 of Mardia, Kent, and Bibby (1979).
  # The approximate F statistic is exercise 3.7.6b.
  x<-as.matrix(x);y<-as.matrix(y)
  n<-dim(x)[1];q1<-dim(x)[2];q2<-dim(y)[2];q<-min(q1,q2)
  S11<-cov(x);S12<-cov(x,y);S21<-t(S12);S22<-cov(y)
  E1<-eigen(solve(S11)%*%S12%*%solve(S22)%*%S21);E2<-eigen(solve(S22)%*%S21%*%solve(S11)%*%S12)
  rsquared<-as.real(E1$values[1:q])
  LR<-NULL;pp<-NULL;qq<-NULL;tt<-NULL
  for (i in 1:q){
    LR<-c(LR,prod(1-rsquared[i:q]))
    pp<-c(pp,q1-i+1)
    qq<-c(qq,q2-i+1)
    tt<-c(tt,n-1-i+1)}
  m<-tt-0.5*(pp+qq+1);lambda<-(1/4)*(pp*qq-2);s<-sqrt((pp^2*qq^2-4)/(pp^2+qq^2-5))
  F<-((m*s-2*lambda)/(pp*qq))*((1-LR^(1/s))/LR^(1/s));df1<-pp*qq;df2<-(m*s-2*lambda);pval<-1-pf(F,df1,df2)
  outmat<-round(cbind(sqrt(rsquared),rsquared,LR,F,df1,df2,pval),dec)
  colnames(outmat)=list("R","RSquared","LR","ApproxF","NumDF","DenDF","pvalue")
  rownames(outmat)=as.character(1:q);xrels<-round(cor(x,x%*%E1$vectors)[,1:q],dec)
  colnames(xrels)<-apply(cbind(rep("U",q),as.character(1:q)),1,paste,collapse="")
  yrels<-round(cor(y,y%*%E2$vectors)[,1:q],dec)
  colnames(yrels)<-apply(cbind(rep("V",q),as.character(1:q)),1,paste,collapse="")
  list(Summary=outmat,a.Coefficients=E1$vectors,b.Coefficients=E2$vectors,
       XUCorrelations=xrels,YVCorrelations=yrels)
} 
## END FUNCTION