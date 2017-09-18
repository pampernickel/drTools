library(RColorBrewer)
library(limma)
library(reshape2)
library(ggplot2)

extractMax <- function(t){
  # function for retrieving x and y elements from a max
  # entry
  as.character(t) -> t
  as.numeric(strsplit(strsplit(t,";")[[1]][2], ",")[[1]]) -> y.fit
  as.numeric(strsplit(strsplit(t,";")[[1]][1], ",")[[1]]) -> x.fit
  list(x.fit, y.fit) -> fit.res
  names(fit.res) <- c("x", "y")
  return(fit.res)
}

rankResponses <- function(ares, assembled, drug.list.all, topK=10, poi=NULL, pos.only=F){
  # given a fit result (ares, which includes the full fit from which all other params)
  # may be derived, compare it with responses of other patients and 
  # rank it accordingly
  if (length(ares) > 1){
    stop("Multi-patient handling for ranking drug response not yet available. Please analyze patients individually.")
  }
  
  t(sapply(ares, function(x) x$logIC50)) -> res.df
  gsub("X6", "6", gsub("\\.", "-", colnames(res.df))) -> colnames(res.df)
  
  # get Emax
  t(sapply(ares[[1]]$max, function(x) extractMax(x)$y[100])) -> emax
  colnames(res.df) -> colnames(emax)
  rownames(res.df) -> rownames(emax)
  
  # get drug mapping -- for Emax, this is more to map the drug names
  .mapResponse(res.df, assembled, drug.list.all) -> assembled.sub
  .mapResponse(emax, assembled, drug.list.all) -> assembled.sub.max
  grep("POI", assembled.sub$id) -> ind
  
  # scale and center assembled sub
  as.matrix(assembled.sub) -> assembled.sub
  as.matrix(assembled.sub[,2:ncol(assembled.sub)]) -> mat
  apply(mat, 2, function(x) as.numeric(as.character(x))) -> mat
  scale(mat, center=T, scale=T) -> csdata
  csdata[ind,] -> poi
  csdata[-ind,] -> csdata
  
  if (length(ind) == 1){
    apply(csdata, 2, function(x) median(x, na.rm=T)) -> med
    med-poi -> diff
    rev(order(diff)) -> ord
    cbind(colnames(csdata)[ord], round(diff[ord], 2), round(med[ord],2)) -> res
    colnames(res) <- c("Drug", "Difference from median IC50", "Median IC50 (centered, scaled)")
    
    # then get other parameters, e.g. AUC and Emax, use them to rank the topK
    # automatically eliminate results where Emax > .15 at tested doses
    sapply(rownames(res), function(x) 
      assembled.sub.max[ind,grep(x, colnames(assembled.sub.max))]) -> me # max.effect
    cbind(res, round(me, 2)) -> res
    colnames(res)[ncol(res)] <- "Emax"
    as.numeric(res[,4])*100 -> res[,4]
    res[which(as.numeric(res[,4]) <= 15),] -> res
    if (nrow(res) > topK){
      res[1:topK,] -> res
    }
    as.data.frame(res) -> res
  } else {
    apply(csdata, 2, function(x) median(x, na.rm=T)) -> med
    apply(poi, 1, function(x){
      med - x -> diff
      rev(order(diff)) -> ord
      cbind(colnames(csdata)[ord], round(diff[ord], 2)) -> res
    }) -> res
  }
  
  if (pos.only == T){
    res[which(as.numeric(as.character(res$`Difference from median IC50`)) > 0),] -> res
  }
  
  res$Emax[which(as.numeric(as.character(res$Emax)) < 0)] <- 0
  rownames(res) <- c(1:nrow(res))
  return(res)
}

formatResponseMatrices <- function(fits.1, fits.2){
  lapply(fits.1, function(x) as.data.frame(x)) -> fits.df.1
  lapply(fits.2, function(x) as.data.frame(x)) -> fits.df.2
  
  # create data frame of containing x positions, y positions, and source
  fits.1.summary <- list()
  for (i in 1:length(fits.df.1)){
    apply(fits.df.1[[i]], 2, function(col) return(length(col$x))) -> x.lengths
    which(x.lengths == 100) -> ind
    fits.df.1.l <- matrix(NA, nrow=100, ncol=ncol(fits.df.1[[1]])+1)
    fits.df.1[[i]][,ind[1]]$x -> fits.df.1.l[,1]
    for (j in 1:ncol(fits.df.1[[i]])){
      if (length(which(is.na(fits.df.1[[i]][,j]$y))) != 100){
        fits.df.1[[i]][,j]$y -> fits.df.1.l[,j+1]
      }
    }
    fits.df.1.l -> fits.1.summary[[i]]
  }
  
  fits.2.summary <- list()
  for (i in 1:length(fits.df.2)){
    apply(fits.df.2[[i]], 2, function(col) return(length(col$x))) -> x.lengths
    which(x.lengths == 100) -> ind
    fits.df.2.l <- matrix(NA, nrow=100, ncol=ncol(fits.df.2[[1]])+1)
    fits.df.2[[i]][,ind[1]]$x -> fits.df.2.l[,1]
    for (j in 1:ncol(fits.df.2[[i]])){
      if (length(which(is.na(fits.df.2[[i]][,j]$y))) != 100){
        fits.df.2[[i]][,j]$y -> fits.df.2.l[,j+1]
      }
    }
    fits.df.2.l -> fits.2.summary[[i]]
  }
  
  names(fits.1.summary) <- names(fits.1)
  names(fits.2.summary) <- names(fits.2)
  list(fits.1.summary, fits.2.summary) -> res.fits.df
  names(res.fits.df) <- c("class1", "class2")
  return(res.fits.df)
}

findPrognosticator <- function(res.fits, significance){
  # easier to do with a for loop than nested applies!
  prognosticator.array <- c()
  dose.array <- c()
  t.stat <- c()
  p.val <- c()
  for (i in 1:length(res.fits$class1)){
    # --- df.1 and df.2 are responses to drugs; take
    # --- point-by-point response (i.e. response at individual dose to
    # --- retrieve dose/drug combinations with prognostic value)
    res.fits$class1[[i]] -> df.1
    res.fits$class2[[i]] -> df.2
    
    for (j in 1:nrow(df.1)){
      as.numeric(df.1[j,2:ncol(df.1)]) -> set.1
      as.numeric(df.2[j,2:ncol(df.2)]) -> set.2
      run.t.test(set.1, set.2) -> test.res
      if (test.res$p <= significance){
        prognosticator.array <- c(prognosticator.array, names(res.fits$class1)[i])
        dose.array <- c(dose.array, df.1[j,1])
        p.val <- c(p.val, test.res$p)
        t.stat <- c(t.stat, test.res$estimate)
      }
    }
  }
  
  cbind(prognosticator.array, dose.array, t.stat, p.val) -> res
  return(res)
}

run.t.test <- function(set.1, set.2){
  if (sd(set.1, na.rm=T) == 0 && sd(set.2, na.rm=T) == 0){
    res.fin <- cbind(1,1)
  } else {
    t.test(set.1, set.2, na.action=na.omit) -> res
    cbind(res$p.value, res$statistic) -> res.fin
  }
  as.data.frame(res.fin) -> res.fin
  colnames(res.fin) <- c("p", "t")
  return(res.fin)
}

getMaxSep <- function(exp.res.mat, factor, class.1, class.2){
  # performs a scan of values in exp.res.mat, including
  # a full Emax scan on the x axis
  # factor: AUC, Emax, logIC50, etc
  # class.1: patients linked to a contrast class.1
  # class.2: patients linked to a contrast class.2
  # in case names are directory names
  
  # now check for matches
  #if (length(which(rownames(exp.res.mat[[1]]) %in% class.1)) == 0 ||
  #      length(which(rownames(exp.res.mat[[1]]) %in% class.2)) == 0){
  #  stop("At least one class contains no matches in your fit results matrix.")
  #}
  
  if (length(grep("fit.res.", rownames(exp.res.mat[[1]]))) > 0){
    which(gsub("fit.res.", "", rownames(exp.res.mat[[1]])) %in% class.1) -> ind.1
    which(gsub("fit.res.", "", rownames(exp.res.mat[[1]])) %in% class.2) -> ind.2
  }
  
  which(names(exp.res.mat) %in% factor) -> list.ind
  
  df <- exp.res.mat[[list.ind]]
  prognosticator.df <- NA
  if (factor %in% "max"){
    # for each drug, retain both column ind (drug ind) and x,y point set
    # that results in a significant difference between classes 1 and 2
    apply(df[ind.1,], 2, function(t) sapply(t, function(x) extractMax(x))) -> fits.1
    apply(df[ind.2,], 2, function(t) sapply(t, function(x) extractMax(x))) -> fits.2
    formatResponseMatrices(fits.1, fits.2) -> res.fits
    findPrognosticator(res.fits, 0.05) -> prognosticator.df
  }
  return(prognosticator.df)
}

euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

imputeNans <- function(exp.res.mat){
  # --- Nans occur in the fits with non-reactive samples. Perform imputations for NaNs
  for (i in 1:length(exp.res.mat)){
    exp.res.mat[[i]] -> df
    for (j in 1:ncol(df)){
      if (length(which(is.nan(df[,j]))) > 0){
        df[which(is.nan(df[,j])),j] <- 4
      }
    }
    df -> exp.res.mat[[i]]
  }
  
  return(exp.res.mat)
}

runPairedTest <- function(matrix, vars, ind.1, ind.2, opt=c("less", "greater")){
  # wrapper for paired t-test
  # input:
  # matrix: matrix of drug responses
  # vars: vector of variables (e.g. IC50, AUC, etc) to be tested
  # ind.1: indices of responses for class 1
  # ind.2: indices of responses for class 2, paired with class.1
  # opt: direction of test
  
  res <- list()
  matrix.orig <- matrix
  for (j in 1:length(vars)){
    grep(vars[j], colnames(matrix.orig)) -> ind
    matrix.orig[,ind] -> matrix
    
    # --- paired test
    # --- paired signif
    paired.resp.res <- c()
    for (i in 1:ncol(matrix)){
      matrix[,i] -> resp
      as.numeric(resp) -> resp
      if (length(which(resp %in% NA)) == 0 &&
            length(which(resp %in% NaN)) == 0 &&
            sd(resp) > 0){
        t.test(resp[ind.1], resp[ind.2], paired=T, alt=opt) -> t
        if (t$p.val < 0.05 && !is.na(t$p.val)){
          paired.resp.res <- c(paired.resp.res, i)
        }
      }
    }
    
    if (length(paired.resp.res) > 0){
      res[[j]] <- paired.resp.res
    }
  }
  
  names(res) <- vars
  return(res)
}

mixModel <- function(exp.res.mat, parameter){
  require(mclust)
  require(mixtools)
  exp.res.mat[[which(names(exp.res.mat) %in% parameter)]] -> mat
  t(mat) -> mat.t
  if (length(which(duplicated(rownames(mat.t))))){
    mat.t[-which(duplicated(rownames(mat.t))),] -> mat.t  
  }
  
  topVar(mat.t) -> mat.top
  t(mat.top) -> mat.top
  #   
  #   # --- generate those with the highest standard deviation
  all.models <- list()
  for (i in 1:ncol(mat.top)){
    # --- remove NaNs, na
    mat.top[,i] -> vect
    if(length(which(is.nan(vect))) > 0){
      vect[-which(is.nan(vect))] -> vect
    }
    as.numeric(vect) -> vect
    mixmdl = normalmixEM(na.omit(vect))
    all.models[[i]] <- mixmdl
  }
  
  return(mat.top)
}

topVar <- function(matrix){
  rownames(matrix)-> probes
  probe.sd <- rep(NA, length(probes))
  
  matrix.df <- data.frame(matrix) #need data frame for apply?
  
  for(i in 1:length(probes)){
    # check, remove NAs, NAns in matrix.df[i,]
    as.numeric(matrix.df[i,]) -> vect
    if (length(which(is.nan(vect))) > 0){
      vect[-which(is.nan(vect))] -> vect
    }
    sd(vect) -> probe.sd[i]
  }
  
  unsorted.probe <- probe.sd #to retain indices
  sorted.probe <- probe.sd[order(-probe.sd)] #Rsorting is ascending by default; minus sign
  top.probes <- sorted.probe[which(sorted.probe >= mean(probe.sd))]
  top.probe.ind <- which(unsorted.probe %in% top.probes)
  matrix.top <- matrix[top.probe.ind, ]
  return(matrix.top)
}

integratedScoreCal <- function(exp.res.mat, vars){
  # nrow= patients, ncol = #drugs
  exp.res.mat[which(names(exp.res.mat) %in% vars)] -> exp.res.mat
  df <- matrix(0,nrow=nrow(exp.res.mat[[1]]), ncol=ncol(exp.res.mat[[1]]))
  for (i in 1:length(exp.res.mat)){
    # take average of all cells across the components
    # of exp.res.mat
    for (j in 1:nrow(df)){
      for (k in 1:ncol(df)){
        df[j,k] <- df[j,k]+exp.res.mat[[i]][j,k]
      } 
    }
  }
  
  df <- df/length(exp.res.mat)
  rownames(df) <- rownames(exp.res.mat[[1]])
  colnames(df) <- colnames(exp.res.mat[[1]])
  return(df)
}

getCorrelations <- function(exp.res.ave, drug.list){
  # --- function to get correlations between paris of variables (e.g. correl bet. IC50 and AUC, etc.)
  all.dfs <- matrix(nrow=0, ncol=length(exp.res.ave[[1]]))
  all.df.labs <- matrix(nrow=0, ncol=2)
  for (i in 1:length(exp.res.ave)){
    df <- as.data.frame(exp.res.ave[[i]])
    cbind(as.character(drug.list), rep(names(exp.res.ave)[i], length(drug.list))) -> df.labs
    rbind(all.dfs, df) -> all.dfs
    rbind(all.df.labs, df.labs) -> all.df.labs
  }
  
  colnames(all.dfs) <- colnames(as.data.frame(exp.res.ave[[1]]))
  
  # --- remove error and ci params
  if (length(which(colnames(all.dfs) %in% c("error.50", "ci.all"))) > 0){
    all.dfs[,-which(colnames(all.dfs) %in% c("error.50", "ci.all"))] -> all.dfs
  }
  
  ind <- c()
  for (i in 1:ncol(all.dfs)){
    c(ind, which(all.dfs[,i] %in% c("NaN", NA))) -> ind
  }
  
  unique(ind) -> ind
  if (length(ind) > 0){
    all.dfs[-ind, ] -> all.dfs
    all.df.labs[-ind,] -> all.df.labs
  }
  
  combn(ncol(all.dfs), 2) -> combs
  # --- create final data frame for visualization
  mat <- matrix(nrow=0, ncol=5)
  for (i in 1:ncol(combs)){
    combs[,i] -> curr.comb
    all.dfs[,curr.comb] -> t
    cbind(t, rep(paste(colnames(all.dfs)[curr.comb[1]],"vs",colnames(all.dfs)[curr.comb[2]],sep=" "),nrow(t)), all.df.labs) -> m
    colnames(m) <- c("", "", "", "", "")
    rbind(mat, m) -> mat
  }
  
  colnames(mat) <- c("v1", "v2", "comparison", "drug", "sample")
  as.data.frame(mat) -> df
  ggplot(df,aes(x=v1,y=v2))+
    geom_point()+
    scale_x_continuous(limits = c(-4, 4))+
    scale_y_continuous(limits = c(-4, 4))+
    #geom_smooth(method = "lm", se=FALSE, color="red")+
    facet_wrap(~comparison) -> p
  ggsave(filename="Fitted.var.correl.png", plot=p)
  return(mat)
}

plotRGL <- function(exp.res.mat, labels){
  # --- create a matrix combining all features from exp.res.mat,
  #as.data.frame(exp.res.mat) -> df
  as.data.frame(exp.res.mat.merged[c(1,6)]) -> df
  # --- remove columns with variance of 0
  df[,as.numeric(which(apply(df, 2, var, na.rm=TRUE) != 0))] -> df
  
  # --- remove columns with NaNs
  removeNans(df) -> df
  # --- data standardization: df has mean 0 and variance 1:
  # --- scaling should take care of this!
  #   for (i in 1:ncol(df)){
  #     df[,i]<- (df[,i] - mean(df[,i]))/sd(df[,i])
  #   }
  for (i in 1:ncol(df)){
    c(which(is.na(df[,i])), which(is.infinite(df[,i]))) -> ind
    if (length(ind) > 0){
      df[-ind,] -> df
    }
  }
  
  # --- pca on patients
  pca = prcomp(t(df),scale=T)
  x= pca$rotation[,1]
  y= pca$rotation[,2]
  z= pca$rotation[,3]
  
  # --- create the different possible colorings based on the labels
  # --- color vector, currently temporary!
  #colors <- list(c("red", "grey50"), c("blue", "grey50"), c()
  for (i in 1:ncol(labels)){
    lab <- rep(NA, length(labels[,i]))
    unique(labels[,i]) -> all.lab
    for (j in 1:length(all.lab)){
      lab[which(labels[,i] %in% all.lab[j])] <- all.lab[j]
    }
    
    if (length(all.lab) >= 3){
      colors <- brewer.pal(length(all.lab),"BrBG")
    } else {
      colors <- c("red", "grey60")
    }
    
    png(paste("PCA",colnames(labels)[i],"png", sep="."), w=800, h=400)
    par(mfrow=c(1,3))
    plot(x,y, col=cbind(c(colors) [as.numeric(as.factor(lab))]), pch=19, cex=1.25)
    text(x,y, labels=lab, col=cbind(c(colors) [as.numeric(as.factor(lab))]))
    plot(y,z, col=cbind(c(colors) [as.numeric(as.factor(lab))]), pch=19, cex=1.25)
    text(y,z, labels=lab, col=cbind(c(colors) [as.numeric(as.factor(lab))]))
    plot(x,z, col=cbind(c(colors) [as.numeric(as.factor(lab))]), pch=19, cex=1.25)
    text(x,z, labels=lab, col=cbind(c(colors) [as.numeric(as.factor(lab))]))
    plot3d(x,y,z,col=cbind(c(colors) [as.numeric(as.factor(all.lab))]), main="PCA by entity")
    #text3d(x,y,z,text=lab,
    #       col=cbind(c(colors) [as.numeric(as.factor(pathology))]),cex=1.0)
    dev.off()
  }
}

runLimma <- function(exp.res.mat, labels, column, p.val, FC){
  if (!column %in% colnames(labels)){
    stop("Selected column must match an entry in 'labels' parameter.")
  } else {
    as.data.frame(exp.res.mat) -> df
    labels[,which(colnames(labels) %in% column)] -> local.labs
    which(local.labs %in% c("", " ", NA)) -> ind
    
    if (length(ind) > 0){
      df[-ind,] -> df
      local.labs[-ind] -> local.labs
    }
    
    as.numeric(as.factor(local.labs)) -> num.labs
    local.labs[order(num.labs)] -> local.labs
    df[order(num.labs),] -> df
    num.labs[order(num.labs)] -> num.labs
    removeNans(df) -> df
    t(df) -> df
    
    # --- perform limma propper
    GRP1 = c(1:length(which(num.labs %in% 1)))
    GRP2 = c(length(GRP1)+1:(length(which(num.labs %in% 2))))
    
    ff <- rep(NA,ncol(df))
    ff[GRP1] <- 0
    ff[GRP2] <- 1
    
    design <- model.matrix(~ -1+factor(ff) )
    colnames(design) <-c("GRP1", "GRP2") 
    fit <- lmFit(df, design) 
    contrast.matrix <- makeContrasts( 
      GRP.1.2 = GRP2 - GRP1,
      levels=design
    )
    
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit2) 
    # --- coef = number of contrasts, in this case, 6
    toptab.GRP = topTable(fit2, coef=c(1), number = dim(df)[1], adjust="BH") 
    toptabGRP = toptab.GRP[which(toptab.GRP [,"adj.P.Val"] <= p.val & abs(toptab.GRP [,"logFC"]) >= FC),]
    responsive = toptabGRP[which(toptabGRP[,"logFC"]>0),]
    resistant = toptabGRP[which(toptabGRP[,"logFC"]<0),]
    list(responsive, resistant, toptabGRP) -> drug.resp
    names(drug.resp) <- c("responsive", "resistant", "toptab")
    return(drug.resp)
  }
}

removeNans <- function(df){
  
  # --- function to remove columns of a dataframe
  # --- with NaNs
  cols <- c()
  for (i in 1:ncol(df)){
    if (length(which(is.nan(df[,i]))) > 0){
      c(cols, i) -> cols
    }
  }
  
  if (length(cols) > 0){
    df[,-cols] -> df
  }
  
  return(df)
}

trimList <- function(exp.res.ave, labels, column, exact.labs){
  # --- function for subsetting exp.res.ave
  # --- exact.labs is a subset of the labels column
  # --- e.g. in the case of the translocation column
  # --- it could correspond to 1;19 or 17;19 or both
  if (!column %in% colnames(labels)){
    stop("Selected column must match an entry in 'labels' parameter.")
  }
  
  labels[,which(colnames(labels) %in% column)] -> local.labs
  #which(!local.labs %in% c("", " ", NA)) -> ind
  which(local.labs %in% exact.labs) -> ind
  exp.res.ave.s <- list()
  if (length(ind) > 0){
    for (i in 1:length(exp.res.ave)){
      if (i %in% ind){
        exp.res.ave[[i]] -> exp.res.ave.s[[(length(exp.res.ave.s)+1)]]
      }
    }
    local.labs[-ind] -> local.labs
    names(exp.res.ave.s) <- names(exp.res.ave)[ind]
  }
  
  return(exp.res.ave.s)
}

formatGSEA <- function(sub.list, param, drug.list){
  as.data.frame(sub.list) -> df
  grep(param, colnames(df)) -> ind
  df[,ind] -> df
  rownames(df) <- drug.list
}

comparePlates <- function(p1, p2){
  # compare a plate/list of plates 2 with a plate 1
  df.sum <- NA
  if (is.list(p2)){
    lapply(p2, function(x){
      if (nrow(p1) != nrow(x)){
        t(p1) -> p1
        if (nrow(p1) != nrow(x)){
          stop("Plates do not have the same dimension.")
        }
      }
      
      melt(p1) -> y
      melt(x) -> x
      cbind(x$value, y$value) -> cor.mat
      return(cor.mat)
    }) -> cor.mats
    names(cor.mats) <- names(p2)
    df.sum <- do.call("rbind", cor.mats)
    unlist(lapply(1:length(cor.mats), function(x) rep(names(cor.mats)[x], nrow(cor.mats[[x]])))) -> labs
    cbind(df.sum, labs) -> df.sum
    as.data.frame(df.sum) -> df.sum
  } else {
    
  }
  return(df.sum) 
}