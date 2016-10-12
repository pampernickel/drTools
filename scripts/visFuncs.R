library(ggplot2)
library(gridExtra)
library(lattice)
library(reshape2)
library(RColorBrewer)
library(NMF)

processMaxCurves <- function(t){
  strsplit(strsplit(t,";")[[1]][2], ",")[[1]] -> y
  strsplit(strsplit(t,";")[[1]][1], ",")[[1]] -> x
  cbind(y, x) -> res
  as.data.frame(res) -> res
  apply(res, 2, function(x) return(as.numeric(x))) -> res
  return(res)
}


visualizePrognosticator <- function(exp.res.mat, class.1, class.2, prognosticator.df, label1, label2){
  unique(prognosticator.df[,1]) -> doi # drugs of interest
  as.numeric(as.character(prognosticator.df[,2])) -> prognosticator.df[,2]
  # find drugs in max array
  exp.res.mat$max -> df
  if (length(grep("fit.res.", rownames(df))) > 0){
    gsub("fit.res.", "", rownames(df)) -> rownames(df)
  }
  df[which(rownames(df) %in% class.1),which(colnames(df) %in% doi)] -> df.1
  df[which(rownames(df) %in% class.2),which(colnames(df) %in% doi)] -> df.2
  
  # then for each prognosticator, create the plots for class.1 and class.1
  for (i in 1:ncol(df.1)){
    as.character(df.1[,i]) -> t
    as.character(df.2[,i]) -> t.1
    lapply(t, function(x) processMaxCurves(x)) -> t.all
    lapply(t.1, function(x) processMaxCurves(x)) -> t.1.all
    
    set.1 <- matrix(nrow=0, ncol=3)
    set.2 <- matrix(nrow=0, ncol=3)
    
    for (j in 1:length(t.all)){
      if (!is.null(nrow(t.all[[j]]))){
        set.1 <- rbind(set.1, cbind(t.all[[j]], rep(label1, nrow(t.all[[j]]))))
      }
    }
    
    for (j in 1:length(t.1.all)){
      if (!is.null(nrow(t.1.all[[j]]))){
        set.2 <- rbind(set.2, cbind(t.1.all[[j]], rep(label2, nrow(t.1.all[[j]]))))
      }
    }
    
    rbind(set.1, set.2) -> set.df
    as.data.frame(set.df) -> set.df
    colnames(set.df) <- c("y", "x", "class")
    as.numeric(as.character(set.df$y)) -> set.df$y
    as.numeric(as.character(set.df$x)) -> set.df$x
    
    # --- subset to match prognosticator for drug for boxplot
    #generate side-by-side plots of response curve and boxplot for prognosticator
    # first, generate plots
    man.scale <- scale_colour_manual(name="Class", values=c(label1="red", label2="blue"))
    man.fill <- scale_fill_manual(name="Class", values=c(label1="red", label2="blue"))
    p <- ggplot(data=set.df, aes(x=x, y=y, color=class))+
      stat_summary(fun.data="mean_sdl", mult=1, geom="smooth",
                   aes(fill=factor(class)))+scale_x_log10()+
      #annotate("rect", xmin=min(as.numeric(prognosticator.df[which(prognosticator.df[,1] %in% doi[i]),2])),
      #         xmax=max(as.numeric(prognosticator.df[which(prognosticator.df[,1] %in% doi[i]),2])),
      #         ymin=min(set.df[,1]),ymax=max(set.df[,1]), alpha=0.15)+
      ggtitle(toupper(doi[i]))+man.fill+man.scale+
      theme(panel.background = element_rect(fill = "grey99"))
    
    # --- then generate a cut through the curve that shows the maximum separation in terms of the t-statistic
    # --- and the p-value
    which(prognosticator.df[,3] %in% min(as.numeric(prognosticator.df[which(prognosticator.df[,1] %in% doi[i]),3]))) -> pt
    prognosticator.df[pt[1],2] -> x.slice.pt
    set.df[which(set.df$x == as.numeric(x.slice.pt)),] -> set.df.sub
    
    p1 <- ggplot(data=set.df.sub, aes(factor(class), y=y))+
      geom_boxplot(outlier.size=0)+geom_point(size=4.5, position=position_jitter(width=0.15))+
      ggtitle(paste(doi[i], "@", round(as.numeric(x.slice.pt),1), "nm", sep=" "))
    
    p2 <- arrangeGrob(p, p1, ncol=2, nrow=1)
    ggsave(paste("main.",doi[i],".pdf",sep=""), p, width=4, height=4)
    ggsave(paste(doi[i],".pdf",sep=""), p2)
  }
}


visualizeResponse <- function(exp.res.mat, class.1, class.2, doi, label1, label2){
  # find drugs in max array
  exp.res.mat$max -> df
  if (length(grep("fit.res.", rownames(df))) > 0){
    gsub("fit.res.", "", rownames(df)) -> rownames(df)
  }
  
  df[which(rownames(df) %in% class.1),which(colnames(df) %in% doi)] -> df.1
  df[which(rownames(df) %in% class.2),which(colnames(df) %in% doi)] -> df.2
  
  # then for each prognosticator, create the plots for class.1 and class.1
  as.character(df.1) -> t
  as.character(df.2) -> t.1
  lapply(t, function(x) processMaxCurves(x)) -> t.all
  lapply(t.1, function(x) processMaxCurves(x)) -> t.1.all
    
  set.1 <- matrix(nrow=0, ncol=3)
  set.2 <- matrix(nrow=0, ncol=3)
    
  for (j in 1:length(t.all)){
    if (!is.null(nrow(t.all[[j]]))){
      set.1 <- rbind(set.1, cbind(t.all[[j]], rep(label1, nrow(t.all[[j]]))))
    }
  }
    
  for (j in 1:length(t.1.all)){
    if (!is.null(nrow(t.1.all[[j]]))){
      set.2 <- rbind(set.2, cbind(t.1.all[[j]], rep(label2, nrow(t.1.all[[j]]))))
    }
  }
    
  rbind(set.1, set.2) -> set.df
  as.data.frame(set.df) -> set.df
  colnames(set.df) <- c("y", "x", "class")
  as.numeric(as.character(set.df$y)) -> set.df$y
  as.numeric(as.character(set.df$x)) -> set.df$x
  
  #print(head(set.df))
  print(label1)
  # --- subset to match prognosticator for drug for boxplot
  #generate side-by-side plots of response curve and boxplot for prognosticator
    # first, generate plots
    man.scale <- scale_colour_manual(name="Class", values=factor(c(label1="red", label2="blue")))
    man.fill <- scale_fill_manual(name="Class", values=factor(c(label1="red", label2="blue")))
    p <- ggplot(data=set.df, aes(x=x, y=y, color=class))+
      stat_summary(fun.data="mean_sdl", mult=1, geom="smooth",
                   aes(fill=factor(class)))+scale_x_log10()+
      ggtitle(toupper(doi))+
    theme(panel.background = element_rect(fill = "grey99"))
    ggsave(paste("main.",doi,".pdf",sep=""), p, width=4, height=4)
}

plotTestRes <- function(matrix, vars, labels, res, drug.list){
  # Function for visualizing paired t-test results;
  # works with runPairedTest()
  # input: 
  # matrix: matrix of response parameters
  # vars: vector of variables (e.g. IC50, AUC, etc)
  # labels: vector of classes (e.g. diagnosis, relapse)
  # res: indices calculated from runPairedTest indicating significant response
  #     per entry in vars
  # drug list: vector of drug names
  
  require(ggplot2)
  
  # --- create plots with following features
  # --- score, label, parameter
  all.res <- matrix(nrow=0, ncol=5)
  colnames(all.res) <- c("response", "drug.name", "labels", "variable", "patient")
  for (j in 1:length(vars)){
    grep(vars[j], colnames(matrix)) -> ind
    matrix[,ind] -> matrix.sub
    
    if (length(res[[j]]) > 0){
      res[[j]] -> resp.ind # --- indices to extract response from mat
      for (i in 1:length(resp.ind)){
        cbind(matrix.sub[,resp.ind[i]],
              as.character(drug.list[resp.ind[i]]),
              labels, 
              rep(vars[j], length(matrix.sub[,resp.ind[i]])),
              rownames(matrix.sub)) -> res.sub
        colnames(res.sub) <- c("response", "drug.name", "labels", "variable", "patient")
        rbind(all.res, res.sub) -> all.res
      }
    }
  }
  
  gsub("a", "", all.res[,5]) -> all.res[,5]
  gsub("b", "", all.res[,5]) -> all.res[,5]
  gsub("c", "", all.res[,5]) -> all.res[,5]
  as.data.frame(all.res) -> all.res
  
  for (i in 1:length(vars)){
    if (length(res[[i]]) > 0){
      all.res.sub <- all.res[which(all.res$variable %in% vars[i]),]
      as.data.frame(all.res.sub) -> all.res.sub
      round(as.numeric(as.character(all.res.sub$response)),2) -> all.res.sub$response
      p <- ggplot(data = all.res.sub, aes(x=labels, y=response))
      p.1 <- p+geom_boxplot(outlier.shape = NA)+facet_wrap(~drug.name)+
        geom_jitter(position=position_jitter(w=0.1, h=0.1), 
                    aes(color=patient),size=5, alpha=0.6)
      ggsave(filename=paste(vars[i],"paired.t.res", "png", sep="."), plot=p.1)
    }
  }
}

plotPerDrug <- function(drug.name, experiment, type=c("raw", "normalized")){
  # --- create directory to put patient plots in
  if (is.list(experiment)){
    # --- find name of drug in the first patient. This is with 
    # --- the assumption that the drug is always placed in the same
    # --- position
    doses <- c(0,1,2,3,4) # log dose
    # final data frame for ggplot
    fin.data <- matrix(NA, 0, 3)
    temp <- as.data.frame(matrix(numeric(0), 5,0))
    if (type=="raw"){
      for (i in 1:5){
        temp <- cbind(temp, experiment[[1]][i][[1]])
      }
    } else {
      temp <- experiment[[1]]
      temp <- temp[,c(2:ncol(temp))]
    }
    
    grep(drug.name, colnames(temp)) -> ind
    # --- create data frame containing all responses for this drug
    fin.data <- as.data.frame(matrix(NA,0,3))
    for (i in 1:length(experiment)){
      temp <- as.data.frame(matrix(numeric(0), 5,0))
      resp <- c()
      dose <- c()
      patient <- c()
      if (type=="raw"){
        for (j in 1:5){
          temp <- cbind(temp, experiment[[i]][j][[1]])
        }
      } else {
        temp <- experiment[[i]]
        temp <- temp[,c(2:ncol(temp))]
      }
      
      temp[,ind] -> temp.0
      # --- take average
      c(resp, rowMeans(temp.0)) -> resp
      c(dose, doses, doses) -> dose
      c(patient, rep(names(experiment)[i], length(doses))) -> patient
      cbind(patient, dose, resp) -> temp.1
      rbind(temp.1, fin.data) -> fin.data
    }
    
    # --- append translocation information to fin.data
    c("1a", "2a", "3a", "4a", "5a") -> t.1
    rep(NA, nrow(fin.data)) -> translocation.type
    for (i in 1:length(t.1)){
      grep(t.1[i], fin.data[,1]) -> ind
      translocation.type[ind] <- "1;19"
    }
    
    translocation.type[which(is.na(translocation.type))] <- "17;19"
    cbind(fin.data, translocation.type) -> fin.data
    as.data.frame(fin.data) -> fin.data
    p <- ggplot(data = fin.data, aes(x = dose, y = as.numeric(as.character(resp)), group=patient, color=translocation.type))
    p.1 <- p+geom_point()+geom_line()+ggtitle(drug.name)+xlab("Log Dose (nm)")+ylab("Live cells")
    if (type=="raw"){
      ggsave(filename=paste(drug.name, "response.summary.raw.png", sep="."), plot=p.1)
    } else if (type=="normalized"){
      ggsave(filename=paste(drug.name, "response.summary.norm.png", sep="."), plot=p.1)
    }
  } else {
    stop("Invalid experiment type.")
  }
}

plotResponses <- function(experiments, exp.res.ave, param, func, cat){
  # --- creates horizontal bar plots similar to those 
  # --- in Garnett et al.
  # --- func param: mean, sd, median, etc
  
  # --- capture local environment to ensure that no scoping problems occur
  .e <- environment()
  if (!param %in% names(exp.res.ave[[1]])){
    stop("Argument does not match any parameter in the experimental results")  
  }
  
  # --- choose exp.res.ave rows based on labels
  # --- get drug.list based on names of experiments, just to be sure
  # --- that any internally added controls, e.g. for 8-pt data, DMSO labels
  # --- are added internally will not cause problems
  getDrugNames.int(experiments) -> drug.list # --- internal function
  resp <- matrix(nrow=length(drug.list),ncol=length(exp.res.ave))
  rownames(resp) <- as.character(drug.list)
  grep(param,names(exp.res.ave[[1]])) -> ind
  for (i in 1:length(exp.res.ave)){
    exp.res.ave[[i]] -> df
    df[[ind]] -> resp[,i]
  }
  
  df <- matrix(nrow=0,ncol=3)
  for (i in 1:nrow(resp)){
    rep(as.character(rownames(resp)[i]), ncol(resp)) -> drug
    rep(as.character(cat[i]), ncol(resp)) -> cat.all
    cbind(drug, resp[i,], cat.all) -> t
    rbind(df, t) -> df
  }
  
  as.data.frame(df) -> df
  
  which(df[,2] %in% c("NaN", NA)) -> ind
  if (length(ind) > 0){
    df[-ind,] -> df
  }
  
  
  # --- plots:
  # 1: count of cases within a certain range
  #range.lims <- c(4, 3.5, 3, 2.5 ,2, 1.5, 1, 0.5, 0)
  #categories <- rep(NA, nrow(df))
  #for (i in 1:length(range.lims)){
  #  paste("cat",i,sep="") -> categories[which(as.numeric(as.character(df[,2])) <= range.lims[i])]
  #}
  
  #cbind(df, categories) -> df
  
  
  colnames(df) <- c("drug", "param", "cat")
  as.numeric(as.character(df[,2])) -> df[,2]
  
  ggplot(df, aes(x=reorder(drug, param, FUN=func),y=param,fill=cat), environment = .e)+
    geom_boxplot()+
    coord_flip()+
    ylab(param)+
    xlab("Drugs") -> p
  
  return(p)
}

getDrugNames.int <-function(experiments){
  drug.list <- c()
  curr.exp <- experiments$normalized[[1]]
  if (is.data.frame(curr.exp)){
    drug.list <- colnames(curr.exp)
  } else if (is.list(curr.exp)){
    for (j in 1:length(curr.exp)){
      drug.list <- c(drug.list, names(curr.exp[[j]]))
    }
  }
  
  if (length(which(drug.list %in% "Concentrations") > 0)){
    drug.list[-which(drug.list %in% "Concentrations")] -> drug.list
  }
  
  temp <- c(1:length(drug.list))
  which(temp %% 2 == 1) -> ind
  drug.list[ind] -> drug.list
  for (i in 1:length(drug.list)){
    substr(drug.list[i], 1, nchar(drug.list[i])-2) -> drug.list[i]
  }
  
  return(drug.list)
}

visSD <- function(df, drug.name){
  # --- create histograms of different parameters per drug
  # --- df is the same type as that used for the pca
  grep(drug.name, colnames(df)) -> cols
  plot.df <- matrix(nrow=0, ncol=2)
  for (i in 1:length(cols)){
    curr.col <- cols[i]
    curr.col.name <- gsub(paste(".",drug.name, sep=""), "", colnames(df)[curr.col])
    cbind(df[,curr.col], rep(curr.col.name, length(df[,curr.col]))) -> t
    rbind(plot.df, t) -> plot.df
  }
  
  as.data.frame(plot.df) -> plot.df
  colnames(plot.df) <- c("val", "parameter")
  ggplot(plot.df, aes(x=round(as.numeric(as.character(val)),2), fill=parameter))+
    geom_histogram(binwidth=0.05, alpha=0.5, position="identity")+
    facet_wrap(~parameter)+
    ggtitle(paste("Distribution of parameters for", drug.name, sep=" "))+
    theme(text = element_text(size=10))+
    guides(fill=guide_legend(title="Parameter"))+
    xlab("Parameter values") -> p
  return(p)
}

visMatrices <- function(exp.res.mat, labels, width, height, annot.width){
  require(nclust)
  color.m <- matrix(data = NA, nrow = nrow(labels), ncol = ncol(labels))
  for (i in 1:ncol(labels)){
    as.numeric(as.factor(labels[,i])) -> temp
    if (length(unique(temp[!is.na(temp)])) == 2){
      # --- binary palette
      color.m[which(temp %in% 1),i] <- "white"
      color.m[which(temp %in% 2),i] <- "black"
    } else {
      # --- use topo colors
      rainbow(length(unique(temp))) -> color.sub
      for (j in 1:length(unique(temp))){
        color.m[which(temp %in% j),i] <- color.sub[j]
      }
      # --- then overwrite everything with blanks
      color.m[c(which(labels[,i] %in% ""),which(is.na(labels[,i]))),i] <- "white"
    }
  }
  colnames(color.m) <- colnames(labels)
  
  for (i in 1:length(exp.res.mat)){
    exp.res.mat[[i]] -> curr.mat
    curr.mat*-1 -> curr.mat
    curr.mat[,-grep("DMSO", colnames(curr.mat))] -> curr.mat
    pat.names <- rownames(curr.mat)
    # --- convert labels to colors
    #list(as.matrix(data.frame("Pathology" = labels[,1], x)))
    nclust.clust <- nclust2(curr.mat, cpwcorr=TRUE)
    # --- save names of patient order (in a cluster) per matrix
    # --- currently needs to be implemented this way since patient
    # --- names are currently in the form of full directory names
    #listPatients(nclust.clust, names(exp.res.mat)[i])
    png(paste(names(exp.res.mat)[i],".png", sep=""), w=width, h=height)
    coldmap(curr.mat, clust = nclust.clust, rannot=list(color.m), clab=colnames(exp.res.mat[[i]]),
            rdw = 5, cdw = 5, rlab=list(pat.names), rannot.width=annot.width, clab.height=7, rlab.width=5, cex.clab=1.5)
    dev.off()
  }
}

createLabels <- function(exp.res.mat, patient.labels){
  toupper(rownames(exp.res.mat[[1]])) -> rows
  toupper(patient.labels[,1]) -> patient.labels[,1]
  labels <- matrix(data = NA, nrow = length(rows), ncol = (ncol(patient.labels)-1))
  
  # --- find correspondence bet. patient labels and rows
  for (i in 1:nrow(patient.labels)){
    grep(patient.labels[i,1],rows, fixed=TRUE) -> ind
    if (length(ind)==1){
      for (j in 2:ncol(patient.labels)){
        labels[ind,(j-1)] <- patient.labels[i,j]
      }
    }
  }
  
  colnames(labels) <- colnames(patient.labels)[2:ncol(patient.labels)]  
  return(labels)
}

# function to extract response point coordinates
# for annotation in a 3d scatterplot
getPatientCoords <- function(df,name,drug,auc,ic50,max){
  # one of original data frames from which data was extracted
  # can also be used for extracting info for a group patients
  # on multiple drugs; only do a grep if an exact match
  # can't be found
  row.ind <- 0
  if (length(grep("fit.res.", rownames(ic50))) > 0){
    gsub("fit.res.", "", rownames(ic50)) -> rownames(ic50)
  }
  
  if (name %in% rownames(ic50)){
    which(rownames(ic50) %in% name) -> row.ind
  } else {
    grep(name,rownames(ic50)) -> row.ind
  }
  
  if (length(row.ind) > 1){
    warning("multiple row indices found.")
    row.ind[1] -> row.ind
  }
  
  which(colnames(df) %in% drug) -> col.ind
  as.numeric(auc[row.ind, col.ind]) -> x
  as.numeric(ic50[row.ind, col.ind]) -> y
  as.numeric(max[row.ind, col.ind]) -> z
  #as.data.frame(cbind(x,y,z)) -> coords
  cbind(x,y,z) -> coords
  return(coords)
}

# http://stackoverflow.com/questions/6681145/how-can-i-arrange-an-arbitrary-number-of-ggplots-using-grid-arrange
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

mapResponse <- function(res.df, assembled, drug.list.all=NULL, no.samples=1){
  # default: create dotplot of res.df results against assembled
  toupper(colnames(res.df)) -> colnames(res.df)
  toupper(colnames(assembled)) -> colnames(assembled)
  if (nrow(res.df) == 1 && length(which(is.na(res.df))) > 0){
    # remove drugs where no fit was performed
    res.df[,-which(is.na(res.df))] -> res.df
    t(as.data.frame(res.df)) -> res.df
  }
  
  sapply(colnames(res.df), function(x) 
    ifelse(length(grep(x, drug.list.all$all.names))==1, 
          drug.list.all$final.name[grep(x, drug.list.all$all.names)],
          x)) -> colnames(res.df)
  
  intersect(colnames(res.df), colnames(assembled)) -> common.drugs
  res.df[,which(colnames(res.df) %in% common.drugs)] -> res.df.sub
  assembled[,which(colnames(assembled) %in% common.drugs)] -> assembled.sub
  
  if (is.numeric(res.df.sub)){
    sapply(colnames(assembled.sub), function(x) 
      which(toupper(names(res.df.sub)) %in% x)) -> c.ind
    max(sapply(c.ind, function(x) length(x))) -> ml
    respMat <- matrix(NA, nrow=ml, ncol=length(c.ind))
    # create 
    #res.df.sub[c.ind] -> res.df.sub
    for (i in 1:length(c.ind)){
      respMat[1:length(c.ind[[i]]),i] <- res.df.sub[c.ind[[i]]]
    }
    colnames(respMat) <- colnames(assembled.sub)
  }
    
  rbind(assembled.sub, respMat) -> assembled.sub
  which(rownames(assembled.sub) %ni% 
          rownames(assembled)) -> r.ind
  rownames(assembled.sub)[which(rownames(assembled.sub) %ni% rownames(assembled))] <- 
    paste("POI",r.ind, sep=",")
  cbind(rownames(assembled.sub), assembled.sub) -> assembled.sub
  colnames(assembled.sub)[1] <- "id"
  melt(assembled.sub, id.vars="id") -> assembled.l
    
  class <- rep("OTHER", nrow(assembled.l))
  class[grep("POI", assembled.l$id)] <- "CHECK"
  cbind(assembled.l, class) -> assembled.l
  assembled.l[which(class %in% "OTHER"),] -> df.o
  assembled.l[which(class %in% "CHECK"),] -> df.i
  
  # order by drug class
  sapply(df.o$variable, function(x) 
    drug.list.all$classification[which(drug.list.all$final.name %in% x)]) -> dc
  sapply(df.o$variable, function(x) 
    drug.list.all$subclass[which(drug.list.all$final.name %in% x)]) -> dc2
  cbind(df.o, dc, dc2) -> df.o
  
  sapply(df.i$variable, function(x) 
    drug.list.all$classification[which(drug.list.all$final.name %in% x)]) -> dc
  sapply(df.i$variable, function(x) 
    drug.list.all$subclass[which(drug.list.all$final.name %in% x)]) -> dc2
  cbind(df.i, dc, dc2) -> df.i
  
  man.col <- scale_colour_manual(values = rep("red", no.samples))
  ggplot(df.o, aes(x=value))+
    stat_density(aes(ymax = ..density..,  ymin = -..density..),
                 fill = "grey50", colour = "grey50",
                 geom = "ribbon", position = "identity", alpha=0.2)+
    #facet_wrap(~variable, nrow=1)+coord_flip()+
    facet_wrap(~dc+variable)+coord_flip()+
    geom_jitter(data = df.o, 
                aes(x=value, y=0), alpha=0.2, width=0.005)+
    geom_jitter(data = df.i, 
                aes(x=value, y=0,shape=class,color=class), alpha=0.7, width=0.005,
                size=4.5)+
    scale_y_continuous(breaks=c(-2, 0, 2))+theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_text(size=14),
          strip.text=element_text(size=11),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          panel.margin = unit(0.0, "lines"))+
    man.col+theme(legend.position="bottom")+xlim(-0.6,4.5)
}

addDrug <- function(curr.df, all.drugs, patient, group){
  vis.df <- matrix(0, nrow=0, ncol=6)
  colnames(vis.df) <- c("Concentration", "Response", "Replicate", "Drug", "Patient", "Group")
  for (i in 1:length(all.drugs)){
    # add to data.frame
    grep(all.drugs[i], colnames(curr.df)) -> m
    for (j in 1:length(m)){
      cbind(curr.df$Concentrations, curr.df[,m[j]], 
            rep(j, nrow(curr.df)), rep(all.drugs[i], nrow(curr.df)),
            rep(patient, nrow(curr.df)),
            rep(group, nrow(curr.df))) -> t
      colnames(t) <- c("Concentration", "Response", "Replicate", "Drug", "Patient", "Group")
      rbind(vis.df, t) -> vis.df
    }
  }
  return(vis.df)
}

plotRaw <- function(all.resp, mode=c("", "r")){
  # plot all raw data based on all.resp
  # vis.df format: Conc, Response, Replicate, Drug
  vis.df <- matrix(0, nrow=0, ncol=6)
  colnames(vis.df) <- c("Concentration", "Response", "Replicate", "Drug", "Patient", "Group")
  for (i in 1:length(all.resp)){
    # find unique drugs, and plot replicates together
    for (j in 1:length(all.resp[[i]])){
      all.resp[[i]][[j]] -> curr.df
      unique(sapply(strsplit(colnames(curr.df)[2:ncol(curr.df)], "_"), function(x)
        x[1])) -> all.drugs
      addDrug(curr.df, all.drugs, names(all.resp)[i], j) -> t
      rbind(vis.df, t) -> vis.df
    }
  }
  as.data.frame(vis.df) -> vis.df
  for (i in 1:2){
    as.numeric(as.character(vis.df[,i])) -> vis.df[,i]
  }
  log10(vis.df$Concentration) -> vis.df$Concentration
  gsub("\\.", "-", vis.df$Drug) -> vis.df$Drug
  toupper(vis.df$Drug) -> vis.df$Drug
  toupper(vis.df$Patient) -> vis.df$Patient
  as.numeric(as.character(vis.df.sub$Replicate)) -> vis.df.sub$Replicate
  unique(vis.df$Drug) -> all.drugs
  for (i in 1:length(all.drugs)){
    # split drugs into separate data frames, then handle per groups
    vis.df[which(vis.df$Drug %in% all.drugs[i]),] -> vis.df.sub
    if (length(unique(vis.df.sub$Group)) > 1){
      # Drug was in more than one plate, more than one concentration
      lapply(unique(vis.df.sub$Group), function(x)
        unique(vis.df.sub$Replicate[which(vis.df.sub$Group %in% x)])) -> all.reps
      unique(sapply(all.reps, function(x) length(x))) -> dups
      if (length(unique(sapply(all.reps, function(x) paste(x, collapse=",")))) < 
          length(unique(vis.df.sub$Group))){
        # for each unique group, change to the next possible number
        ctr <- 1
        for (j in 1:length(unique(vis.df.sub$Group))){
          as.numeric(as.character(vis.df.sub$Replicate[which(vis.df.sub$Group %in% vis.df.sub$Group[j])])) -> curr.reps
          for (k in 1:length(unique(curr.reps))){
            curr.reps[which(curr.reps %in% unique(curr.reps)[k])] <- ctr
            ctr <- ctr+1
          }
          vis.df.sub$Replicate[which(vis.df.sub$Group %in% vis.df.sub$Group[j])] <- curr.reps
        }
      }
    }
    vis.df[which(vis.df$Drug %in% all.drugs[i]),] <- vis.df.sub
  }
  

  if (mode == ""){
    for (i in 1:length(unique(vis.df$Patient))){
      print(paste("Plotting response for patient", unique(vis.df$Patient)[i], "...", sep=" "))
      vis.df[which(vis.df$Patient %in% unique(vis.df$Patient)[i]),] -> vis.df.sub
      ggplot(vis.df.sub, aes(x=Concentration, y=Response, color=factor(Replicate), 
                 group=factor(Replicate)))+geom_point()+geom_line()+
                 facet_wrap(~Drug)+theme_bw()+
                 ggtitle(unique(vis.df$Patient)[i])-> p
      print(p)
    }
  } else {
    return(vis.df)
  }
}

visPlates <- function(dir){
  # get all .csv/.txt files and create a heatmap of each
  list.dirs(dir, recursive=T) -> sd
  lapply(sd, function(x) getFiles(x)) -> files
  names(files) <- as.character(sapply(sd, function(x) split_path(x)[1]))
  
  unique(unlist(files)) -> ff
  for (k in 1:length(ff)){
    if (length(grep(".csv", as.character(ff[k]))) > 0){
      t(as.matrix(read.csv(as.character(ff[k]), stringsAsFactors=T, head=F))) -> p
      aheatmap(p, Rowv=NA, Colv=NA)
    }
  }
}

plotFits <- function(exp.res){
  plotRaw(exp.res$experiments, mode="r") -> vis.df
}