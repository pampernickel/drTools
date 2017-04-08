# PDobay for AG Bourquin, KISPI

#@...
# Routines for reading and processing combination experiments
#@...

library(reshape2)

readCombos <- function(dir, res.dir, mode=c("", "normalized")){
  # no.pat: number of patients per plate
  
  print(paste("Processing directory ", dir, "...", sep=""))
  unlist(strsplit(dir, "/")) -> dir.name
  
  c(list.files(path = dir, pattern = ".txt", recursive=TRUE, full.names = TRUE), 
    list.files(path = dir, pattern = ".csv", recursive=TRUE, full.names = TRUE)) -> files  
  c(list.files(path = res.dir, pattern = ".txt", recursive=TRUE, full.names = TRUE), 
    list.files(path = res.dir, pattern = ".csv", recursive=TRUE, full.names = TRUE)) -> res.files
  
  if (length(files) == 0){
    stop(paste("No files found in ", dir, ".", sep=""))
  }
  
  if (length(res.files) == 0){
    stop(paste("No files found in ", res.dir, ".", sep=""))
  }
  
  if (length(grep(".txt", files)) > 0 && mode=="normalized"){ # 
    lapply(files, function(x) suppressWarnings(readLines(x))) -> contents
    lapply(contents, function(x) c(grep("mM", x),grep("mL", x))) -> infoLines
    lapply(infoLines, function(x) sort(x)) -> infoLines
    lapply(c(1:length(contents)), function(x)
      strsplit(contents[[x]][infoLines[[x]]], "\t")) -> meta
    
    if (length(infoLines[[1]]) == 2){
      # case of one pair of combos per plate (can be with one or two
      # patients)
      sapply(infoLines[[1]], function(x)
        sapply(strsplit(contents[[1]][x], "\t"), function(y) y[1])) -> drug.names
      getPlate(infoLines, contents, 1, 1) -> p1
      getPlate(infoLines, contents, 1, 2) -> p2
      which(p1 != "", arr.ind = T) -> p1.ind
      which(p2 != "", arr.ind = T) -> p2.ind
      
      # check area of overlap bet p1.ind and p2.ind
      apply(p1.ind, 1, function(x) 
        ifelse(x[1] %in% p2.ind[,1] && x[2] %in% p2.ind[,2], T, F)) -> o.1
      apply(p2.ind, 1, function(x) 
        ifelse(x[1] %in% p1.ind[,1] && x[2] %in% p1.ind[,2], T, F)) -> o.2
      
      p1.ind[which(o.1 %in% T),] -> o.coords
      p1.ind[which(o.1 %in% F),] -> d1.only
      p2.ind[which(o.2 %in% F),] -> d2.only
      
      combo.mats <- NA
      if (length(unique(d1.only[,1])) > length(unique(d1.only[,2])) &&
          length(unique(d2.only[,1])) < length(unique(d2.only[,2]))){
        lapply(res.files, function(x){
          read.csv(x, header=F) -> content
          t(content) -> content
          combo.mat <- list()
          
          length(unique(d1.only[,2])) -> no.combos
          c(nrow(o.coords)/no.combos, nrow(o.coords)) -> lims
          start <- start.1 <- 1
          for (j in 1:no.combos){
            as.numeric(as.character(p1[d1.only[which(d1.only[,2] %in% 
                                                       unique(d1.only[,2])[j]),1],d1.only[j,2]])) -> d1.doses
            as.numeric(as.character(p2[unique(d2.only[1:length(d1.doses),1]),
                                       d2.only[1:length(d1.doses),2]])) -> d2.doses
            
            mat <- matrix(NA, nrow=length(d1.doses)+2, ncol=length(d2.doses)+2)
            mat[1,1] <- "XX"
            mat[2:nrow(mat),1] <- c(d1.doses,0)
            mat[1,2:ncol(mat)] <- c(0,d2.doses)
            
            for (k in start:lims[j]){
              mat[which(mat[,1] %in% p1[o.coords[k,1],o.coords[k,2]]),
                  which(mat[1,] %in% p2[o.coords[k,1],o.coords[k,2]])] <-
                content[o.coords[k,1],o.coords[k,2]]
            }
            start <- lims[j]+1
            
            # fill content for single drug treatment
            c(nrow(d1.only)/length(unique(d1.only[,2])),
              nrow(d1.only)) -> lims.1
            mat[2:(nrow(mat)-1),2] <- as.numeric(content[d1.only[c(start.1:lims.1[j]),1],
                                                         unique(d1.only[c(start.1:lims.1[j]),2])])
            mat[nrow(mat), 3:ncol(mat)] <- as.numeric(content[unique(d2.only[c(start.1:lims.1[j]),1]),
                                                              d2.only[c(start.1:lims.1[j]),2]])
            start.1 <- lims.1[j]+1
            
            # then get all contents that are NOT part of o.coords, d1.only or d2.only,
            # which correspond with the DMSO-containing wells; the value here would be in the 0,0, and
            # used in the normalization of the plate
            rbind(o.coords, d1.only, d2.only) -> all.coords
            which(content != 0, arr.ind = T) -> plate.coords
            plate.coords[which(apply(plate.coords, 1, function(x)
              x[1] %in% all.coords[,1] && x[2] %in% all.coords[,2]) %in% F),] -> dmso.coords
            dmso.coords[order(dmso.coords[,1]),] -> dmso.coords
            
            # split dmso.coords as well
            dmso.coords[which(dmso.coords[,1] %in% unique(dmso.coords[,1])[j]),] -> dmso.coords.sub
            mean(apply(dmso.coords.sub, 1, function(x)
              content[x[1],x[2]])) -> dmso.mean
            mat[nrow(mat),2] <- dmso.mean
            apply(mat[2:nrow(mat),2:ncol(mat)],2, function(y) as.numeric(as.character(y))) -> mat_n
            100*mat_n/dmso.mean -> mat_n
            rownames(mat_n) <- mat[2:nrow(mat),1]
            colnames(mat_n) <- mat[1, 2:ncol(mat)]
            mat_n -> combo.mat[[j]]
          }
          names(combo.mat) <- rep(paste(drug.names, collapse="_"), length(combo.mat))
          return(combo.mat)
        }) -> combo.mats
        names(combo.mats) <- res.files
      }
    }
  }
  return(combo.mats)
}

processCombos <- function(combos, additivity=c("HSA", "Loewe", "Bliss")){
  # prepare combos for visualization (double-plot)
  # flatten list, i.e. have all combo data frames in a single structure
  unlist(combos, recursive = FALSE) -> cl
  
  all.combos <- list()
  for (i in 1:length(cl)){
    temp <- matrix(NA, nrow=0, ncol=4)
    colnames(temp) <- c("x", "y", "combo", "patient")
    
    cl[[i]] -> main
    getComboProperties(cl, i) -> meta

    drug1.doses <- as.numeric(rownames(cl[[i]]))
    drug2.doses <- as.numeric(colnames(cl[[i]]))
    
    # as this is a double-axis plot, need to keep the x-axis uniform (i.e. just intervals, no
    # actual doses, especially if the doses of drugs1 and 2 are not the same)
    drug1.doses <- drug1.doses[1:length(drug1.doses)-1]
    drug2.doses <- drug2.doses[2:length(drug2.doses)]
    
    d1.doses.proxy <- 1:length(drug1.doses)
    d2.doses.proxy <- 1:length(drug2.doses)
    if (length(which(diff(drug1.doses) < 0)) == length(drug1.doses)-1){
      # drug1 doses are descending
      d1.doses.proxy <- length(drug1.doses):1
    }
    
    if (length(which(diff(drug2.doses) < 0)) == length(drug2.doses)-1){
      # drug1 doses are descending
      d2.doses.proxy <- length(drug2.doses):1
    }
    
    main[nrow(main),2:ncol(main)] -> drug2.resp
    main[2:nrow(main)] -> drug1.resp
    
    main[c(2:nrow(main)),c(2:ncol(main))] -> resp.matrix
    
    t <- cbind(d1.doses.proxy, smooth(drug1.resp, kind = "3R"), 
               rep(meta$drug1, length(drug1.resp)), 
               rep(meta$pat, length(drug1.resp)))
    colnames(t) <- colnames(temp)
    rbind(temp, t) -> temp
    
    t <- cbind(d2.doses.proxy, smooth(drug2.resp, kind = "3R"), 
               rep(meta$drug2, length(drug2.resp)), 
               rep(meta$pat, length(drug2.resp))) # replace test with var fin.name
    colnames(t) <- colnames(temp)
    rbind(temp, t) -> temp
    
    # drug names
    if (nchar(meta$drug1) > 4){
      d1 <- substr(meta$drug1, 1, 4)
    } else {
      d1 <- meta$drug1
    }
    
    if (nchar(meta$drug2) > 4){
      d2 <- substr(meta$drug2, 1, 4)
    } else {
      d2 <- meta$drug2
    }
    
    # TO DO: add additivity line; select from options
    # Loewe additivity: effect of a drug if it were combined with itself ye = y1(x1+x2), x1 & 2 are doses
    # HSA: min(y1, y2), i.e. highest monotherapy effect
    # Bliss: ye=y1+y2-y1y2
    
    additivity.line <- NA
    if (additivity == "HSA"){
      apply(rbind(smooth(drug1.resp),
                  smooth(drug2.resp)), 2, function(x) min(x)) -> additivity.line
    } else if (additivity == "Loewe"){
      # fit individual drug responses to approximate the effect of doubling the dose
      addFit(drug1.doses, as.numeric(smooth(drug1.resp)), max(drug1.doses)) -> f1
      addFit(drug2.doses, as.numeric(smooth(drug2.resp)), max(drug2.doses)) -> f2
      
      # cross section based on drug with lower IC50 -- or do a default cross section
      # through drug2 if the dr is with respect to drug1
      ff <- extractMax(f2$max)
      if (f1$logIC50 <= f2$logIC50){
        extractMax(f1$max) -> ff
        # choose cross-sections based on the more active drug:
        # cross-section through columns
        for (k in 1:ncol(resp.matrix)){
          as.numeric(resp.matrix[,k]) -> combo.resp
          t <- cbind(d1.doses.proxy, as.numeric(smooth(combo.resp, kind = "3R")), 
                     rep(paste(round(drug2.doses[k],2), " ", d2, "; ",
                               round(drug1.doses[k],2), " ", d1, sep=""), 
                         length(combo.resp)), 
                     rep(meta$pat, length(combo.resp)))
          colnames(t) <- c("x", "y", "combo", "patient")
          rbind(temp, t) -> temp
        }
      } else {
        # choose cross-sections based on the more active drug:
        # cross-section through rows
        for (k in 1:nrow(resp.matrix)){
          as.numeric(resp.matrix[k,]) -> combo.resp
          t <- cbind(d2.doses.proxy, as.numeric(smooth(combo.resp, kind = "3R")), 
                     rep(paste(round(drug1.doses[k],2), " ", d1, "; ",
                               round(drug2.doses[k],2), " ", d2, sep=""), 
                         length(combo.resp)), 
                     rep(meta$pat, length(combo.resp)))
          colnames(t) <- c("x", "y", "combo", "patient")
          rbind(temp, t) -> temp
        }
      }
      
      max(c(drug2.doses, drug1.doses))/2^seq(1,9,by=1) -> sq
      
      
      
      # for each sq, check closest dose in fit
      sapply(sq, function(x){
        abs(ff$x-x) -> diff
        which(diff < 5) -> ind
        ff$x[which(diff %in% min(diff[ind]))][1] -> diff
      }) -> vals
      
      sapply(sq*2, function(x){
        abs(ff$x-x) -> diff
        which(diff < 5) -> ind
        ff$x[which(diff %in% min(diff[ind]))][1] -> diff
      }) -> vals.1
      ff$y[which(ff$x %in% vals.1)] -> additivity.line #x values for these will be vals, and not vals.1
      
      # check temp$x for orientation of additivity line
      rbind(seq(length(vals)-(length(vals)-1),nrow(temp),by=length(vals)),
            seq(length(vals),nrow(temp),by=length(vals))) -> stops
      additivity.line.fin <- rep(0, nrow(temp))
      for (k in 1:ncol(stops)){
        x <- stops[,k]
        if ((isDescending(as.numeric(temp[c(x[1]:x[2]),1])) && isDescending(additivity.line)) ||
            (!isDescending(as.numeric(temp[c(x[1]:x[2]),1])) && !isDescending(additivity.line))){
          additivity.line.fin[c(x[1]:x[2])] <- additivity.line[length(additivity.line):1]
        } else if ((!isDescending(as.numeric(temp[c(x[1]:x[2]),1])) && isDescending(additivity.line)) ||
                   isDescending(as.numeric(temp[c(x[1]:x[2]),1])) && !isDescending(additivity.line)){
          additivity.line.fin[c(x[1]:x[2])] <- additivity.line
        }
      }
      
      cbind(temp, additivity.line.fin) -> temp
      as.data.frame(temp) -> temp
      as.numeric(as.character(temp$x)) -> temp$x
      as.numeric(as.character(temp$y)) -> temp$y
      as.numeric(as.character(temp$additivity.line.fin)) -> temp$additivity.line.fin
      cols <- colorRampPalette(brewer.pal(8, "RdBu"))(11)[11:1]
      temp$combo <- factor(temp$combo, levels=unique(temp$combo)) 
      ggplot(temp, aes(x=x, y=y, group = combo, colour = combo, 
                       ymin = 0, ymax = additivity.line.fin))+
        geom_ribbon(alpha=0.05, color="grey", fill="grey")+
        geom_line(size=1)+
        scale_color_manual(values=c(cols[1], cols[length(cols)], cols[2:(length(cols)-1)]))+
        scale_x_log10()+
        ylim(0, 200)
      
    } else if (additivity == "Bliss"){
      
    }
    
    temp -> all.combos[[i]]
  }
  return(all.combos)
}

calcCI <- function(combos){
  unlist(combos, recursive = FALSE) -> cl
  df <- matrix(0, nrow=0, ncol=4)
  colnames(df) <- c("patient", "drug1", "drug2", "CI")
  for (i in 1:length(cl)){
    cl[[i]] -> main  
    getComboProperties(cl, i) -> meta
    cbind(rownames(main), main) -> main
    colnames(main)[1] <- meta$drug2
    apply(main, 2, function(x) as.numeric(as.character(x))) -> main
    shapeA(as.data.frame(main), drug1 = meta$drug1, drug2 = meta$drug2) -> drMatrix
    as.numeric(IC50(drMatrix)[3]) -> f
    cbind(meta$pat, meta$drug1, meta$drug2, f) -> t
    colnames(t) <- colnames(df)
    rbind(df, t) -> df
  }
  as.data.frame(df) -> df
  as.numeric(as.character(df$CI)) -> df$CI
  
  # check cases tending asymptotically to infinity or to 0
  df$CI[which(df$CI > 3)] <- 3
  df$CI[which(df$CI < 0.01)] <- 0.01
  return(df)
}

getComboProperties <- function(cl,i){
  unlist(strsplit(sapply(strsplit(names(cl)[i], "\\."), function(x) x[length(x)]), "_"))[1] -> d1
  unlist(strsplit(sapply(strsplit(names(cl)[i], "\\."), function(x) x[length(x)]), "_"))[2] -> d2
  unlist(strsplit(sapply(strsplit(names(cl)[i], "\\."), function(x) x[2]), "/")) -> pn
  gsub("_cumul_1_v_1_2", "", pn) -> pn
  paste(pn[length(pn)], i, sep="_") -> pat
  list(d1, d2, pat) -> meta
  names(meta) <- c("drug1", "drug2", "pat")
  return(meta)
}

# ::: synergyFinder
calcSynergy <- function(responses){
  if(!is.loaded("synergyfinder")) library(synergyfinder)
  # responses is a list containing 
  # various dose response matrices with the same
  # drug pairs and concentration ranges
  # the last element of the list is called drug.pairs;
  # matrices are nested lists `dose.response.mats'
  
}