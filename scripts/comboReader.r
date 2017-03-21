# PDobay for AG Bourquin, KISPI

#@...
# Routines for reading and processing combination experiments
#@...

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
      getPlate(infoLines, contents, i, 1) -> p1
      getPlate(infoLines, contents, i, 2) -> p2
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
            mat_n/dmso.mean -> mat_n
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

processCombos <- function(combos){
  # flatten list, i.e. have all combo data frames in a single structure
  unlist(combos, recursive = FALSE) -> cl
  
  for (i in 1:length(cl)){
    temp <- matrix(NA, nrow=0, ncol=4)
    colnames(temp) <- c("x", "y", "combo", "patient")
    
    cl[[i]] -> main
    unlist(strsplit(sapply(strsplit(names(cl)[i], "\\."), function(x) x[2]), "/")) -> pn
    gsub("_cumul_1_v_1_2", "", pn) -> pn
    paste(pn[length(pn)], i, sep="_") -> pat
    
    unlist(strsplit(sapply(strsplit(names(cl)[i], "\\."), function(x) x[length(x)]), "_"))[1] -> d1
    unlist(strsplit(sapply(strsplit(names(cl)[i], "\\."), function(x) x[length(x)]), "_"))[2] -> d2
    
    drug1.doses <- as.numeric(rownames(cl[[i]]))
    drug2.doses <- as.numeric(colnames(cl[[i]]))
    
    drug1.doses <- drug1.doses[1:length(drug1.doses)-1]
    drug2.doses <- drug2.doses[2:length(drug2.doses)]
    
    main[nrow(main),2:ncol(main)] -> drug2.resp
    main[2:nrow(main)] -> drug1.resp
    
    main[c(2:nrow(main)),c(2:ncol(main))] -> resp.matrix
    
    t <- cbind(drug1.doses, smooth(drug1.resp, kind = "3R"), rep(d1, length(drug1.resp)), 
               rep(pat, length(drug1.resp)))
    colnames(t) <- colnames(temp)
    rbind(temp, t) -> temp
    
    t <- cbind(drug2.doses, smooth(drug2.resp, kind = "3R"), 
               rep(d2, length(drug2.resp)), 
               rep(pat, length(drug2.resp))) # replace test with var fin.name
    colnames(t) <- colnames(temp)
    rbind(temp, t) -> temp
    
    # cross-section through drug 1
    for (k in 1:nrow(resp.matrix)){
      as.numeric(resp.matrix[k,]) -> combo.resp
      t <- cbind(drug2.doses, smooth(combo.resp, kind = "3R"), 
                 rep(paste(d1, d2, round(drug1.doses[k],2),sep="."), 
                     length(combo.resp)), 
                 rep(pat, length(combo.resp)))
      colnames(t) <- c("x", "y", "combo", "patient")
      rbind(temp, t) -> temp
    }
  }
  
  # add additivity line
  
}