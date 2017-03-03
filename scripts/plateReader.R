# PDobay for AG Bourquin, KISPI

#@...
# Layout processing functions
#@...


readFormat <- function(dir, replicates=2, dups, dup.mode, dilution=12.5){
  # function that returns a list of length n, where
  # n corresponds to the number of drugs on the plate + one slot for the control
  # for each drug, returns a structure
  # that contains the doses, location (rows and columns) for each
  # set of replicates; in all cases, dose rows and column rows MUST
  # be the same length as the number of replicates
  
  # dups and dup mode are vectors the same length as the number of plates/
  # layouts; dups can be implicit or explicit,
  # dup.mode can be rows or columns
  print(paste("Processing directory ", dir, "...", sep=""))
  unlist(strsplit(dir, "/")) -> dir.name
    
  c(list.files(path = dir, pattern = ".txt", recursive=TRUE, full.names = TRUE), 
  list.files(path = dir, pattern = ".csv", recursive=TRUE, full.names = TRUE)) -> files  
  if (length(files) == 0){
    stop(paste("No files found in ", dir, ".", sep=""))
  }
  
  if (length(dups) < length(files) |
      length(dup.mode) < length(files)){
    stop("Length of duplicate modes not the same as number of layout templates.")
  }
  
  if (length(grep(".txt", files)) > 0){ # 
    lapply(files, function(x) suppressWarnings(readLines(x))) -> contents
    lapply(contents, function(x) c(grep("mM", x),grep("mL", x))) -> infoLines
    lapply(infoLines, function(x) sort(x)) -> infoLines
    lapply(c(1:length(contents)), function(x)
      strsplit(contents[[x]][infoLines[[x]]], "\t")) -> meta
    meta1 <- list()
    for (x in 1:length(contents)){
      all.res <- list()
      for (j in 1:length(infoLines[[x]])){
        getPlate(infoLines, contents, x, j) -> plate
        # case of all drugs, where there is a serial dilution
        if (length(grep("control", tolower(meta[[x]][[j]]))) == 0){
            getDoses(plate, dups[x]) -> doses
            if (dups[x] == "ei"){ #explicit, irregular
              all.res[[j]] <- doses
            } else if (replicates==2 && dups[x] == "i") {
              as.numeric(plate[doses$dose.rows[1],]) -> t # used for finding the columns/rows of dups
              all.res[[j]] <- handleImplicit(dup.mode[x], t, doses, dilution)
            } else if (dups[x] == "e"){
              # replicates are explicitly indicated
              as.numeric(plate[doses$dose.rows[1],]) -> t # used for finding the columns/rows of dups
              all.res[[j]] <- handleExplicit(dup.mode[x], t, doses, dilution, replicates)
            }
        } else {
            # control: single dose only; problem is to localize it
            # in particular, there's a need to solve problems related to
            # tabs filled in automatically by R when the plate shape is uneven
            which(apply(plate, 1, function(y)
              length(which(y %in% ""))) < ncol(plate)) -> dose.rows
            # for each row, find the columns
            lapply(dose.rows, function(y)
              which(plate[y,] %ni% "")) -> control.cols
            
            # t(apply(plate[dose.rows,], 1, 
            #         function(x) suppressWarnings(as.numeric(x)))) -> controls
            # apply(controls, 1, function(x) length(which(is.na(x)))) -> r 
            # 
            # controls[which(r < ncol(plate)),] -> controls
            # dose.rows[which(r < ncol(plate))] -> dose.rows
            # t(apply(controls, 1, function(x) which(!is.na(x)))) -> control.cols
            # unique(as.numeric(dose.rows)) -> dose.rows
            # lapply(1:nrow(control.cols), function(x) unlist(list(control.cols[x,]))) -> control.cols
            # 
            # # get measurement for control
            list(NA, as.numeric(dose.rows), 
                  control.cols) -> all.res[[j]]
            names(all.res[[j]]) <- c("doses", "rows", "cols")
          }
      }
      
      meta1[[x]] <- all.res
      names(meta1[[x]]) <- sapply(meta[[x]], function(z) return(z[1]))
    }
  }
  names(meta1) <- paste("plate", c(1:length(meta1)), sep="")
  return(meta1)
}


getFiles <- function(dir){
  #print(paste("Checking directory ", dir, "...", sep=""))
  
  #list.files(path=dir, full.names=T) -> f
  c(list.files(path = dir, pattern = ".txt", recursive=TRUE, full.names = TRUE), 
    list.files(path = dir, pattern = ".csv", recursive=TRUE, full.names = TRUE)) -> files
  # draft option to keep subdirectory structure
  # lapply(f, function(x){
  #   l <- x
  #   if (length(grep(".csv", x)) > 0){
  #     return(l)
  #   } else {
  #     c(list.files(path = x, pattern = ".txt", recursive=TRUE, full.names = TRUE), 
  #       list.files(path = x, pattern = ".csv", recursive=TRUE, full.names = TRUE)) -> l
  #   }
  #   return(l)
  # }) -> files
  return(files)
}

getPlate <- function(infoLines, contents, x, j){
  plate <- NA
  infoLines[[x]][j]+2 -> s
  if (j < length(infoLines[[x]])){
    infoLines[[x]][j+1]-1 -> e
  } else {
    if (length(which(diff(infoLines[[x]][c(2:length(infoLines[[x]]))]-
                          infoLines[[x]][c(1:(length(infoLines[[x]])-1))]) != 0)) == 0){
      e <-  length(contents[[x]])
    }
  }

  contents[[x]][c(s:e)] -> plate
  if (length(c(which(plate %in% ""),which(is.na(plate)))) > 0){
    plate[-c(which(plate %in% ""),
             which(is.na(plate)))] -> plate
  }
  
  # in the case of DMSO, results in repetition of "P" in the last line
  # case where input has problems i.e. different no .of rows
  sapply(plate, function(y) strsplit(y, "\t")) -> splits
  unique(sapply(splits, function(x) length(x))) -> ncols
  if (length(ncols) > 1){
    warning("Some entries do not have the same number of columns. Check the tab-delimited file. Attempting to pad...")
    # get longest n of cols
    which(sapply(splits, function(x) length(x)) %in% max(ncols))[1] -> max.col

    # create matrix based on max.col and nrows
    plate.sub <- matrix(NA, nrow=length(splits), ncol=max(ncols))
    for (i in 1:length(splits)){
      plate.sub[i,1:length(splits[[i]])] <- splits[[i]]
    }
    plate.sub -> plate
  } else {
    t(as.data.frame(sapply(plate, function(y) strsplit(y, "\t")))) -> plate
    rownames(plate) <- plate[,1]
    plate[,-1] -> plate  
  }
  return(plate)
}

getDoses <- function(plate, dups){
  if (dups == "i" | dups == "e"){
    which(apply(plate, 1, function(y)
      length(which(y %in% ""))) < ncol(plate)) -> dose.rows
    # also check if the doses are found in a single column or not
    # if the doses are NOT found in a single column for multi-row
    # cases, these are priming spots
    which(apply(plate, 2, function(y)
      length(which(y %in% ""))) < nrow(plate)) -> dose.cols
    
    which(plate[dose.rows,] %in% 
            max(as.numeric(plate[dose.rows,]),na.rm=T)) -> max.loc
    which(plate[dose.rows,] %in% 
            min(as.numeric(plate[dose.rows,]),na.rm=T)) -> min.loc
    
    # check if there are replicates or not
    as.numeric(plate[dose.rows[1],which(!is.na(as.numeric(plate[dose.rows[1],])))]) -> doses
    #unique(doses)/dilution -> doses
    as.numeric(unique(doses)) -> doses
    list(doses, dose.rows) -> res
    names(res) <- c("doses", "dose.rows")
  } else if (dups == "ei"){ #ei: explicit, irregular
    # test on other "e" cases
    which(apply(plate, 1, function(y)
      length(which(y %in% ""))) < ncol(plate)) -> dose.rows
    which(apply(plate, 2, function(y)
      length(which(y %in% ""))) < nrow(plate)) -> dose.cols
    
    plate[dose.rows,] -> dr
    if (length(dose.rows) > 1 && length(dose.cols) > 1){
      as.vector(dr) -> dr
      as.numeric(dr[which(dr %ni% "")]) -> dr
      
      # check intervals: detect and restrict to serial dilutions
      dr[rev(order(dr))] -> dr
      table(round((dr[1:(length(dr)-1)]/dr[2:length(dr)]), 3)) -> int.counts
      dr[which(round((dr[1:(length(dr)-1)]/dr[2:length(dr)]), 3) %in%
              names(int.counts)[which(int.counts %in% max(int.counts))])] -> doses
      
      if (length(dose.rows) < length(dose.cols)/2){
        # indicates that the duplicates are spread across rows; treat like 
        # explicit cases; check dose cols associated with each dose.row
        lapply(dose.rows, function(x){
          which(plate[x,] %ni% "")
        }) -> cols
        as.numeric(plate[dose.rows[1], cols[[1]]]) -> doses #make sure that doses are in a correct orientation
        list(doses, dose.rows, cols) -> res
        names(res) <- c("doses", "rows", "cols")
      } else if (length(dose.rows) == length(dose.cols) ||
                 length(dose.rows) != length(dose.cols)) {
        # replicates distributed over rows and columns
        # dose.rows: columns
        # dose.cols: rows
        # find these again as contigs
        find.cols <- find.rows <- list()
        for (i in 1:length(doses)){
          # find those in a row
          which(sapply(apply(plate, 1, function(x) 
            grep(doses[i], x)), function(y) length(y)) != 0) -> find.rows[[i]]
          
          which(sapply(apply(plate, 2, function(x) 
            grep(doses[i], x)), function(y) length(y)) != 0) -> find.cols[[i]]
        }
        
        col <- names(table(unlist(find.cols)))[which(table(unlist(find.cols)) %in% 
                                                       max(table(unlist(find.cols))))]
        row <- names(table(unlist(find.rows)))[which(table(unlist(find.rows)) %in% 
                                                       max(table(unlist(find.rows))))]
        
        # then for both row and col, find cols and rows, respectively
        which(plate[as.numeric(row),] %in% doses) -> row.cols
        which(plate[,as.numeric(col)] %in% doses) -> col.rows
        
        row.cols[c(which(diff(row.cols) %in% 1),
                   (length(which(diff(row.cols) %in% 1))+1))] -> row.cols
        col.rows[c(which(diff(col.rows) %in% 1),
                   (length(which(diff(col.rows) %in% 1))+1))] -> col.rows
        as.numeric(plate[col.rows,as.numeric(col)]) -> doses #make sure that doses are in a correct orientation
        list(doses, 
             list(col.rows, as.numeric(row)), 
             list(as.numeric(col), row.cols)) -> res
        names(res) <- c("doses", "rows", "cols")
      } else {
        stop("Tough luck -- this part is not yet handled. Contact me at maria.pamela.david@gmail.com.")
      }
      
    } else {
      as.vector(dr) -> dr
      as.numeric(dr[which(dr %ni% "")]) -> doses
  
      # for each dose, get row and column coords
      as.numeric(sapply(doses, function(x)
        unlist(apply(plate, 1, function(y) which(y %in% x))))) -> cols
      
      as.numeric(sapply(doses, function(x)
        unlist(apply(plate, 2, function(y) which(y %in% x))))) -> rows
      
      
      # check if there are replicates or not
      rbind(rows, cols) -> coords
      list(doses, coords) -> res
      names(res) <- c("doses", "coords")
    }
  }
  return(res)
}

handleImplicit <- function(mode, t, doses, dilution){
  res <- list()
  rep.ind <- NA
  if (length(grep("c", mode)) > 0){
    # get duplicate columns
    which(t %in% doses$doses) -> col.ind
    col.ind+1 -> rep.ind
    list(doses$doses/dilution, c(as.numeric(doses$dose.rows), 
                                 as.numeric(doses$dose.rows)), 
         list(col.ind, rep.ind)) -> res
    names(res) <- c("doses", "rows", "cols")
  }
  return(res)  
}

handleExplicit <- function(mode, t, doses, dilution, replicates){
  res <- list()
  if (length(grep("r", mode)) > 0 | length(grep("n", mode)) > 0){
    sapply(doses$doses, function(y)
      which(t %in% y)) -> all.inds
    all.cols <- all.rows <- list()
    for (i in 1:length(replicates)){
      all.cols[[i]] <- all.inds
    }
    list(doses$doses/dilution, as.numeric(doses$dose.rows), 
         all.cols) -> res
    names(res) <- c("doses", "rows", "cols")
  } else if (length(grep("c", mode)) > 0){
    # cols processing, also handle case
    # when a drug appears twice on a plate
    sapply(doses$doses, function(y)
      which(t %in% y)) -> all.inds
    t(all.inds) -> cols
    
    # check if the number of rows and columns match
    # as.numeric(rep(doses$dose.rows, replicates))  -> rows
    rows <- c()
    for (i in 1:length(doses$dose.rows)){
      c(rows, rep(doses$dose.rows[i], replicates)) -> rows
    }
    
    all.inds <- list()
    for (i in 1:ncol(cols)){
      cols[,i] -> all.inds[[i]]
    }
    
    if (length(rows) > ncol(cols)){
      all.inds <- list()
      ctr <- 1
      for (j in 1:replicates){
        for (i in 1:ncol(cols)){
          cols[,i] -> all.inds[[ctr]]
          ctr <- ctr+1
        }
      }
    }
    
    list(doses$doses/dilution, 
         rows,
         all.inds) -> res
    names(res) <- c("doses", "rows", "cols")
  }
  return(res)
}

#@...
# File compliance checks
#@...
checkFiles <- function(files, layout){
  print("Checking file:layout match...")
  if (length(layout) > 1){
    sapply(files, function(x) ifelse(length(x)==length(layout),T,F)) -> l
    if (length(which(l %in% F)) > 0){
      stop("Folder", which(l %in% F), " does not have the same number 
                    of files as layouts.", sep=" ")
    }
      # then check if we can match the layout names to each file name
    lapply(files, function(y) 
      as.numeric(sapply(names(layout), function(x)
        grep(paste(x, "_", sep=""), y, ignore.case=T)))) -> matches
    sapply(matches, function(x)
      ifelse(length(which(is.na(x))==0),T,F)) -> p
    
    if (length(which(p %in% T)) > 0){
      stop("Folder", which(p %in% T), " does not have the same number 
                    of files as layouts.", sep=" ")
    }
    
    lapply(files, function(y){
      y[as.numeric(unlist(sapply(names(layout), function(x){
        grep(paste(x, "_", sep=""), y, ignore.case=T)
      })))] -> y
    }) -> files
  }
  
  print("All files matched to layout.")
  return(files)
}

#@...
# Experiment processing functions
#@...
readExperiment <- function(files, layout, mode="", pos.control="",
                           resolve.warnings=F, historical.data="", drug.list.all=""){
  if (resolve.warnings %in% T && is.character(historical.data)){
    stop("Enter historical data in the form of a patient x drug response matrix to resolve response inconsistencies.")  
  }
  
  if (resolve.warnings %in% T && is.character(drug.list.all)){
    stop("Include a matrix matching the list of your drug names against the drug names in historical data.\nThe matrix should have at least two columns:
         all.names and final.name.")  
  }
  
  lapply(files, function(f){
    lapply(f, function(y){
      if (length(grep(".csv", y)) == 1){
        data <- read.csv(as.character(y), header = FALSE)
      } else if (length(grep(".csv", y)) == 1) {
        data <- read.table(as.character(y), header = FALSE, sep = "\t")
      }
    }) -> all.dat
    
    # start with layout
    fin.resp <- list()
    controls <- list()
    for (i in 1:length(layout)){
      if (mode==""){
        t(all.dat[[i]]) -> curr.plate      
      } else {
        all.dat[[i]] -> curr.plate
      }
      
      layout[[i]] -> curr.layout
      
      # check column batches; here, check if there are a bunch
      # of drugs that have the same columns
      ## treat controls separately
      curr.layout[grep("control", tolower(names(curr.layout)))] -> control
      processControl(control, curr.plate)$mean -> c.mean # dmso mean for plate
      as.numeric(processControl(control, curr.plate)$all.controls) -> controls[[i]]
      
      if (mode == "indirect"){
        curr.layout[which(names(curr.layout) %in% pos.control)] -> pos_control
        processPosControl(pos_control, curr.plate) -> pos
      }
      
      # get control, and return as a single value (i.e. mean of all
      # control values + sd)
      curr.layout[-grep("control", tolower(names(curr.layout)))] -> curr.layout
      
      # change if duplicates are implicit (e.g. in "ei" cases)
      names(curr.layout) -> drug.names
      processPlate(curr.layout, curr.plate, mode) -> fin.resp[[i]]
    }

    lapply(fin.resp, function(x) 
      lapply(x, function(y) y$resp)) -> all.resp
    lapply(fin.resp, function(x)
      lapply(x, function(y) y[2:3])) -> warnings.list
    # consolidate into plates that have the same concentrations
    lapply(all.resp, function(x) lapply(x, function(y) 
      return(y[,1]))) -> all.conc
    findUniqueConc(all.conc) -> unique.conc
    
    # for each class, group responses
    groupResponses(all.resp, unique.conc) -> all.resp.fin
    
    # then assign drug names based on unique conc; this is mainly
    # for single-rep instances
    if (length(files)==1 && length(grep("coords", names(layout[[1]][[1]]))) > 0){
      # explicit coords
      unique(unlist(unique.conc)) -> conc
      if (length(drug.names) == length(unlist(unique.conc))){
        unlist(unique.conc) -> uconc
        sapply(conc, function(x) drug.names[which(uconc %in% x)]) -> cnames
        for (i in 1:length(all.resp.fin)){
          colnames(all.resp.fin[[i]])[1] <- "Concentrations"
          colnames(all.resp.fin[[i]])[2:ncol(all.resp.fin[[i]])] <- cnames[[i]]
        }
      }
    }
    
    # visualize controls
    list(all.resp.fin, warnings.list, controls) -> all.resp.fin
    names(all.resp.fin) <- c("resp", "warnings", "controls")
    
    ## visualizeControls(all.resp.fin$controls)
    return(all.resp.fin)
  }) -> all.resp
  
  # process warnings within all.resp; add tag if warnings are to be resolved
  if (resolve.warnings %in% T){
    lapply(all.resp, function(x) x$warnings) -> all.warnings
    lapply(all.warnings, function(x){
      unlist(x) -> xu
      cbind(xu[seq(1,length(xu),by=2)], xu[seq(2,length(xu),by=2)]) -> mat
    }) -> warning.mat
    lapply(all.resp, function(x) x$controls) -> all.controls
    lapply(all.resp, function(x) x$resp) -> all.resp
    resolveWarnings(all.resp, warning.mat, historical.data, drug.list.all) -> all.resp
  } else {
    lapply(all.resp, function(x) x$resp) -> all.resp
  }
  
  return(all.resp)
}

resolveWarnings <- function(all.resp, warning.mat, historical.data, drug.list.all){
  # for each response with a warning, fit the data
  for (i in 1:length(all.resp)){
    all.resp[[i]][[1]] -> curr.resp
    gsub("([.-])|[[:punct:]]", "\\.", colnames(curr.resp)) -> colnames(curr.resp)
    as.character(warning.mat[[i]][,2]) -> curr.warnings
    gsub(" ", "\\.", gsub("([.-])|[[:punct:]]", "\\.", curr.warnings)) -> curr.warnings
    lapply(curr.warnings, function(x) c(1,grep(x, colnames(curr.resp), ignore.case=T))) -> cols
    
    # do a check to make sure that everything  is matched
    if (length(unique(sapply(cols, function(x) length(x)))) == 1){
      # fit with respect to position in concentration scale rather
      # than absolute concentration to approximate doses used in
      # historical data, i.e. scale as a function of range; note that 
      # in historical data, max conc was log10 = 4, and min was 
      # log10 ~ -0.5
      lapply(cols, function(x){
        apply(curr.resp[,x[-1]], 2, function(y) 
          removeOutliers(y, 10^seq(1:5))) -> curr.resp[,x[-1]]
        10^seq(1:5) -> curr.resp[,1]
        fit(curr.resp[,x], 1, 10000) -> f
        return(f)
      }) -> res
      
      lapply(res, function(x) x$logIC50) -> fits
      names(fits) <- curr.warnings
      
      # find closest drug in historical data
      gsub("\\.", "-", names(fits)) -> sub.name
      findClosestMatch(sub.name, drug.list.all) -> name.match
      lapply(1:length(fits), function(x){
        fits[[x]] -> y
        median(historical.data[,which(colnames(historical.data) %in% 
                                     name.match[x])], na.rm = T) -> med.act
        
        # check only cases where at least one is reported active;
        # high discordance can be expected anyway if the compound
        # is highly inactive
        if (length(which(y == 4)) != length(y)){
          # take the one that's closest to the historical median
          # select the first element with the minimum difference,
          # which(diff %in% min(diff))[1]
          abs(y-med.act) -> diff
          cols[[x]][2:length(cols[[x]])][which(diff %in% min(diff))[1]] ->
            cols[[x]][2:length(cols[[x]])]
        }
        return(cols[[x]])
      }) -> fin.fits
      
      # then change all.resp where fin.fits is different from cols
      which(mapply(function(x,y) 
        length(unique(x)) ==length(unique(y)), fin.fits, cols) %in% F) -> c.ind
      for (j in 1:length(c.ind)){
        cols[[c.ind[j]]][-1] -> orig.cols
        fin.fits[[c.ind[j]]][-1] -> rep.cols
        curr.resp[,orig.cols] <- curr.resp[,rep.cols]
        curr.resp[,1] <- all.resp[[i]][[1]]$Concentrations
      }
    }
    
    gsub("\\.2", "_2", gsub("\\.1", "_1", colnames(curr.resp))) -> colnames(curr.resp)
    curr.resp -> all.resp[[i]][[1]]
  }
  return(all.resp)
}

processPlate <- function(curr.layout, curr.plate, mode){
  lapply(1:length(curr.layout), function(x){
    # check which format of the layout is available
    resp <- NULL
    curr.drug <- names(curr.layout)[x]
    f.warn <- c()
    d.warn <- c()
    if (length(grep("rows", names(curr.layout[[x]])))>0){
      curr.layout[[x]]$rows -> rows
      curr.layout[[x]]$cols -> cols
      
      ## handle different row/col combinations
      resp <- matrix(NA, nrow=length(curr.layout[[x]]$doses),
                     ncol=(length(rows)+1))
      resp[,1] <- curr.layout[[x]]$doses
      for (i in 1:length(rows)){
        if (is.list(rows) && is.list(cols)){
          curr.plate[rows[[i]],cols[[i]]] -> resp[,(i+1)]
        } else if (is.numeric(rows) && is.numeric(cols)){
          curr.plate[rows[i], cols[i]] -> resp[,(i+1)]
        } else if (is.numeric(rows) && !is.numeric(cols)){
          as.numeric(curr.plate[rows[i], cols[[i]]]) -> resp[,(i+1)]
        } else if (!is.numeric(rows) && is.numeric(cols)){
          as.numeric(curr.plate[rows[i], cols[[i]]]) -> resp[,(i+1)]
        }
      }
      
      # check if resp elements are lists; if this is the case, convert
      # to numeric
      apply(resp, 2, function(x) as.numeric(as.character(x))) -> resp
      length(curr.layout[[x]]$rows) -> reps
      colnames(resp) <- c("Concentrations", sapply(1:reps, function(y) 
        paste(names(curr.layout)[x],y,sep="_")))
      
      if (length(1:reps) > 1){
        combn(reps, 2) -> all.combs
        apply(all.combs, 2, function(y)
          cor.test(resp[,y[1]+1], resp[,y[2]+1], method="s")$estimate) -> cors
        if (cors < 0.65){
          warning(paste("Replicates of ", names(curr.layout)[x], " in ",
                        f[1], " have cor < 0.65.", sep=""))
          print(f[1])
          c(f.warn, f[1]) -> f.warn
          c(d.warn, names(curr.layout)[x]) -> d.warn
        }
      }
      
      as.data.frame(resp) -> resp
      apply(resp, 2, function(x) as.numeric(as.character(x))) -> resp
      if (mode == "indirect"){
        (resp[,2:ncol(resp)]-pos)/(c.mean-pos) -> resp[,2:ncol(resp)]
      } else {
        resp[,2:ncol(resp)]/c.mean -> resp[,2:ncol(resp)]
      }
    } else if (length(grep("coords", names(curr.layout[[x]])))>0){
      as.numeric(curr.layout[[x]]$coords[1,]) -> rows
      as.numeric(curr.layout[[x]]$coords[2,]) -> cols
      cbind(curr.layout[[x]]$doses, sapply(1:length(rows), function(y){
        curr.plate[rows[y], cols[[y]]]
      })) -> resp
      apply(resp, 2, function(x) as.numeric(as.character(x))) -> resp
      if (mode == "indirect"){
        (resp[,2:ncol(resp)]-pos)/(c.mean-pos) -> resp[,2:ncol(resp)]
      } else {
        resp[,2:ncol(resp)]/c.mean -> resp[,2:ncol(resp)]
      }
      colnames(resp) <- c("Concentrations", rep(curr.drug, length(2:ncol(resp))))
    }
    list(resp, f.warn, d.warn) -> fin.resp
    names(fin.resp) <- c("resp", "f.warn", "d.warn")
    return(fin.resp)
  }) -> res
  return(res)
}

processControl <- function(control, curr.plate){
  lapply(1:length(control[[1]]$rows), function(x){
    curr.plate[as.numeric(unlist(control[[1]]$rows[x])), 
               as.numeric(unlist(control[[1]]$cols[x]))] -> c
    return(c)
  }) -> all.controls
 mean(unlist(all.controls)) -> control.mean
 
 # visualize control distibution
 list(unlist(all.controls), control.mean) -> res
 names(res) <- c("all.controls", "mean")
 return(res)
}

processPosControl <- function(control, curr.plate){
  lapply(1:length(control[[1]]$rows), function(x){
    curr.plate[as.numeric(unlist(control[[1]]$rows[x])), 
               as.numeric(unlist(control[[1]]$cols[x]))] -> c
    return(c)
  }) -> all.controls
  
  mean(unlist(lapply(all.controls, function(x)
    x[which(control[[1]]$doses %in% max(control[[1]]$doses))]))) -> res
  return(res)
}

findUniqueConc <- function(all.conc){
  lapply(all.conc, function(x) 
    lapply(x, function(y) paste(y, collapse=" "))) -> conc.col
  unique(unlist(conc.col)) -> conc.classes
  as.numeric(as.factor(conc.classes)) -> class
  
  lapply(conc.col, function(x)
    unlist(lapply(x, function(y) 
      class[which(conc.classes %in% y)]))) -> class.labels
  return(class.labels)
}

groupResponses <- function(all.resp, unique.conc){
  # groups responses based on concentrations
  unique(unlist(unique.conc)) -> all.conc
  all.df <- list()
  for (i in 1:length(all.conc)){
    lapply(1:length(unique.conc), function(x){
      which(unique.conc[[x]] %in% all.conc[i]) -> ind
      return(ind)
    }) -> m
    
    loc.df <- matrix(0, nrow=0, ncol=0)
    for (j in 1:length(m)){
      if (length(m[[j]]) > 0 & nrow(loc.df) == 0){
        # first entry
        as.data.frame(all.resp[[j]][m[[j]]]) -> loc.df
      } else if (length(m[[j]]) > 0 & nrow(loc.df) > 0){
        # append
        cbind(loc.df, as.data.frame(all.resp[[j]][m[[j]]])) -> loc.df
      }
    }
  
    if (ncol(loc.df) > 0 && length(grep("Concentrations", colnames(loc.df), ignore.case=T))>1){
      grep("Concentrations", colnames(loc.df), ignore.case=T)[2:length(grep("Concentrations", colnames(loc.df), ignore.case=T))] -> rm.ind
      rm.ind[which(!is.na(rm.ind))] -> rm.ind
      if (length(rm.ind) > 0){
        loc.df[,-rm.ind] -> all.df[[i]]
      }
    } else if (ncol(loc.df) > 0 && length(grep("Concentrations", colnames(loc.df), ignore.case=T))==1){
      loc.df -> all.df[[i]]
    }
  }
  
  return(all.df)
}

#@...
# Other routines
#@...
split_path <- function(path) {
  if (dirname(path) %in% c(".", path)) return(basename(path))
  return(c(basename(path), split_path(dirname(path))))
}