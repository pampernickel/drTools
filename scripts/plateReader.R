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
            as.numeric(plate[doses$dose.rows[1],]) -> t # used for finding the columns/rows of dups
            if (replicates == 1 && dups == "ei"){ #explicit, irregular
              all.res[[j]] <- doses
            } else if (replicates==2 && dups[x] == "i") {
              all.res[[j]] <- handleImplicit(dup.mode[x], t, doses, dilution)
            } else if (replicates==2 && dups[x] == "e"){
              # replicates are explicitly indicated
              all.res[[j]] <- handleExplicit(dup.mode[x], t, doses, dilution)
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
  print(paste("Checking directory ", dir, "...", sep=""))
  c(list.files(path = dir, pattern = ".txt", recursive=TRUE, full.names = TRUE), 
    list.files(path = dir, pattern = ".csv", recursive=TRUE, full.names = TRUE)) -> files
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
    plate[dose.rows,] -> dr
    if (length(dose.rows) > 1){
      as.vector(dr) -> dr
      as.numeric(dr[which(dr %ni% "")]) -> dr
      
      # check intervals: detect and restrict to serial dilutions
      dr[rev(order(dr))] -> dr
      table(round((dr[1:(length(dr)-1)]/dr[2:length(dr)]), 3)) -> int.counts
      dr[which(round((dr[1:(length(dr)-1)]/dr[2:length(dr)]), 3) %in%
              names(int.counts)[which(int.counts %in% max(int.counts))])] -> doses
    } else {
      as.vector(dr) -> dr
      as.numeric(dr[which(dr %ni% "")]) -> doses
    }
    
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

handleExplicit <- function(mode, t, doses, dilution){
  res <- list()
  if (length(grep("r", mode)) > 0){
    sapply(doses$doses, function(y)
      which(t %in% y)) -> all.inds
    list(doses$doses/dilution, as.numeric(doses$dose.rows), 
         list(all.inds, all.inds)) -> res
    names(res) <- c("doses", "rows", "cols")
  } else if (length(grep("c", mode)) > 0){
    # cols processing
    sapply(doses$doses, function(y)
      which(t %in% y)) -> all.inds
    t(all.inds) -> cols
    all.inds <- list()
    for (i in 1:ncol(cols)){
      cols[,i] -> all.inds[[i]]
    }
    list(doses$doses/dilution, 
         c(as.numeric(doses$dose.rows),
          as.numeric(doses$dose.rows)),
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
readExperiment <- function(files, layout){  
  lapply(files, function(f){
    lapply(f, function(y){
      if (length(grep(".csv", y)) == 1){
        data <- read.csv(as.character(y), header = FALSE)
      } else if (length(grep(".csv", y)) == 1) {
        data <- read.table(as.character(y), header = FALSE, sep = "\t")
      }
    }) -> all.dat
    
    # start with layout
    all.resp <- list()
    for (i in 1:length(layout)){
      t(all.dat[[i]]) -> curr.plate      
      layout[[i]] -> curr.layout
      # check column batches; here, check if there are a bunch
      # of drugs that have the same columns
      ## treat controls separately
      curr.layout[grep("control", tolower(names(curr.layout)))] -> control
      processControl(control, curr.plate) -> c.mean # dmso mean for plate
      
      # get control, and return as a single value (i.e. mean of all
      # control values + sd)
      curr.layout[-grep("control", tolower(names(curr.layout)))] -> curr.layout
      names(curr.layout) -> drug.names
      lapply(1:length(curr.layout), function(x){
        # check which format of the layout is available
        resp <- NULL
        curr.drug <- names(curr.layout)[x]
        if (length(grep("rows", names(curr.layout[[x]])))>0){
          curr.layout[[x]]$rows -> rows
          curr.layout[[x]]$cols -> cols
          cbind(curr.layout[[x]]$doses, sapply(1:length(rows), function(y){
            curr.plate[rows[y], cols[[y]]]
          })) -> resp
          
          length(curr.layout[[x]]$rows) -> reps
          colnames(resp) <- c("Concentrations", sapply(1:reps, function(y) paste(names(curr.layout)[x],y,sep="_")))
          
          if (length(1:reps) > 1){
            combn(reps, 2) -> all.combs
            apply(all.combs, 2, function(y)
              cor.test(resp[,y[1]+1], resp[,y[2]+1], method="s")$estimate) -> cors
            if (cors < 0.65){
              warning(paste("Replicates of ", names(curr.layout)[x], " in ",
                          f[1], " have cor < 0.65.", sep=""))
            }
          }
          
          as.data.frame(resp) -> resp
          apply(resp, 2, function(x) as.numeric(as.character(x))) -> resp
          resp[,2:ncol(resp)]/c.mean -> resp[,2:ncol(resp)]
        } else if (length(grep("coords", names(curr.layout[[x]])))>0){
          as.numeric(curr.layout[[x]]$coords[1,]) -> rows
          as.numeric(curr.layout[[x]]$coords[2,]) -> cols
          cbind(curr.layout[[x]]$doses, sapply(1:length(rows), function(y){
            curr.plate[rows[y], cols[[y]]]
          })) -> resp
          resp[,2:ncol(resp)]/c.mean -> resp[,2:ncol(resp)]
          colnames(resp) <- c("Concentrations", rep(curr.drug, length(2:ncol(resp))))
        }
        return(resp)
      }) -> all.resp[[i]]
    }
    
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
  
    return(all.resp.fin)
  }) -> all.resp
  return(all.resp)
}

processControl <- function(control, curr.plate){
  lapply(1:length(control[[1]]$rows), function(x){
    curr.plate[as.numeric(unlist(control[[1]]$rows[x])), 
               as.numeric(unlist(control[[1]]$cols[x]))] -> c
    return(c)
  }) -> all.controls
 mean(unlist(all.controls)) -> control.mean
 return(control.mean)
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