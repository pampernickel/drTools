# PDobay for AG Bourquin, KISPI

#@...
# Routines for reading and processing combination experiments
#@...

library(reshape2)

readCombos <- function(dir, res.dir, mode=c("", "normalized"), dil.factor=1,
                       singleLayout=F){
  # singleLayout: applicable to the XML case, where individual plates are detected
  # this option 
  print(paste("Processing directory ", dir, "...", sep=""))
  unlist(strsplit(dir, "/")) -> dir.name
  
  c(list.files(path = dir, pattern = ".txt", recursive=TRUE, full.names = TRUE), 
    list.files(path = dir, pattern = ".csv", recursive=TRUE, full.names = TRUE),
    list.files(path = dir, pattern = ".xml", recursive=TRUE, full.names = TRUE)) -> files  
  c(list.files(path = res.dir, pattern = ".txt", recursive=TRUE, full.names = TRUE), 
    list.files(path = res.dir, pattern = ".csv", recursive=TRUE, full.names = TRUE)) -> res.files
  
  if (length(files) == 0){
    stop(paste("No files found in ", dir, ".", sep=""))
  }
  
  if (length(res.files) == 0){
    stop(paste("No files found in ", res.dir, ".", sep=""))
  }
  
  combo.mats <- NA
  if (length(grep(".xml", files)) > 0){
    # pass to an xml reader internally
    readXML(files) -> coords
    readFileXML(coords, res.files, dil.factor, singleLayout) -> combo.mats
  } else if (length(grep(".txt", files)) > 0 && mode=="normalized"){ # 
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
          start <- start.1 <- d.coords.s <- 1
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
            
            # split dmso.coords as well
            nrow(dmso.coords)/2 -> d.coords.e
            dmso.coords[d.coords.s:(d.coords.s+d.coords.e),] -> dmso.coords.sub
            d.coords.s <- d.coords.e
            mean(apply(dmso.coords.sub, 1, function(x)
              content[x[1],x[2]])) -> dmso.mean
            mat[nrow(mat),2] <- dmso.mean
            print(dmso.mean)
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
        paste(names(combo.mats), sapply(combo.mats, function(x) names(x)), sep=":") -> nn
        unlist(combo.mats, recursive=F) -> combo.mats
        names(combo.mats) <- nn
      }
    } else {
      # multiple combos per plate
      # infoLines, meta, contents
      sapply(infoLines[[1]], function(x)
        sapply(strsplit(contents[[1]][x], "\t"), function(y) y[1])) -> drug.names
      lapply(1:length(drug.names), function(x) getPlate(infoLines, contents, 1, x)) -> plates
      
      # check areas of overlaps between various plates; plate that overlaps with all other plates
      # is the drug that's used in combination with everything else
      sapply(plates, function(x) length(which(as.vector(x) %ni% ""))) -> cc
      which(cc %in% max(cc)) -> mp # master plate
      plates[[mp]] -> mpl; plates[-mp] -> plates; drug.names[mp] -> md; drug.names[-mp] -> od
      
      # get coords for mp
      which(mpl != "", arr.ind = T) -> p1.ind
      
      # then for each of the n other plates, check areas of overlap bet p1.ind and the nth
      # plate, and create all.combos from these
      lapply(res.files, function(y){
        read.csv(y, header=F) -> content
        t(content) -> content
        lapply(1:length(plates), function(x){
          plates[[x]] -> cp
          which(cp != "", arr.ind = T) -> cp.ind
          apply(p1.ind, 1, function(x) 
            ifelse(x[1] %in% cp.ind[,1] && x[2] %in% cp.ind[,2], T, F)) -> o.1
          apply(cp.ind, 1, function(x) 
            ifelse(x[1] %in% p1.ind[,1] && x[2] %in% p1.ind[,2], T, F)) -> o.2
          cp.ind[which(o.2 %in% T),] -> o.coords
          
          # now find where the content are d1.only and d2.only; in the case of d1
          # find the borders of o.2
          cp.ind[which(o.2 %in% F),] -> d2.only
          
          p1.ind[which(o.1 %in% F),] -> d1.only
          cp.ind[which(o.2 %in% T)[1],] -> tl
          cp.ind[which(o.2 %in% T)[length(which(o.2 %in% T))],] -> br
          c(tl[1], br[2]) -> tr
          
          # then subtract one from both tl, tr
          cbind(rep(c(tl[1]-1), length(tl[2]:tr[2])), 
                tl[2]:tr[2]) -> temp
          apply(d1.only, 1, function(x) 
            ifelse(x[1] %in% temp[,1] && x[2] %in% temp[,2], T, F)) -> d1.c
          d1.only[which(d1.c %in% T),] -> d1.only
          
          as.numeric(as.character(mpl[d1.only])) -> d1.doses
          as.numeric(as.character(cp[d2.only])) -> d2.doses
          mat <- matrix(NA, nrow=length(d1.doses)+2, ncol=length(d2.doses)+2)
          mat[1,1] <- "XX"
          mat[2:nrow(mat),1] <- c(d1.doses,0)
          mat[1,2:ncol(mat)] <- c(0,d2.doses)
          for (k in 1:nrow(o.coords)){
            mat[which(mat[,1] %in% mpl[o.coords[k,1],o.coords[k,2]]),
                which(mat[1,] %in% cp[o.coords[k,1],o.coords[k,2]])] <-
              content[o.coords[k,1],o.coords[k,2]]
          }
          
          # fill content for single drug treatment
          mat[2:(nrow(mat)-1),2] <- as.numeric(content[d1.only])
          mat[nrow(mat), 3:ncol(mat)] <- as.numeric(content[d2.only])
          
          # then get the dmso-containing well
          tl-1 -> dmso
          content[dmso[1],dmso[2]] -> mat[which(mat[,1] %in% 0),which(mat[1,] %in% 0)]
          apply(mat[2:nrow(mat),2:ncol(mat)],2, function(y) as.numeric(as.character(y))) -> mat_n
          100*mat_n/as.numeric(content[dmso[1],dmso[2]]) -> mat_n
          rownames(mat_n) <- mat[2:nrow(mat),1]
          colnames(mat_n) <- mat[1, 2:ncol(mat)]
          return(mat_n)
        }) -> res
        paste(md, od, sep="_") -> names(res)  
        return(res)
      }) -> combo.mats
      names(combo.mats) <- sapply(strsplit(sapply(strsplit(res.files, "/"), 
                                           function(x) x[length(x)]), "_"), function(y) y[1])
    }
  }
  return(combo.mats)
}

readFileXML <- function(coords, res.files, dil.factor, singleLayout){
  # create function to read a results file specified in res.files into
  # all.combos format; the readFileXML function is based on the default
  # where there is a one-to-one xml:result folder
  layout <- NA
  all.combos <- list()
  if (singleLayout == T){
    # case where only one of the layouts out of the coords file will be
    # taken; otherwise, there would be one layout per plate, and there should
    # be a concordance between the plate name and the layout (i.e. it's possible to match
    # the plate names to the layout)
    if (length(coords) == 1){
      coords[[1]][[1]] -> layout
    } else {
      lapply(coords, function(x) x[[1]]) -> layout
    }
   
    # go through res.files
    lapply(res.files, function(x){
      t(read.csv(x, header=F)) -> f  
    }) -> plates
    names(plates) <- res.files
    gsub("_", "-", gsub(".csv", "", sapply(strsplit(names(plates), "\\/"), function(x) x[length(x)]))) -> names(plates)
    for (j in 1:length(plates)){
      plates[[j]] -> curr.plate
      combos <- list()
      for (i in 1:length(layout[[1]][[1]])){
        layout[[1]][[1]][[i]] -> curr.layout
        curr.plate[curr.layout$combo_coords$r,curr.layout$combo_coords$c] -> p
       
        # check if the plate is not a blank
        if (length(unique(as.vector(p))) > 1){
          # check orientation of drugs
          if (length(unique(curr.layout$d1c$c)) == 1){
            # d1 is on a single column; add 0,0
            if (isDescending(curr.layout$doses[[1]])){
              rownames(p) <- c(curr.layout$doses[[1]], 0)
            } else if (!isDescending(curr.layout$doses[[1]])){
              rownames(p) <- c(0, curr.layout$doses[[1]])
            }
            
            if (isDescending(curr.layout$doses[[2]])){
              colnames(p) <- c(curr.layout$doses[[2]], 0)
            } else if (!isDescending(curr.layout$doses[[2]])){
              colnames(p) <- c(0, curr.layout$doses[[2]])
            }
          }
          
          # normalize against dmsos; also check if some dmsos have values = 0, which
          # could be the case when dmso is dispensed, but there are not enough cells
          # on the plate
          p[which(rownames(p) %in% 0),which(colnames(p) %in% 0)] -> zz # zero zero
          curr.plate[unique(curr.layout$dmso_coords$r),unique(curr.layout$dmso_coords$c)] -> dmsos
          if (length(which(dmsos %in% 0)) > 0){
            dmsos[-which(dmsos %in% 0)] -> dmsos
            warning("Removed potentially blank DMSO wells.")
          }
          
          mean(dmsos, na.rm=T) -> dmso.mean
          if (zz/dmso.mean >= 1.25){
            dmso.mean <- zz
          }
          p/dmso.mean -> combos[[i]]
        } else {
          # blank plate
          NA -> combos[[i]]
        }
      }
      names(combos) <- paste(curr.layout$drug1, curr.layout$drug2, 1:length(layout[[1]][[1]]), sep="_")
      combos -> all.combos[[j]]
    }
    names(all.combos) <- names(plates)
    
    # get all combo names
    unlist(all.combos, recursive = F) -> all.combos
    all.combos[which(sapply(all.combos, function(x) ifelse(is.matrix(x), T, F)) %in% T)] -> all.combos
  } else {
    # check if length of layouts match length of res.files
    if (length(coords) == 1){
      if (length(coords[[1]]) != length(res.files)){
        stop("Number of results files not equal to number of plates detected in the .xml layout file.")
      } else if (length(coords[[1]]) == length(res.files)) {
        coords[[1]] -> layout
      }
    }
  }
  return(all.combos)
}

readFile <- function(coords, res.files, dil.factor){
  t(read.csv(res.files, header=F)) -> f
  # create mat based on d2.list, coords[[i]]$d1.mat[,i]
  df <- list()
  for (i in 1:length(coords)){
    df.sub <- matrix(0, nrow=(length(coords[[i]]$rows)),
                     ncol=(length(coords[[i]]$columns[[i]][[1]])))
    for (j in 1:length(coords[[i]]$columns[[1]])){
      f[coords[[i]]$rows,coords[[i]]$columns[[i]][[j]]] -> df.sub
      f[coords[[i]]$dmso.row,coords[[i]]$columns[[i]][[j]]] -> dmsos
      
      # append rownames as doses and columns
      if (isDescending(coords[[i]]$d1.mat[,i])){
       c(coords[[i]]$d1.mat[,i], 0) -> rn 
      } else {
        c(0, coords[[i]]$d1.mat[,i]) -> rn 
      }
      
      if (isDescending(coords[[i]]$d2.list[[j]])){
        c(coords[[i]]$d2.list[[j]], 0) -> cn 
      } else {
        c(0, coords[[i]]$d2.list[[j]]) -> cn 
      }
      
      rn/dil.factor -> rownames(df.sub)
      cn/dil.factor -> colnames(df.sub)
      df.sub/mean(dmsos) -> df.sub
      df.sub -> df[[j]]
    }
  }
  names(df) <- paste(res.files, c(1:length(df)), sep="_")
  return(df)
}

calcCI <- function(combos){
  if (is.list(combos[[1]])){
    unlist(combos, recursive = FALSE) -> cl
  } else {
    combos -> cl
  }
  
  df <- matrix(0, nrow=0, ncol=4)
  colnames(df) <- c("patient", "drug1", "drug2", "CI")
  for (i in 7:length(cl)){
    print(i)
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

checkCombos <- function(all.combos){
  # go through combos and check if some might have fitting issues
  # specifically, check for cases where contents for single drugs are missing
  checkMissingDrugs(all.combos) -> cl
}

checkMissingDrugs <- function(all.combos){
  which(sapply(all.combos, function(x){
    ifelse(length(unique(x[,which(colnames(x) %in% 0)])) == 1 &&
             unique(x[,which(colnames(x) %in% 0)]) == 0, T, F)
    }) %in% T) -> rm.ind1
  which(sapply(all.combos, function(x){
    ifelse(length(unique(as.vector(x[which(rownames(x) %in% 0),]))) == 1 &&
             unique(as.vector(x[which(rownames(x) %in% 0),])) == 0, T, F)
  }) %in% T) -> rm.ind2
  setdiff(1:length(all.combos), unique(c(rm.ind1, rm.ind2))) -> ind
  all.combos[ind] -> all.combos
  return(all.combos)
}

getComboProperties <- function(cl,i){
  # based on brice's combos
  unlist(strsplit(sapply(strsplit(names(cl)[i], "\\:"), function(x) x[length(x)]), "_"))[1] -> d1
  unlist(strsplit(sapply(strsplit(names(cl)[i], "\\:"), function(x) x[length(x)]), "_"))[2] -> d2
  unlist(strsplit(sapply(strsplit(names(cl)[i], "\\:"), function(x) x[1]), "/")) -> d
  d[length(d)] -> pn
  gsub("_cumul_1_v_1_2", "", pn) -> pn
  gsub("-cumul-1-v-1-2", "", pn) -> pn
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