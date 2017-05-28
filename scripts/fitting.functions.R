library(drc)
library(pracma)

'%ni%' <- Negate('%in%')

fitData <- function(experiments, min=NA, max.x=NA){
  # --- update function based on current experiment structure
  # --- use normalized component for fitting,
  # --- unless one wants to obtain the ratio of live cells to a
  # --- second parameter
  # --- defaults for min and max based on ALL data
  res <- list()
  for (i in 1:length(experiments)){
    curr.exp <- experiments[[i]]
    if (is.data.frame(curr.exp) || is.matrix(curr.exp)){
      # --- step 1: remove all outliers in a data set
      #apply(curr.exp[,2:ncol(curr.exp)], 2, function(x) removeOutliers(x, curr.exp[,1])) -> curr.exp[,2:ncol(curr.exp)]
      if (is.na(min)){
        min <- min(curr.exp[,1])
      }
      
      if (is.na(max.x)){
        max.x <- max(curr.exp[,1])
      }
      
      fit(curr.exp, min, max.x) -> res[[i]]
      names(res)[[i]] <- names(experiments)[i]
      curr.exp -> experiments$normalized[[i]] # --- replace values in matrix with version without outliers
    } else if (is.list(curr.exp)){
      # --- res is a list of 5 lists, corresponding to different fitted
      # --- parameters
      print(paste("Processing ", names(experiments)[i], "...", sep=""))
      orig.min <- min; orig.max <- max.x
      logIC50 <- logEC10 <- logEC50 <- logEC90 <- inf.pt <- max <- AUC <- Error.50 <- ci.all <- c()
      all.res.list <- list()
      for (j in 1:length(curr.exp)){
        print(paste("Plate", j, sep=" "))
        curr.exp[[j]] -> df
        apply(as.data.frame(df[,2:ncol(df)]), 2, function(x) 
          removeOutliers(x, df[,1])) -> df[,2:ncol(df)]
        if (is.na(min)){
          min <- min(df[,1])
        }
        
        if (is.na(max.x)){
          max.x <- max(df[,1])
        }
        
        fit(df, min, max.x) -> all.res.list[[j]]
        df -> curr.exp[[j]]
        min <- orig.min; max.x <- orig.max #reset
      }
      curr.exp -> experiments[[i]]
      res.loc <- list(logIC50, logEC10, logEC50, logEC90, inf.pt, AUC, max, Error.50, ci.all)
      
      for (k in 1:length(all.res.list)){
        for (l in 1:length(res.loc)){
          c(res.loc[[l]], all.res.list[[k]][[l]]) -> res.loc[[l]]
        }  
      }
      
      sapply(experiments[[i]], function(x) colnames(x)[2:ncol(x)]) -> nn
      as.character(unlist(nn)) -> nn
      
      for (k in 1:length(res.loc)){
        res.loc[[k]] -> x
        names(x) <- nn
        x -> res.loc[[k]]
      }
      
      names(res.loc) <- c("logIC50", "logEC10", "logEC50", "logEC90", "inf.pt", "AUC", "max", "error.50", "ci.all")
      res[[i]] <- res.loc
      names(res)[i] <- names(experiments)[i]
    }
  }
  
  list(res, experiments) -> res.fin
  names(res.fin) <- c("res", "experiments")
  return(res.fin)
}

getResponseClass <- function(y.dat, x.dat, curr.exp=NA){
  # make curr.exp optional, double-check if variable is needed
  # at all
  fit.class <- NA
  if (length(which(y.dat == 0)) < length(x.dat) &&
        length(which(is.na(y.dat))) <= (length(x.dat)/2) && 
        length(which(!is.na(y.dat))) >= 4){
    if (length(which(is.na(y.dat))) > 0){
      x.dat[-which(is.na(y.dat))] -> x.dat
      y.dat[-which(is.na(y.dat))] -> y.dat
    }
    
    # for creating the pts matrix, ensure that y.dat and x.dat are in descending
    # order; smoothen y.dat already
    order(x.dat) -> ind
    x.dat[ind] -> x.dat
    y.dat[ind] -> y.dat
    getSlopes(y.dat) -> slopes.u
    
    y.dat.u <- y.dat # unsmoothed y.dat
    smooth(y.dat) -> y.dat
    getSlopes(y.dat) -> slopes
    slopes[c(floor((length(slopes)/2)):length(slopes))] -> last
    
    # lm(formula=y.dat~log10(x.dat), na.action = na.omit)$coefficients[2] < -0.07 &&
    #  mean(slopes) <= -0.1 
    if (mean(y.dat, na.rm=T) < 0.76 ||
        length(which(slopes.u < 0)) == length(slopes.u) ||
        length(which(last < 0)) == length(last)){
      # --- added case where all slopes are downwards;
      # --- also included a check on the slope of the last two points
      # --- class to be fitted
      fit.class <- 1
    } else if (length(which(is.na(y.dat))) >= length(y.dat)/2 ||
                 length(which(!is.na(y.dat))) < 4){
      # all outliers, skip DMSO!
      fit.class <- 2 #max
    } else {
      if (lm(formula=y.dat~log10(x.dat), na.action = na.omit)$coefficients[2] > 0.06 ||
            (lm(formula=y.dat~log10(x.dat), na.action = na.omit)$coefficients[1] <= 2.0 &&
               lm(formula=y.dat~log10(x.dat), na.action = na.omit)$coefficients[1] >= 0.75 &&
               mean(y.dat, na.rm=T) >= 0.7 && round(lm(formula=y.dat~log10(x.dat), na.action = na.omit)$coefficients[2], 2) > -0.2)){
        fit.class <- 2 # another max
      } else if (lm(formula=y.dat~log10(x.dat), na.action=na.omit)$coefficients[2] <= 0.1 &&
                  lm(formula=y.dat~log10(x.dat), na.action=na.omit)$coefficients[2] >= -0.1 &&
                  lm(formula=y.dat~log10(x.dat), na.action=na.omit)$coefficients[1] <= 0.55 &&
                  lm(formula=y.dat.u~log10(x.dat), na.action=na.omit)$coefficients[2] <= 0.1 &&
                  lm(formula=y.dat.u~log10(x.dat), na.action=na.omit)$coefficients[2] >= -0.1 &&
                  lm(formula=y.dat.u~log10(x.dat), na.action=na.omit)$coefficients[1] <= 0.55) {
        # --- check if even the linear fit with log10(x) is flat; condition has to be met for both y.dat.u, y.dat
        fit.class <- 3 # min
      } else if (lm(formula=y.dat~log10(x.dat), na.action=na.omit)$coefficients[2] <= 0.1 &&
                   lm(formula=y.dat~log10(x.dat), na.action=na.omit)$coefficients[2] >= -0.2 &&
                   (lm(formula=y.dat~log10(x.dat), na.action=na.omit)$coefficients[1] > 0.55 ||
                      mean(y.dat, na.rm=T) > 0.8)){
        fit.class <- 2 # yet another max
      } else {
        fit.class <- 1
      }
    }
  }
  return(fit.class)
}

getMatchIndex <- function(result.name, all.names){
  which(all.names %in% result.name) -> ind
  return(ind)
}

formatResList <- function(results, df, response.class, min.x, max.x){
  res.loc <- list()
  logIC50 <- rep(NA, (ncol(df)-1)); logEC10 <- rep(NA, (ncol(df)-1))
  logEC50 <- rep(NA, (ncol(df)-1)); logEC90 <- rep(NA, (ncol(df)-1))
  inf.pt <- rep(NA, (ncol(df)-1)); max <- rep(NA, (ncol(df)-1))
  AUC <- rep(NA, (ncol(df)-1)); Error.50 <- rep(NA, (ncol(df)-1))
  ci.all <- rep(NA, (ncol(df)-1))
  res.loc <- list(logIC50, logEC10, logEC50, logEC90, inf.pt, AUC, max, Error.50, ci.all)
  names(res.loc) <- c("logIC50", "logEC10", "logEC50", "logEC90", "inf.pt", "AUC", "max", "error.50", "ci.all")
  
  if (length(results) > 0){
    colnames(df)[2:ncol(df)] -> all.names
    unlist(sapply(names(results), function(x) getMatchIndex(x,all.names))) -> ind
    for (i in 1:length(ind)){
      which(names(results) %in% all.names[ind[i]]) -> list.ind
      results[[list.ind]] -> curr.result
      res.loc$logIC50[ind[i]] <- as.numeric(curr.result$logIC50)
      res.loc$logEC10[ind[i]] <- as.numeric(curr.result$logEC10); res.loc$logEC50[ind[i]] <- as.numeric(curr.result$logEC50)
      res.loc$logEC90[ind[i]] <- as.numeric(curr.result$logEC90); res.loc$inf.pt[ind[i]] <- as.numeric(curr.result$inf.pt)
      res.loc$AUC[ind[i]] <- as.numeric(curr.result$AUC); res.loc$max[ind[i]] <- curr.result$max #as.numeric(curr.result$max)
      res.loc$error.50[ind[i]] <- as.numeric(curr.result$error.50); res.loc$ci.all[ind[i]] <- as.numeric(curr.result$ci.all)
    }
  }
  
  # --- class 2
  #max(df[,1])-min(df[,1]) -> range; range/99 -> int  
  #seq(min(df[,1]), max(df[,1]), by=int) -> interval
  
  #max.x-min.x -> range; range/99 -> int  
  which(response.class %in% 2) -> ind
  exp(seq(log(min.x), log(max.x), length = 100)) -> interval
  
  # --- project IC50 based on current slope *if* max(df[,1]) < max.x
  # --- check if the last two points have a negative slope
  if (max(df[,1]) < max.x){
    ind+1 -> df.ind
    sapply(df.ind, function(x) getSlope(nrow(df)-1, df[,x], df[,1])) -> slopes
    df.ind[which(slopes <= -0.5)] -> ind.2
    if (length(ind.2) > 0){
      ind[-which(ind %in% (ind.2-1))] -> ind
    }
  }
  
  # --- put placeholders for probably real flat liners
  res.loc$logIC50[ind] <- log10(max.x); res.loc$logEC10[ind] <- log10(max.x); res.loc$logEC50[ind] <- log10(max.x)
  res.loc$logEC90[ind] <- log10(max.x); res.loc$inf.pt[ind] <- log10(max.x);
  res.loc$error.50[ind] <- log10(max.x); res.loc$ci.all[ind] <- NA;
  res.loc$max[ind] <- paste(c(paste(interval, collapse=","), paste(rep(1, length(interval)), collapse=",")), collapse=";")
  
  # --- calculate AUC from linear fit
  if (length(ind) >= 1){
    ind+1 -> df.ind
    if (length(ind) > 1){
      apply(df[,df.ind], 2, function(x) lm(x~df[,1],na.action=na.omit)) -> lm.fits
    } else {
      lm.fits <- list(lm(df[,df.ind]~df[,1],na.action=na.omit))
    }
    
    lapply(lm.fits, function(x) return(x$coefficients[2]*interval+x$coefficients[1])) -> res
    lapply(res, function(x) abs(trapz(interval,x))) -> auc.calc
    res.loc$AUC[ind] <- unlist(auc.calc)
  }
  
  # --- then process those in ind.2 differently: fit with LL.4(), extrapolate
  # --- until max.x *if* max(df[,1]) < max.x
  if (max(df[,1]) < max.x){
    if (length(ind.2) > 0){
      if (length(ind.2) > 1){
        new.x <- c(df[,1], max.x)
        apply(df[,ind.2], 2, function(x) return(c(x, 0))) -> df.p
        apply(df.p, 2, function(x) addFit(new.x, x, max.x)) -> curr.fit
      } else {
        new.x <- c(df[,1], max.x)
        new.y <- c(df[,ind.2], 0)
        list(addFit(new.x, new.y, max.x)) -> curr.fit
      }
      
      ind.2-1 -> ind.2 # revert to list index, rather than column index
      for (i in 1:length(curr.fit)){
        res.loc$logIC50[ind.2[i]] <- curr.fit[[i]]$logIC50; res.loc$logEC10[ind.2[i]] <- curr.fit[[i]]$logIC50 
        res.loc$logEC50[ind.2[i]] <- curr.fit[[i]]$logEC50; res.loc$logEC90[ind.2[i]] <- curr.fit[[i]]$logEC90 
        res.loc$inf.pt[ind.2[i]] <- curr.fit[[i]]$inf.pt; res.loc$error.50[ind.2[i]] <- curr.fit[[i]]$error.50
        res.loc$ci.all[ind.2[i]] <- curr.fit[[i]]$ci.all; res.loc$max[ind.2[i]] <- curr.fit[[i]]$max
        res.loc$AUC[ind.2[i]] <- curr.fit[[i]]$AUC
      }
    }
  }
  
  which(response.class %in% 3) -> ind
  res.loc$logIC50[ind] <- log10(min.x); res.loc$logEC10[ind] <- log10(min.x); res.loc$logEC50[ind] <- log10(min.x)
  res.loc$logEC90[ind] <- log10(min.x); res.loc$inf.pt[ind] <- log10(min.x);
  res.loc$error.50[ind] <- log10(min.x); res.loc$ci.all[ind] <- paste(c(paste(interval, collapse=","), paste(rep(0, length(interval)), collapse=",")), collapse=";")
  
  if (length(ind) >= 1){
    ind+1 -> df.ind
    # --- get all points between minimum and maximum df[,df.ind]
    lm.fits <- list()
    auc.calc <- list()
    if (length(ind) > 1){
      for (i in 1:length(df.ind)){
        which(df[,df.ind[i]] %in% max(df[,df.ind[i]], na.rm=TRUE)) -> start
        which(df[,df.ind[i]] %in% min(df[,df.ind[i]], na.rm=TRUE)) -> end
        lm(df[c(start:end),df.ind[i]]~df[c(start:end),1],na.action=na.omit) -> lm.fits[[i]]
        # get min max of interval and fit from there
        min(df[c(start:end),1]) -> new.min
        max(df[c(start:end),1]) -> new.max
        new.max-new.min -> new.range; new.range/99 -> new.int 
        seq(new.min, new.max, by=new.int) -> new.interval
        lapply(lm.fits, function(x) return(x$coefficients[2]*new.interval+x$coefficients[1])) -> res
        lapply(res, function(x) abs(trapz(new.interval,x))) -> auc.calc
      }
    } else {
      which(df[,df.ind] %in% max(df[,df.ind], na.rm=TRUE)) -> start
      which(df[,df.ind] %in% min(df[,df.ind], na.rm=TRUE)) -> end
      lm.fits <- lm(df[c(start:end),df.ind]~df[c(start:end),1],na.action=na.omit)
      min(df[c(start:end),1]) -> new.min
      max(df[c(start:end),1]) -> new.max
      new.max-new.min -> new.range; new.range/99 -> new.int 
      seq(new.min, new.max, by=new.int) -> new.interval
      lm.fits$coefficients[2]*new.interval+lm.fits$coefficients[1] -> res
      abs(trapz(new.interval,res)) -> auc.calc
    }
    res.loc$AUC[ind] <- unlist(auc.calc)
    res.loc$max[ind] <- paste(c(paste(interval, collapse=","), paste(rep(0, length(interval)), collapse=",")), collapse=";")
  }
  
  return(res.loc)
}

fit <- function(df, min.x, max.x){
  # --- restructured in the following form:
  # --- removal of all outliers that could cause non-convergent fits
  # --- retention of switches, i.e. if statements to identify
  # --- non-convergent cases that are not due to outliers
  apply(as.data.frame(df[,2:ncol(df)]), 2, function(x) 
    getResponseClass(x, df[,1], df)) -> response.class
  names(response.class) <- colnames(df)[c(2:ncol(df))] # --- to be sure that the name is associated with the var
  
  # --- override fit class for empty cells, DMSO
  if (length(grep("Empty", colnames(df)) > 0)){
    response.class[grep("Empty", colnames(df))-1] <- NA
  } else if (length(grep("DMSO", colnames(df)) > 0)){
    response.class[grep("Empty", colnames(df))-1] <- 2
  }
  
  # --- then apply fitting only to those with a class response of 1
  which(response.class %in% 1)+1 -> ind
  results <- list()
  if (length(ind) > 0){ # addition in the event of non-reactives across a plate (ext.data)
    apply(as.data.frame(df[,ind]), 2, function(x) addFit(df[,1], x, max.x)) -> results
    names(results) <- colnames(df)[ind]
  }
  formatResList(results, df, response.class, min.x, max.x) -> res.loc
  return(res.loc)
}

addFit <- function(x.dat, y.dat, max.x){
  require(drc)
  require(pracma)
  # --- if max(x.dat) < max.x, estimate y.dat at max.x
  # --- with a linear fit
  exp.model <- NA
  if (max(x.dat) < max.x){
    # get linear fits of last three points
    getSlope(length(y.dat)-1, y.dat, x.dat) -> slope
    if (slope >= -0.1){
      # extrapolate max.x as last point of y.dat
      # need to check direction (format of cgp is different from JP's data; easiest
      # is to flip the data frame on the horizontal axis, but we could also do the
      # extrapolation here)
      if (all(x.dat == cummin(x.dat))){
        # decreasing, add to start
        c(max.x, x.dat) -> x.dat
        if (all(y.dat[which(!is.na(y.dat))] == 
                  cummin(y.dat[which(!is.na(y.dat))]))){
          # y also decreasing, add to start
          c(y.dat[length(y.dat)], y.dat) -> y.dat
        } else {
          # add to end
          c(y.dat, y.dat[length(y.dat)]) -> y.dat
        }
      } else if (all(x.dat == cummax(x.dat))){
      #   # increasing, add to end: for garnett data only!!!
      #   if (lm(y.dat.new[c((length(y.dat.new)-2):length(y.dat.new))]~
      #            x.dat.new[c((length(y.dat.new)-2):length(y.dat.new))], 
      #          na.action = na.omit)$coeff[2] < 0){
      #     # -- if slope of last three points is positive, extrapolate
      #     c(y.dat,y.dat[length(y.dat)]) -> y.dat
      #     c(x.dat, max.x) -> x.dat        
      #   }
      }
      exp.model <- drm(y.dat~x.dat, fct = LL.4(), na.action = na.omit, control=drmc(errorm = F))
    } else if (slope < -0.1){
      # estimate y at max.x
      lm(y.dat~x.dat, na.action = na.omit) -> fit
      
      fit$coefficients[2]*max.x+fit$coefficients[1] -> max.y
      if (max.y < 0){ 
        max.y <- 0
      }
      c(x.dat, max.x) -> x.dat
      c(y.dat, max.y) -> y.dat
      exp.model <- drm(y.dat~x.dat, fct = LL.4(), na.action = na.omit, control=drmc(errorm = F))
    }
  } else {
    exp.model <- drm(y.dat~x.dat, fct = LL.4(), na.action = na.omit, control=drmc(errorm = F))
    # control=drmc(errorm = F, rmNA=T, relTol = 1e-10)
  }
  
  if (exp.model$fit$convergence){
    exp.model$fit$par -> params
    params[1] -> b # -- Hill coeff
    params[2] -> c # -- min
    params[3] -> d # -- max
    params[4] -> e # -- x_1/2
    
    if (d < 0.5 || c == 0.5){
      # d < 0.5 can occur with any steep drop, whether it's at the start or the
      # end of the response curve; trick would be to find where the drop occurs
      logIC50 <- log10(e)
    } else if (b < 0 || c > 0.5) {
      logIC50 <- log10(max(x.dat))
    } else {
      logIC50 <- log10(e*(((d-c)/(abs(0.5-c)))-1)^(1/b))
    }
    
    logEC10 <- log10(ED(exp.model,10)[[1]])
    logEC50 <- log10(ED(exp.model,50)[[1]])
    logEC90 <- log10(ED(exp.model,90)[[1]])
    inf.pt <- log10(ED(exp.model,50)[[1]])+(1/b)
    Error.50 <- log10(ED(exp.model,50)[[2]])
    ci.all <- mean(c(abs(b-confint(exp.model, "b")),
                     abs(c-confint(exp.model, "c")),
                     abs(d-confint(exp.model, "d")),
                     abs(confint(exp.model, "e"))))
    # --- get points of fit for curve fitting
    dataList <- exp.model[["dataList"]]
    dose <- dataList[["dose"]]
    #xLimits <- c(min(dose), max(dose)) # use full range of 5/8 pt data for fit approx
    x.fit <- exp(seq(log(min(x.dat)), log(max.x), length = 100))
    y.fit <- (exp.model$"curve")[[1]](x.fit)
    
    # --- get max from y.fit at x.fit = x
    max <- paste(c(paste(x.fit, collapse=","), paste(y.fit, collapse=",")), collapse = ";")
    AUC <- abs(trapz(x.fit,y.fit[,1]))
  }
  
  res.loc <- list(logIC50, logEC10, logEC50, logEC90, inf.pt, AUC, max, Error.50, ci.all)
  names(res.loc) <- c("logIC50", "logEC10", "logEC50", "logEC90", "inf.pt", "AUC", "max", "error.50", "ci.all")
  return(res.loc)
}

respFromFunc <- function(x.dat, y.dat, new.x){
  # creates y values given fit with x.dat, y.dat
  new.y <- NA
  exp.model <- drm(y.dat~x.dat, fct = LL.4(), na.action = na.omit, control=drmc(errorm = F))
  if (exp.model$fit$convergence){
    # form in drc:
    # f(x) = c + (d − c)/(1 + exp(b(log(x) − log(e))))
    
    exp.model$fit$par -> params
    params[1] -> b # -- Hill coeff
    params[2] -> c # -- min
    params[3] -> d # -- max
    params[4] -> e # -- x_1/2
    new.y <- (c + (d-c))/(1+exp(b*(log(new.x)-log(e))))
  }
  return(new.y)
}

getSlope <- function(ind, y.dat, x.dat){
  # -- check 
  if (length(which(is.na(y.dat[c(ind:(ind+1))])))== 0){
    slope <- lm(y.dat[c(ind:(ind+1))]~log10(x.dat[c(ind:(ind+1))]))$coefficients[2]
  } else {
    slope <- 1
  }
  return(slope)
}

removeOutliers <- function(y.dat, x.dat, limit=0.35){  
  # --- first, remove values post-normalization that are significantly above 1
  if (length(which(y.dat > 1.25)) > 0){
    y.dat[intersect(which(y.dat > 1.25), which(y.dat < 2.0))] <- 1.0
    if (length(which(y.dat > 2.0)) > 0){
      y.dat[which(y.dat > 2.0)] <- NA # real outlier
    }
  }
  
  # --- handle some patterns
  # --- go through all points and check if the removal of ANY point increases the negative slope
  # --- of a drug response curve
  c(1:(length(y.dat)-1)) -> ind
  sapply(ind, function(x) getSlope(x, y.dat, x.dat)) -> slopes
  y.dat[which(slopes >= limit)] <- NA # for external data: CGP
  if (slopes[length(slopes)] > limit){ #check slope between last two points
    y.dat[length(y.dat)] <- NA
  }
  
  if (length(which(is.na(y.dat))) > 0 && 
        length(which(!is.na(y.dat))) >= 4){
    y.dat -> y.dat.old
    x.dat -> x.dat.old
    x.dat[-which(is.na(y.dat))] -> x.dat
    y.dat[-which(is.na(y.dat))] -> y.dat
    c(1:(length(y.dat)-1)) -> ind
    sapply(ind, function(x) getSlope(x, y.dat, x.dat)) -> slopes
    y.dat[which(slopes >= 0.5)] <- NA
    if (length(which(y.dat.old %ni% y.dat)) > 0){
      y.dat.old[which(y.dat.old %ni% y.dat)] <- NA
    }
    y.dat.old -> y.dat
  }
  return(y.dat)
}

averageRes <- function(exp.res){
  exp.res.ave <- list()
  for (k in 1:length(exp.res)){
    exp.res[[k]] -> curr.exp
    
    # using the first feature, ic50, determine how many replicates there are;
    # if the number of replicates are not even, i.e. some have two, while some have four
    unique(sapply(strsplit(names(curr.exp$logIC50), "_"), function(x) x[1])) -> nn
    sapply(strsplit(names(curr.exp$logIC50), "_"), function(x) x[1]) -> nno
    logIC50 <- logEC90 <- logEC50 <- logEC10 <- inf.pt <- AUC <- fit <- max <- c()
    res.loc <- list(logIC50, logEC10, logEC50, logEC90, inf.pt, AUC, fit, max)
    names(res.loc) <- c("logIC50", "logEC10", "logEC50", "logEC90", "inf.pt", "AUC", "fit", "max")
    vars <- names(curr.exp)
    for (i in 1:length(nn)){
      which(nno %in% nn[i]) -> ind
        for (j in 1:length(curr.exp)){
          if (vars[j] %in% c("logIC50", "logEC10", "logEC50",
                          "logEC90", "inf.pt", "AUC")){
            c(res.loc[[which(names(res.loc) %in% vars[j])]],
              mean(as.numeric(curr.exp[[j]][ind]), na.rm=T)) -> 
              res.loc[[which(names(res.loc) %in% vars[j])]]
          } else if (vars[j] %in% "max") {
            # nrow is constant, given that fits return 100 points
            y.fits <- matrix(NA, nrow=100, ncol=length(ind))
            x.fits <- matrix(NA, nrow=100, ncol=length(ind))
            for (l in 1:length(ind)){
              if (!is.na(curr.exp$max[ind[l]])){
                as.numeric(strsplit(strsplit(curr.exp$max[ind[l]], ";")[[1]], ",")[[1]]) -> x.fits[,l]
                as.numeric(strsplit(strsplit(curr.exp$max[ind[l]], ";")[[1]], ",")[[2]]) -> y.fits[,l]
              }
            }
            
            # then check if all columns are populated by NAs
            apply(y.fits, 2, function(x) length(which(is.na(x)))) -> check
            if (length(which(check %in% 100)) < ncol(y.fits)){
              apply(y.fits, 1, function(x) mean(x, na.rm=T)) -> mean.y
              apply(x.fits, 1, function(x) mean(x, na.rm=T)) -> mean.x
              c(res.loc[[which(names(res.loc) %in% vars[j])]], paste(c(paste(mean.x, collapse=","), 
                  paste(mean.y, collapse=",")), collapse = ";")) -> res.loc[[which(names(res.loc) %in% vars[j])]]
            } else {
              c(res.loc[[which(names(res.loc) %in% vars[j])]], NA) -> res.loc[[which(names(res.loc) %in% vars[j])]]
            }
          }
        }
    }
    
    # make all elements of res.loc a named a vector
    for (i in 1:length(res.loc)){
      if (!is.null(res.loc[[i]])){
        names(res.loc[[i]]) <- nn
      }
    }
    
    exp.res.ave[[k]] <- res.loc
  }
  names(exp.res.ave) <- names(exp.res)
  return(exp.res.ave)
}

toMatrix <- function(exp.res.ave, drug.list){
  # --- current form assumes that all the drugs are represented on all the plates
  # --- create patient x drug matrix for each of the parameters in exp.res.ave
  logIC50 <- logEC90 <- logEC50 <- logEC10 <- matrix(data = NA, nrow = length(exp.res.ave), ncol = length(drug.list))
  inf.pt <- AUC <- Emax <- matrix(data = NA, nrow = length(exp.res.ave), ncol = length(drug.list))
  #Error.50 <- matrix(data = NA, nrow = length(exp.res.ave), ncol = length(drug.list))
  
  for (i in 1:length(exp.res.ave)){
    # --- get appropriate entry from exp.res.ave
    for (j in 1:length(drug.list)){
      logIC50[i,j] <- exp.res.ave[[i]]$logIC50[j]
      logEC90[i,j] <- exp.res.ave[[i]]$logEC90[j]
      logEC50[i,j] <- exp.res.ave[[i]]$logEC50[j]
      logEC10[i,j] <- exp.res.ave[[i]]$logEC10[j]
      inf.pt[i,j] <- exp.res.ave[[i]]$inf.pt[j]
      AUC[i,j] <- exp.res.ave[[i]]$AUC[j]
      Emax[i,j] <- exp.res.ave[[i]]$Emax[j]
      #Error.50[i,j] <- exp.res.ave[[i]]$error.50[j]
    }
  }
  
  res <- list(logIC50, logEC10, logEC50, logEC90, inf.pt, AUC, Emax)
  names(res) <- c("logIC50", "logEC10", "logEC50", "logEC90", "inf.pt", "AUC", "Emax")
  for (i in 1:length(res)){
    colnames(res[[i]]) <- drug.list
    rownames(res[[i]]) <- names(exp.res.ave)
  }
  
  return(res)
}

toMatrix.asymmetric <- function(experiments, res, drug.list, mode="md"){
  # --- generates an object similar to exp.res.mat;
  # --- params: experiments object, res (from fitData) and drug.list(vector)
  # --- first, extract drug list from each experiment
  if (mode %in% "md"){
    lapply(experiments$normalized, function(x) unlist(lapply(x, function(y) 
      return(colnames(y)[-which(colnames(y) %in% c("Concentrations", "doses"))])))) -> drug.names
  } else if (mode %in% "sd") {
    # single dose for all plates
    lapply(experiments$normalized, function(x) 
      return(colnames(x)[-which(colnames(x) %in% c("Concentrations", "doses"))])) -> drug.names
  }
  
  # --- for each parameter, create a data frame with all the drugs and 
  # --- entries for all the experiments
  params <- names(res[[1]][[1]])
  matrices <- list()
  for (i in 1:length(params)){
    curr.mat <- matrix(NA, nrow=length(experiments$normalized), ncol=length(drug.list))
    colnames(curr.mat) <- drug.list
    rownames(curr.mat) <- names(res$res)
    # --- then for each result, populate curr.mat
    lapply(drug.names, function(x) which(drug.list %in% x)) -> col.ind
    lapply(res$res, function(x) return(x[[i]])) -> content
    #mapply(function(x,y) return(x[[i]][y]), res$res, col.ind) -> content
    for (j in 1:length(content)){
      curr.mat[j,col.ind[[j]]] <- content[[j]]
    }
    matrices[[i]] <- curr.mat
  }
  names(matrices) <- params
  return(matrices)
}

trimMat <- function(exp.res.mat, perc){
  # --- removes drugs for which all patients do not react
  all.cols <- c()
  for (i in 1:length(exp.res.mat)){
    exp.res.mat[[i]] -> curr.mat
    # --- for each drug, check if the no. of non-reactive cases are greater than
    # --- the specified percentage
    for (j in 1:ncol(curr.mat)){
      round(curr.mat[,j], 0) -> estimate
      if (length(which(estimate == 4)) >= perc*nrow(curr.mat)){
        all.cols <- c(all.cols, j)
      }
    }
  }
  unique(all.cols) -> all.cols
  
  for (i in 1:length(exp.res.mat)){
    exp.res.mat[[i]] -> curr.mat
    curr.mat[,-all.cols] -> curr.mat
    
    # --- also remove duplicated columns
    curr.mat[,-which(duplicated(colnames(curr.mat)))] -> curr.mat
    curr.mat -> exp.res.mat[[i]]
  }
  
  return(exp.res.mat)
}

getErrors <- function(exp.res){
  all.errors <- c()
  for (i in 1:length(exp.res)){
    c(all.errors, exp.res[[i]]$error.50) -> all.errors
  }
  # --- remove NAs
  if (length(is.na(all.errors)) > 0){
    all.errors[-which(is.na(all.errors))] -> all.errors
  }
  return(all.errors)
}

getCIs <- function(exp.res){
  all.cis <- c()
  for (i in 1:length(exp.res)){
    c(all.cis, exp.res[[i]]$ci.all) -> all.cis
  }
  # --- remove NAs
  if (length(is.na(all.cis)) > 0){
    all.cis[-which(is.na(all.cis))] -> all.cis
  }
  
  abs(log10(all.cis)) -> all.cis
  return(all.cis)
}