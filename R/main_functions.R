#' Nonparametric Estimation of a Smoothed Turnbull Survival Curve
#'
#' This function estimates survival curves (and time-to-event curves) from interval censored data
#' using the method of Turnbull (1976) and subsequently finds an optimal smoothing bandwidth which
#' minimizes the a penalized log-likelihood function (sBIC) as described in our manuscript.
#'
#' The function takes a matrix or data frame as input, where each row represents a subject.
#' The first column should be the left interval bounds, i.e. the last time which the subject was
#' observed to be event-free, with possible NA if a subject is left-censored. Similarly, the
#' second column are the right interval bounds, i.e. the first time which the subject was observed
#' to have experienced an event, with possible NA if a subject is right-censored.
#'
#' The output is a list containing the original and smoothed Turnbull survival and time-to-event
#' distributions among other sample and algorithm characteristics.
#'
#' @param dat A data.frame or matrix where rows are subjects and columns are left and right interval bounds.
#' @param n.obs The number of observations per subject. Used for calculation of effective N. Defaults to 2.
#' @param left.bound The earliest possible time which an event can occur. Defaults to 0.
#' @param penalty The penalty/penalties to use when calculating the sBIC. Possible values are "logNe",
#' "logNm", or "logN". Default is "logNe".
#' @param n.dec The number of decimal places in the observed data.
#' @param tolerance The tolerance for change in bandwidth when performing local optimization of the sBIC.
#' @param inflection.threshold Threshold used when counting the number of turning points in the time
#' to event density curve. Note that deviations from the default value have not been extensively tested.
#'
#' @export
smoothTB <- function(dat,
                     n.obs = 2,
                     left.bound = 0,
                     penalty = "logNe",
                     n.dec = 2,
                     tolerance = NA,
                     inflection.threshold = 1e-2) {
  
  #BIC constant
  N <- nrow(dat)
  
  #time intervals
  time.int <- (1/(10^(n.dec)))
  
  #vector of times under consideration
  time <- round(seq(from = left.bound, to = max(dat[,1:2], na.rm = T), by = time.int), n.dec)
  
  #round the data to n.dec
  dat[,1:2] <- round(dat[,1:2], n.dec)
  
  #apply TB
  icfit.mod <- interval::icfit(survival::Surv(dat[,1], dat[,2], type = "interval2") ~ 1)
  
  #get survival
  tbe.surv <- interval::getsurv(time, icfit.mod, nonUMLE.method = "interpolation")
  
  TB.s <- tbe.surv[[1]]$S
  
  #get onset
  TB.o <- diff(1-TB.s)
  TB.o[TB.o < 0] <- 0 #rounding error in diff
  
  #Effective Sample Size
  Ne.vec <- vector()
  censor.type <- vector(length = nrow(dat))
  censor.type[is.na(dat[,2])] <- "R"
  censor.type[is.na(dat[,1])] <- "L"
  censor.type[!is.na(dat[,1]) & !is.na(dat[,2])] <- "I"
  
  for(i in 1:nrow(dat)){
    if(censor.type[i] == "I") {
      Ne.vec[i] <- 1-(TB.s[match(dat[i,1], time)]-TB.s[match(dat[i,2], time)])
    }
    if(censor.type[i] == "L") {
      Ne.vec[i] <- 1-(1-TB.s[match(dat[i,2], time)])
    }
    if(censor.type[i] == "R") {
      Ne.vec[i] <- 1-(TB.s[match(dat[i,1], time)])
    }
  }
  
  Ne <- sum(Ne.vec, na.rm = T)
  
  #set optimizer parameters
  if(is.na(tolerance)){tolerance = time.int}
  #global optimization
  opts <- list("algorithm" = "NLOPT_GN_DIRECT_L",
               "xtol_abs" = tolerance*500,
               "maxeval" = 1000
  )
  
  #local optimization
  opts2 <- list("algorithm" = "NLOPT_LN_BOBYQA",
                "xtol_abs" = tolerance,
                "maxeval" = 5000
  )
  
  #optimize
  k.penal <- list("logN" = log(N),"logNm" = log(N*n.obs),"logNe" = log(Ne))
  bic.result <- list()
  
  for(i in 1:length(penalty)){
    
    kk <- as.numeric(k.penal[penalty][i])
    
    bic <- nloptr::nloptr(0, bic.func, lb = 0, ub = (length(TB.o)*time.int), opts = opts,
                          dat = dat, uo = TB.o, time = time, inflection.threshold = inflection.threshold,
                          k = kk)
    bic <- nloptr::nloptr(bic$solution, bic.func, lb = 0, ub = (length(TB.o)*time.int), opts = opts2,
                          dat = dat, uo = TB.o, time = time, inflection.threshold = inflection.threshold,
                          k = kk)
    
    bic.result[[i]] <- bic
  }
  
  #number of censoring types
  Nl <- sum(censor.type == "L")
  Ni <- sum(censor.type == "I")
  Nr <- sum(censor.type == "R")
  
  #write output
  output <- list()
  output[1:7] <- list(time,time[-1],N,Ne,Nl,Ni,Nr)
  names(output)[1:7] <- c("time.s","time.o","N","Ne","Nl","Ni","Nr")
  
  #number of iterations
  for(i in 1:length(bic.result)){
    output[[7+i]] <- bic.result[[i]]$iterations
    names(output)[7+i] <- paste0(penalty[i],"_iterations")
  }
  
  #bic output
  OL <- length(output)
  for(i in 1:length(bic.result)){
    output[[OL+i]] <- bic.result[[i]]$objective
    names(output)[OL+i] <- paste0(penalty[i],"_bic")
  }
  
  #smoothed curves
  OL <- length(output)
  output[(OL+1):(OL+2)] <- list(TB.o,TB.s)
  names(output)[(OL+1):(OL+2)] <- c("TB.o","TB.s")
  OL <- length(output)
  for(i in 1:length(bic.result)){
    smooth.tmp=get.smooth.o(bic.result[[i]]$solution, TB.o, time)
    output[[OL+2*i-1]] <- smooth.tmp
    names(output)[OL+2*i-1] <- paste0("TB.o.smooth.",penalty[i])
    output[[OL+2*i]] <- get.surv(smooth.tmp)
    names(output)[OL+2*i] <- paste0("TB.s.smooth.",penalty[i])
  }
  
  #chosen bandwidth
  OL <- length(output)
  for(i in 1:length(bic.result)){
    output[[OL+i]] <- bic.result[[i]]$solution
    names(output)[OL+i] <- paste0(penalty[i],"_bw")
  }
  
  output[[length(output)+1]] <- dat
  names(output)[length(output)] <- "dat"
  
  return(output)
}

#' Nonparametric Estimation of a Smoothed Kaplan-Meier Survival Curve
#'
#' This function estimates survival curves (and time-to-event curves) from interval censored data
#' using the method of Kaplan & Meier (1958) and subsequently finds an optimal smoothing bandwidth which
#' minimizes the a penalized log-likelihood function (sBIC) as described in our manuscript.
#'
#' The function takes a matrix or data frame as input, where each row represents a subject.
#' The first column should be either the time at event or the time at last follow up (if the
#' subject is right-censored). The second column is a binary variable indicating whether the subject
#' was observed to experience an event (1) or not (0).
#'
#' The output is a list containing the original and smoothed Kaplan-Meier survival and time-to-event
#' distributions among other sample and algorithm characteristics.
#'
#' @param dat A data.frame or matrix where rows are subjects and columns are the event/censoring time and event indicator.
#' @param n.obs The number of observations per subject. Used for calculation of effective N. Defaults to 2.
#' @param left.bound The earliest possible time which an event can occur. Defaults to 0.
#' @param penalty The penalty/penalties to use when calculating the sBIC. Possible values are "logNe",
#' "logNm", or "logN". Default is "logNe".
#' @param n.dec The number of decimal places in the observed data.
#' @param tolerance The tolerance for change in bandwidth when performing local optimization of the sBIC.
#' @param inflection.threshold Threshold used when counting the number of turning points in the time
#' to event density curve. Note that deviations from the default value have not been extensively tested.
#'
#' @export
smoothKM <- function(dat,
                     n.obs = 2,
                     left.bound = 0,
                     penalty = "logNe",
                     n.dec = 2,
                     tolerance = NA,
                     inflection.threshold = 1e-2) {
  
  #BIC constant
  N <- nrow(dat)
  
  n.dec <- n.dec
  time.int <- (1/(10^(n.dec)))
  
  #vector of times under consideration
  time <- round(seq(from = left.bound, to = max(dat[,1:2], na.rm = T), by = time.int), n.dec)
  
  #round the data to n.dec
  dat[,1] <- round(dat[,1], n.dec)
  
  #apply kaplan-meier
  survfit.mod <- survival::survfit(survival::Surv(time = dat[,1], event = dat[,2]) ~ 1)
  
  #get survival
  #survival::summary.survfit
  KM.s <- summary(survfit.mod,time = time)$surv
  KM.s[(length(KM.s)+1):length(time)] <- KM.s[length(KM.s)]
  
  #get onset
  KM.o <- diff(1-KM.s)
  KM.o[KM.o < 0] <- 0 #rounding error in diff
  
  #Effective Sample Size
  Ne.vec <- vector()
  censor.type <- vector(length = nrow(dat))
  censor.type[(dat[,2]) == 0] <- "R"
  censor.type[(dat[,2]) == 1] <- "E"
  
  for(i in 1:nrow(dat)){
    if(censor.type[i] == "E") {
      Ne.vec[i] <- 1
    }
    if(censor.type[i] == "R") {
      Ne.vec[i] <- 1-(KM.s[match(dat[i,1], time)])
    }
  }
  
  Ne <- sum(Ne.vec)
  
  dat.tmp <- rbind(dat[,1], dat[,1])
  dat.tmp <- rbind(dat.tmp, dat[,2])
  dat.tmp[dat.tmp[,3] == 0, 2] <- NA
  
  #set optimizer parameters
  if(is.na(tolerance)){tolerance = time.int}
  #global optimization
  opts <- list("algorithm" = "NLOPT_GN_DIRECT_L",
               "xtol_abs" = tolerance*500,
               "maxeval" = 1000
  )
  
  #local optimization
  opts2 <- list("algorithm" = "NLOPT_LN_BOBYQA",
                "xtol_abs" = tolerance,
                "maxeval" = 5000
  )
  
  #optimize
  k.penal <- list("logN" = log(N), "logNm" = log(N*n.obs), "logNe" =  log(Ne))
  bic.result <- list()
  
  for(i in 1:length(penalty)){
    
    kk <- as.numeric(k.penal[penalty][i])
    
    bic <- nloptr::nloptr(0, bic.func, lb = 0, ub = (length(KM.o)*time.int), opts = opts,
                          dat = dat.tmp, uo = KM.o, time = time, inflection.threshold = inflection.threshold,
                          k = kk)
    bic <- nloptr::nloptr(bic$solution, bic.func, lb = 0, ub = (length(KM.o)*time.int+time.int), opts = opts2,
                          dat = dat.tmp, uo = KM.o, time = time, inflection.threshold = inflection.threshold,
                          k = kk)
    
    bic.result[[i]] <- bic
  }
  
  #number of censoring types
  Ne <- sum(censor.type == "E")
  Nr <- sum(censor.type == "R")
  
  #write output
  output <- list()
  output[1:5] <- list(time,time[-1],N,Ne,Nr)
  names(output)[1:5] <- c("time.s","time.o","N","Ne","Nr")
  
  #number of iterations
  for(i in 1:length(bic.result)){
    output[[5+i]] <- bic.result[[i]]$iterations
    names(output)[5+i] <- paste0(penalty[i],"_iterations")
  }
  
  #bic output
  OL <- length(output)
  for(i in 1:length(bic.result)){
    output[[OL+i]] <- bic.result[[i]]$objective
    names(output)[OL+i] <- paste0(penalty[i],"_bic")
  }
  
  #smoothed curves
  OL <- length(output)
  output[(OL+1):(OL+2)] <- list(KM.o,KM.s)
  names(output)[(OL+1):(OL+2)] <- c("KM.o","KM.s")
  OL <- length(output)
  for(i in 1:length(bic.result)){
    smooth.tmp=get.smooth.o(bic.result[[i]]$solution, KM.o, time)
    output[[OL+2*i-1]] <- smooth.tmp
    names(output)[OL+2*i-1] <- paste0("KM.o.smooth.",penalty[i])
    output[[OL+2*i]] <- get.surv(smooth.tmp)
    names(output)[OL+2*i] <- paste0("KM.s.smooth.",penalty[i])
  }
  
  #chosen bandwidth
  OL <-length(output)
  for(i in 1:length(bic.result)){
    output[[OL+i]] <- bic.result[[i]]$solution
    names(output)[OL+i] <- paste0(penalty[i],"_bw")
  }
  
  output[[length(output)+1]] <- dat
  names(output)[length(output)] <- "dat"
  
  return(output)
}

#' Imputation of Event Time
#'
#' Given a data frame or matrix of subject event intervals, this function imputes an estimate
#' of the exact onset time as the expectation given an assumed time-to-event distribution.
#'
#' The output is a vector of predicted event times for the appropriate subjects.
#'
#' @param dat A data.frame or matrix where rows are subjects and columns are left and right interval bounds.
#' @param event.distribution A vector of event probabilities, e.g. the TB.o.smooth.logNe vector output from
#' running the smoothTB function.
#' @param time A vector of event times corresponding to the probabilities in the event.distribution variable,
#' e.g. the time.o vector output from running the smoothTB function.
#' @param type The type of observation(s) you would like to predict. Mutually exclusive options are "LI" (left-
#' and interval-censored observations), "I" (interval-censored observations), "LIR" (left-, interval-, and right-
#' censored observations). Default is "I". Note that depending on the underlying event mechanism, it may be
#' inappropriate to attempt imputation for left- or right-censored subjects.
#' @param n.dec The number of decimal places in the observed data.
#'
#' @export
impute.time <- function(dat, 
                          event.distribution, 
                          time, 
                          type = "I", 
                          n.dec = 2) {
  
  out.pred <- vector()
  
  for(i in 1:nrow(dat)) {
    
    dat.i <- as.numeric(unname(dat[i,]))
    
    if(type == "LI") {
      if(is.na(dat.i[2])) {
        pred <- NA
      } else {
        if(is.na(dat.i[1])) dat.i[1] <- 1
        times <- round(seq(dat.i[1],dat.i[2], by = (1/(10^(n.dec)))), n.dec)
        times.p <- (match(times, time))
        event.distribution[event.distribution == 0] <- .Machine$double.eps
        pred <- sum(times*event.distribution[times.p]/sum(event.distribution[times.p], na.rm = T), na.rm = T)
      }
    }
    
    if(type == "LIR") {
      if(is.na(dat.i[2])) {
        times <- round(seq(dat.i[1],max(time), by = (1/(10^(n.dec)))), n.dec)
        times.p <- (match(times, time))
        event.distribution[event.distribution == 0] <- .Machine$double.eps
        pred <- sum(times*event.distribution[times.p]/sum(event.distribution[times.p], na.rm = T), na.rm = T)
      } else {
        if(is.na(dat.i[1])) dat.i[1] <- 1
        times <- round(seq(dat.i[1],dat.i[2], by = (1/(10^(n.dec)))), n.dec)
        times.p <- (match(times, time))
        event.distribution[event.distribution == 0] <- .Machine$double.eps
        pred <- sum(times*event.distribution[times.p]/sum(event.distribution[times.p], na.rm = T), na.rm = T)
      }
    }
    
    if(type == "I") {
      if(is.na(dat.i[2]) | is.na(dat.i[1])) {
        pred <- NA
      } else {
        times <- round(seq(dat.i[1],dat.i[2], by = (1/(10^(n.dec)))), n.dec)
        times.p <- (match(times, time))
        event.distribution[event.distribution == 0] <- .Machine$double.eps
        pred <- sum(times*event.distribution[times.p]/sum(event.distribution[times.p], na.rm = T), na.rm = T)
      }
    }
    
    out.pred[i] <- pred
  }
  
  return(out.pred)
}