
#' Internal function for calculating the sBIC
#'
bic.func <- function(z, dat, uo, time, inflection.threshold, k) {
  
  if(z < (time[2]-time[1])) {uo.smooth <- uo
  } else {
    uo.smooth <- get.smooth.o(z, uo, time)
  }
  
  uo.smooth[is.na(uo.smooth)] <- 0
  
  likelihood <- get.likelihood(dat, uo.smooth, time)
  
  #changes in gradient
  gc.tmp <- diff(uo.smooth)
  
  #set below threshold to 0
  gc.tmp[(abs(gc.tmp)) < inflection.threshold*mean(abs(gc.tmp), na.rm = T) & (abs(gc.tmp)) != 0] <- NA
  gc.tmp <- zoo::na.locf(gc.tmp)
  
  #inflection points
  ip.tmp <- sign(gc.tmp)
  
  ip <- sum(diff(ip.tmp) != 0) + 1 #for slope
  
  #adds one if the curve doesn't start at 0
  if(uo.smooth[1] != 0) ip <- ip+1
  
  #loglik
  bic <- -2*(sum(log(likelihood), na.rm = T)) + k*(ip)
  
  return(bic)
}

#' Internal function for calculating the data likelihood given an assumed time-to-event density. Used by bic.func().
#'
get.likelihood <- function(dat, uo, time) {
  
  #set empty storage vectors
  likelihood <- vector()
  
  dat[is.na(dat[,1]),1] <- 0
  
  #get numerator and denominator contributions
  for (i in 1:nrow(dat)) {
    if(!is.na(dat[i,2])) {
      #likelihood
      likelihood[i] <- caTools::sumexact(uo[(match(dat[i,1],time)):(match(dat[i,2],time)-1)], na.rm = T)
    } else {
      likelihood[i] <- 1-(caTools::cumsumexact(uo)[match(dat[i,1], time)-1])
    }
  }
  return(likelihood)
}


#' Internal function for counting the number of turning points in a time-to-event curve. Used by bic.func().
#'
get.tp <- function(uo.smooth, inflection.threshold) {
  #changes in gradient
  gc.tmp <- diff(uo.smooth)
  #set below threshold to 0
  gc.tmp[(abs(gc.tmp)) < inflection.threshold*mean(abs(gc.tmp), na.rm = T) & (abs(gc.tmp)) != 0] <- NA
  gc.tmp <- na.locf(gc.tmp)
  #inflection points
  ip.tmp <- sign(gc.tmp)
  
  ip <- sum(diff(ip.tmp) != 0) + 1 #for slope
  
  #adds one if the curve doesn't start at 0
  if(uo.smooth[1] != 0) ip <- ip+1
  
  return(ip)
}

#' Internal function for converting a time-to-event distribution to a survival distribution.
#'
get.surv <- function(uo) {
  surv.u <- vector()
  surv.u[1] <- 1
  for(i in 2:(length(uo)+1)) {
    surv.u[i] <- surv.u[i-1] - uo[i-1]
  }
  return(surv.u)
}

#' Internal function for applying the Nadaraya–Watson kernel regression smoother to a time-to-event distribution.
#'
get.smooth.o <- function(z, turnbull.o, time) {
  uo <- turnbull.o
  
  if(z < (time[2]-time[1])) {
    uo.smooth <- turnbull.o
  } else {
    
    uo.smooth <- stats::ksmooth(x = time[-1], y = turnbull.o, kernel = "normal", bandwidth = z, n.points = length(turnbull.o))
    uo.smooth <- uo.smooth$y
  }
  
  if(sum(uo.smooth, na.rm = T) > 1) {
    uo.smooth <- uo.smooth/sum(uo.smooth, na.rm = T)
  }
  
  return(uo.smooth)
}