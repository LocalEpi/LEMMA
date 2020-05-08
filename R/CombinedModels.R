#' @import ggplot2
#' @import data.table
#' @import matrixStats

#TODO: RunSim needs some error checking that params (or maybe other data.tables) have the right names - I had a bug at one point because I was passing intervention1_multiplier instead of intervention1.multiplier and it didn't throw an error, it just ignored intervention1.multiplier

GetBetaMatrix <- function(beta, N, k) {
  stopifnot(length(beta) == 2, length(N) == 2, length(k) == 2)
  matrix(c(k[1], 1 - k[1], 1 - k[2], k[2]), nrow = 2, ncol = 2) * 
    matrix(sum(N) / N * beta, nrow = 2, ncol = 2, byrow = T)
}

#initial.new.exposures (num.param.sets x num.init.exp) matrix or NULL - sets new.exposures at start.date
#population num.pops x 1 - assumes this is S, all others 0
#start.date scalar Date
#end.date scalar Date
#params data.table (num.param.sets x num.params)
SEIR <- function(initial.new.exposures, population, start.date, end.date, params) {
  stopifnot(is.null(params$use.hosp.rate)) #use.hosp.rate removed, assumed TRUE

  base.compartment.names <- c("S", "E", "IH", "IR", "R", "HP", "DC")

  num.pops <- length(population)
  stopifnot(num.pops == 2)
  initial.conditions <- sapply(base.compartment.names, function (z) rep(0, num.pops), simplify = F) # q: the base compartments
  initial.conditions$S <- population #length = num.pops
  
  N <- sum(population)
  
  #support mpfr
  zero <- 0 * N
  na <- NA_real_ * N
  
  num.param.sets <- nrow(params) 
  if (is.null(initial.new.exposures)) {
    num.init.exp <- 1
  } else {
    stopifnot(nrow(initial.new.exposures) == num.param.sets)
    num.init.exp <- ncol(initial.new.exposures)
  }
  
  p <- params #for readability
  
  sigma <- 1 / p$latent.period
  
  gamma.r <- 1 / p$illness.length.given.nonhosp
  gamma.h <- 1 / (p$exposed.to.hospital - p$latent.period) #infected.to.hospital = exposed.to.hospital - latent.period
  stopifnot(p$exposed.to.hospital > p$latent.period) #otherwise gamma.h is nonsense
  
  has.E <- p$latent.period > 0 
  
  stopifnot("prop.hospitalized.1" %in% names(p), "prop.hospitalized.2" %in% names(p), !("prop.hospitalized" %in% names(p)))
  prop.hospitalized <- cbind(p$prop.hospitalized.1, p$prop.hospitalized.2)
  # prop.hospitalized <- c(prop.hospitalized) #FIXME is this right?? keeping as vector avoids non-conformable arrays error in loop, but need to check that the recycling works (if so may as well be here instead of in loop?)
  
  stopifnot("hosp.length.of.stay.1" %in% names(p), "hosp.length.of.stay.2" %in% names(p), !("hosp.length.of.stay" %in% names(p)))
  hosp.length.of.stay <- cbind(p$hosp.length.of.stay.1, p$hosp.length.of.stay.2)
  psi <- 1 / hosp.length.of.stay
  
  
  dates <- seq(start.date, end.date, by = "day")
  num.days <- length(dates)
  stopifnot(num.days > 1) #single day may cause problems
  
  overall.gamma <- gamma.r  #TODO: figure out what this should be (weighted average of gamma.r and gamma.h?)
  beta <- GetBeta(dates, params, overall.gamma) #num.days x num.param.sets
  
  # S = susceptible
  # E = exposed (but not yet infectious)
  # IR = infectious but not severe (will not be hospitalized)
  # IH = infectious severe (will be hospitalized, although currently not in hospital)
  # HP = hospital population
  # DC = discharged (or died) from hospital
  # R = recovered (did not go to hospital)
  
  empty.compartment <- array(na, dim = c(num.days, num.param.sets, num.pops, num.init.exp), dimnames = list(as.character(dates), NULL, NULL, NULL))
  q <- sapply(base.compartment.names, function (z) empty.compartment, simplify = F) # q: the base compartments
  d <- sapply(base.compartment.names, function (z) NULL) #d: change in base compartments 
  
  total.infected <- new.exposures <- empty.compartment
  
  #date.index <- GetDateIndex(initial.conditions$date, dates) 
  #stopifnot(identical(which(date.index), 1L), initial.conditions$date == start.date) #needs work if a range of dates (need to at least change the starting tt)
  date.index <- 1 #FIXME
  
  
  for (i in base.compartment.names) { #TODO: improve speed here? this is taking a decent amout of time
    for (k in 1:num.pops) {
      q[[i]][date.index, , k, ] <- initial.conditions[[i]][k] 
    }
  }
  #N <- sum(sapply(base.compartment.names, function (z) initial.conditions[[z]][1])) #TODO: clean this up  #FIXME: should be num.pop x 1
  
  d$E <- new.infections <- new.admits <- new.dischanges <- array(na, dim = c(num.param.sets, num.pops, num.init.exp))
  # run SEIR model
  for (tt in 1:num.days) { #most of this loop is 1:(num.days - 1), but we need total.infected also on num.days
    total.infected[tt, , , ] <- q$IR[tt, , , ] + q$IH[tt, , , ] + (p$patients.in.hosp.are.infectious) * q$HP[tt, , , ]
    stopifnot(!is.na(total.infected[tt, , , ])) #temp
    if (tt == num.days) break
    
    if (tt == 1 && !is.null(initial.new.exposures)) {
      #new.exposures[tt, , , ] <- 0
      #new.exposures[tt, , 1, ] <- initial.new.exposures #initial.new.exposures is only in pop1
      # new.exposures[tt, , 2, ] <- initial.new.exposures; warning("temp!")
      #FIXME: super slow, just for dev
      for (j in 1:num.param.sets) {
        for (k in 1:num.init.exp) {
          new.exposures[tt, j, , k] <- initial.new.exposures[j, k] * (population / N)
        }
      }
    } else {
      #FIXME: super slow, just for dev
      for (j in 1:num.param.sets) {
        for (k in 1:num.init.exp) {
          beta.mat <-  beta[tt, j, , ] #num.pops x num.pops
          S.N <- q$S[tt, j, , k] / N
          new.exposures[tt, j, , k] <- S.N * (beta.mat %*% total.infected[tt, j, , k])
        }
      }
      cat("new code: tt = ", tt, "\n")
      print(p[1])
      cat("new.exposures=\n")
      print(new.exposures)
      if (tt == 3) {
        print(beta[tt, 1, , ])
        print(q$S[tt, 1, , ])
        print(total.infected[tt, 1, , ])
        print(N)
      }
      if (F) {
        # warning("temp - truncation disabled")
      } else {
        new.exposures[tt, , , ] <- pmax(new.exposures[tt, , , ], 0)
        new.exposures[tt, , , ] <- pmin(new.exposures[tt, , , ], q$S[tt, , , ])
      }
      
    }
    new.nonhosp.recovered <- gamma.r * adrop(q$IR[tt, , , , drop = F], 1) 
    
    d$S <- adrop(-new.exposures[tt, , , , drop = F], 1)
    d$E[has.E, , ] <- new.exposures[tt, has.E, , ] - sigma[has.E] * q$E[tt, has.E, , ]
    d$E[!has.E, , ] <- 0
    new.infections[has.E, , ] <- sigma[has.E] * q$E[tt, has.E, , ]
    new.infections[!has.E, , ] <- new.exposures[tt, !has.E, , ]
    
    new.admits <- gamma.h * adrop(q$IH[tt, , , , drop = F], 1)
    new.dischanges <- c(psi) * adrop(q$HP[tt, , , , drop = F], 1) #TODO: check this
    
    d$IR <- new.infections * (1 - c(prop.hospitalized)) - c(new.nonhosp.recovered) #TODO: check this
    d$IH <- new.infections * c(prop.hospitalized) - c(new.admits) #TODO: check this
    
    d$HP <- new.admits - new.dischanges
    d$DC <- new.dischanges
    d$R <- new.nonhosp.recovered
    
    total.d <- array(zero, dim = c(num.param.sets, num.pops, num.init.exp)) #just for error checking, could remove
    for (i in base.compartment.names) {
      q[[i]][tt + 1, , , ] <- adrop(q[[i]][tt, , , , drop = F], 1) + d[[i]] 
      total.d <- total.d + c(d[[i]]) #just for error checking, could remove
    }
    
    
    stopifnot(abs(total.d) < 0.001) #d adds to zero in each compartment for each param.set and init.exp - just for error checking, could remove
    for (i in base.compartment.names) {
      stopifnot(q[[i]][tt + 1, ,, ] >= 0)
    }
    
  }
  
  if (class(new.exposures) == "mpfrArray") {
    q <- lapply(q, asNumeric) #if using mpfr
    total.infected <- asNumeric(total.infected)
    new.exposures <- asNumeric(new.exposures)
  }
  
  
  icu <- c(matrix(p$prop.icu, nrow = num.days, ncol = num.param.sets, byrow = T)) * q$HP  
  vent <-  c(matrix(p$prop.vent, nrow = num.days, ncol = num.param.sets, byrow = T)) * icu
  list.return <- c(q, list(I = total.infected, 
                           icu = icu, 
                           vent = vent,
                           R.total = q$R + q$DC, #R is recovered (did not go to hospital)
                           active.cases = q$E + q$IR + q$IH + q$HP, #active.cases includes those in hospital whether or not infectious 
                           total.cases = q$E + q$IR + q$IH + q$R + q$HP + q$DC)) #everyone ever exposed (implicitly includes dead) (everyone except S)
  names(list.return)[names(list.return) == "R"] <- "R.nonhosp"
  list.return$hosp <- adrop(list.return$HP[, , 1, , drop = F], 3)
  
  list.return$new.exposures <- new.exposures
  return(list.return) #list of num.days x num.param x num.init.exp  S E ...
}


#utility wrapper
Seir <- function(initial.new.exposures, initial.conditions, start.date, end.date, params) {
  if (!is.null(initial.new.exposures)) {
    if (is.matrix(initial.new.exposures)) {
      stopifnot(nrow(initial.new.exposures) == nrow(params))
    } else {
      stopifnot(length(initial.new.exposures) == nrow(params))
      initial.new.exposures <- as.matrix(initial.new.exposures)
    }
  }
  stopifnot(end.date > start.date)
  list.return <- SEIR(initial.new.exposures, initial.conditions, start.date, end.date, params)
  num.param.sets <- dim(list.return$hosp)[2]
  num.init.exp <- dim(list.return$hosp)[3]
  if (num.init.exp == 1) {
    list.return <- lapply(list.return, drop)
    if (num.param.sets == 1) {
      return(cbind(date = seq(start.date, end.date, by = "day"), as.data.table(list.return))) #data.table date S E ...
    } else {
      return(list.return) #list of num.days x num.params S E ...
    }
  } else {
    return(list.return) #list of num.days x num.param x num.init.exp  S E ...
  }
}

#initial.new.exposures (num.param.sets x num.init.exp) matrix
FitSEIR <- function(initial.new.exposures, total.population, start.date, observed.data, params) {
  num.param.sets <- nrow(params)
  num.init.exp <- ncol(initial.new.exposures)
  
  seir <- SEIR(initial.new.exposures, population = total.population, start.date, end.date = max(observed.data$date), params)
  
  # TODO: add weights if we're fitting data other than hospital (eg. ICU, total cases)
  # weights <- lapply(observed.data[, -"date"], function (obs.data.col) min(1, 1 / mean(obs.data.col)))
  weights <- lapply(observed.data[, -"date"], function (obs.data.col) 1)
  error <- CalcError(seir, observed.data, weights)
  
  #a large fixed penalty caused numeric precision problems - this relative penalty works but means badness can only be compared within one param and set of initial.new.exposures (don't compare across params or within params across calls to FitSEIR)
  penalty <- 1.01 * rowMaxs(error$sum.error.sq) #num.param.sets x 1
  
  #sum.error.sq and peak.before.first.observed: num.param.sets x num.init.exp
  badness <- error$sum.error.sq + error$peak.before.first.observed * penalty #num.param.sets x num.init.exp
  #TODO: figure out if ties.method = "last" always works - theoretically there should be no ties, but there can be because of precision - consider using Rmpfr package
  best.fit.index <- max.col(-badness, ties.method = "last")
  optimal.initial.new.exposures <- rowCollapse(error$optimal.initial.new.exposures.scale * initial.new.exposures, best.fit.index)
  
  best <- list(index = best.fit.index, initial.new.exposures = initial.new.exposures[cbind(seq_len(num.param.sets), best.fit.index)], optimal.initial.new.exposures = optimal.initial.new.exposures, badness = badness[cbind(seq_len(num.param.sets), best.fit.index)])
  
  hosp.at.last.obs.date <- matrix(seir$hosp[nrow(seir$hosp), , ], nrow = num.param.sets, ncol = num.init.exp) #fix if dropped to vector
  #TODO: rename allz
  allz <- list(initial.new.exposures=initial.new.exposures, badness=badness, hosp=hosp.at.last.obs.date, peak=error$peak.before.first.observed) #badness and peak are just for debugging
  return(list(best = best, allz = allz))
}


# all.dates %in% dates.to.find with error checking
GetDateIndex <- function(dates.to.find, all.dates) {
  is.sorted <- function(x) !is.unsorted(x, strictly = T)  #x is strictly sorted
  
  stopifnot(dates.to.find %in% all.dates, is.sorted(dates.to.find), is.sorted(all.dates))
  return(all.dates %in% dates.to.find)
}

#simulated.data list: num.days x num.param.sets x num.init.exp    with rownames = dates
#observed.data  data.table: date and subset of column names in simulated.data (e.g. S, I, IH, active.cases, etc)
#weights: list with same names as observed.data
CalcError <- function(simulated.data, observed.data, weights) {
  hosp <- simulated.data$hosp # num.obs.days x num.param.sets x num.init.exp
  num.param.sets <- dim(hosp)[2]
  num.init.exp <- dim(hosp)[3]
  
  sim.date <- as.Date(rownames(hosp)) #rownames are same for any element
  obs.date <- observed.data$date
  sim.data.index <- GetDateIndex(obs.date, sim.date)
  
  sum.error.sq <- sum.sim.sq <- sum.sum.obs <- matrix(0, nrow = num.param.sets, ncol = num.init.exp)
  for (data.name in setdiff(names(observed.data), "date")) {
    sim.data <- simulated.data[[data.name]][sim.data.index, , , drop = F] # num.obs.days x num.param.sets x num.init.exp
    obs.data <- observed.data[[data.name]] # num.obs.days x 1
    
    weight <- weights[[data.name]]
    sum.error.sq <- sum.error.sq + colSums(weight^2 * (sim.data - obs.data)^2, na.rm = T)
    sum.sim.sq <- sum.sim.sq + colSums(weight^2 * sim.data^2, na.rm = T)
    sum.sum.obs <- sum.sum.obs + colSums(weight^2 * (sim.data * obs.data), na.rm = T)
  }
  optimal.initial.new.exposures.scale <- sum.sum.obs / sum.sim.sq  #multiply by initial.new.exposures to get best guess for optimal initial.new.exposures
  
  #assumes obs.date > min(sim.date)
  peak.before.first.observed <- hosp[sim.date == obs.date[1] - 1, , ] > hosp[sim.date == obs.date[1], , ]
  dim(peak.before.first.observed) <- c(num.param.sets, num.init.exp) #in case a dim was dropped
  
  return(list(sum.error.sq = sum.error.sq, peak.before.first.observed = peak.before.first.observed, optimal.initial.new.exposures.scale = optimal.initial.new.exposures.scale))
}

GetBeta.noint <- function(dates, params, overall.gamma) {
  #num.days x num.param.sets x num.pops (in: S) x num.pops (by: I)
  num.days <- length(dates)
  num.param.sets <- nrow(params)
  num.pops <- 2 #FIXME
  prior.beta <- array(NA_real_, dim = c(num.days, num.param.sets, num.pops, num.pops))
  prior.beta[, , 1, 1] <- matrix(params$r0.initial.11 * overall.gamma, nrow = num.days, ncol = num.param.sets, byrow = T)
  #affects new infections in 1 caused by 2
  prior.beta[, , 1, 2] <- matrix(params$r0.initial.12 * overall.gamma, nrow = num.days, ncol = num.param.sets, byrow = T)
  prior.beta[, , 2, 1] <- matrix(params$r0.initial.21 * overall.gamma, nrow = num.days, ncol = num.param.sets, byrow = T)
  prior.beta[, , 2, 2] <- matrix(params$r0.initial.22 * overall.gamma, nrow = num.days, ncol = num.param.sets, byrow = T)
  
  #FIXME: doesn't do anything with intervention yet
  cum.multiplier <- 1
  beta <- prior.beta * cum.multiplier
  return(beta)
}

#for each parameter set, we take num.popsxnum.pops matrix of (list of initialR0, mult, date, smooth) 
#=> output beta: num.days x num.param.sets x num.pops x num.pops

GetBeta <- function(dates, params, overall.gamma) {
  ExpandToDays <- function(param) {
    matrix(param, nrow = num.days, ncol = num.param.sets, byrow = T)
  }
  GetSingleBeta <- function(pop.index1, pop.index2) {
    num.interventions <- length(grep("^intervention[[:digit:]]+\\.date$", names(params)))
    multiplier <- matrix(1, nrow = num.days, ncol = num.param.sets)
    for (intervention.num in seq_len(num.interventions)) {
      int.str <- paste0("intervention", intervention.num, ".")
      intervention.date <- params[[paste0(int.str, "date")]]
      intervention.smooth.days <- params[[paste0(int.str, "smooth.days")]]
      intervention.multiplier <- params[[paste0(int.str, "multiplier.", pop.index1, pop.index2)]]
      stopifnot(!is.null(intervention.smooth.days), !is.null(intervention.date), !is.null(intervention.multiplier))
      
      multiplier.mat <- ExpandToDays(intervention.multiplier ^ (1 / intervention.smooth.days))
      index <- date.mat >= ExpandToDays(intervention.date) & date.mat <= ExpandToDays(intervention.date + intervention.smooth.days - 1)
      multiplier[index] <- multiplier[index] * multiplier.mat[index]
    }
    cum.multiplier <- colCumprods(multiplier) #num.days x num.param.sets
    r0.initial <- params[[paste0("r0.initial.", pop.index1, pop.index2)]]
    stopifnot(!is.null(r0.initial))
    prior.beta <- matrix(r0.initial * overall.gamma, nrow = num.days, ncol = num.param.sets, byrow = T)
    return(prior.beta * cum.multiplier)
  }
  
  #scale beta by int_mult
  num.days <- length(dates)
  num.param.sets <- nrow(params)
  date.mat <- matrix(dates, nrow = num.days, ncol = num.param.sets)
  
  num.pops <- 2
  beta <- array(NA_real_, dim = c(num.days, num.param.sets, num.pops, num.pops))
  for (pop.index1 in 1:num.pops) {
    for (pop.index2 in 1:num.pops) {
      beta[, , pop.index1, pop.index2] <- GetSingleBeta(pop.index1, pop.index2)
    }
  }
  return(beta)
}


GetBeta.orig <- function(dates, params, overall.gamma) {
  #scale beta by int_mult
  num.days <- length(dates)
  num.param.sets <- nrow(params)
  multiplier <- matrix(1, nrow = num.days, ncol = num.param.sets)
  intervention.multiplier.str <- grep("^intervention[[:digit:]]+\\.multiplier$", names(params), value = T)
  #intervention1.smooth.days, intervention2.smooth.days, etc.
  #intervention1.date, intervention2.date, etc.
  #intervention1.multiplier, intervention2.multiplier, etc.
  date.mat <- matrix(dates, nrow = num.days, ncol = num.param.sets)
  ExpandToDays <- function(param) matrix(param, nrow = num.days, ncol = num.param.sets, byrow = T)
  for (int.mult.str in intervention.multiplier.str) {
    intervention.multiplier <- params[[int.mult.str]]
    intervention.smooth.days <- params[[sub("multiplier", "smooth.days", int.mult.str)]]
    intervention.date <- params[[sub("multiplier", "date", int.mult.str)]]
    stopifnot(!is.null(intervention.multiplier), !is.null(intervention.smooth.days), !is.null(intervention.date))
    
    multiplier.mat <- ExpandToDays(intervention.multiplier ^ (1 / intervention.smooth.days))
    index <- date.mat >= ExpandToDays(intervention.date) & date.mat <= ExpandToDays(intervention.date + intervention.smooth.days - 1)
    multiplier[index] <- multiplier[index] * multiplier.mat[index]
  }
  cum.multiplier <- colCumprods(multiplier) #num.days x num.param.sets
  
  prior.beta <- matrix(params$r0.initial * overall.gamma, nrow = num.days, ncol = num.param.sets, byrow = T)
  beta <- prior.beta * cum.multiplier
  return(beta)
}

#probably a better way to do this
SeqAround <- function(lo, hi, est, n) {
  s <- (1:n) ^ 10
  d1 <- (est - lo) / sum(s)
  d2 <- (hi - est) / sum(s)
  c(est - d1 * rev(s), est, est + d2 * s)
}

GetNextX <- function(fit, first.iter, expander, num.init.exp) {
  #TODO: vectorize? (I think the time saved vs complexity is not worth it)
  #TODO: make a smarter guess of x.min/x.max - inner bucket should be estimate of the width of initial.exposures needed to establish convergence, outer buckets should try to establish an interior solution
  
  #if expander is too small, optimal.initial.new.exposures won't fit in x.min to x.max during (for iter > 1); if expander is too big, won't detect convergence
  #if you see a lot of "expanded max" or "opt.outside.minmax" -> make expander bigger
  #if convergence is low on iter = 2 => try expander smaller or bigger or increase num.init.exp (this is harder to say because it can detect convergence between different indexes)
  
  #num.init.exp is number to use on all iterations except the initial fit (which uses one initial exposure)
  stopifnot(num.init.exp %% 2 == 1) #num.init.exp must be odd
  
  num.params <- length(fit$best$index) #this is params that have not yet converged
  x.set.next <- matrix(NA_real_, nrow = num.params, ncol = num.init.exp)
  hosp.span <- rep(NA, num.params)
  num.expand.min <- num.expand.max <- num.opt.outside.minmax <- 0
  
  for (j in 1:num.params) {
    if (first.iter) {
      x.set.prev <- fit$best$optimal.initial.new.exposures[j]  #TODO: make this less confusing to read - it's a bit hacky
      hosp.prev <- NULL #won't be used
    } else {
      x.set.prev <- fit$allz$initial.new.exposures[j, ]
      hosp.prev <- fit$allz$hosp[j, ]
    }
    
    index <- fit$best$index[j]
    if (index == 1) {
      if (!first.iter) num.expand.min <- num.expand.min + 1
      min.x <- x.set.prev[index] / expander
      min.hosp <- -Inf
    } else {
      min.x <- x.set.prev[index - 1]
      min.hosp <- hosp.prev[index - 1]
    }
    if (index == length(x.set.prev)) { #note: length(x.set.prev) is either 1 or num.init.exp
      if (!first.iter) num.expand.max <- num.expand.max + 1
      max.x <- x.set.prev[index] * expander
      max.hosp <- Inf
    } else {
      max.x <- x.set.prev[index + 1]
      max.hosp <- hosp.prev[index + 1]
    }
    
    hosp.span[j] <- max.hosp - min.hosp
    #converged if interior solution and diff between hosp [-1] and hosp [+1] < 0.5
    
    if (fit$best$optimal.initial.new.exposures[j] > min.x & fit$best$optimal.initial.new.exposures[j] < max.x) {
      x.set.next[j, ] <- SeqAround(min.x, max.x, fit$best$optimal.initial.new.exposures[j], num.init.exp %/% 2)
    } else {
      num.opt.outside.minmax <- num.opt.outside.minmax + 1
      x.set.next[j, ] <- seq(min.x, max.x, length.out = num.init.exp)
    }
  }
  if (num.expand.min > 0) cat("expanded min ", num.expand.min, "/", num.params, "\n")
  if (num.expand.max > 0) cat("expanded max ", num.expand.max, "/", num.params, "\n")
  if (num.opt.outside.minmax > 0) cat("opt.outside.minmax ", num.opt.outside.minmax, "/", num.params, "\n")
  
  
  stopifnot(!anyNA(x.set.next))
  converged.out <- abs(hosp.span) < 0.5
  return(list(x.set.next = x.set.next[!converged.out, , drop = F], converged = converged.out, hosp.span = hosp.span)) #note: x.set.next is among non-converged, others are not
}


RunSim <- function(total.population, observed.data, start.date, end.date, params, search.args) {
  ff <- function(xmin, xmax) {
    #temp for debugging - plot initial.new.exposure vs badness for single param
    stopifnot(num.param.sets == 1)
    temp.x.set <- t(seq(xmin, xmax, length.out = 1000))
    temp.fit <- FitSEIR(temp.x.set, total.population, start.date, observed.data, params)
    plot.dt <- rbind(data.table(x = temp.x.set[1, ], badness = temp.fit$allz$badness[1, ], type = "extra"), data.table(x = fit$allz$initial.new.exposures[1, ], badness = fit$allz$badness[1, ], type = "x.set"))
    plot.dt <- plot.dt[x > xmin & x < xmax]
    print(ggplot(plot.dt, aes(x, badness)) + geom_point(aes(color=factor(type))) + geom_vline(xintercept = fit$best$optimal.initial.new.exposures) + ggtitle(paste0("iter = ", iter)) + coord_cartesian(ylim=plot.dt[, range(badness)]) + scale_y_log10())
  }
  
  num.param.sets <- nrow(params)
  best.dt <- data.table(converged = rep(F, num.param.sets))
  fit <- FitSEIR(matrix(1e-30, nrow = num.param.sets, ncol = 1), total.population, start.date, observed.data, params)
  for (iter in 1:search.args$max.iter) {
    cat("---------- iter = ", iter, " -------------\n")
    
    #get next x.set (num.not.converged x num.init.exp) -- fit only has num.not.converged.
    next.list <- GetNextX(fit, first.iter = iter == 1, expander = search.args$expander, num.init.exp = search.args$num.init.exp)
    
    best.dt[converged == F, initial.new.exposures := fit$best$initial.new.exposures]
    best.dt[converged == F, objective := fit$best$badness]
    best.dt[converged == F, hosp.span := next.list$hosp.span]
    best.dt[converged == F, converged := next.list$converged]
    
    if (F) {
      print(fit$best)
      print(next.list$x.set.next)
      print(best.dt)
    }
    print(best.dt[, .(.N, mean.initial.new.exposures = mean(initial.new.exposures), mean.hosp.span = mean(hosp.span), max.hosp.span = max(hosp.span), mean.objective = mean(objective)), keyby=converged])
    if (all(best.dt$converged)) break
    
    fit <- FitSEIR(next.list$x.set.next, total.population, start.date, observed.data, params[!best.dt$converged])
  }
  if (!all(best.dt$converged)) {
    if (mean(!best.dt$converged) > search.args$max.nonconverge) {
      cat("did not converge:\n")
      print(which(!best.dt$converged))
      stop("failed to converge")
    }
  }
  
  #TODO: return seir at best.dt$index and use as initial condition for final projection in RunSim? (saves a little time but adds complexity)
  seir <- Seir(best.dt$initial.new.exposures, initial.conditions = total.population, start.date, end.date, params)
  return(seir)
}

#just for debugging
RunSim.slow <- function(total.population, observed.data, start.date, end.date, params) {
  f <- function(x, p, check.peak) {
    seir <- SEIR(matrix(x, 1, 1), initial.conditions = total.population, start.date, end.date = max(observed.data$date), p)
    
    weights <- lapply(observed.data[, -"date"], function (obs.data.col) 1)
    error <- CalcError(seir, observed.data, weights)
    
    badness <- error$sum.error.sq
    if (check.peak) {
      return(c(badness=badness, peak = error$peak.before.first.observed))
    } else {
      cat("x = ", x, "f(x) = ", badness, "\n")
      badness
    }
  }
  best.fit <- rep(NA_real_, nrow(params))
  for (j in 1:nrow(params)) {
    
    #find upper bound
    x <- c(10^(-30:0), observed.data[1, hosp], length.out = 100)
    fx <- sapply(x, f, p = params[j], check.peak = T)
    if (!all(fx["peak", ]  == 1)) {
      fx["badness", ][fx["peak", ] == 1] <- Inf
    }
    upper <- x[which.min(fx["badness", ]) + 1]
    stopifnot(is.finite(upper))
    
    opt <- optimize(f, interval = c(0, upper), p = params[j], check.peak = F,  tol=1e-30)
    best.fit[j] <- opt$minimum
    if (opt$objective > 10) {
      print(opt)
    }
  }
  
  
  seir <- Seir(best.fit, initial.conditions = total.population, start.date, end.date, params)
  
}

RunExample <- function() {
  if (T) {
    set.seed(NULL); sseed <- round(runif(1, 1, 10000)); set.seed(sseed); cat("sseed = ", sseed, "\n")
  } else {
    set.seed(2267)
  }
  
  start.date <- as.Date("2020/1/23")
  hosp.date <- as.Date("2020/3/25")
  
  num.params <- 20
  latent.period <- RandInt(num.params, 0, 6)
  params <- data.table(
    latent.period = latent.period,
    illness.length.given.nonhosp = RandInt(num.params, 2, 10),
    exposed.to.hospital = RandInt(num.params, latent.period + 1, 12),
    hosp.length.of.stay = RandInt(num.params, 5, 15),
    prop.hospitalized = RandInt(num.params, 4, 6) / 100,
    prop.icu = 0,
    prop.vent = 0,
    patients.in.hosp.are.infectious = T,
    use.hosp.rate = T,
    intervention1_smooth_days = 7,
    intervention1_date = as.Date("2020/3/6"),
    intervention1_multiplier = sample(c(.6, .4, .5, .8, 1), size = num.params, replace = T),
    intervention2_smooth_days = 7,
    intervention2_date = as.Date("2020/3/16"),
    intervention2_multiplier = sample(c(.45, .3, .4, .8, 1), size = num.params, replace = T),
    r0.initial = runif(num.params, 2.5, 4.5)
  )
  
  observed.data <- data.table(date = hosp.date + 0:5, hosp = c(20, 22, 26, 27, 30, 35))
  end.date <- max(observed.data$date)
  
  #RunSim example
  if (T) {
    sim <- RunSim.new(params, total.population = 880000, observed.data = observed.data, start.date = start.date, end.date = end.date)
  }
  
  #Seir examples
  if (F) {
    #single params, single initial.new.exposures
    print(Seir(initial.new.exposures = 0.1, initial.conditions = 880000, start.date, end.date, params[1]))
    #multiple params, single initial.new.exposures
    print(Seir(initial.new.exposures = 5, initial.conditions = 880000, start.date, end.date, params[1]))
    
    #multiple params, multiple initial.new.exposures
    print(Seir(initial.new.exposures = runif(nrow(params)), initial.conditions = 880000, start.date, end.date, params))
  }
}

FindSearchBug <- function(use.slow) {
  
  
  start.date <- as.Date("2020/1/23")
  hosp.date <- as.Date("2020/3/25")
  
  num.params <- 20
  latent.period <- RandInt(num.params, 0, 6)
  params <- data.table(
    latent.period = latent.period,
    illness.length.given.nonhosp = RandInt(num.params, 2, 10),
    exposed.to.hospital = RandInt(num.params, latent.period + 1, 12),
    hosp.length.of.stay = RandInt(num.params, 5, 15),
    prop.hospitalized = RandInt(num.params, 4, 6) / 100,
    prop.icu = 0,
    prop.vent = 0,
    patients.in.hosp.are.infectious = T,
    use.hosp.rate = T,
    intervention1_smooth_days = 7,
    intervention1_date = as.Date("2020/3/6"),
    intervention1_multiplier = sample(c(.6, .4, .5, .8, 1), size = num.params, replace = T),
    intervention1_smooth_days = 7,
    intervention1_date = as.Date("2020/3/16"),
    intervention1_multiplier = sample(c(.45, .3, .4, .8, 1), size = num.params, replace = T),
    r0.initial = runif(num.params, 2.5, 4.5)
  )
  
  # params <- params[c(12)]
  
  # observed.data <- data.table(date = hosp.date, hosp = 20)
  observed.data <- data.table(date = hosp.date + 0:5, hosp = c(20, 22, 26, 27, 30, 35))
  end.date <- max(observed.data$date)
  
  if (use.slow) {
    sim <- RunSim.slow(params, total.population = 880000, observed.data = observed.data, start.date = start.date, end.date = end.date)
  } else {
    sim <- RunSim.new(params, total.population = 880000, observed.data = observed.data, start.date = start.date, end.date = end.date)
  }
  
  cat("\n\n")
  
  zz <- sim$hosp[as.character(observed.data$date), ]
  if (is.vector(zz) == 1) {
    print(summary(zz))
  } else {
    print(cbind(observed.data, rowQuantiles(zz)))
  }
  sq.err <- colSums((zz - observed.data$hosp)^2)
  print(summary(sq.err))
  cat()
}

#adrop as written doesn't work for mpfr
adrop <- function (x, drop = TRUE, named.vector = TRUE, one.d.array = FALSE, 
                   ...) 
{
  if (is.null(dim(x))) 
    stop("require an object with a dim attribute")
  if (length(list(...))) 
    if (length(names(list(...)))) 
      stop("have unrecognized ... arguments for adrop.default: ", 
           paste(names(list(...)), collapse = ", "))
  else stop("have unrecognized unnamed ... arguments for adrop.default")
  x.dim <- dim(x)
  if (is.logical(drop)) {
    if (length(drop) != length(x.dim)) 
      stop("length of drop is not equal length of dim(x)")
    drop <- which(drop)
  }
  else if (is.character(drop)) {
    if (any(is.na(i <- match(drop, names(x.dim))))) 
      stop("dimension names ", paste("'", drop[is.na(i)], 
                                     "'", sep = "", collapse = " "), " not found in x")
    drop <- i
  }
  else if (is.null(drop)) {
    drop <- numeric(0)
  }
  if (!is.numeric(drop) || any(is.na(drop)) || any(drop < 1 | 
                                                   drop > length(x.dim))) 
    stop("drop must contain dimension numbers")
  if (!all(x.dim[drop] == 1)) 
    stop("dimensions to drop (", paste(drop, collapse = ", "), 
         ") do not have length 1")
  x.dimnames <- dimnames(x)
  #dimnames(x) <- NULL
  dimnames(x) <- vector("list", length(dimnames(x)))
  dim(x) <- NULL
  keep <- setdiff(seq(len = length(x.dim)), drop)
  if (length(x.dim[keep]) > 1 || (length(x.dim[keep]) == 1 && 
                                  one.d.array)) {
    dim(x) <- x.dim[keep]
    if (!is.null(x.dimnames)) 
      dimnames(x) <- x.dimnames[keep]
  }
  else if (length(x.dim[keep]) == 1 && named.vector) {
    names(x) <- x.dimnames[keep][[1]]
  }
  else {
  }
  x
}

sweep1 <- function(x, MARGIN, STATS, FUN) {
  sweep(X, MARGIN, as.array(STATS), FUN)
}
