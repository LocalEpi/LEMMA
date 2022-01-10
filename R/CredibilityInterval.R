#' @import data.table
#' @import matrixStats

#` Compute credibility interval without writing output
CredibilityIntervalData <- function(inputs, fit.to.data = NULL) {
  inputs$all.inputs.str <- ToString(inputs)
  inputs.copy <- copy(inputs)

  if (is.null(fit.to.data)) {
    fit.to.data <- RunSim(inputs)
  }
  fit.extended <- ExtendSim(list(inputs = inputs, fit.to.data = fit.to.data))
  #projection <- GetProjection(fit.extended, inputs)
  projection <- NULL #temp
  return(
    list(
      fit.to.data = fit.to.data,
      fit.extended = fit.extended,
      projection = projection,
      inputs = inputs.copy
    )
  )
}

#` Main function to calculate credibility interval
CredibilityInterval <- function(inputs, fit.to.data = NULL) {
  TestOutputFile(inputs$internal.args$output.filestr)

  ci_data <- CredibilityIntervalData(inputs, fit.to.data)

  return(ci_data)

  excel.output <- GetExcelOutput(ci_data$projection, ci_data$fit.to.data, ci_data$inputs)
  gplot <- GetPdfOutput(ci_data$fit.extended, ci_data$projection, ci_data$inputs)
  invisible(
    list(
      fit.to.data = ci_data$fit.to.data,
      fit.extended = ci_data$fit.extended,
      projection = ci_data$projection,
      gplot = gplot,
      excel.output = excel.output,
      inputs = ci_data$inputs
    )
  )
}

ProjectScenario <- function(lemma.object, new.inputs) {
  orig.inputs <- lemma.object$inputs
  new.inputs$all.inputs.str <- ToString(new.inputs)
  #check that new.inputs only changes values after end of observed data
  CompareInputs(orig.inputs, new.inputs)

  fit.to.data <- lemma.object$fit.to.data
  invisible(CredibilityInterval(new.inputs, lemma.object$fit.to.data))
}

#order needs to match LEMMA.stan:
# int obs_hosp_census = 1;
# int obs_cases = 2;
DataTypes <- function() c("hosp", "cases")

GetStanInputs <- function(inputs) {
  dt <- melt(inputs$params, id.vars = "name")
  setkey(dt, "name")

  # mu and sigma from inputs$params
  seir_inputs <- as.list(dt[, value])
  names(seir_inputs) <- dt[, paste0(variable, "_", name)]

  # the start date at which to start the simulation
  day0 <- inputs$internal.args$simulation.start.date

  # the number of days to run the model for
  nt <- as.numeric(inputs$model.inputs$end.date - day0)
  seir_inputs[['nt']] = nt

  seir_inputs[['nobs_types']] <- length(DataTypes())

  obs.data <- copy(inputs$obs.data)

  # Observed Data
  num.data.types <- length(DataTypes())
  nobs <- rep(NA, num.data.types)
  obs.mat <- tobs.mat <- matrix(-1, nrow(inputs$obs.data), num.data.types)
  for (i in 1:num.data.types) {
    data.type <- DataTypes()[i]
    date <- inputs$obs.data$date
    conf <- inputs$obs.data[, get(paste0(data.type, ".conf"))]
    pui <- inputs$obs.data[, get(paste0(data.type, ".pui"))]
    if (!all(is.na(pui))) {
      if (any(is.na(conf) != is.na(pui))) {
        stop(data.type, ": If some of Data PUI is not NA (blank), then all dates where Confirmed is NA should be have PUI is NA also and vice versa")
      }
    } else {
      pui <- rep(0, length(pui))
    }
    frac_pui <- inputs$frac_pui[name == data.type, mu]
    combined <- conf + frac_pui * pui
    tobs <- as.numeric(date - day0)

    index <- !is.na(combined)
    nobs[i] <- sum(index)
    if (nobs[i] > 0) {
      obs.mat[1:nobs[i], i] <- combined[index]
      tobs.mat[1:nobs[i], i] <- tobs[index]
    }
  }
  nobs_max <- max(nobs)
  seir_inputs[['nobs']] <- nobs
  seir_inputs[['nobs_max']] <- nobs_max
  seir_inputs[['obs_data']] <- t(obs.mat[1:nobs_max, ])
  seir_inputs[['tobs']] <- t(tobs.mat[1:nobs_max, ])

  # total population
  seir_inputs[['npop']] = inputs$model.inputs$total.population

  Loess <- function(values, span = 0.5) {
    dt <- data.table(values, index = 1:length(values))
    m <- loess(values ~ index, data = dt, degree = 1, span = span)
    pmax(0, predict(m, newdata = dt))
  }
  EstSigmaObs <- function(data.type) {
    j <- which(data.type == DataTypes())
    if (nobs[j] > 10) {
      y <- obs.mat[1:nobs[j], j]
      yhat <- Loess(y)
      return(pmax(0.001, sd(y - yhat))) #pmax to avoid problems with sd = 0
    } else {
      return(1)
    }
  }
  sigma_obs <- sapply(DataTypes(), EstSigmaObs)
  stopifnot(length(inputs$internal.args$weights) == num.data.types)
  sigma_obs <- sigma_obs / inputs$internal.args$weights #e.g. weight 0.5 = make prior on sigma_obs twice as large - less exact fit to that data type
  seir_inputs[['sigma_obs_est_inv']] <- 1 / sigma_obs

  # real<lower=0.0> frac_hosp_lemma;
  # real<lower=0.0> VE_infection;
  # real<lower=0.0> VE_infection_delta;
  # int<lower=0> holiday_start;
  # int<lower=0> holiday_end;
  # real<lower=0.0> holiday_multiplier;
  # real<lower=0.0> omicron_recovered_booster_scale;
  # real<lower=0.0> num_boosters[nt];
  # real<lower=0.0> booster_VE_infection;
  # real<lower=0.0> booster_VE_severe;
  seir_inputs <- c(seir_inputs, inputs$omicron) #temp
  seir_inputs$holiday_start <- as.numeric(seir_inputs$holiday_start - day0)
  seir_inputs$holiday_end <- as.numeric(seir_inputs$holiday_end - day0)
  seir_inputs$num_boosters <- seir_inputs$num_boosters[1:nt]

  # lambda parameter for initial conditions of infected
  seir_inputs[['lambda_initial_infected']] = 1 / inputs$internal.args$intial_infected

  return(seir_inputs)
}

RunSim <- function(inputs) {
  inputs$model.inputs$end.date <- max(inputs$obs.data$date)
  seir_inputs <- GetStanInputs(inputs)
  internal.args <- inputs$internal.args

  GetInit <- function(chain_id) {
    init.names <- grep("^mu_", names(seir_inputs), value = T)
    init <- seir_inputs[init.names]
    names(init) <- sub("mu_", "", init.names)
    init <- c(init, list(sigma_obs = 1 / seir_inputs$sigma_obs_est_inv, ini_exposed = 1 / seir_inputs$lambda_ini_exposed))
    return(init)
  }
  # message('NOTE: You may see an error message (non-finite gradient, validate transformed params, model is leaking).\nThat is fine - LEMMA is working properly as long as it says "Optimization terminated normally"')
  fit <- rstan::sampling(stanmodels$LEMMA,
                         data = seir_inputs,
                         seed = inputs$internal.args$random.seed,
                         init = GetInit,
                         # iter = 100, #temp
                         verbose = F,
                         refresh = F,
                         cores = 4,
  )
  # if (fit$return_code != 0) {
  #   warning("Stan code did not converge! Results are not reliable. return_code = ", fit$return_code)
  # }
  return(fit)
}

ExtendSim_old <- function(inputs, fit.to.data) {
  GetInit <- function(chain_id) {
    init <- fit.to.data$par
    return(init)
  }

  seir_inputs <- GetStanInputs(inputs)
  internal.args <- inputs$internal.args

  stan_seir_fit <- rstan::sampling(stanmodels$LEMMA,
                                   data = seir_inputs,
                                   seed = internal.args$random.seed,
                                   iter = 1,
                                   algorithm = "Fixed_param",
                                   cores = 1,
                                   chains = 1,
                                   refresh = 0,
                                   init = GetInit)

  #convert to format used by optimizing
  fit <- list(par = lapply(rstan::extract(stan_seir_fit), drop))
  return(fit)
}

ExtendSim <- function(lemma.object, new.interventions=NULL, extend.iter=NULL) {
  ExtractIter <- function(fit_element, chain_id) {
    ndim <- length(dim(fit_element))
    if (ndim == 1) {
      fit_element[chain_id]
    } else if (ndim == 2) {
      fit_element[chain_id, ]
    } else if (ndim == 3) {
      fit_element[chain_id, , ]
    } else {
      stop("unexpected ndim")
    }
  }
  GetInit <- function(chain_id) {
    init <- lapply(params, ExtractIter, chain_id)
    if (!is.null(new.interventions)) {
      n <- nrow(new.interventions)
      new_beta_multiplier <- pmax(0.01, rnorm(n, new.interventions$mu_beta_inter, new.interventions$sigma_beta_inter))
      new_t_inter <- pmax(1.01, rnorm(n, new.interventions$mu_t_inter - day0, new.interventions$sigma_t_inter))
      new_len_inter <- pmax(1.01, rnorm(n, new.interventions$mu_len_inter, new.interventions$sigma_len_inter))
      init$beta_multiplier <- as.array(c(init$beta_multiplier, new_beta_multiplier))
      init$t_inter <- as.array(c(init$t_inter, new_t_inter))
      init$len_inter <- as.array(c(init$len_inter, new_len_inter))
    }
    return(init)
  }
  inputs <- lemma.object$inputs
  if (!is.null(new.interventions)) {
    max.obs.data.date <- max(inputs$obs.data$date)
    if (any(new.interventions$mu_t_inter <= max.obs.data.date)) {
      stop("dates in new.interventions must be after last observed data")
    }
    inputs$interventions <- rbind(inputs$interventions, new.interventions)
  }

  day0 <- inputs$internal.args$simulation.start.date

  fit.to.data <- lemma.object$fit.to.data
  params <- rstan::extract(fit.to.data, pars = c("x", "sim_data", "new_cases", "soon_positive", "sim_data_with_error", "lp__"), include = F)
  seir_inputs <- GetStanInputs(inputs)
  seir_inputs$extend <- 1L
  internal.args <- inputs$internal.args
  total.chains <- nrow(as.matrix(fit.to.data))
  if (is.null(extend.iter)) {
    extend.iter <- total.chains
    chain.id <- 1:total.chains
  } else {
    if (extend.iter > total.chains) {
      stop("extend.iter cannot be greater than total.chains")
    }
    chain.id <- sample.int(total.chains, extend.iter)
  }


  out <- capture.output(
    run_time <- system.time({
      stan_seir_fit <- rstan::sampling(stanmodels$LEMMA,
                                       data = seir_inputs,
                                       seed = internal.args$random.seed,
                                       iter = 1,
                                       algorithm = "Fixed_param",
                                       chains = extend.iter,
                                       chain_id = chain.id,
                                       cores = 1,
                                       refresh = 0,
                                       init = GetInit
      )
    })
  )
  if (any(grepl("error", out, ignore.case = T))) {
    print(out)
  }
  print(run_time)
  return(stan_seir_fit)
}


GetProjection <- function(fit, inputs) {
  date <- seq(inputs$internal.args$simulation.start.date + 1, inputs$model.inputs$end.date, by = "day")
  projection <- sapply(DataTypes(), function (i) {
    sim.data.index <- switch(i, hosp = 1, cases = 2, stop("unexpected bounds name"))
    return(fit$par$sim_data[sim.data.index, ])
  }, simplify = FALSE)

  if (F) {
    projection$rt <- fit$par$Rt

    # int Su = 1;
    # int Sv = 2;
    # int Eu = 3;
    # int Ev = 4;
    # int Imildu = 5; //note: on compartment diagram this is labelled "Inonhospu"
    # int Imildv = 6; //note: on compartment diagram this is labelled "Inonhospv"
    # int Iprehu = 7;
    # int Iprehv = 8;
    # int Hmodu  = 9;
    # int Hmodv  = 10;
    # int Hicuu  = 11;
    # int Hicuv  = 12;
    # int Rliveu = 13;
    # int Rlivev = 14;
    # int Rpremort_nonhospu = 15;
    # int Rpremort_nonhospv = 16;
    # int Rmortu = 17;
    # int Rmortv = 18;
    # int Sv_fail = 19;
    x <- fit$par$x

    projection$exposed <- colSums(x[3:4, ])
    projection$infected <- colSums(x[5:8, ])
    projection$activeCases <- colSums(x[3:12, ])
    projection$totalCases <- fit$par$total_cases
    projection$susceptibleUnvax <- x[1, ]
    projection$vaccinated <- colSums(x[c(2, 4, 6, 8, 10, 12, 14, 16, 19), ])
    projection$vaccinatedSuccessfully <- colSums(x[c(2, 4, 6, 8, 10, 12, 14, 16), ])
    #excludes vaccine failed

    projection$deathsU <- x[17, ]
    projection$deathsV <- x[18, ]
    projection$admitsU <- fit$par$new_admitsu
    projection$admitsV <- fit$par$new_admitsv
    projection$totalCasesU <- fit$par$total_casesu
    projection$totalCasesV <- fit$par$total_casesv

    projection$relativeEffectiveContactRate <- fit$par$beta / fit$par$beta[1]
    projection$effectiveContactRate <- fit$par$beta
  }

  return(data.table(date, as.data.table(projection)))
}

ToString <- function(inputs.orig) {
  #Make a human readable string from the inputs

  inputs <- copy(inputs.orig)

  inputs$time.of.run <- as.character(Sys.time())
  inputs$LEMMA.version <- getNamespaceVersion("LEMMA")

  prev.width <- getOption("width")
  prev.print.nrows <- getOption("datatable.print.nrows")
  prev.max.print <- getOption("max.print")
  options(width = 300, datatable.print.nrows = 10000, max.print = 10000)
  all.inputs.str <- utils::capture.output(print(inputs))
  options(width = prev.width, datatable.print.nrows = prev.print.nrows, max.print = prev.max.print)
  all.inputs.str <- c("NOTE: set font to Courier to read", all.inputs.str)
  return(all.inputs.str)
}

CompareInputs <- function(orig.inputs, new.inputs) {
  orig.inputs <- copy(orig.inputs)
  new.inputs <- copy(new.inputs)

  match.not.needed <- c("all.inputs.str", "internal.args", "model.inputs")
  match.needed <- c("params", "frac_pui", "obs.data")
  stopifnot(setequal(names(orig.inputs), names(new.inputs)))
  stopifnot(setequal(names(orig.inputs), c(match.not.needed, match.needed)))
  for (i in match.needed) {
    eq <- all.equal(orig.inputs[[i]], new.inputs[[i]], tolerance = 1e-4, check.attributes = F)
    if (!isTRUE(eq)) {
      cat(i, " does not match\n")
      print(eq)
      stop("in ProjectScenario, lemma.object$inputs must be the same as new.inputs, other than changes that are after the last observed data")
    }
  }
  invisible(NULL)
}

