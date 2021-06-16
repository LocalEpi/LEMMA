#' @import data.table
#' @import matrixStats


#` Main function to calculate credibility interval
CredibilityInterval <- function(inputs, fit.to.data = NULL) {
  TestOutputFile(inputs$internal.args$output.filestr)
  inputs$all.inputs.str <- ToString(inputs)
  inputs.copy <- copy(inputs)

  new.interventions <- inputs$interventions[mu_t_inter > max(inputs$obs.data$date)]
  inputs$interventions <- inputs$interventions[mu_t_inter <= max(inputs$obs.data$date)]

  if (is.null(fit.to.data)) {
    fit.to.data <- RunSim(inputs)
  }
  fit.extended <- ExtendSim(inputs, fit.to.data, new.interventions)
  projection <- GetProjection(fit.extended, inputs)

  excel.output <- GetExcelOutput(projection, fit.to.data, inputs)
  gplot <- GetPdfOutput(fit.extended, projection, inputs)
  invisible(list(fit.to.data = fit.to.data, fit.extended = fit.extended, projection = projection, gplot = gplot, excel.output = excel.output, inputs = inputs.copy))
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
# int obs_icu_census = 2;
# int obs_cum_deaths = 3;
# int obs_admits = 4;
# int obs_cases = 5;
# int obs_seroprev = 6;
DataTypes <- function() c("hosp", "icu", "deaths", "admits", "cases", "seroprev")

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

  if ("cum.admits.conf" %in% names(inputs$obs.data)) {
    stop("cum.admits has been replaced by (new) admits. Please update your input file.")
  }

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

  seir_inputs <- c(seir_inputs, inputs$initial.state)

  # total population
  seir_inputs[['npop']] = inputs$model.inputs$total.population

  # lambda parameter for initial conditions of "exposed"
  mean.ini <- 1e-5 * seir_inputs[['npop']]
  seir_inputs[['lambda_ini_exposed']] = 1 / mean.ini

  # interventions
  inputs$interventions$mu_t_inter <- as.numeric(inputs$interventions$mu_t_inter - day0)
  seir_inputs <- c(seir_inputs, lapply(inputs$interventions, as.array)) #as.array fixes problems if only one intervention

  # number of interventions
  seir_inputs[['ninter']] = nrow(inputs$interventions)

  dates <- seq(day0 + 1, day0 + nt, by = "day")
  vaccines <- inputs$vaccines[date %in% dates]
  stopifnot(isTRUE(all.equal(vaccines$date, dates)))
  stopifnot(setequal(names(vaccines), c("date", "vaccinated_per_day", "vaccine_efficacy_for_susceptibility", "duration_vaccinated",
                                        "duration_natural", "frac_hosp_multiplier_vaccinated", "frac_hosp_multiplier_unvaccinated", "frac_icu_multiplier_vaccinated", "frac_icu_multiplier_unvaccinated",
                                        "frac_mort_multiplier_vaccinated", "frac_mort_multiplier_unvaccinated", "transmission_variant_multiplier")))
  vaccines$date <- NULL
  seir_inputs <- c(seir_inputs, vaccines)
  seir_inputs[['prior_infection_vaccine_scale']] <- inputs$internal$prior.infection.vaccine.scale

  stopifnot(seir_inputs[['vaccine_efficacy_for_susceptibility']] <= seir_inputs[['vaccine_efficacy_against_progression']])

  #duration_natural and duration_vaccinated are years in Excel input but should be days here
  stopifnot(seir_inputs[["duration_natural"]] > 90)
  stopifnot(seir_inputs[["duration_vaccinated"]] > 90)

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
    names(init) <- sub("beta_inter", "beta_multiplier", names(init)) #beta_multiplier is inconsistently named
    init <- c(init, list(sigma_obs = 1 / seir_inputs$sigma_obs_est_inv, ini_exposed = 1 / seir_inputs$lambda_ini_exposed))
    if (!is.null(inputs$internal.args$init_frac_mort_nonhosp)) {
      #some small counties need a smaller initial value of frac_mort_nonhosp to converge
      init$frac_mort_nonhosp <- inputs$internal.args$init_frac_mort_nonhosp
    }
    return(init)
  }
  message('NOTE: You may see an error message (non-finite gradient, validate transformed params, model is leaking).\nThat is fine - LEMMA is working properly as long as it says "Optimization terminated normally"')
  fit <- rstan::optimizing(stanmodels$LEMMA,
                    data = seir_inputs,
                    seed = inputs$internal.args$random.seed,
                    init = GetInit,
                    iter = inputs$internal.args$optimize.iter,
                    verbose = T,
                    as_vector = F
  )
  if (fit$return_code != 0) {
    warning("Stan code did not converge! Results are not reliable. return_code = ", fit$return_code)
  }
  return(fit)
}

ExtendSim <- function(inputs, fit.to.data, new.interventions) {
  GetInit <- function(chain_id) {
    init <- fit.to.data$par
    if (!is.null(new.interventions)) {
      n <- nrow(new.interventions)
      init$beta_multiplier <- as.array(c(init$beta_multiplier, new.interventions$mu_beta_inter))
      init$t_inter <- as.array(c(init$t_inter, new.interventions$mu_t_inter - day0))
      init$len_inter <- as.array(c(init$len_inter, new.interventions$mu_len_inter))
    }
    return(init)
  }
  if (!is.null(new.interventions)) {
    max.obs.data.date <- max(inputs$obs.data$date)
    if (any(new.interventions$mu_t_inter <= max.obs.data.date)) {
      stop("dates in new.interventions must be after last observed data")
    }
    inputs$interventions <- rbind(inputs$interventions, new.interventions)
  }

  day0 <- inputs$internal.args$simulation.start.date
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


GetProjection <- function(fit, inputs) {
  date <- seq(inputs$internal.args$simulation.start.date + 1, inputs$model.inputs$end.date, by = "day")
  projection <- sapply(DataTypes(), function (i) {
    sim.data.index <- switch(i, hosp = 1, icu = 2, deaths = 3, admits = 4, cases = 5, seroprev = 6, stop("unexpected bounds name"))
    return(fit$par$sim_data[sim.data.index, ])
  }, simplify = FALSE)

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
  x <- fit$par$x

  projection$exposed <- colSums(x[3:4, ])
  projection$infected <- colSums(x[5:8, ])
  projection$activeCases <- colSums(x[3:12, ])
  projection$totalCases <- fit$par$total_cases
  projection$susceptibleUnvax <- x[1, ]
  projection$vaccinated <- colSums(x[c(2, 4, 6, 8, 10, 12, 14, 16), ])

  projection$deathsU <- x[17, ]
  projection$deathsV <- x[18, ]
  projection$admitsU <- fit$par$new_admitsu
  projection$admitsV <- fit$par$new_admitsv
  projection$totalCasesU <- fit$par$total_casesu
  projection$totalCasesV <- fit$par$total_casesv

  return(data.table(date, as.data.table(projection)))
}

ToString <- function(inputs.orig) {
  #Make a human readable string from the inputs

  inputs <- copy(inputs.orig)
  if (inputs$internal.args$hide.nonpublic.data) {
    inputs$obs.data$seroprev.conf <- rep("nonpublic", nrow(inputs$obs.data))
  }

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

  match.not.needed <- c("vaccines_nonstan", "all.inputs.str", "internal.args", "model.inputs")
  match.needed <- c("params", "frac_pui", "obs.data")
  partial.match.needed <- c("vaccines", "interventions")
  stopifnot(setequal(names(orig.inputs), names(new.inputs)))
  stopifnot(setequal(names(orig.inputs), c(match.not.needed, match.needed, partial.match.needed)))
  for (i in match.needed) {
    eq <- all.equal(orig.inputs[[i]], new.inputs[[i]], tolerance = 1e-4, check.attributes = F)
    if (!isTRUE(eq)) {
      cat(i, " does not match\n")
      print(eq)
      stop("in ProjectScenario, lemma.object$inputs must be the same as new.inputs, other than changes that are after the last observed data")
    }
  }
  last.obs.date <-  max(orig.inputs$obs.data$date)
  for (i in partial.match.needed) {
    orig <- orig.inputs[[i]]
    new <- new.inputs[[i]]
    if (i == "vaccines") {
      orig <- orig[date <= last.obs.date]
      new <- new[date <= last.obs.date]
    } else if (i == "interventions") {
      orig <- orig[mu_t_inter <= last.obs.date]
      new <- new[mu_t_inter <= last.obs.date]
    } else {
      stop("not expected")
    }
    eq <- all.equal(orig, new, tolerance = 1e-4, check.attributes = F)
    if (!isTRUE(eq)) {
      cat("last.obs.date = ", as.character(last.obs.date), "\n")
      cat(i, " does not match up to last.obs.date\n")
      print(eq)
      stop("in ProjectScenario, lemma.object$inputs must be the same as new.inputs, other than changes that are after the last observed data")
    }
  }
  invisible(NULL)
}

