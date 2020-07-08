#' @import data.table
#' @import matrixStats

#` Main function to calculate credibility interval
CredibilityInterval <- function(inputs) {
  TestOutputFile(inputs$internal.args$output.filestr)

  new.interventions <- inputs$interventions[mu_t_inter > max(inputs$obs.data$date)]
  inputs$interventions <- inputs$interventions[mu_t_inter <= max(inputs$obs.data$date)]

  cat("Fitting to observed data\n")
  fit.to.data <- RunSim(inputs)
  cat("Projecting\n")

  fit.extended <- ExtendSim(list(inputs = inputs, fit.to.data = fit.to.data), new.interventions, extend.iter = NULL)
  posterior.quantiles <- GetQuantiles(fit.extended, inputs)
  excel.output <- GetExcelOutput(posterior.quantiles, inputs)
  gplot <- GetPdfOutput(fit.extended, posterior.quantiles, inputs)
  invisible(list(fit.to.data = fit.to.data, fit.extended = fit.extended, posterior.quantiles = posterior.quantiles, gplot = gplot, excel.output = excel.output, inputs = inputs))
}


#Example:
# z <- LEMMA::CredibilityIntervalFromExcel("~/Dropbox/LEMMA_shared/JS code branch/lemma input and output/SF June 13/SFJune13sd0.1.xlsx")
# #new.interventions is the same structure as sheets$Interventions (it can have more than one row)
# new.int <- structure(list(mu_t_inter = structure(18437, class = "Date"),
#                sigma_t_inter = 2, mu_beta_inter = 1.5, sigma_beta_inter = 1e-04,
#                mu_len_inter = 7, sigma_len_inter = 2), row.names = c(NA, -1L), class = "data.frame")
# ProjectScenario(z, new.int, "~/Dropbox/LEMMA_shared/JS code branch/lemma input and output/SF June 13/SFJune13sd0.1-new1.5")
ProjectScenario <- function(lemma.object, new.interventions, new.output.filestr = NULL, extend.iter = NULL) {
  inputs <- lemma.object$inputs
  fit.to.data <- lemma.object$fit.to.data
  if (!is.null(new.output.filestr)) {
    inputs$internal.args$output.filestr <- new.output.filestr
  }
  TestOutputFile(inputs$internal.args$output.filestr)
  fit.extended <- ExtendSim(lemma.object, new.interventions, extend.iter)
  posterior.quantiles <- GetQuantiles(fit.extended, inputs)
  excel.output <- GetExcelOutput(posterior.quantiles, inputs)
  gplot <- GetPdfOutput(fit.extended, posterior.quantiles, inputs)
  invisible(list(fit.to.data = fit.to.data, fit.extended = fit.extended, posterior.quantiles = posterior.quantiles, gplot = gplot, excel.output = excel.output, inputs = inputs))
}

#order needs to match LEMMA.stan:
# int obs_hosp_census = 1;
# int obs_icu_census = 2;
# int obs_cum_deaths = 3;
# int obs_cum_admits = 4;
DataTypes <- function() c("hosp", "icu", "deaths", "cum.admits")

ConvertNa <- function(dt) {
  t(data.matrix(dt[, lapply(.SD, function (z) ifelse(is.na(z), -1, z))]))
}

GetStanInputs <- function(inputs) {
  data.types <- DataTypes()
  dt <- melt(inputs$params, id.vars = "name")
  setkey(dt, "name")

  # mu and sigma from inputs$params
  seir_inputs <- as.list(dt[, value])
  names(seir_inputs) <- dt[, paste0(variable, "_", name)]

  # the start date at which to start the simulation
  day0 <- inputs$internal.args$simulation.start.date

  # the number of days to run the model for
  seir_inputs[['nt']] = as.numeric(inputs$model.inputs$end.date - day0)

  seir_inputs[['nobs_types']] <- length(data.types)
  seir_inputs[['nobs']] <- nrow(inputs$obs.data)
  seir_inputs[['tobs']] <- as.numeric(inputs$obs.data$date - day0)

  # Observed Data - Confirmed
  seir_inputs[['obs_data_conf']] <- ConvertNa(inputs$obs.data[, paste0(DataTypes(), ".", "conf")])
  # Observed Data - PUI
  seir_inputs[['obs_data_pui']] <- ConvertNa(inputs$obs.data[, paste0(DataTypes(), ".", "pui")])

  # total population
  seir_inputs[['npop']] = inputs$model.inputs$total.population


  # lambda parameter for initial conditions of "exposed"
  seir_inputs[['lambda_ini_exposed']] = inputs$internal.args$lambda_ini_exposed

  # interventions
  inputs$interventions$mu_t_inter <- as.numeric(inputs$interventions$mu_t_inter - day0)
  seir_inputs <- c(seir_inputs, inputs$interventions)

  # number of interventions
  seir_inputs[['ninter']] = nrow(inputs$interventions)

  seir_inputs$len_inter_icu <- 60
  seir_inputs$t_inter_icu <- 59
  seir_inputs$mu_frac_icu_multiplier <- 1
  seir_inputs$sigma_frac_icu_multiplier <- 0.2
  seir_inputs$mu_frac_mort_multiplier <- 1
  seir_inputs$sigma_frac_mort_multiplier <- 0.2

  # fraction of PUI that are true positive
  stopifnot(identical(inputs$frac_pui$name, data.types))
  frac_pui <- list(mu_frac_pui = inputs$frac_pui$mu, sigma_frac_pui = inputs$frac_pui$sigma)
  seir_inputs <- c(seir_inputs, frac_pui)
  return(seir_inputs)
}

RunSim <- function(inputs) {
  inputs$model.inputs$end.date <- max(inputs$obs.data$date)
  seir_inputs <- GetStanInputs(inputs)
  seir_inputs$extend <- 0L
  internal.args <- inputs$internal.args

  GetInit <- function(chain_id) {
    init.names <- grep("^mu_", names(seir_inputs), value = T)
    init <- seir_inputs[init.names]
    names(init) <- sub("mu_", "", init.names)
    names(init) <- sub("beta_inter", "beta_multiplier", names(init)) #beta_multiplier is inconsistently named
    names(init) <- sub("frac_pui", "frac_PUI", names(init)) #frac_PUI is inconsistently named
    init <- c(init, list(sigma_obs = rep(1, length(init$frac_PUI)), ini_exposed = 1 / seir_inputs$lambda_ini_exposed))
    return(init)
  }
  if (is.null(inputs$internal.args$warmup) || is.na(inputs$internal.args$warmup)) {
    warmup <- floor(internal.args$iter / 2)
  } else {
    warmup <- inputs$internal.args$warmup
  }
  run_time <- system.time({
    stan_seir_fit <- rstan::sampling(stanmodels$LEMMA,
  # stan_seir_fit <- rstan::stan("inst/stan/LEMMA.stan",
                                     data = seir_inputs,
                                     seed = internal.args$random.seed,
                                     iter = internal.args$iter,
                                     warmup = warmup,
                                     cores = internal.args$cores,
                                     refresh = internal.args$refresh,
                                     control = list(max_treedepth = internal.args$max_treedepth, adapt_delta = internal.args$adapt_delta),
                                     pars = c("error", "x"),
                                     init = GetInit,
                                     include = FALSE
    )
  })
  print(run_time)
  return(stan_seir_fit)
}

ExtendSim <- function(lemma.object, new.interventions, extend.iter) {
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
      init$beta_multiplier <- c(init$beta_multiplier, new_beta_multiplier)
      init$t_inter <- c(init$t_inter, new_t_inter)
      init$len_inter <- c(init$len_inter, new_len_inter)
    }
    return(init)
  }
  inputs <- lemma.object$inputs
  if (!is.null(new.interventions)) {
    max.obs.data.date <- max(inputs$obs.data$date)
    if (any(new.interventions$mu_t_inter <= max.obs.data.date)) {
      stop("dates in new.interventions must be after last observed data")
    }
  }
  inputs$interventions <- rbind(inputs$interventions, new.interventions)
  day0 <- inputs$internal.args$simulation.start.date

  fit.to.data <- lemma.object$fit.to.data
  params <- rstan::extract(fit.to.data)
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
                                       refresh = internal.args$refresh,
                                       pars = c("error", "x"),
                                       include = FALSE,
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


GetQuantiles <- function(fit, inputs) {
  dates <- seq(inputs$internal.args$simulation.start.date + 1, inputs$model.inputs$end.date, by = "day")

  sim.data <- rstan::extract(fit, pars = "sim_data")[[1]]

  quantiles <- sapply(DataTypes(), function (i) {
    sim.data.index <- switch(i, hosp = 1, icu = 2, deaths = 3, cum.admits = 4, stop("unexpected bounds name"))
    q <- colQuantiles(sim.data[, sim.data.index, ], probs = seq(0, 1, by = 0.05))
    rownames(q) <- as.character(dates)
    return(q)
  }, simplify = FALSE)

  rt.date <- max(inputs$obs.data$date) - 14
  rt.all <- rstan::extract(fit, pars = "Rt")[[1]]
  rt.quantiles <- colQuantiles(rt.all, probs = seq(0, 1, by = 0.05))
  rownames(rt.quantiles) <- as.character(dates)
  rt.quantiles <- rt.quantiles[dates <= rt.date, ]
  quantiles <- c(quantiles, list(rt = rt.quantiles))
  return(quantiles)
}




