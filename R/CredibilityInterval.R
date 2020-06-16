#' @import data.table
#' @import matrixStats

#` Main function to calculate credibility interval
CredibilityInterval <- function(inputs) {
  TestOutputFile(inputs$internal.args$output.filestr)
  fit <- RunSim(inputs)
  posterior.quantiles <- GetQuantiles(fit, inputs)
  excel.output <- GetExcelOutput(posterior.quantiles, inputs)
  gplot <- GetPdfOutput(fit, posterior.quantiles, inputs)
  return(list(fit = fit, gplot = gplot, excel.output = excel.output, inputs = inputs))
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

  # number of interventions
  seir_inputs[['ninter']] = nrow(inputs$interventions)

  # interventions
  inputs$interventions[, mu_t_inter := as.numeric(mu_t_inter - day0)]
  seir_inputs <- c(seir_inputs, inputs$interventions)

  # fraction of PUI that are true positive
  stopifnot(identical(inputs$frac_pui$name, data.types))
  frac_pui <- list(mu_frac_pui = inputs$frac_pui$mu, sigma_frac_pui = inputs$frac_pui$sigma)
  seir_inputs <- c(seir_inputs, frac_pui)
  return(seir_inputs)
}

RunSim <- function(inputs) {
  seir_inputs <- GetStanInputs(inputs)
  internal.args <- inputs$internal.args
  cat("Starting Stan\n")
  run_time <- system.time({
    stan_seir_fit <- rstan::sampling(stanmodels$LEMMA,
      data = seir_inputs,
      seed = internal.args$random.seed,
      iter = internal.args$iter,
      cores = internal.args$cores,
      refresh = internal.args$refresh,
      control = list(max_treedepth = internal.args$max_treedepth, adapt_delta = internal.args$adapt_delta),
      pars = c("error", "beta", "ini_exposed", "sigma_obs", "x"),
      include = FALSE
    )
  })
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
  return(quantiles)
}




