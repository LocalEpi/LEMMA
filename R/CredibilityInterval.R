#' @import data.table
#' @import matrixStats

SimDataTypes <- function() c("hosp", "icu", "deaths", "cum.admits")

ConvertNa <- function(dt) {
  t(data.matrix(dt[, lapply(.SD, function (z) ifelse(is.na(z), -1, z))]))
}

GetInputs <- function(upp.params, model.inputs, bounds.list, internal.args, obs.data) {
  #List to fill with inputs for stan
  seir_inputs <- list()
  
  # the start date at which to start the simulation
  startdate = internal.args$simulation.start.date #this is day 0
  
  # the number of days to run the model for
  seir_inputs[['nt']] = as.numeric(model.inputs$end.date - startdate)
  
  seir_inputs[['nobs_types']] <- 4  
  seir_inputs[['nobs']] <- nrow(obs.data)
  seir_inputs[['tobs']] <- as.numeric(obs.data$date - startdate)
  
  # int obs_hosp_census = 1;
  # int obs_icu_census = 2;
  # int obs_cum_deaths = 3;
  # int obs_cum_admits = 4; 
  seir_inputs[['obs_data_conf']] <- ConvertNa(obs.data[, .(hosp.lower, icu.lower, deaths.lower, cum.admits.lower)])
  seir_inputs[['obs_data_pui']] <- ConvertNa(obs.data[, .(hosp.upper, icu.upper, deaths.upper, cum.admits.upper)])
  
  # the population for each age group
  seir_inputs[['npop']] = model.inputs$total.population
  
  # the fraction of hospitalized cases (ICU + non-ICU)
  seir_inputs[['mu_frac_hosp']] = upp.params$prop.hospitalized 
  seir_inputs[['sigma_frac_hosp']] = 0.02 
  
  # the fraction of ICU cases of hosp
  seir_inputs[['mu_frac_icu']] = upp.params$prop.icu
  seir_inputs[['sigma_frac_icu']] = 0.02 
  
  # the death rate of ICU
  seir_inputs[['mu_frac_mort']] = upp.params$prop.death
  seir_inputs[['sigma_frac_mort']] = 0.04
  
  # mean duration in "exposed" stage
  seir_inputs[['mu_duration_lat']] = upp.params$latent.period 
  # standard deviation (sd) of duration in "exposed" stage
  seir_inputs[['sigma_duration_lat']] = 1.5 
  
  # mean duration in "infectious" stage for mild cases
  seir_inputs[['mu_duration_rec_mild']] = upp.params$illness.length.given.nonhosp 
  # sd of duration in "infectious" stage for mild cases
  seir_inputs[['sigma_duration_rec_mild']] = 2.0 
  
  # mean duration in "infectious" stage for hospitalized cases
  seir_inputs[['mu_duration_pre_hosp']] = upp.params$exposed.to.hospital - upp.params$latent.period  
  # sd of duration in "infectious" stage for hospitalized cases
  seir_inputs[['sigma_duration_pre_hosp']] = 2.0 
  
  # mean duration in hospital for non-ICU cases
  seir_inputs[['mu_duration_hosp_mod']] = upp.params$hosp.length.of.stay
  # sd of duration in hospital for non-ICU cases
  seir_inputs[['sigma_duration_hosp_mod']] = 2.0 
  
  # mean duration in hospital for ICU cases
  seir_inputs[['mu_duration_hosp_icu']] = upp.params$hosp.length.of.stay  #TEMP
  # sd of duration in hospital for ICU cases
  seir_inputs[['sigma_duration_hosp_icu']] = 4.0
  
  # lambda parameter for initial conditions of "exposed"
  seir_inputs[['lambda_ini_exposed']] = 0.3
  
  # mean initial beta estimate
  seir_inputs[['mu_R0']] = upp.params$r0.initial
  # sd initial beta estimate
  seir_inputs[['sigma_R0']] = 0.5 
  
  # number of interventions
  seir_inputs[['ninter']] = 3
  
  # Note: The inputs below must be lists or arrays, even if one intervention is specified 
  
  # start time of each interventions
  intervention_date <- c(upp.params$intervention1.date, upp.params$intervention2.date, upp.params$intervention3.date)
  seir_inputs[['mu_t_inter']] = array(as.numeric(intervention_date - startdate + 1))
  # length of each intervention
  seir_inputs[['mu_len_inter']] = array(c(upp.params$intervention1.smooth.days, upp.params$intervention2.smooth.days, upp.params$intervention3.smooth.days))
  
  # mean change in beta through intervention 
  seir_inputs[['mu_beta_inter']] = array(c(upp.params$intervention1.multiplier, upp.params$intervention2.multiplier, upp.params$intervention3.multiplier))
  # sd change in beta through intervention
  seir_inputs[['sigma_beta_inter']] = array(rep(0.2, seir_inputs[['ninter']]))
  return(seir_inputs)
}

RunSim <- function(upp.params, model.inputs, bounds.list, internal.args, obs.data) {
  seir_inputs <- GetInputs(upp.params, model.inputs, bounds.list, internal.args, obs.data)
  
  quick_test <- F
  if (quick_test) {
    chains <- 100 
    cores <- 1
    iter <- 1
    control <- NULL
    refresh <- 0
    algorithm <- "Fixed_param"
    init = function(chain_id) list(beta_multiplier = seir_inputs[['mu_beta_inter']], R0 = chain_id/30)
  } else {
    extra_iter <- F
    chains <- cores <- 4
    iter <- if (extra_iter) 2000 else 1000
    control <- if (extra_iter) list(max_treedepth = 15, adapt_delta = 0.9) else NULL
    init <- "random"
    algorithm <- "NUTS"
    refresh <-  250
  }
  
  cat("Starting Stan\n")
  run_time <- system.time({
    stan_seir_fit <- stan(
      file = "LEMMA.stan",
      data = seir_inputs,
      #seed = 75624,
      chains = chains,
      iter = iter,
      cores = cores,
      refresh = refresh,
      algorithm = algorithm,
      init = init,
      control = control,
      pars = c("error", "beta", "ini_exposed", "sigma_obs", "x"), 
      include = FALSE
    )
  })
  print(run_time)
  
  return(stan_seir_fit)
}

GetQuantiles <- function(fit, dates) {
  sim.data <- extract(fit, pars = "sim_data")[[1]]
  quantiles <- list()
  for (i in SimDataTypes()) {
    sim.data.index <- switch(i, hosp = 1, icu = 2, deaths = 3, cum.admits = 4, stop("unexpected bounds name"))
    quantiles[[i]] <- colQuantiles(sim.data[, sim.data.index, ], probs = seq(0, 1, by = 0.05))
    rownames(quantiles[[i]]) <- as.character(dates)
  }
  return(quantiles)
}

#` Main function to calculate credibility interval
CredibilityInterval <- function(upp.params, model.inputs, bounds.list, internal.args, extras, obs.data) {
  options("openxlsx.numFmt" = "0.0")
  devlist <- grDevices::dev.list()
  sapply(devlist[names(devlist) == "pdf"], grDevices::dev.off) #shuts down any old pdf (if there was a crash part way)
  
  all.inputs.str <- "fixme" #the old code is pretty slow
  #all.inputs.str <- utils::capture.output(print(sapply(ls(), function(z) get(z)))) #I'm sure there's a better way to do this
  rm(extras) #extra is only used to save extra information to output file
  
  filestr <- paste0(internal.args$output.filestr, if (internal.args$add.timestamp.to.filestr) date() else "")
  TestOutputFile(filestr)
  
  fit <- RunSim(upp.params, model.inputs, bounds.list, internal.args, obs.data)
  dates <- seq(internal.args$simulation.start.date + 1, model.inputs$end.date, by = "day")
  posterior.quantiles <- GetQuantiles(fit, dates)
  excel.output <- GetExcelOutput(posterior.quantiles, model.inputs, filestr, all.inputs.str)
  
  gplot <- GetPdfOutput(fit, posterior.quantiles, filestr, bounds.list, internal.args, model.inputs)
  return(list(fit = fit, gplot = gplot, excel.output = excel.output, filestr = filestr, all.inputs.str = all.inputs.str))
}


