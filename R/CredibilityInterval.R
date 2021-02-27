#' @import data.table
#' @import matrixStats

#no longer used in xlsx input - sigma_len_inter, sigma_frac_pui
branchname <- function() "CIbugFixed" #fixme - remove this

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
  nt <- as.numeric(inputs$model.inputs$end.date - day0)
  seir_inputs[['nt']] = nt

  seir_inputs[['nobs_types']] <- length(data.types)

  obs.data <- copy(inputs$obs.data)
  if (IsValidInput(inputs$internal.args$initial.deaths)) {
    obs.data[, deaths.conf := deaths.conf - inputs$internal.args$initial.deaths]
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
        stop(i, ": If some of Data PUI is not NA (blank), then all dates where Confirmed is NA should be have PUI is NA also and vice versa")
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

  # real<lower=0.0> mu_iniE;
  # real<lower=0.0> sigma_iniE;
  # real<lower=0.0> mu_ini_Imild;
  # real<lower=0.0> sigma_ini_Imild;
  # real<lower=0.0> mu_ini_Ipreh;
  # real<lower=0.0> sigma_ini_Ipreh;
  # real<lower=0.0> mu_ini_Rlive;
  # real<lower=0.0> sigma_ini_Rlive;
  # real<lower=0.0> mu_ini_cases;
  # real<lower=0.0> sigma_ini_cases;

  stopifnot(inputs$initial.state$from_beginning %in% 0:1)
  if (inputs$initial.state$from_beginning) {
    na <- 999 #not actually used in LEMMA.stan
    inputs$initial.state <- c(inputs$initial.state, sigma_iniE = na, mu_ini_Imild = na, sigma_ini_Imild = na, mu_ini_Ipreh = na, sigma_ini_Ipreh = na, mu_ini_Rlive = na, sigma_ini_Rlive = na, mu_ini_cases = na, sigma_ini_cases = na)
  }
  seir_inputs <- c(seir_inputs, inputs$initial.state)

  # total population
  seir_inputs[['npop']] = inputs$model.inputs$total.population


  # lambda parameter for initial conditions of "exposed"
  seir_inputs[['lambda_ini_exposed']] = inputs$internal.args$lambda_ini_exposed

  # interventions
  inputs$interventions$mu_t_inter <- as.numeric(inputs$interventions$mu_t_inter - day0)
  setnames(inputs$interventions, "mu_len_inter", "len_inter")
  seir_inputs <- c(seir_inputs, lapply(inputs$interventions, as.array)) #as.array fixes problems if only one intervention

  # number of interventions
  seir_inputs[['ninter']] = nrow(inputs$interventions)


  #each of these is vector length nt - may be passed in as up to end of simulation but in RunSim we only run up to end of observed data

  seir_inputs[['vaccinated_per_day']] <- inputs$vaccines$vaccinated_per_day[1:nt]
  seir_inputs[['vaccine_efficacy_for_susceptibility']] <- inputs$vaccines$vaccine_efficacy_for_susceptibility[1:nt]
  seir_inputs[['vaccine_efficacy_against_progression']] <- inputs$vaccines$vaccine_efficacy_against_progression[1:nt]

  stopifnot(seir_inputs[['vaccine_efficacy_for_susceptibility']] <= seir_inputs[['vaccine_efficacy_against_progression']])

  seir_inputs[['duration_vaccinated']] <- inputs$vaccines$duration_vaccinated[1:nt]
  seir_inputs[['duration_natural']] <- inputs$vaccines$duration_natural[1:nt]
  seir_inputs[['frac_hosp_multiplier']] <- inputs$vaccines$frac_hosp_multiplier[1:nt]
  seir_inputs[['frac_icu_multiplier']] <- inputs$vaccines$frac_icu_multiplier[1:nt]
  seir_inputs[['frac_mort_multiplier']] <- inputs$vaccines$frac_mort_multiplier[1:nt]

  # fraction of PUI that are true positive
  stopifnot(identical(inputs$frac_pui$name, data.types))
  frac_pui <- list(mu_frac_pui = inputs$frac_pui$mu, sigma_frac_pui = inputs$frac_pui$sigma)
  seir_inputs <- c(seir_inputs, frac_pui)
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
      return(pmax(1, sd(y - yhat))) #pmax(1, ) to avoid problems with sd = 0
    } else {
      return(1)
    }
  }
  sigma_obs <- sapply(data.types, EstSigmaObs)
  seir_inputs[['sigma_obs_est_inv']] <- 1 / sigma_obs

  return(seir_inputs)
}

RunSim <- function(inputs) {
  inputs$model.inputs$end.date <- max(inputs$obs.data$date)
  seir_inputs <- GetStanInputs(inputs)
  seir_inputs$extend <- 0L
  internal.args <- inputs$internal.args

  GetInit <- function(chain_id) {
    init.names <- grep("^mu_", names(seir_inputs), value = T)
    if (inputs$ini$from_beginning) {
      init.names <- setdiff(init.names, c("mu_ini_Imild", "mu_ini_Ipreh", "mu_ini_Rlive"))
    }
    init <- seir_inputs[init.names]
    names(init) <- sub("mu_", "", init.names)
    names(init) <- sub("beta_inter", "beta_multiplier", names(init)) #beta_multiplier is inconsistently named
    init <- c(init, list(sigma_obs = 1 / seir_inputs$sigma_obs_est_inv))
    if (!inputs$ini$from_beginning) {
      init$ini_Imild <- as.array(init$ini_Imild)
      init$ini_Ipreh <- as.array(init$ini_Ipreh)
      init$ini_Rlive <- as.array(init$ini_Rlive)
    }
    return(init)
  }
  if (IsValidInput(inputs$internal.args$warmup)) {
    warmup <- inputs$internal.args$warmup
  } else {
    warmup <- floor(internal.args$iter / 2)
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
                                     #pars = c("error"),
  init = GetInit #init = GetInit,
                                     #include = FALSE
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
    if (!lemma.object$inputs$ini$from_beginning) {
      init$ini_Imild <- as.array(init$ini_Imild)
      init$ini_Ipreh <- as.array(init$ini_Ipreh)
      init$ini_Rlive <- as.array(init$ini_Rlive)
    }
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
                                       pars = c("error"),
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
  GetQuant <- function(mat) {
    q <- colQuantiles(mat, probs = seq(0, 1, by = 0.05))
    rownames(q) <- as.character(dates)
    return(q)
  }

  dates <- seq(inputs$internal.args$simulation.start.date + 1, inputs$model.inputs$end.date, by = "day")
  sim.data <- rstan::extract(fit, pars = "sim_data")[[1]]
  sigma.obs <- rstan::extract(fit, pars = "sigma_obs")[[1]]


  quantiles <- sapply(DataTypes(), function (i) {
    sim.data.index <- switch(i, hosp = 1, icu = 2, deaths = 3, cum.admits = 4, stop("unexpected bounds name"))
    nrep <- 100
    sim.data.without.error <- sim.data[, sim.data.index, ]
    niter <- nrow(sim.data.without.error)
    num.days <- ncol(sim.data.without.error)
    sim.data.without.error.rep <- matrix(rep(sim.data.without.error, each=nrep), niter * nrep, num.days)
    error.sd <- sigma.obs[, sim.data.index]
    error.term.rep <- matrix(rnorm(num.days * niter * nrep) * error.sd, niter * nrep, num.days) #recycles error.sd
    sim.data.with.error <- sim.data.without.error.rep + error.term.rep
    q <- GetQuant(sim.data.with.error)
    q[q < 0] <- 0
    return(q)
  }, simplify = FALSE)


  rt.date <- max(inputs$obs.data$date) + 5 #output up to end of observed data here (add 5 to make sure LEMMA Rt is included in Ensemble), but cut off last 14 days in pdf output (keep these last 14 in xlsx output for CalCAT)
  rt.all <- rstan::extract(fit, pars = "Rt")[[1]]
  rt.quantiles <- GetQuant(rt.all)
  rt.quantiles <- rt.quantiles[dates <= rt.date, ]

  # int Su = 1;
  # int Sv = 2;
  # int Eu = 3;
  # int Ev = 4;
  # int Imildu = 5;
  # int Imildv = 6;
  # int Iprehu = 7;
  # int Iprehv = 8;
  # int Hmodu  = 9;
  # int Hmodv  = 10;
  # int Hicuu  = 11;
  # int Hicuv  = 12;
  # int Rliveu = 13;
  # int Rlivev = 14;
  # int Rmort = 15;

  #these don't have a sigma_obs, would need to add it if we had an observed quantity for any of these (not likely)
  x <- rstan::extract(fit, pars = "x")[[1]]

  exposed <- GetQuant(x[, 3, ] + x[, 4, ])
  infected <- GetQuant(x[, 5, ] + x[, 6, ] + x[, 7, ] + x[, 8, ])
  active.cases <- GetQuant(x[, 3, ] + x[, 4, ] + x[, 5, ] + x[, 6, ] + x[, 7, ] + x[, 8, ] + x[, 9, ] + x[, 10, ]+ x[, 11, ] + x[, 12, ])
  total.cases <- GetQuant(rstan::extract(fit, pars = "total_cases")[[1]])
  Su <- GetQuant(x[, 1, ])

  quantiles <- c(quantiles, list(rt = rt.quantiles, exposed = exposed, infected = infected, activeCases = active.cases, totalCases = total.cases, Su = Su))
  if (IsValidInput(inputs$internal.args$initial.deaths)) {
    quantiles$deaths <- quantiles$deaths + inputs$internal.args$initial.deaths
  }
  return(quantiles)
}

IsValidInput <- function(x) {
  stopifnot(length(x) %in% 0:1)
  if (is.null(x)) return(FALSE)
  return(is.finite(x))
}




