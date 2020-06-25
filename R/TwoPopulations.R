TwoPop <- function(input.file1, input.file2) {
  inputs1 <- LEMMA:::GetInputs(input.file1)
  inputs2 <- LEMMA:::GetInputs(input.file2)

  seir_inputs1 <- LEMMA:::GetStanInputs(inputs1)
  seir_inputs2 <- LEMMA:::GetStanInputs(inputs2)


  same.params <- c("mu_duration_latent", "sigma_duration_latent", "mu_duration_rec_mild", "sigma_duration_rec_mild", "nt", "nobs_types", "nobs", "tobs", "ninter", "mu_t_inter", "sigma_t_inter", "mu_len_inter", "sigma_len_inter", "mu_frac_pui", "sigma_frac_pui")
  cat("The following are assumed equal in both populations:\n")
  print(same.params)
  cat("They are read from input.file1 and are ignored in input.file2\n\n")

  diff.params <- setdiff(names(seir_inputs1), same.params)
  seir_inputs <- seir_inputs1[same.params]
  for (i in diff.params) {
    if (length(seir_inputs1[[i]]) == 1) {
      seir_inputs[[i]] <- c(seir_inputs1[[i]], seir_inputs2[[i]])
    } else {
      ndim <- length(dim(seir_inputs1[[i]]))
      if (ndim == 0) {
        seir_inputs[[i]] <- cbind(seir_inputs1[[i]], seir_inputs2[[i]])
      } else if (ndim == 2) {
        seir_inputs[[i]] <- abind::abind(seir_inputs1[[i]], seir_inputs2[[i]], along = 3)
      } else {
        stop("??")
      }
    }
  }
  seir_inputs$npops <- 2
  seir_inputs$population <- seir_inputs$npop
  seir_inputs$npop <- NULL #this was renamed

  internal.args <- inputs1$internal.args

  run_time <- system.time({
    stan_seir_fit <- rstan::sampling(stanmodels$LEMMA,
                                 data = seir_inputs,
                                 seed = internal.args$random.seed,
                                 iter = internal.args$iter,
                                 chains = 4,
                                 cores = internal.args$cores,
                                 refresh = internal.args$refresh,
                                 control = list(max_treedepth = internal.args$max_treedepth, adapt_delta = internal.args$adapt_delta),
                                 pars = c("error", "beta", "x"),
                                 include = FALSE
    )
  })
  print(run_time)
  quantiles1 <- GetQuantiles2(stan_seir_fit, inputs1, 1)
  quantiles2 <- GetQuantiles2(stan_seir_fit, inputs1, 2)
  g1 <- GetPlots(quantiles1, inputs1, "white")
  g2 <- GetPlots(quantiles2, inputs2, "African American")
  return(list(fit = stan_seir_fit, g1 = g1, g2 = g2))
}


GetInputs <- function(input.file) {
  use.data <- "census"
  cat("only census data is currently used\n")
  sheets <- LEMMA:::ReadInputs(input.file)
  if (use.data == "both") {
    #nothing
  } else if (use.data == "census") {
    sheets$Data$cum.admits.conf <- NA_integer_
  } else if (use.data == "admits") {
    sheets$Data$hosp.conf <- NA_integer_
  } else {
    stop("??")
  }

  inputs <- LEMMA:::ProcessSheets(sheets, input.file)
  return(inputs)
}


GetQuantiles2 <- function(fit, inputs, ipop) {
  dates <- seq(inputs$internal.args$simulation.start.date + 1, inputs$model.inputs$end.date, by = "day")

  sim.data <- rstan::extract(fit, pars = "sim_data")[[1]]

  quantiles <- sapply(DataTypes(), function (i) {
    sim.data.index <- switch(i, hosp = 1, icu = 2, deaths = 3, cum.admits = 4, stop("unexpected bounds name"))
    q <- colQuantiles(sim.data[, sim.data.index, , ipop], probs = seq(0, 1, by = 0.05))
    rownames(q) <- as.character(dates)
    return(q)
  }, simplify = FALSE)
  return(quantiles)
}

GetPlots <- function(quantiles, inputs, subtitl) {
  filestr.out <- paste0(inputs$internal.args$output.filestr, ".pdf")
  short.term <- long.term <- list()
  for (i in names(quantiles)) {
    short.term[[i]] <- GetProjectionPlot(short.term = T, quantiles = quantiles, data.type = i, inputs = inputs)
    long.term[[i]] <- GetProjectionPlot(short.term = F, quantiles = quantiles, data.type = i, inputs = inputs)
  }
  gg <- lapply(c(short.term, long.term), function (z) z + labs(subtitle = subtitl))
  grDevices::pdf(file = filestr.out, width = 9.350, height = 7.225)
  print(gg)
  dev.off()
  cat("\nPDF output: ", filestr.out, "\n")
  return(gg)
}


