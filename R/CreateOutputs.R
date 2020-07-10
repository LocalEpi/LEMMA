#' @import ggplot2

GetExcelOutput <- function(quantile.list, inputs) {
  options("openxlsx.numFmt" = "0.00")

  output.list <- lapply(quantile.list, function (quant) {
    display.date.range <- as.character(seq(inputs$model.inputs$start.display.date, inputs$model.inputs$end.date, by = "day"))
    index <- rownames(quant) %in% display.date.range
    output <- data.table(date = rownames(quant)[index], quant[index, ])
    return(output)
  })
  output.list$all.inputs = inputs$all.inputs.str
  filestr.out <- paste0(inputs$internal.args$output.filestr, ".xlsx")
  openxlsx::write.xlsx(output.list, file = filestr.out)
  cat("\nExcel output: ", filestr.out, "\n")
  return(output.list)
}

#modification of rstan::stan_hist to deal with dates (also rstan::stan_hist crashes when using ExtendSim)
StanHist <- function (object, pars, base.date, ...) {
  dots <- rstan:::.add_aesthetics(list(...), c("fill", "color"))
  samp <- melt(as.data.table(as.matrix(object, pars = pars)), id.vars = integer(0))

  if (!is.null(base.date)) {
    samp$value <- base.date + samp$value
  }
  thm <- rstan:::rstanvis_hist_theme()
  base <- ggplot2::ggplot(samp, ggplot2::aes_string(x = "value", y = "..density.."))
  graph <- base + do.call(ggplot2::geom_histogram, dots) +
    ggplot2::xlab("") + thm
  if (uniqueN(samp$variable) == 1) {
    graph + ggplot2::xlab(unique(samp$variable))
  } else {
    graph + ggplot2::facet_wrap(~variable, scales = "free")
  }
}

PlotHist <- function(fit, pars, base.date) {
  base.date <- if (pars == "t_inter") base.date else NULL
  return(StanHist(fit, pars = pars, fill="steelblue3", alpha = 0.2, bins=30, base.date = base.date))
}

GetYLabel <- function(data.type, long.name) {
  switch(data.type,
         hosp = "Number of COVID19 Patients in Hospital",
         icu = "Number of COVID19 Patients in ICU",
         deaths = "Number of COVID19 Deaths",
         cum.admits = "Number of Cumulative COVID19 Admissions to Hospital",
         stop("unexpected data.type"))
}

GetTitleLabel <- function(data.type) {
  switch(data.type,
         hosp = "Hospitalization",
         icu = "ICU",
         deaths = "Death",
         cum.admits = "Cumulative Admission",
         stop("unexpected data.type"))
}

GetProjectionPlot <- function(short.term, quantiles, data.type, inputs) {
  obs.data <- inputs$obs.data[, .(date, conf = get(paste0(data.type, ".conf")), pui = get(paste0(data.type, ".pui")))]
  if (all(is.na(obs.data$conf))) return(NULL)
  quantiles.dt <- data.table(date = as.Date(rownames(quantiles[[data.type]])), quantiles[[data.type]])
  dt.plot <- merge(obs.data, quantiles.dt, all = T, by = "date")

  if (short.term) {
    date.breaks <- "1 week"
    max.date <- obs.data[!is.na(conf), max(date)] + 3
    title1 <- paste("Short Term", GetTitleLabel(data.type), "Projection")
    plot.observed.data <- inputs$internal.args$plot.observed.data.short.term
  } else {
    date.breaks <- "1 month"
    max.date <- max(dt.plot$date) - 3 #makes it look a little nicer when ending on first of the month
    title1 <- paste("Long Term", GetTitleLabel(data.type), "Projection")
    plot.observed.data <- inputs$internal.args$plot.observed.data.long.term
  }

  if (plot.observed.data) {
    min.date <- obs.data[!is.na(conf), min(date)] - 7
  } else {
    min.date <- inputs$model.inputs$start.display.date - 7
  }

  obs.size <- 1.5
  frac.pui <- inputs$frac_pui[name == data.type, mu]
  lb <- "Confirmed"
  ub <- paste0("Confirmed + ", 100 * frac.pui, "%PUI")
  dt.plot[, upper := conf + frac.pui * pui]

  dt.plot <- dt.plot[date >= min.date & date <= max.date]
  gg <- ggplot(dt.plot, aes(x=date)) +
    theme_light() +
    geom_line(aes(y = `50%`, color = "Median"))


  if (plot.observed.data) {
    gg <- gg + geom_point(aes(y=conf, color = lb), size = obs.size, shape = 4, na.rm = T)
    if (!all(is.na(dt.plot$upper))) {
      gg <- gg + geom_point(aes(y=upper, color = ub), size = obs.size, shape = 4, na.rm = T)
      over.aes <- list(shape = c(NA, 4, 4), linetype = c(1, 0, 0))
    } else {
      over.aes <- list(shape = c(NA, 4), linetype = c(1, 0))
    }
  } else {
    over.aes <- list()
  }

  gg <- gg + geom_ribbon(aes(ymin=`25%`, ymax=`75%`, alpha = "25%-75%"), fill = "blue") +
    geom_ribbon(aes(ymin=`15%`, ymax=`85%`, alpha = "15%-85%"), fill = "blue") +
    geom_ribbon(aes(ymin=`5%`, ymax=`95%`, alpha = "5%-95%"), fill = "blue")

  gg <- gg +
    xlab("") +
    ylab(GetYLabel(data.type)) +
    labs(title = title1, caption = "localepi.github.io/LEMMA") +
    scale_color_manual("", values = c("blue", "palegreen4", "red4"), breaks = c("Median", lb, ub)) +
    scale_alpha_manual("", values = c(0.2, 0.3, 0.4), breaks = c("5%-95%", "15%-85%", "25%-75%")) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    scale_x_date(date_breaks = date.breaks, date_labels = "%b %d", expand = expansion()) +
    guides(color = guide_legend("", order = 1, override.aes = over.aes), alpha = guide_legend("", order = 2)) +
    theme(legend.position = "bottom", plot.margin = margin(0.1, 0.2, 0, 0.1, "in"),
          axis.text.x=element_text(hjust = 0.7))
  print(gg)
  return(gg)
}

GetRtPlot <- function(quantiles, inputs) {
  sim.dates <- inputs$internal.args$simulation.start.date + 1:nrow(quantiles)
  min.date <- as.Date("2020/4/1")
  dt.plot <- data.table(date = sim.dates, quantiles)[date > min.date]

  gg <- ggplot(dt.plot, aes(x=date)) +
    theme_light() +
    geom_line(aes(y = `50%`, color = "Median"))
  gg <- gg + geom_ribbon(aes(ymin=`25%`, ymax=`75%`, alpha = "25%-75%"), fill = "blue") +
    geom_ribbon(aes(ymin=`15%`, ymax=`85%`, alpha = "15%-85%"), fill = "blue") +
    geom_ribbon(aes(ymin=`5%`, ymax=`95%`, alpha = "5%-95%"), fill = "blue")

  gg <- gg +
    xlab("") +
    ylab("Re") +
    labs(title = "Effective Reproduction Number") +
    scale_color_manual("", values = "blue") +
    scale_alpha_manual("", values = c(0.2, 0.3, 0.4), breaks = c("5%-95%", "15%-85%", "25%-75%")) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    scale_x_date(date_breaks = "1 week", date_labels = "%b %d", expand = expansion()) +
    guides(color = guide_legend("", order = 1), alpha = guide_legend("", order = 2)) +
    theme(legend.position = "bottom", plot.margin = margin(0.1, 0.2, 0, 0.1, "in"),
          axis.text.x=element_text(hjust = 0.7)) +
    geom_hline(yintercept = 1, lty = "dashed")
  print(gg)
  return(gg)
}


GetPdfOutput <- function(fit, quantiles, inputs) {
  devlist <- grDevices::dev.list()
  sapply(devlist[names(devlist) == "pdf"], grDevices::dev.off) #shuts down any old pdf (if there was a crash part way)

  filestr.out <- paste0(inputs$internal.args$output.filestr, ".pdf")
  grDevices::pdf(file = filestr.out, width = 9.350, height = 7.225)

  short.term <- long.term <- list()
  for (i in DataTypes()) {
    short.term[[i]] <- GetProjectionPlot(short.term = T, quantiles = quantiles, data.type = i, inputs = inputs)
    long.term[[i]] <- GetProjectionPlot(short.term = F, quantiles = quantiles, data.type = i, inputs = inputs)
  }

  rt <- quantiles$rt[1:(nrow(quantiles$rt) - 14), ] #cutoff Rt plot 14 days before last observed data
  rt.plot <- GetRtPlot(rt, inputs)
  rt.date <- rownames(rt)[nrow(rt)]
  date.index <- as.numeric(as.Date(rt.date) - inputs$internal.args$simulation.start.date)
  rt.pars <- paste0("Rt[", date.index, "]")
  rt.quantiles <- round(quantile(rstan::extract(fit, pars = rt.pars)[[1]], c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)), 2)
  print(rt.quantiles)
  subtitl <- paste0(capture.output(print(rt.quantiles)), collapse = "\n")
  g <- PlotHist(fit, rt.pars) + xlab(NULL) + labs(title = paste("Rt as of", rt.date), subtitle = subtitl) + theme(plot.subtitle=element_text(family="Courier"))
  print(g)

  pars <- c("r0", "duration_latent", "duration_rec_mild", "duration_pre_hosp", "duration_hosp_mod",
            "duration_hosp_icu", "frac_hosp", "frac_icu", "frac_mort",
            "beta_multiplier", "t_inter", "len_inter", "frac_PUI", "ini_exposed")
  lapply(pars, function (p) print(PlotHist(fit, p, base.date = inputs$internal.args$simulation.start.date)))

  grDevices::dev.off()
  cat("\nPDF output: ", filestr.out, "\n")
  return(list(short.term = short.term, long.term = long.term, rt = rt.plot))
}

TestOutputFile <- function(filestr) {
  success <- file.create(filestr) #filestr has no extension, it's not the actual .pdf or .xlsx output
  if (!success) {
    stop("Unable to write to output file: ", filestr, "\nCheck your working directory and the value of internal.args$output.filestr")
  }
  file.remove(filestr)
}
