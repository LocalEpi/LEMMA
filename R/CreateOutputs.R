#' @import ggplot2

GetExcelOutput <- function(projection, fit, inputs) {
  options("openxlsx.numFmt" = "0.00")
  display.date.range <- seq(inputs$model.inputs$start.display.date, inputs$model.inputs$end.date, by = "day")
  output.list <- list(projection = projection[date %in% display.date.range])

  pars <- c("r0", "duration_latent", "duration_rec_mild", "duration_pre_hosp", "duration_hosp_mod",
            "duration_hosp_icu", "duration_mort_nonhosp",  "frac_hosp", "frac_icu", "frac_mort", "frac_mort_nonhosp" , "frac_tested",
            "ini_exposed")
  output.list$posteriorParams <- data.table(param = pars, posterior = unlist(fit$par[pars]))
  output.list$posteriorInterventions <- data.table(beta_multiplier = as.vector(fit$par$beta_multiplier), t_inter = as.character(fit$par$t_inter + inputs$internal.args$simulation.start.date), len_inter = as.vector(fit$par$len_inter))

  output.list$all.inputs = inputs$all.inputs.str
  output.list$all.inputs.asrun <- inputs$all.inputs.str.asrun
  filestr.out <- paste0(inputs$internal.args$output.filestr, ".xlsx")
  openxlsx::write.xlsx(output.list, file = filestr.out)
  cat("\nExcel output: ", filestr.out, "\n")
  return(output.list)
}


GetYLabel <- function(data.type, long.name) {
  switch(data.type,
         hosp = "Number of COVID19 Patients in Hospital",
         icu = "Number of COVID19 Patients in ICU",
         deaths = "Number of COVID19 Deaths",
         admits = "Number of New COVID19 Admissions to Hospital",
         cases = "Number of COVID19 Cases",
         seroprev = "Seroprevalence (natural or vaccine) Fraction",
         stop("unexpected data.type"))
}

GetTitleLabel <- function(data.type) {
  switch(data.type,
         hosp = "Hospitalization",
         icu = "ICU",
         deaths = "Death",
         admits = "Admissions",
         cases = "Cases",
         seroprev = "Seroprevalence",
         stop("unexpected data.type"))
}

GetProjectionPlot <- function(short.term, projection, data.type, inputs) {
  obs.data <- inputs$obs.data[, .(date, conf = get(paste0(data.type, ".conf")), pui = get(paste0(data.type, ".pui")))]
  if (all(is.na(obs.data$conf))) return(NULL)
  projection.dt <- projection[, c("date", data.type), with = F]
  names(projection.dt)[2] <- "proj"
  dt.plot <- merge(obs.data, projection.dt, all = T, by = "date")

  if (short.term) {
    max.date <- obs.data[!is.na(conf), max(date)] + 3
    title1 <- paste("Short Term", GetTitleLabel(data.type), "Projection")
    plot.observed.data <- inputs$internal.args$plot.observed.data.short.term
  } else {
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
    geom_line(aes(y = proj, color = "Median"))


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

  gg <- gg +
    xlab("") +
    ylab(GetYLabel(data.type)) +
    labs(title = title1, caption = "localepi.github.io/LEMMA") +
    scale_color_manual("", values = c("blue", "palegreen4", "red4"), breaks = c("Median", lb, ub)) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %d", expand = expansion()) +
    guides(color = guide_legend("", order = 1, override.aes = over.aes), alpha = guide_legend("", order = 2)) +
    theme(legend.position = "bottom", plot.margin = margin(0.1, 0.2, 0, 0.1, "in"),
          axis.text.x=element_text(hjust = 0.7))
  return(gg)
}

GetRtPlot <- function(projection, inputs) {
  dt.plot <- projection[date <= (max(inputs$obs.data$date) - 14), .(date, rt)]

  min.date <- as.Date("2020/4/1")
  if (dt.plot[, max(date)] > min.date) {
    dt.plot <- dt.plot[date >= min.date]
  }

  gg <- ggplot(dt.plot, aes(x=date)) +
    theme_light() +
    geom_line(aes(y = rt, color = "Median"))

  gg <- gg +
    xlab("") +
    ylab("Re") +
    labs(title = "Effective Reproduction Number", subtitle = paste0("Rt as of ", as.character(dt.plot[, max(date)]), " = ", dt.plot[.N, sprintf("%2.2f", rt)])) +
    scale_color_manual("", values = "blue") +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %d", expand = expansion()) +
    guides(color = guide_legend("", order = 1), alpha = guide_legend("", order = 2)) +
    theme(legend.position = "bottom", plot.margin = margin(0.1, 0.2, 0, 0.1, "in"),
          axis.text.x=element_text(hjust = 0.7)) +
    geom_hline(yintercept = 1, lty = "dashed")
  print(gg)
  return(gg)
}

expansion <- function() c(0, 0, 0, 0) #fixes a problem if old ggplot version

GetPdfOutput <- function(fit, projection, inputs) {
  devlist <- grDevices::dev.list()
  sapply(devlist[names(devlist) == "pdf"], grDevices::dev.off) #shuts down any old pdf (if there was a crash part way)

  filestr.out <- paste0(inputs$internal.args$output.filestr, ".pdf")
  grDevices::pdf(file = filestr.out, width = 9.350, height = 7.225)

  short.term <- long.term <- list()
  for (i in DataTypes()) {
    if (i == "seroprev" & inputs$internal.args$hide.nonpublic.data) {
      #do not show seroprev (nonpublic)
    } else {
      short.term[[i]] <- GetProjectionPlot(short.term = T, projection = projection, data.type = i, inputs = inputs)
      if (!is.null(short.term[[i]])) print(short.term[[i]])
    }
    long.term[[i]] <- GetProjectionPlot(short.term = F, projection = projection, data.type = i, inputs = inputs)
    if (!is.null(long.term[[i]])) print(long.term[[i]])
  }

  rt.plot <- GetRtPlot(projection, inputs)
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
