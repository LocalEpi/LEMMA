#' @import ggplot2

GetExcelOutput <- function(quantile.list, model.inputs, filestr, all.inputs.str) {
  display.date.range <- as.character(seq(model.inputs$start.display.date, model.inputs$end.date, by = "day"))
  output.list <- lapply(quantile.list, function (quant) {
    output <- data.table(date = display.date.range, quant[display.date.range, ])
    return(output)
  })
  output.list$all.inputs = all.inputs.str
  filestr.out <- paste0(filestr, ".xlsx")
  openxlsx::write.xlsx(output.list, file = filestr.out)
  cat("\nExcel output: ", filestr.out, "\n")
  return(output.list)
} 
  
PlotHist <- function(fit, pars) {
  g <- stan_hist(fit, pars = pars, fill="steelblue3", alpha = 0.2, bins=30)
  return(g)
}

GetYLabel <- function(compartment.name, long.name) {
  switch(compartment.name, hosp = "Number of COVID19 Patients in Hospital", 
         deaths = "Number of COVID19 Deaths",
         long.name)
}

GetProjectionPlot <- function(short.term, quantiles, bounds.list, plot.observed.data, compartment.name, model.inputs) {
  date <- as.Date(rownames(quantiles))
  dt.plot <- merge(data.table(date), bounds.list$bounds, all.x = T, by = "date")
  if (!is.null(quantiles)) {
    dt.plot <- merge(dt.plot, data.table(date = as.Date(rownames(quantiles)), quantiles), by = "date")
  }
  
  if (short.term) {
    date.breaks <- "1 week"
    max.date <- bounds.list$bounds[!is.na(lower), max(date)] + 3
    title1 <- paste("Short Term", bounds.list$long.name, "Projection")
  } else {
    date.breaks <- "1 month"
    max.date <- max(dt.plot$date) - 3 #makes it look a little nicer when ending on first of the month
    title1 <- paste("Long Term", bounds.list$long.name, "Projection")
  }
 
  if (plot.observed.data) {
    min.date <- bounds.list$bounds[!is.na(lower), min(date)] - 7 
  } else {
    min.date <- model.inputs$start.display.date - 7
  }
  
  obs.size <- 1.5
  lb <- bounds.list$lower.bound.label
  ub <- bounds.list$upper.bound.label
  
  dt.plot <- dt.plot[date >= min.date & date <= max.date]
  gg <- ggplot(dt.plot, aes(x=date)) + 
    theme_light() +
    geom_line(aes(y = `50%`, color = "Median"))
 
 
  if (plot.observed.data) {
    gg <- gg + geom_point(aes(y=lower, color = lb), size = obs.size, shape = 4, na.rm = T)
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
    ylab(GetYLabel(compartment.name, bounds.list$long.name)) +
    labs(title = title1) + 
    scale_color_manual("", values = c("blue", "palegreen4", "red4"), breaks = c("Median", lb, ub)) +
    scale_alpha_manual("", values = c(0.2, 0.3, 0.4), breaks = c("5%-95%", "15%-85%", "25%-75%")) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + 
    scale_x_date(date_breaks = date.breaks, date_labels = "%b %d", expand = expansion()) +
    guides(color = guide_legend("", order = 1, override.aes = over.aes), alpha = guide_legend("", order = 2)) +
    theme(legend.position = "bottom", plot.margin = margin(0.1, 0.2, 0, 0.1, "in"),
          axis.text.x=element_text(hjust = 0.7))
  return(gg)
}

GetPdfOutput <- function(fit, quantiles, filestr, bounds, internal.args, model.inputs) {
  filestr.out <- paste0(filestr, ".pdf")
  grDevices::pdf(file = filestr.out, width = 9.350, height = 7.225)
  
  short.term <- long.term <- list()
  for (i in names(bounds)) {
    b <- copy(bounds[[i]])
    b$bounds[, upper := lower + 0.3 * upper] #temp hack (upper is PUI for now)
    print(short.term[[i]] <- GetProjectionPlot(short.term = T, quantiles = quantiles[[i]], bounds.list = b, plot.observed.data = internal.args$plot.observed.data.short.term, compartment.name = i, model.inputs = model.inputs))
    print(long.term[[i]] <- GetProjectionPlot(short.term = F, quantiles = quantiles[[i]], bounds.list = b, plot.observed.data = internal.args$plot.observed.data.long.term, compartment.name = i, model.inputs = model.inputs))
  }
  
  rt.date <- Sys.Date() - 14
  date.index <- as.numeric(rt.date - internal.args$simulation.start.date)
  rt.pars <- paste0("Rt[", date.index, "]")
  rt <- extract(fit, pars = rt.pars)[[1]]
  rt.quantiles <- round(quantile(rt, c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)), 2)
  print(rt.quantiles)
  subtitl <- paste0(capture.output(print(rt.quantiles)), collapse = "\n")
  g <- PlotHist(fit, rt.pars) + xlab(NULL) + labs(title = paste("Rt as of", rt.date), subtitle = subtitl) + theme(plot.subtitle=element_text(family="Courier"))
  print(g)
  
  pars <- c("R0", "duration_lat", "duration_rec_mild", "duration_pre_hosp", "duration_hosp_mod", 
            "duration_hosp_icu", "frac_hosp", "frac_icu", "frac_mort", 
            "beta_multiplier", "t_inter", "len_inter", "frac_PUI")
  lapply(pars, function (p) print(PlotHist(fit, p)))
  
  grDevices::dev.off()
  cat("\nPDF output: ", filestr.out, "\n")
  return(list(short.term = short.term, long.term = long.term))
}

TestOutputFile <- function(filestr) {
  z <- file.create(filestr) #filestr has no extension, it's not the actual .pdf or .xlsx output
  if (!z) {
    stop("Unable to write to output file: ", filestr, "\nCheck your working directory and the value of internal.args$output.filestr")
  }
  file.remove(filestr) 
}
