#' @import ggplot2

# upp = User's Prior Projection

GetModelName <- function(dt) {
  model.name <- dt[, factor(paste0(as.integer(hasE), as.integer(hospInf), as.integer(hospRate)))]
  return(model.name)
}

GetPlotTitle <- function(posterior.niter) {
  if (posterior.niter == 0) {
    warn.str <- "\nNo credibility interval available"
  } else if (posterior.niter < 100) {
    warn.str <- "\nVery low number of iterations, do not use for inference"
  } else if (posterior.niter < 500) {
    warn.str <- "\nLow number of iterations, use with caution"
  } else {
    warn.str <- ""
  }
  paste0("Posterior Distribution, niter = ", posterior.niter, warn.str)
}

GetExcelOutput <- function(sim, upp, in.bounds, upp.in.bounds, date.range, filestr, all.inputs.str) {
  output.list <- list(hosp = NULL, icu = NULL, vent = NULL, active.cases = NULL, total.cases = NULL)
  output.names <- names(output.list)
  
  probs2 <- c(0.95, 1, 0.15, 0.25, seq(0.55, 0.9, by = 0.05))
  for (j in output.names) {
    sim.accepted <- sim[[j]][, in.bounds, drop = F]
    if (is.null(sim)) {
      quant1 <- quant2 <- NULL
    } else {
      quant1 <- rowQuantiles(sim.accepted, probs = c(0, 0.05, 0.5))
      quant2 <- rowQuantiles(sim.accepted, probs = probs2)
    }
    output <- data.table(date = date.range, quant1, upp = upp[[j]], quant2)
    output[, notes := ""]
    output[1, notes := GetPlotTitle(posterior.niter = sum(in.bounds))]
    output[2, notes := paste0("User's Prior Projection ", ifelse(upp.in.bounds, "accepted", "rejected"))]
    output.list[[j]] <- cbind(output.list[[j]], output)
  }
  output.list$all.inputs = all.inputs.str
  filestr.out <- paste0(filestr, ".xlsx")
  openxlsx::write.xlsx(output.list, file = filestr.out)
  cat("\nExcel output: ", filestr.out, "\n")
  return(output.list)
}

PlotHist <- function(x, posterior.title, cap, xlab, in.bounds) {
  if (uniqueN(x) > 1) {
    g.list <- list()
    for (prior in c(T, F)) {
      if (prior) {
        d <- data.table(x)
        breaks <- NULL
        title1 <- paste0("Prior Distribution, niter = ", length(x))
      } else {
        d <- data.table(x = x[in.bounds])
        breaks <- unlist(ggplot_build(g.list[[1]])$data[[1]][, c("xmin", 'xmax')])
        breaks <- sort(unique(round(breaks, digits = 8))) #round to fix problems with nearly identical breaks
        title1 <- posterior.title
      }
      g <- ggplot(d, aes(x=x))
      if (class(x) %in%c("logical", "integer", "Date", "factor") || all(x == round(x))) {
        g <- g + geom_bar()
      } else {
        g <- g + geom_histogram(aes(y=..density..),     
                                #   binwidth=.1,
                                bins = 20,
                                breaks = breaks,
                                color="black", fill="white") +
          geom_density(alpha=.2, fill="steelblue3") 
      }
      g <- g + labs(title = title1, caption = cap) + xlab(xlab)
      g.list[[length(g.list) + 1]] <- g
    }
  } else {
    g.list <- NULL
    cat("No histogram for ", xlab, ", constant value = ", as.character(unique(x)), "\n", sep = "")
  }
  return(g.list)
}

GetProjectionPlot <- function(short.term, niter, hosp.quantiles, bounds.without.multiplier, bounds.labels, plot.observed.data) {
  if (short.term) {
    date.breaks <- "1 week"
    max.date <- bounds.without.multiplier[!is.na(lower), max(date)] + 3
    title1 <- "Short Term Hospitalization Projections"
  } else {
    date.breaks <- "1 month"
    max.date <- max(hosp.quantiles$date) - 3 #makes it look a little nicer when ending on first of the month
    title1 <- "Long Term Hospitalization Projections"
  }
  
  if (niter < 500) {
    subtitle1 <- GetPlotTitle(niter)
  } else {
    subtitle1 <- NULL
  }
  
  plot.cred.int <- niter > 1
  plot.upp <- niter < 100
 
  min.date <- min(bounds.without.multiplier$date) - 7 
  obs.size <- 3
  lb <- bounds.labels[1]
  ub <- bounds.labels[2]
  
  dt.plot <- merge(hosp.quantiles, bounds.without.multiplier, all.x = T, by = "date")
  dt.plot <- dt.plot[date >= min.date & date <= max.date]
  gg <- ggplot(dt.plot, aes(x=date)) +    
    theme_light()
  
  if (niter > 0) {
    gg <- gg + geom_line(aes(y = `50%`, color = "Median"))
  }
  
  if (plot.observed.data) {
    gg <- gg + 
      geom_point(aes(y=upper, color = ub), size = obs.size, shape = 4, na.rm = T) +
      geom_point(aes(y=lower, color = lb), size = obs.size, shape = 4, na.rm = T)
    over.aes <- list(shape = c(if (niter > 0) NA, 4, 4), linetype = c(if (niter > 0) 1, 0, 0))
  } else {
    over.aes <- list()
  }
  
  if (plot.cred.int) {
    gg <- gg + geom_ribbon(aes(ymin=`25%`, ymax=`75%`, alpha = "25%-75%"), fill = "blue") +
      geom_ribbon(aes(ymin=`15%`, ymax=`85%`, alpha = "15%-85%"), fill = "blue") +
      geom_ribbon(aes(ymin=`5%`, ymax=`95%`, alpha = "5%-95%"), fill = "blue") 
  }
  
  if (plot.upp) {
    gg <- gg + geom_line(aes(y = upp), color = "yellow1", size = 1, linetype = "dashed")  + 
      theme(panel.background = element_rect(fill = "grey80")) + 
      annotate("text", x = median(dt.plot$date), y = 10, label = "Yellow line is User's Prior Projection (column E of input)")
  } 
  
  gg <- gg + 
    xlab("") + 
    ylab("Number of COVID19 Patients in Hospital") +
    labs(title = title1, subtitle = subtitle1) + 
    scale_color_manual("", values = c("blue", "palegreen4", "red4"), breaks = c("Median", lb, ub)) +
    scale_alpha_manual("", values = c(0.2, 0.3, 0.4), breaks = c("5%-95%", "15%-85%", "25%-75%")) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + 
    scale_x_date(date_breaks = date.breaks, date_labels = "%b %d", expand = expansion()) +
    guides(color = guide_legend("", order = 1, override.aes = over.aes), alpha = guide_legend("", order = 2)) +
    theme(legend.position = "bottom", plot.margin = margin(0.1, 0.2, 0, 0.1, "in"),
          axis.text.x=element_text(hjust = 0.7))
  return(gg)
}

GetPdfOutput <- function(hosp, in.bounds, all.params, filestr, bounds.without.multiplier, internal.args) {
  post.niter <- sum(in.bounds)
  posterior.title <- GetPlotTitle(post.niter)
  filestr.out <- paste0(filestr, ".pdf")
  grDevices::pdf(file = filestr.out, width = 9.350, height = 7.225)
  
  bounds.labels <- c(internal.args$lower.bound.label, internal.args$upper.bound.label)
  print(short.term <- GetProjectionPlot(short.term = T, niter = post.niter, hosp.quantiles = hosp, bounds.without.multiplier = bounds.without.multiplier, bounds.labels = bounds.labels, plot.observed.data = internal.args$plot.observed.data.short.term))
  print(long.term <- GetProjectionPlot(short.term = F, niter = post.niter, hosp.quantiles = hosp, bounds.without.multiplier = bounds.without.multiplier, bounds.labels = bounds.labels, plot.observed.data = internal.args$plot.observed.data.long.term))
  
  intervention.multiplier.str <- grep("^intervention[[:digit:]]+\\.multiplier$", names(all.params), value = T)
  captions <- list()
  for (int.num in 0:length(intervention.multiplier.str)) {
    cumRe <- all.params[, r0.initial]
    cumRe.name <- paste0(".Re.", int.num)
    captions[[cumRe.name]] <- paste0(cumRe.name, " = r0.initial")
    for (i in seq_len(int.num)) {
      cur.multiplier <- all.params[[paste0("intervention", i, ".multiplier")]]
      stopifnot(!is.null(cur.multiplier))
      cumRe <- cumRe * cur.multiplier
      captions[[cumRe.name]] <- paste0(captions[[cumRe.name]], " * intervention", i, ".multiplier")
    }
    all.params[[cumRe.name]] <- cumRe
  }
  suppressWarnings(all.params[, model := factor(paste0(as.integer(latent.period > 0), as.integer(patients.in.hosp.are.infectious), as.integer(use.hosp.rate)))]) #suppress "Invalid .internal.selfref detected and fixed by ..." (not important)
  captions$model <- "(HasE  InfectInHosp   RateToHosp)"
  if (post.niter >= 1) {
    for (param.name in sort(names(all.params[,-"r0.initial"]))) {
      g.list <- PlotHist(all.params[[param.name]], posterior.title, captions[[param.name]], param.name, in.bounds)
      sapply(g.list, print)
    }
  }
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
