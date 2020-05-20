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

GetExcelOutput <- function(quantile.list, model.inputs, filestr, all.inputs.str) {
  display.date.range <- as.character(seq(model.inputs$start.display.date, model.inputs$end.date, by = "day"))
  output.list <- lapply(quantile.list, function (quant) {
    output <- data.table(date = display.date.range, quant[display.date.range, ], notes = "")
    output[1, notes := GetPlotTitle(posterior.niter = attr(quantile.list, "posterior.niter"))]
    return(output)
  })
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
                                bins = 20,
                                breaks = breaks,
                                color="black", fill="white") +
          geom_density(alpha=.2, fill="steelblue3") 
      }
      
      if (grepl(".Re.", xlab, fixed = T) && !prior) {
        cat("\n", cap, ":\n")
        re.quantiles <- round(quantile(d$x, c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)), 2)
        print(re.quantiles)
        subtitl <- paste0(capture.output(print(re.quantiles)), collapse = "\n")
      } else {
        subtitl <- NULL
      }
         
      g <- g + 
        labs(title = title1, subtitle = subtitl, caption = cap) + 
        xlab(xlab) + 
        theme(plot.subtitle=element_text(family="Courier"))
      g.list[[length(g.list) + 1]] <- g
    }
  } else {
    g.list <- NULL
    cat("No histogram for ", xlab, ", constant value = ", as.character(unique(x)), "\n", sep = "")
  }
  return(g.list)
}

GetYLabel <- function(compartment.name, long.name) {
  switch(compartment.name, hosp = "Number of COVID19 Patients in Hospital", 
         deaths = "Number of COVID19 Deaths",
         long.name)
}

GetProjectionPlot <- function(short.term, niter, quantiles, bounds.list, plot.observed.data, compartment.name, model.inputs, upp.sim) {
  dt.plot <- merge(data.table(date = as.Date(names(upp.sim)), upp = upp.sim), bounds.list$bounds, all.x = T, by = "date")
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
  
  if (niter < 500) {
    subtitle1 <- GetPlotTitle(niter)
  } else {
    subtitle1 <- NULL
  }
  
  plot.cred.int <- niter > 1
  plot.upp <- niter < 100
 
  if (plot.observed.data) {
    min.date <- bounds.list$bounds[!is.na(lower), min(date)] - 7 
  } else {
    min.date <- model.inputs$start.display.date - 7
  }
  
  obs.size <- 3
  lb <- bounds.list$lower.bound.label
  ub <- bounds.list$upper.bound.label
  
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
      annotate("text", x = median(dt.plot$date), y = 1, label = "Yellow line is User's Prior Projection (column E of input)")
  } 
  
  gg <- gg + 
    xlab("") + 
    ylab(GetYLabel(compartment.name, bounds.list$long.name)) +
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

GetPdfOutput <- function(quantiles, in.bounds, all.params, filestr, bounds, internal.args, model.inputs, upp.sim) {
  post.niter <- sum(in.bounds)
  posterior.title <- GetPlotTitle(post.niter)
  filestr.out <- paste0(filestr, ".pdf")
  grDevices::pdf(file = filestr.out, width = 9.350, height = 7.225)
  print(ggplot() + labs(title = paste("LEMMA", getNamespaceVersion("LEMMA")), subtitle = "This page is left blank to facilitate viewing in Mac Preview 'Two Page' mode"))
  
  short.term <- long.term <- list()
  for (i in names(bounds)) {
    print(short.term[[i]] <- GetProjectionPlot(short.term = T, niter = post.niter, quantiles = quantiles[[i]], bounds.list = bounds[[i]], plot.observed.data = internal.args$plot.observed.data.short.term, compartment.name = i, model.inputs = model.inputs, upp.sim = upp.sim[[i]]))
    print(long.term[[i]] <- GetProjectionPlot(short.term = F, niter = post.niter, quantiles = quantiles[[i]], bounds = bounds[[i]], plot.observed.data = internal.args$plot.observed.data.long.term, compartment.name = i, model.inputs = model.inputs, upp.sim = upp.sim[[i]]))
  }
  
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
