#' @import ggplot2
#' @import data.table
#' @import matrixStats
#' @importFrom stats median optimize runif
#' @importFrom graphics barplot

#TODO: show observed even if only fitting initial to some

GetModelName <- function(dt) {
  model.name <- dt[, factor(paste0(as.integer(hasE), as.integer(hospInf), as.integer(hospRate)))]
  return(model.name)
}

#TRUE if at least required.in.bounds fraction of vector x is in bounds
#if both lower and upper are NA, ignore that bound
#if one of lower/upper is NA but not the other, error
InBounds <- function(x, bounds, required.in.bounds) {
  stopifnot(bounds[xor(is.na(lower), is.na(upper)), .N] == 0)
  if (is.vector(x)) {
    x <- as.matrix(x)
  }
  in.bounds <- x >= bounds$lower & x <= bounds$upper
  return(colMeans(in.bounds, na.rm = T) >= required.in.bounds)
}

RunSim1 <- function(params1, model.inputs, observed.data, internal.args, date.range) {
  if (!internal.args$show.progress) {
    sink("CredibilityInterval progress log.txt")
  }
  sim <- RunSim(total.population = model.inputs$total.population, observed.data = observed.data, start.date = internal.args$simulation.start.date, end.date = model.inputs$end.date, params = params1, search.args = list(max.iter = internal.args$search.max.iter, expander = internal.args$search.expander, num.init.exp = internal.args$search.num.init.exp, max.nonconverge = internal.args$max.nonconverge))
  if (!internal.args$show.progress) {
    sink()
  }
  if (nrow(params1) == 1) {
    sim <- sim[date %in% date.range]
  } else {
    sim <- lapply(sim, function (z) z[as.character(date.range), ])
  }
  return(sim)
}

GetPlotTitle <- function(posterior.niter) {
  if (posterior.niter < 100) {
    warn.str <- "\nVery low number of iterations, DO NOT use for inference"
  } else if (posterior.niter < 1000) {
    warn.str <- "\nLow number of iterations, use with caution"
  } else {
    warn.str <- ""
  }
  paste0("Posterior Distribution, niter = ", posterior.niter, warn.str)
}

GetExcelOutput <- function(sim, best.guess, in.bounds, best.guess.in.bounds, date.range, filestr, all.inputs.str) {
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
    output <- data.table(date = date.range, quant1, bestguess = best.guess[[j]], quant2)
    output[, notes := ""]
    output[1, notes := GetPlotTitle(posterior.niter = sum(in.bounds))]
    output[2, notes := paste0("bestguess ", ifelse(best.guess.in.bounds, "accepted", "rejected"))]
    output.list[[j]] <- cbind(output.list[[j]], output)
  }
  output.list$all.inputs = all.inputs.str
  filestr.out <- paste0(filestr, ".xlsx")
  openxlsx::write.xlsx(output.list, file = filestr.out)
  cat("\nExcel output: ", filestr.out, "\n")
  return(output.list)
}

PlotHist <- function(x, posterior.title, sub, xlab, in.bounds) {
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
      g <- g + labs(title = title1, caption = sub) + xlab(xlab)
      print(g)
      g.list[[length(g.list) + 1]] <- g
    }
  } else {
    cat("No histogram for ", xlab, ", constant value = ", as.character(unique(x)), "\n", sep = "")
  }
  return(NULL)
}

GetPdfOutput <- function(hosp, in.bounds, all.params, filestr, bounds.without.multiplier, internal.args) {
  posterior.title <- GetPlotTitle(posterior.niter = sum(in.bounds))
  filestr.out <- paste0(filestr, ".pdf")
  grDevices::pdf(file = filestr.out, width = 9.350, height = 7.225)
  dt.plot <- merge(hosp, bounds.without.multiplier, all.x = T, by = "date")
  gg <- ggplot(dt.plot, aes(x=date)) +
    xlab("Date") + 
    ylab("Hospitalizations") +
    geom_line(aes(y = bestguess, color = "Best Guess")) +
    labs(title = posterior.title, caption = if (internal.args$include.plot.caption) 'Upper Bound and Lower Bound are from "Hospitalization Data" sheet in Excel input') +
    scale_color_manual("Projections", values = c("red", "yellow"), breaks = c("Median", "Best Guess")) +
    scale_alpha_manual("Range", values = c(0.2, 0.3, 0.4), breaks = c("5%-95%", "15%-85%", "25%-75%"))
  
  if (internal.args$plot.observed.data) {
    gg <- gg + geom_point(aes(y=upper, shape = "Upper Bound"), fill = "black", na.rm = T) +
      geom_point(aes(y=lower, shape = "Lower Bound"), fill = "black", na.rm = T) +
      scale_shape_manual("Data", values = c("triangle filled", "triangle down filled"), breaks = c( "Upper Bound", "Lower Bound")) 
  }
    
  if (sum(in.bounds) >= 1) {
    gg <- gg + geom_ribbon(aes(ymin=`25%`, ymax=`75%`, alpha = "25%-75%")) +
      geom_ribbon(aes(ymin=`15%`, ymax=`85%`, alpha = "15%-85%")) +
      geom_ribbon(aes(ymin=`5%`, ymax=`95%`, alpha = "5%-95%")) +
      geom_line(aes(y = `50%`, color = "Median"))
  }
  print(gg)
  
  if (sum(in.bounds) >= 1) {
    for (param.name in c("currentRe", names(all.params), "model")) {
      sub <- NULL
      if (param.name == "model") {
        cur.param <- GetModelName(all.params[, .(hasE = latent.period > 0, hospInf = patients.in.hosp.are.infectious, hospRate = use.hosp.rate)])
        sub <- "(HasE  InfectInHosp   RateToHosp)"
      } else if (param.name == "currentRe") {
        #TODO: this should be based on the current date vs the dates of the multipliers
        cur.param <- all.params[, r0.initial * intervention1.multiplier * intervention2.multiplier] #note: doesn't include int_mult3 
      } else {
        cur.param <- all.params[[param.name]]
      }
      
      PlotHist(cur.param, posterior.title, sub, param.name, in.bounds)
    }
  }
  grDevices::dev.off()
  cat("\nPDF output: ", filestr.out, "\n")
  return(gg)
}

TestOutputFile <- function(filestr) {
  z <- file.create(filestr) #filestr has no extension, it's not the actual .pdf or .xlsx output
  if (!z) {
    stop("Unable to write to output file: ", filestr, "\nCheck your working directory and the value of internal.args$output.filestr")
  }
  file.remove(filestr) 
}

#` Main function to calculate credibility interval
CredibilityInterval <- function(all.params, model.inputs, hosp.bounds, best.guess.params, observed.data, internal.args, extras) {
  options("openxlsx.numFmt" = "0.0")
  devlist <- grDevices::dev.list()
  sapply(devlist[names(devlist) == "pdf"], grDevices::dev.off) #shuts down any old pdf (if there was a crash part way)
  sapply(seq_len(sink.number()), sink, file=NULL) #same for sink
  
  all.inputs.str <- utils::capture.output(print(sapply(ls(), function(z) get(z)))) #I'm sure there's a better way to do this
  rm(extras) #extra is only used to save extra information to output file
  
  filestr <- paste0(internal.args$output.filestr, if (internal.args$add.timestamp.to.filestr) date() else "")
  TestOutputFile(filestr)
  
  date.range <- seq(model.inputs$start.display.date, model.inputs$end.date, by = "day")
  
  bounds.without.multiplier <- merge(data.table(date = date.range), hosp.bounds, all.x = T)
  bounds.with.multiplier <- copy(bounds.without.multiplier)
  bounds.with.multiplier[, lower := internal.args$lower.bound.multiplier * lower]
  bounds.with.multiplier[, upper := internal.args$upper.bound.multiplier * upper]
  
  cat("\nProjecting single Best Guess scenario:\n")
  best.guess.sim <- RunSim1(params1 = best.guess.params, model.inputs = model.inputs, observed.data = observed.data, internal.args = internal.args, date.range = date.range)
  
  best.guess.in.bounds <- InBounds(best.guess.sim$hosp, bounds.with.multiplier, required.in.bounds = internal.args$required.in.bounds)
  if (!best.guess.in.bounds) {
    cat("best.guess$hosp is not compatible with bounds\n")
    dt.print <- cbind(best.guess.sim[, .(best.guess.hosp = round(hosp, 1))], bounds.with.multiplier)
    dt.print[, OK := best.guess.hosp >= lower & best.guess.hosp <= upper]
    print(dt.print[!is.na(lower) & !is.na(upper)])
  }
  
  if (nrow(all.params) > 1) {
    cat("\nProjecting", nrow(all.params), "scenarios based on priors:\n")
    sim <- RunSim1(params1 = all.params, model.inputs = model.inputs, observed.data = observed.data, internal.args = internal.args, date.range = date.range)
    in.bounds <- InBounds(sim$hosp, bounds.with.multiplier, required.in.bounds = internal.args$required.in.bounds)
  } else {
    sim <- NULL
    in.bounds <- FALSE
  }
  
  output.list <- GetExcelOutput(sim, best.guess.sim, in.bounds, best.guess.in.bounds, date.range, filestr, all.inputs.str)
  gplot <- GetPdfOutput(hosp = output.list$hosp, in.bounds, all.params, filestr, bounds.without.multiplier, internal.args)
  return(list(sim = sim, gplot = gplot, output.list = output.list, best.guess.sim = best.guess.sim, in.bounds = in.bounds, best.guess.in.bounds = best.guess.in.bounds, date.range = date.range, filestr = filestr, all.inputs.str = all.inputs.str))
}


