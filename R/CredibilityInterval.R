#' @import data.table
#' @import matrixStats

# upp = User's Prior Projection


#TRUE if at least required.in.bounds fraction of matrix x (dates x params) is in bounds
#if both lower and upper are NA, ignore that bound
#if one of lower/upper is NA but not the other, error
#used for a single data type (hosp, active.cases, etc)
InBounds1 <- function(x, bounds.obj) {
  stopifnot(bounds.obj$bounds[xor(is.na(lower), is.na(upper)), .N] == 0)
  if (is.vector(x)) {
    x <- as.matrix(x)
  }
  stopifnot(!is.null(rownames(x)))
  x <- x[as.character(bounds.obj$bounds$date), ]
  in.bounds <- x >= (bounds.obj$bounds$lower * bounds.obj$lower.bound.multiplier) & 
    x <= (bounds.obj$bounds$upper * bounds.obj$upper.bound.multiplier)
  return(colMeans(in.bounds, na.rm = T) >= bounds.obj$required.in.bounds)
}

InBounds <- function(sim, bounds.list) {
  stopifnot(is.matrix(sim[[1]]))
  num.param.sets <- ncol(sim[[1]]) #should be the same for all sim
  in.bounds <- rep(T, num.param.sets)
  for (i in names(bounds.list)) {
    in.bounds <- in.bounds & InBounds1(sim[[i]], bounds.list[[i]])
  }
  return(in.bounds)
}

RunSim1 <- function(params1, model.inputs, internal.args, bounds.list) {
  if (!internal.args$show.progress) {
    sink("CredibilityInterval progress log.txt")
  }
  observed.data <- bounds.list[[1]]$bounds[, .(date)]
  for (i in names(bounds.list)) {
    temp.dt <- bounds.list[[i]]$bounds[, .(temp = (lower + upper) / 2)]
    setnames(temp.dt, "temp", i)
    observed.data <- cbind(observed.data, temp.dt)
  }
  if (!is.na(internal.args$min.obs.date.to.fit)) {
    observed.data <- observed.data[date >= internal.args$min.obs.date.to.fit]
  }
  if (!is.na(internal.args$max.obs.date.to.fit)) {
    observed.data <- observed.data[date <= internal.args$max.obs.date.to.fit]
  }
  
  sim <- RunSim(total.population = model.inputs$total.population, observed.data = observed.data, start.date = internal.args$simulation.start.date, end.date = model.inputs$end.date, params = params1, search.args = list(max.iter = internal.args$search.max.iter, expander = internal.args$search.expander, num.init.exp = internal.args$search.num.init.exp, max.nonconverge = internal.args$max.nonconverge))
  if (!internal.args$show.progress) {
    sink()
  }
  return(sim)
}

GetQuantiles <- function(sim, in.bounds) {
  quantile.list <- lapply(sim, function (z) {
    sim.accepted <- z[, in.bounds, drop = F]
    rowQuantiles(sim.accepted, probs = seq(0, 1, by = 0.05))
  })
  attr(quantile.list, "posterior.niter") <- sum(in.bounds)
  return(quantile.list)
}

SmoothBounds <- function(bounds.list) {
  for (i in names(bounds.list)) {
    span <- bounds.list[[i]]$loess.span
    
    bounds <- copy(bounds.list[[i]]$bounds)
    bounds[, orig.lower := lower]
    bounds[, orig.upper := upper]
    
    if (span > 0) {
      bounds[, date.index := 1:.N]
      bounds$lower <- predict(loess(lower ~ date.index, data = bounds, span = span, na.action = na.exclude), newdata = bounds)
      bounds$upper <- predict(loess(upper ~ date.index, data = bounds, span = span, na.action = na.exclude), newdata = bounds)
    }
    
    gg <- ggplot(bounds, aes(x = date)) +
      geom_point(aes(y=orig.upper), color = "red4", shape = 4, na.rm = T) +
      geom_point(aes(y=orig.lower), color = "palegreen4", shape = 4, na.rm = T) +
      ylab(bounds.list[[i]]$long.name) +
      ggtitle(bounds.list[[i]]$long.name)
    
    if (span > 0) {
      gg <- gg +  
        geom_line(aes(y = lower), na.rm = T) +
        geom_line(aes(y = upper), na.rm = T) +
        labs(subtitle = paste("Pre and Post Data Smoothing (loess.span = ", span, ")"))
      bounds.list[[i]]$bounds <- bounds[, .(date, lower, upper)]
    } else {
      gg <- gg + labs(subtitle = "No Smoothing (loess.span = 0)")
    }
    print(gg)
  }

  return(bounds.list)
} 

#` Main function to calculate credibility interval
CredibilityInterval <- function(all.params, model.inputs, bounds.list, upp.params, internal.args, extras) {
  options("openxlsx.numFmt" = "0.0")
  devlist <- grDevices::dev.list()
  sapply(devlist[names(devlist) == "pdf"], grDevices::dev.off) #shuts down any old pdf (if there was a crash part way)
  sapply(seq_len(sink.number()), sink, file=NULL) #same for sink
  
  all.inputs.str <- utils::capture.output(print(sapply(ls(), function(z) get(z)))) #I'm sure there's a better way to do this
  rm(extras) #extra is only used to save extra information to output file
  
  filestr <- paste0(internal.args$output.filestr, if (internal.args$add.timestamp.to.filestr) date() else "")
  TestOutputFile(filestr)
  
  bounds.list <- SmoothBounds(bounds.list)
  
  cat("\nProjecting single User's Prior Projection scenario:\n")
  upp.sim <- RunSim1(params1 = upp.params, model.inputs = model.inputs, internal.args = internal.args, bounds.list = bounds.list)
  if (nrow(all.params) > 1) {
    cat("\nProjecting", nrow(all.params), "scenarios based on priors:\n")
    sim <- RunSim1(params1 = all.params, model.inputs = model.inputs, internal.args = internal.args, bounds.list = bounds.list)
    in.bounds <- InBounds(sim, bounds.list)
    posterior.quantiles <- GetQuantiles(sim[union(c("hosp", "icu", "vent", "active.cases", "total.cases"), names(bounds.list))], in.bounds)
    excel.output <- GetExcelOutput(posterior.quantiles, model.inputs, filestr, all.inputs.str)
  } else {
    sim <- posterior.quantiles <- excel.output <- NULL
    in.bounds <- FALSE
  }
  gplot <- GetPdfOutput(posterior.quantiles, in.bounds, all.params, filestr, bounds.list, internal.args, model.inputs, upp.sim)
  return(list(sim = sim, gplot = gplot, excel.output = excel.output, upp.sim = upp.sim, in.bounds = in.bounds, filestr = filestr, all.inputs.str = all.inputs.str))
}


