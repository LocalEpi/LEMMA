#' @import data.table
#' @import matrixStats

# upp = User's Prior Projection

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


#` Main function to calculate credibility interval
CredibilityInterval <- function(all.params, model.inputs, hosp.bounds, upp.params, observed.data, internal.args, extras) {
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
  
  cat("\nProjecting single User's Prior Projection scenario:\n")
  upp.sim <- RunSim1(params1 = upp.params, model.inputs = model.inputs, observed.data = observed.data, internal.args = internal.args, date.range = date.range)
  
  upp.in.bounds <- InBounds(upp.sim$hosp, bounds.with.multiplier, required.in.bounds = internal.args$required.in.bounds)
  if (!upp.in.bounds) {
    cat("\nnote: User's Prior Projection is not compatible with bounds\n")
    dt.print <- cbind(upp.sim[, .(upp.hosp = round(hosp, 1))], bounds.with.multiplier)
    dt.print[, OK := upp.hosp >= lower & upp.hosp <= upper]
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
  
  output.list <- GetExcelOutput(sim, upp.sim, in.bounds, upp.in.bounds, date.range, filestr, all.inputs.str)
  gplot <- GetPdfOutput(hosp = output.list$hosp, in.bounds, all.params, filestr, bounds.without.multiplier, internal.args)
  return(list(sim = sim, gplot = gplot, output.list = output.list, upp.sim = upp.sim, in.bounds = in.bounds, upp.in.bounds = upp.in.bounds, date.range = date.range, filestr = filestr, all.inputs.str = all.inputs.str))
}


