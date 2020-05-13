#' @import ggplot2
#' @import data.table
#' @import matrixStats
#' @importFrom stats median optimize runif
#' @importFrom graphics barplot

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
    sim$HP2 <- sim$HP.V2 
  } else {
    sim <- lapply(sim, function (z) {
      if (length(dim(z)) == 3) {
        z[as.character(date.range), , ]
      } else {
        z[as.character(date.range), ] #hack for hosp
      }
    })
    sim$HP2 <- sim$HP[, , 2] 
  }
  
  return(sim)
}

#` Main function to calculate credibility interval
CredibilityInterval <- function(all.params, model.inputs, hosp.bounds, HP2.bounds, best.guess.params, observed.data, internal.args, extras) {
  options("openxlsx.numFmt" = "0.0")
  sapply(grDevices::dev.list(), grDevices::dev.off) #shuts down any old pdf (if there was a crash part way)
  sapply(seq_len(sink.number()), sink, file=NULL) #same for sink

  all.inputs.str <- utils::capture.output(print(sapply(ls(), function(z) get(z)))) #I'm sure there's a better way to do this
  rm(extras) #extra is only used to save extra information to output file

  date.range <- seq(model.inputs$start.display.date, model.inputs$end.date, by = "day") 
  bounds <- rbind(merge(data.table(date = date.range), hosp.bounds, all.x = T),
                  merge(data.table(date = date.range), HP2.bounds, all.x = T))

  best.guess.sim <- RunSim1(params1 = best.guess.params, model.inputs = model.inputs, observed.data = observed.data, internal.args = internal.args, date.range = date.range)

  best.guess.in.bounds <- InBounds(c(best.guess.sim$hosp, best.guess.sim$HP2), bounds, required.in.bounds = internal.args$required.in.bounds)
  if (!best.guess.in.bounds) {
    cat("best.guess$hosp is not compatible with bounds\n")
    dt.print <- cbind(best.guess.sim[, .(best.guess.hosp = c(round(hosp, 1), round(HP2, 1)))], bounds)
    dt.print[, OK := best.guess.hosp >= lower & best.guess.hosp <= upper]
    print(dt.print[!is.na(lower) & !is.na(upper)])
  }

  if (nrow(all.params) > 2) {
    sim <- RunSim1(params1 = all.params, model.inputs = model.inputs, observed.data = observed.data, internal.args = internal.args, date.range = date.range)
    in.bounds <- InBounds(rbind(sim$hosp, sim$HP2), bounds, required.in.bounds = internal.args$required.in.bounds)
    
    filestr <- paste0(internal.args$output.filestr, if (internal.args$add.timestamp.to.filestr) date() else "")
    
    output.list <- GetExcelOutput(sim, best.guess.sim, in.bounds, best.guess.in.bounds, date.range, filestr, all.inputs.str)
    #this needs to get updated to master - return gplot, run GetPdfOutput even if sum(in.bounds) <=1, maybe other stuff
    if (sum(in.bounds) <= 1) {
      cat("niter = ", sum(in.bounds), " / ", length(in.bounds), "in bounds. No pdf output written.\n")
    } else {
      GetPdfOutput(hosp = output.list$hosp, in.bounds = in.bounds, all.params = all.params, filestr = filestr, bounds.without.multiplier = NULL) #need to fix bounds.without.multiplier
    }
  } else {
    sim <- in.bounds <- filestr <- output.list <- NULL
  }
  return(list(sim = sim, output.list = output.list, best.guess.sim = best.guess.sim, in.bounds = in.bounds, best.guess.in.bounds = best.guess.in.bounds, date.range = date.range, filestr = filestr, all.inputs.str = all.inputs.str))
}


