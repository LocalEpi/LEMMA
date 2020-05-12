ReadExcel <- function(path, col_types, sheet, range = NULL) {
  x <- as.data.table(readxl::read_excel(path = path, col_types = col_types, sheet = sheet, range = range))
  attr(x, "sheetname") <- sheet
  return(x)
}

TableToList <- function(x) {
  names1 <- x$internal.name
  values <- x$value
  names(values) <- names1
  return(as.list(values))
}

DistToList <- function(x) {
  names1 <- x$internal.name
  values <- list()
  for (i in 1:nrow(x)) {
    xx <- x[i]
    values[[i]] <- lapply(x[i][, c("low", "midlow", "mid", "midhigh", "high")], function (z) unlist(z[[1]]))
  }
  names(values) <- names1
  lst <- as.list(values)
  return(lst)
}

Unlist <- function(zz) {
  #regular 'unlist' causes problems with dates
  do.call(c, zz)
}

SampleParam <- function(p, probs, niter) {
  p <- Unlist(p)
  if (class(p) == "logical") {
    return(sample(p, size = niter, replace = T, prob = probs))
  }
  discrete <- sample(p, size = 1000, replace = T, prob = probs)
  cont <- rnorm(niter, mean = mean(discrete), sd = sd(discrete))
  if (class(p) == "Date") {
    x <- as.Date(round(cont), origin = "1970-01-01")
  } else if (class(p) == "numeric") {
    x <- pmax(0, cont) #all parameters are >= 0
    if (all(p == round(p))) {
      x <- round(x) #if inputs are integers, return integers
      if (all(p > 0)) {
        x <- pmax(1, x) #don't return 0 unless it was an input
      }
    } else {
      x <- pmax(min(p) / 10, x) #don't return 0 unless it was an input
    }
  } else {
    stop("unexpected class in SampleParam")
  }
  return(x)
}

GetParams <- function(param.dist, niter, get.best.guess, N, make.pops.same) {
  probs <- unlist(param.dist$parameter.weights)
  stopifnot(sum(probs) == 1)

  weight.labels <- param.dist$weight.labels
  param.dist1 <- param.dist
  param.dist1$parameter.weights <- param.dist1$weight.labels <-  NULL
  if (get.best.guess) {
    stopifnot(weight.labels$mid == "Best Guess")
    params <- lapply(param.dist1, function (z) z$mid)
  } else {
    params <- lapply(param.dist1, SampleParam, niter = niter, probs = probs)
  }
  params <- as.data.table(params)
  params[, exposed.to.hospital := latent.period + infectious.to.hospital]
  params$infectious.to.hospital <- NULL
  
  if (make.pops.same) {
    #for testing that same pops with uniforming mixing = one pop 
    params[, r0.initial.2 := r0.initial.1]
    params[, prop.hospitalized.2 := prop.hospitalized.1]
    params[, hosp.length.of.stay.2 := hosp.length.of.stay.1]
    params[, intervention1.multiplier.2 := intervention1.multiplier.1]
    params[, intervention2.multiplier.2 := intervention2.multiplier.1]
    params[, m1 := 0]
    params[, m2 := 0]
  }
  
  GetCrossTermCoefs <- function(niter) {
    x <- matrix(runif(4 * niter), nrow = niter, ncol = 4)
    x <- x / rowSums(x)
    return(x)
  }
  
  cross.term.coefs <- as.data.table(cbind(GetCrossTermCoefs(niter), GetCrossTermCoefs(niter)))
  setnames(cross.term.coefs, paste0("cross.term.coefs.", c("a", "b", "c", "d"), rep(1:2, each=4)))
  params <- cbind(params, cross.term.coefs)

  return(params)
}

ReadInputs <- function(path, make.pops.same = F) {
  sheets <- list(ReadExcel(path, col_types = c("text", "text", "list", "list", "list", "list", "list", "skip"), sheet = "Parameters with Distributions"), #assumes text under last row has been moved
                 ReadExcel(path, col_types = c("text", "text", "list"), sheet = "Model Inputs"),
                 ReadExcel(path, col_types = c("date", "numeric", "numeric", "skip"), sheet = "Hospitilization Data"),
                 ReadExcel(path, col_types = c("text", "list", "skip"), sheet = "Internal"))
  names(sheets) <- sapply(sheets, function (z) attr(z, "sheetname"))

  sheets <- rapply(sheets, as.Date, classes = "POSIXt", how = "replace") #convert dates

  param.dist <- DistToList(sheets$`Parameters with Distributions`)
 
  # for (i in names(temp.scale)) {
  #   stopifnot(i %in% names(param.dist))
  #   for (j in names(param.dist[[i]])) {
  #     param.dist[[i]][[j]] <- param.dist[[i]][[j]] * temp.scale[[i]]
  #   }
  # }
  
  model.inputs <- TableToList(sheets$`Model Inputs`)
  stopifnot(!is.null(model.inputs$total.population1), !is.null(model.inputs$total.population2))
  model.inputs$total.population <- c(model.inputs$total.population1, model.inputs$total.population2)
  model.inputs$total.population1 <- model.inputs$total.population2 <- NULL
  
  hosp.data <- sheets$`Hospitilization Data`
  internal <- TableToList(sheets$Internal)

  observed.data <- hosp.data[, .(date = Date, hosp = (LowerBound + UpperBound) / 2)] #TODO: make this more flexible?
  if (!is.na(internal$min.obs.date.to.fit)) {
    observed.data <- observed.data[date >= internal$min.obs.date.to.fit]
  }
  if (!is.na(internal$max.obs.date.to.fit)) {
    observed.data <- observed.data[date <= internal$max.obs.date.to.fit]
  }

  hosp.bounds <- hosp.data[, .(date = Date, lower = internal$lower.bound.multiplier * LowerBound, upper = internal$upper.bound.multiplier * UpperBound)]

  set.seed(internal$random.seed)
  if (internal$random.seed == 63613) set.seed(NULL)
  
  params <- GetParams(param.dist, internal$main.iterations, get.best.guess = F, N = model.inputs$total.population, make.pops.same)
  best.guess <- GetParams(param.dist, internal$main.iterations, get.best.guess = T, N = model.inputs$total.population, make.pops.same)
  if (is.na(internal$output.filestr)) {
    internal$output.filestr <- sub(".xlsx", " output", path, fixed = T)
  }
  sheets$time.of.run <- as.character(Sys.time())
  return(list(all.params = params, model.inputs = model.inputs, hosp.bounds = hosp.bounds, observed.data = observed.data, internal.args = internal, best.guess.params = best.guess, excel.input = sheets))
}

#' Run Credibility Interval based on Excel inputs
#'
#' @param input.file A .xlsx file
#' @return NULL
#' @export
CredibilityIntervalFromExcel <- function(input.file) {
  inputs <- ReadInputs(input.file)

  CredibilityInterval(all.params = inputs$all.params, model.inputs = inputs$model.inputs, hosp.bounds = inputs$hosp.bounds, best.guess.params = inputs$best.guess.params, observed.data = inputs$observed.data, internal.args = inputs$internal.args, extras = inputs$excel.input)
  cat("\n\nDone")
}
