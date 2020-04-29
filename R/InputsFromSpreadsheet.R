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

GetParams <- function(param.dist, niter, get.best.guess) {
  probs <- unlist(param.dist$parameter.weights)
  stopifnot(sum(probs) == 1)

  weight.labels <- param.dist$weight.labels
  param.dist1 <- param.dist
  param.dist1$parameter.weights <- param.dist1$weight.labels <-  NULL
  if (get.best.guess) {
    stopifnot(weight.labels$mid == "Best Guess")
    params <- lapply(param.dist1, function (z) z$mid)
  } else {
    params <- lapply(param.dist1, function (z) sample(Unlist(z), size = niter, replace = T, prob = probs))
  }
  params <- as.data.table(params)
  params[, exposed.to.hospital := latent.period + infectious.to.hospital]
  params$infectious.to.hospital <- NULL
  return(params)
}

ReadInputs <- function(path) {
  sheets <- list(ReadExcel(path, col_types = c("text", "text", "list", "list", "list", "list", "list", "skip"), sheet = "Parameters with Distributions", range = "A1:H22"), #don't read the whole sheet - there are hidden data validation lists below
                 ReadExcel(path, col_types = c("text", "text", "list"), sheet = "Model Inputs"),
                 ReadExcel(path, col_types = c("date", "numeric", "numeric", "skip"), sheet = "Hospitilization Data"),
                 ReadExcel(path, col_types = c("text", "list", "skip"), sheet = "Internal"))
  names(sheets) <- sapply(sheets, function (z) attr(z, "sheetname"))

  sheets <- rapply(sheets, as.Date, classes = "POSIXt", how = "replace") #convert dates

  param.dist <- DistToList(sheets$`Parameters with Distributions`)
  model.inputs <- TableToList(sheets$`Model Inputs`)
  if (!("start.display.date" %in% names(model.inputs))) {
    model.inputs$start.display.date <- as.Date("2020/3/1")
  }
  hosp.data <- sheets$`Hospitilization Data`
  internal <- TableToList(sheets$Internal)

  observed.data <- hosp.data[, .(date = Date, hosp = (LowerBound + UpperBound) / 2)] #TODO: make this more flexible?
  if (!is.na(internal$min.obs.date.to.fit)) {
    observed.data <- observed.data[date >= internal$min.obs.date.to.fit]
  }
  if (!is.na(internal$max.obs.date.to.fit)) {
    observed.data <- observed.data[date <= internal$max.obs.date.to.fit]
  }

  hosp.bounds <- hosp.data[, .(date = Date, lower = LowerBound, upper = UpperBound)]

  set.seed(internal$random.seed)
  params <- GetParams(param.dist, internal$main.iterations, get.best.guess = F)
  best.guess <- GetParams(param.dist, internal$main.iterations, get.best.guess = T)

  if (is.na(internal$output.filestr)) {
    internal$output.filestr <- sub(".xlsx", " output", path, fixed = T)
  }
  sheets$time.of.run <- as.character(Sys.time())
  sheets$LEMMA.version <- getNamespaceVersion("LEMMA")
 
  return(list(all.params = params, model.inputs = model.inputs, hosp.bounds = hosp.bounds, observed.data = observed.data, internal.args = internal, best.guess.params = best.guess, excel.input = sheets))
}

#' Run Credibility Interval based on Excel inputs
#'
#' @param input.file A .xlsx file
#' @return list with all outputs (invisible)
#' @export
CredibilityIntervalFromExcel <- function(input.file) {
  inputs <- ReadInputs(input.file)

  cred.int <- CredibilityInterval(all.params = inputs$all.params, model.inputs = inputs$model.inputs, hosp.bounds = inputs$hosp.bounds, best.guess.params = inputs$best.guess.params, observed.data = inputs$observed.data, internal.args = inputs$internal.args, extras = inputs$excel.input)
  cat("\nDone\n\n")
  cat("Current LEMMA version: ", inputs$excel.input$LEMMA.version, "\n")
  cat("LEMMA is in early development. Please reinstall from github daily.\n")
  invisible(cred.int)
}
