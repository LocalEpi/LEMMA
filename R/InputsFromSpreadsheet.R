ReadExcel <- function(path, col_types, sheet, ...) {
  x <- as.data.table(readxl::read_excel(path = path, col_types = col_types, sheet = sheet, ...))
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

GetParams <- function(param.dist, niter, get.upp) {
  probs <- unlist(param.dist$parameter.weights)
  stopifnot(sum(probs) == 1)

  weight.labels <- param.dist$weight.labels
  param.dist1 <- param.dist
  param.dist1$parameter.weights <- param.dist1$weight.labels <-  NULL
  params <- lapply(param.dist1, function (z) z$mid)
  params <- as.data.table(params)
  params[latent.period > 0 & latent.period < 1, latent.period := round(latent.period)]
  params[illness.length.given.nonhosp < 1, illness.length.given.nonhosp := 1]
  params[infectious.to.hospital < 1, infectious.to.hospital := 1]
  params[hosp.length.of.stay < 1, hosp.length.of.stay := 1]
  params[, exposed.to.hospital := latent.period + infectious.to.hospital]
  params$infectious.to.hospital <- NULL
  return(params)
}

ReadInputs <- function(path) {
  sheets <- list(ReadExcel(path, col_types = c("text", "text", "list", "list", "list", "list", "list", "skip"), sheet = "Parameters with Distributions"), 
                 ReadExcel(path, col_types = c("text", "text", "list"), sheet = "Model Inputs"),
                 # ReadExcel(path, col_types = c("date", "numeric", "numeric", "skip"), sheet = "Hospitilization Data"),
                 
                 ReadExcel(path, sheet = "Data", col_types = "guess", skip = 5),
                 ReadExcel(path, col_types = c("text", "list", "skip"), sheet = "Internal"))
  names(sheets) <- sapply(sheets, function (z) attr(z, "sheetname"))
  sheets$bounds.args <- ReadExcel(path, sheet = "Data", col_types = "guess", n_max = 4)
  sheets$`Parameters with Distributions` <- sheets$`Parameters with Distributions`[1:(which.max(is.na(internal.name)) - 1)]
  
  sheets <- rapply(sheets, as.Date, classes = "POSIXt", how = "replace") #convert dates
  return(sheets)
}

ProcessSheets <- function(sheets, path, generate.params = TRUE) {
  param.dist <- DistToList(sheets$`Parameters with Distributions`)
  model.inputs <- TableToList(sheets$`Model Inputs`)
  if (!("start.display.date" %in% names(model.inputs))) {
    model.inputs$start.display.date <- as.Date("2020/3/1")
  }
  #hosp.data <- sheets$`Hospitilization Data`
  #hosp.bounds <- hosp.data[, .(date = Date, lower = LowerBound, upper = UpperBound)]
  
  internal <- TableToList(sheets$Internal)
  if (!("plot.observed.data.long.term" %in% names(internal))) {
    internal$plot.observed.data.long.term <- FALSE
  }
  if (!("plot.observed.data.short.term" %in% names(internal))) {
    internal$plot.observed.data.short.term <- TRUE
  }
  if (!("lower.bound.label" %in% names(internal))) {
    internal$lower.bound.label <- "Confirmed COVID19"
  }
  if (!("upper.bound.label" %in% names(internal))) {
    internal$upper.bound.label <- "Probable COVID19"
  }
  
  
  valid.bounds.names <- SimDataTypes()
  
  data.columns <- setdiff(names(sheets$Data), "date")
  for (i in data.columns) {
    if (!sub(".lower|.upper", "", i) %in% valid.bounds.names) {
      stop("unexpected data column: ", i)
    }
  }
  bounds.args.columns <- setdiff(names(sheets$bounds.args), c("internal.name", "not.used"))
  for (i in data.columns) {
    if (!sub(".lower|.upper", "", i) %in% valid.bounds.names) {
      stop("unexpected bounds.args column: ", i)
    }
  }
  bounds.names <- sub(".lower|.upper", "", data.columns)
  bounds.list <- list()
  for (i in bounds.names) { 
    lb <- paste0(i, ".lower")
    ub <- paste0(i, ".upper")
    bounds <- sheets$Data[, .(date, lower = get(lb), upper = get(ub))]
    lower.mult <- sheets$bounds.args[internal.name == "bound.multiplier", as.numeric(get(lb))]
    upper.mult <- sheets$bounds.args[internal.name == "bound.multiplier", as.numeric(get(ub))]
    loess.span <- sheets$bounds.args[internal.name == "loess.span", as.numeric(get(lb))] #stored in the "lower" column only
    req.in.bounds <- sheets$bounds.args[internal.name == "required.in.bounds", as.numeric(get(lb))] #stored in the "lower" column only
    long.name <- sheets$bounds.args[internal.name == "long.name", get(lb)] #stored in the "lower" column only
    lower.label <- sheets$bounds.args[internal.name == "bound.label", get(lb)]
    upper.label <- sheets$bounds.args[internal.name == "bound.label", get(ub)]
    if (!(all(is.na(bounds$lower)) && all(is.na(bounds$upper)))) {
      bounds.arguments <- c(req.in.bounds, lower.mult, upper.mult, loess.span)
      # if (anyNA(bounds.arguments) || length(bounds.arguments) != 4) {
      #   print(bounds.arguments)
      #   stop("required.in.bounds, lower.bound.multiplier, upper.bound.multiplier, loess.span, must all be specified for ", i)
      # }
      bounds.list[[i]] <- list(required.in.bounds = req.in.bounds, 
                               lower.bound.multiplier = lower.mult, 
                               upper.bound.multiplier = upper.mult,
                               loess.span = loess.span,
                               long.name = long.name, 
                               lower.bound.label = lower.label,
                               upper.bound.label = upper.label,
                               bounds = bounds)
    }
  }

  set.seed(internal$random.seed)
  upp <- GetParams(param.dist, internal$main.iterations, get.upp = T)

  if (is.na(internal$output.filestr)) {
    internal$output.filestr <- sub(".xlsx", " output", path, fixed = T)
  }
  sheets$time.of.run <- as.character(Sys.time())
  sheets$LEMMA.version <- getNamespaceVersion("LEMMA")
 
  return(list(all.params = NA, model.inputs = model.inputs, bounds.list = bounds.list, internal.args = internal, upp.params = upp, excel.input = sheets, param.dist = param.dist))
}

#' Run Credibility Interval based on Excel inputs
#'
#' @param input.file A .xlsx file
#' @return list with all outputs (invisible)
#' @export
CredibilityIntervalFromExcel <- function(input.file) {
  sheets <- ReadInputs(input.file)
  inputs <- ProcessSheets(sheets, input.file)

  cred.int <- CredibilityInterval(model.inputs = inputs$model.inputs, bounds.list = inputs$bounds.list, upp.params = inputs$upp.params, internal.args = inputs$internal.args, extras = inputs$excel.input, obs.data = sheets$Data)
  cat("\nDone\n\n")
  cat("Current LEMMA version: ", inputs$excel.input$LEMMA.version, "\n")
  cat("LEMMA is in early development. Please reinstall from github daily.\n")
  invisible(c(cred.int, inputs = list(inputs)))
}

#' Run Credibility Interval based on Excel inputs
#'
#' @param input.file A .xlsx file
#' @return a ggplot object
#' @export
VaryOneParameter <- function(input.file, parameter.name = "") {
  sheets <- ReadInputs(input.file)
  inputs <- ProcessSheets(sheets, input.file, generate.params = F)
  params <- inputs$upp.params
  if (parameter.name %in% names(params)) {
    parameter.range <- Unlist(inputs$param.dist[[parameter.name]])
    params <- params[rep(1, length(parameter.range))]
    params[[parameter.name]] <- parameter.range
  } else {
    cat("parameter.name needs to be one of the following:\n")
    print(names(params))
    stop("unrecognized parameter.name")
  }
  date.range <- seq(inputs$model.inputs$start.display.date, inputs$model.inputs$end.date, by = "day")
  sim <- RunSim1(params1 = params, model.inputs = inputs$model.inputs, bounds.list = inputs$bounds.list, internal.args = inputs$internal.args)
  
  colnames(sim$hosp) <- parameter.range
  dt.plot <- melt(data.table(date = as.Date(rownames(sim$hosp)), sim$hosp), id = "date", value.name = "Hospitalizations", variable.name = "parameter")
  gplot <- ggplot(dt.plot, aes(x=date, y=Hospitalizations, group = parameter)) +
    xlab("Date") + 
    geom_line(aes(color = parameter)) +
    scale_color_hue(name = parameter.name)
   return(list(gplot = gplot, sim = sim)) 
}
