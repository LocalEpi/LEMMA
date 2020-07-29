#' Run Credibility Interval based on Excel inputs
#'
#' @param input.file1 A .xlsx file for population1
#' @param input.file2 A .xlsx file for population2
#' @return list with all outputs (invisible)
#' @export
CredibilityIntervalFromExcel <- function(input.file1, input.file2) {
  inputs1 <- GetInputs(input.file1)
  inputs2 <- GetInputs(input.file2)
  cred.int <- TwoPop(inputs1, inputs2)
  invisible(cred.int)
}

ReadExcel <- function(path, col_types = "guess", sheet, ...) {
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

ToString <- function(sheets) {
  #Make a human readable string from the raw Excel input

  sheets$time.of.run <- as.character(Sys.time())
  sheets$LEMMA.version <- getNamespaceVersion("LEMMA")

  prev.width <- getOption("width")
  prev.print.nrows <- getOption("datatable.print.nrows")
  options(width = 300, datatable.print.nrows = 300)
  all.inputs.str <- utils::capture.output(print(sheets))
  options(width = prev.width, datatable.print.nrows = prev.print.nrows)
  all.inputs.str <- c("NOTE: set font to Courier to read", all.inputs.str)
  return(all.inputs.str)
}

ReadInputs <- function(path) {
  sheets <- list(ReadExcel(path, sheet = "Parameters with Distributions"),
                 ReadExcel(path, col_types = c("text", "text", "list"), sheet = "Model Inputs"),
                 ReadExcel(path, sheet = "Interventions", skip = 2),
                 ReadExcel(path, sheet = "Data", skip = 3),
                 ReadExcel(path, sheet = "PUI Details", skip = 4),
                 ReadExcel(path, col_types = c("text", "list", "skip", "skip"), sheet = "Internal"))
  names(sheets) <- sapply(sheets, function (z) attr(z, "sheetname"))
  sheets <- rapply(sheets, as.Date, classes = "POSIXt", how = "replace") #convert dates
  return(sheets)
}

AddInterventions <- function(interventions, max.date) {
  interval <- 14

  min.date <- min(interventions$mu_t_inter)
  d.set <- seq(max.date, min.date, by = "-1 day")
  for (i in seq_along(d.set)) {
    d <- d.set[i]
    if (min(abs(as.numeric(interventions$mu_t_inter - d))) > interval) {
      if (i == 1) {
        sd1 <- 0.1
      } else {
        sd1 <- 0.3
      }
      new.int <- data.table(mu_t_inter = d, sigma_t_inter = 2, mu_beta_inter = 1, sigma_beta_inter = sd1, mu_len_inter = 7, sigma_len_inter = 2)
      interventions <- rbind(interventions, new.int)
    }
  }
  setkey(interventions, "mu_t_inter")
  return(interventions)
}

ProcessSheets <- function(sheets, path) {
  # seir_inputs <- list()
  params <- sheets$`Parameters with Distributions`[, .(name = internal.name, mu = Mean, sigma = `Standard Deviation`)]
  params[, sigma := pmax(sigma, mu / 100)] #Stan crashes if sigma = 0
  frac_pui <- sheets$`PUI Details`[, .(name = internal.name, mu = Mean, sigma = `Standard Deviation`)]
  frac_pui[, sigma := pmax(sigma, mu / 100)]

  model.inputs <- TableToList(sheets$`Model Inputs`)
  internal.args <- TableToList(sheets$Internal)
  interventions <- sheets$Interventions
  obs.data <- sheets$Data

  all.na <- rowAlls(is.na(as.matrix(obs.data[, -"date"])))
  obs.data <- obs.data[all.na == F]

  if (is.null(internal.args$automatic.interventions)) {
    internal.args$automatic.interventions <- T #this was left off early versions
  }
  if ("skip1" %in% names(interventions)) {
    interventions$skip1 <- interventions$skip2 <- NULL #these are just notes, but not in all versions
  }
  interventions[, sigma_t_inter := pmax(sigma_t_inter, 0.01)]
  interventions[, sigma_beta_inter := pmax(sigma_beta_inter, 0.001)]
  interventions[, sigma_len_inter := pmax(sigma_len_inter, 0.01)]
  if (internal.args$automatic.interventions) {
    interventions <- AddInterventions(interventions, max.date = obs.data[, max(date)])
  }

  if (is.na(internal.args$output.filestr)) {
    internal.args$output.filestr <- sub(".xlsx", " output", path, fixed = T)
  }
  if (internal.args$add.timestamp.to.filestr) {
    internal.args$output.filestr <- paste0(internal.args$output.filestr, date())
  }
  if (is.na(internal.args$cores)) {
    internal.args$cores <- parallel::detectCores()
  }

  all.inputs.str <- ToString(sheets)
  return(list(params = params, frac_pui = frac_pui, model.inputs = model.inputs, internal.args = internal.args, interventions = interventions, obs.data = obs.data, all.inputs.str = all.inputs.str))
}



