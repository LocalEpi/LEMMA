#' Run Credibility Interval based on Excel inputs
#'
#' @param input.file1 A .xlsx file for population1
#' @param input.file2 A .xlsx file for population2
#' @return list with all outputs (invisible)
#' @export
CredibilityIntervalFromExcel <- function(input.file1, input.file2) {
  cred.int <- TwoPop(input.file1, input.file2)
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
  options(width = 300)
  all.inputs.str <- utils::capture.output(print(sheets))
  options(width = prev.width)
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

ProcessSheets <- function(sheets, path) {
  # seir_inputs <- list()
  params <- sheets$`Parameters with Distributions`[, .(name = internal.name, mu = Mean, sigma = `Standard Deviation`)]
  frac_pui <- sheets$`PUI Details`[, .(name = internal.name, mu = Mean, sigma = `Standard Deviation`)]

  model.inputs <- TableToList(sheets$`Model Inputs`)
  internal.args <- TableToList(sheets$Internal)
  interventions <- sheets$Interventions
  obs.data <- sheets$Data

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



