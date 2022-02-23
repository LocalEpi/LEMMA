# --------------------------------------------------------------------------------
#   Functions:
#   1. CredibilityIntervalFromExcel
#   2. ReadExcel
#   3. TableToList
#   4. ReadInputs
#   5. AddInterventions
#   6. ProcessSheets
# --------------------------------------------------------------------------------

#' Run Credibility Interval based on Excel inputs
#'
#' @param input.file A .xlsx file
#' @return list with all outputs (invisible)
#' @export
CredibilityIntervalFromExcel <- function(input.file) {
  sheets <- ReadInputs(input.file)
  inputs <- ProcessSheets(sheets)

  cred.int <- CredibilityInterval(inputs)
  cat("\nDone\n\n")
  cat("Current LEMMA version: ", getNamespaceVersion("LEMMA"), "\n")
  invisible(cred.int)
}

ReadExcel <- function(path, col_types = "guess", sheet, ...) {
  x <- as.data.table(readxl::read_excel(path = path, col_types = col_types, sheet = sheet, ...))
  attr(x, "sheetname") <- sheet
  return(x)
}

TableToList <- function(x, value.name = NULL) {
  names1 <- x$internal.name
  if (is.null(value.name)) {
    values <- x$value
    names(values) <- names1
  } else {
    values <- x[[value.name]]
    names(values) <- paste0(names1, "_", value.name)
  }
  return(as.list(values))
}

ReadInputs <- function(path) {
  #suppress warning in Vaccine Doses about date
  sheets <- list(ReadExcel(path, sheet = "Parameters with Distributions"),
                 ReadExcel(path, col_types = c("text", "text", "list"), sheet = "Model Inputs"),
                 ReadExcel(path, sheet = "Interventions", skip = 2),
                 ReadExcel(path, sheet = "Data", skip = 3),
                 ReadExcel(path, sheet = "PUI Details", skip = 4),
                 ReadExcel(path, col_types = c("text", "list", "skip", "skip"), sheet = "Internal"))
  names(sheets) <- sapply(sheets, function (z) attr(z, "sheetname"))
  sheets <- rapply(sheets, as.Date, classes = "POSIXt", how = "replace") #convert dates

  sheets$Internal <- rbind(sheets$Internal, data.table(internal.name = "input.file", value = path))
  return(sheets)
}

ProcessSheets <- function(sheets) {
  params <- sheets$`Parameters with Distributions`[, .(name = internal.name, mu = Mean, sigma = `Standard Deviation`)]
  params[, sigma := pmax(sigma, mu / 100)] #Stan crashes if sigma = 0
  frac_pui <- sheets$`PUI Details`[, .(name = internal.name, mu = Mean)]

  model.inputs <- TableToList(sheets$`Model Inputs`)
  internal.args <- TableToList(sheets$Internal)
  interventions <- sheets$Interventions
  obs.data <- sheets$Data

  all.na <- rowAlls(is.na(as.matrix(obs.data[, -"date"])))
  obs.data <- obs.data[all.na == F]

  if (is.na(internal.args$output.filestr)) {
    internal.args$output.filestr <- sub(".xlsx", " output", internal.args$input.file, fixed = T)
  }
  if (internal.args$add.timestamp.to.filestr) {
    internal.args$output.filestr <- paste0(internal.args$output.filestr, gsub(":", "-", as.character(date()), fixed = T))
  }

  internal.args$weights <- rep(1, length(DataTypes())) #TODO - make this an input?

  return(list(params = params, interventions = interventions, frac_pui = frac_pui, model.inputs = model.inputs, internal.args = internal.args, obs.data = obs.data))
}



