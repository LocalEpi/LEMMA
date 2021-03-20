#' Run Credibility Interval based on Excel inputs
#'
#' @param input.file A .xlsx file
#' @return list with all outputs (invisible)
#' @export
CredibilityIntervalFromExcel <- function(input.file) {
  sheets <- ReadInputs(input.file)
  inputs <- ProcessSheets(sheets, input.file)

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

TableToList <- function(x) {
  names1 <- x$internal.name
  values <- x$value
  names(values) <- names1
  return(as.list(values))
}

ReadInputs <- function(path) {
  sheets <- list(ReadExcel(path, sheet = "Parameters with Distributions"),
                 ReadExcel(path, col_types = c("text", "text", "list"), sheet = "Model Inputs"),
                 ReadExcel(path, sheet = "Interventions", skip = 2),
                 ReadExcel(path, sheet = "Data", skip = 3),
                 ReadExcel(path, sheet = "PUI Details", skip = 4),
                 ReadExcel(path, col_types = c("text", "list", "skip", "skip"), sheet = "Internal"),
                 ReadExcel(path, sheet = "VaccinesDoses", skip = 1),
                 ReadExcel(path, sheet = "VaccinesVariants", skip = 1),
                 ReadExcel(path, sheet = "VaccinesPopulation", skip = 1),
                 ReadExcel(path, col_types = c("text", "list", "skip"), sheet = "VaccinesMisc"))

  names(sheets) <- sapply(sheets, function (z) attr(z, "sheetname"))
  sheets <- rapply(sheets, as.Date, classes = "POSIXt", how = "replace") #convert dates
  return(sheets)
}

AddInterventions <- function(interventions, min.date, max.date) {
  d.set <- seq(max.date - 10, min.date, by = "-1 day")
  for (i in seq_along(d.set)) {
    d <- d.set[i]
    if (d < (max.date - 60)) {
      interval <- 30
      sigma_t_inter <- 5
      mu_len_inter <- 15
      sigma_len_inter <- 5
    } else {
      interval <- 14
      sigma_t_inter <- 2
      mu_len_inter <- 7
      sigma_len_inter <- 2
    }
    if ((length(interventions$mu_t_inter) == 0) || (min(abs(as.numeric(interventions$mu_t_inter - d))) >= interval)) {
      new.int <- data.table(mu_t_inter = d, sigma_t_inter, mu_beta_inter = 1, sigma_beta_inter = 0.1, mu_len_inter, sigma_len_inter)
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
  frac_pui <- sheets$`PUI Details`[, .(name = internal.name, mu = Mean)]

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
    if (nrow(interventions) == 0) {
      min.date <- internal.args$simulation.start.date
    } else {
      min.date <- min(interventions$mu_t_inter)
    }
    interventions <- AddInterventions(interventions, min.date = min.date, max.date = obs.data[, max(date)])
  }

  if (is.na(internal.args$output.filestr)) {
    internal.args$output.filestr <- sub(".xlsx", " output", path, fixed = T)
  }
  if (internal.args$add.timestamp.to.filestr) {
    internal.args$output.filestr <- paste0(internal.args$output.filestr, gsub(":", "-", as.character(date()), fixed = T))
  }
  internal.args$weights <- rep(1, length(DataTypes())) #TODO - make this an input?

  sheets$VaccinesMisc <- TableToList(sheets$VaccinesMisc)
  start_date <- internal.args$simulation.start.date + 1
  end_date <- model.inputs$end.date
  vaccines.list <- GetVaccineParams(doses_actual = sheets$VaccinesDoses,
                    doses_per_day_base = sheets$VaccinesMisc$doses_per_day_base,
                    doses_per_day_increase  = sheets$VaccinesMisc$doses_per_day_increase,
                    doses_per_day_maximum = sheets$VaccinesMisc$doses_per_day_maximum,
                    start_increase_day = sheets$VaccinesMisc$start_increase_day,
                    start_date = start_date,
                    end_date = end_date,
                    population = sheets$VaccinesPopulation,
                    vax_uptake = sheets$VaccinesMisc$vax_uptake,
                    max_second_dose_frac = rep(sheets$VaccinesMisc$max_second_dose_frac, length(start_date:end_date)),
                    variants = sheets$VaccinesVariants,
                    variant_day0 = sheets$VaccinesMisc$variant_day0)
  return(list(params = params, frac_pui = frac_pui, model.inputs = model.inputs, internal.args = internal.args, interventions = interventions, obs.data = obs.data, vaccines = vaccines.list$vaccines, vaccines_nonstan = vaccines.list$vaccines_nonstan))
}



