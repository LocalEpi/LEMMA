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
                 ReadExcel(path, sheet = "Vaccine Distribution", skip = 1),
                 ReadExcel(path, sheet = "Vaccine Doses - Observed", skip = 1),
                 ReadExcel(path, sheet = "Vaccine Doses - Future", skip = 1, col_types = c("text", "skip", "list", "list")),
                 ReadExcel(path, sheet = "Variants", skip = 2),
                 ReadExcel(path, sheet = "PUI Details", skip = 4),
                 ReadExcel(path, col_types = c("text", "list", "skip", "skip"), sheet = "Internal"))
  names(sheets) <- sapply(sheets, function (z) attr(z, "sheetname"))
  sheets <- rapply(sheets, as.Date, classes = "POSIXt", how = "replace") #convert dates

  sheets$Internal <- rbind(sheets$Internal, data.table(internal.name = "input.file", value = path))
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
  setkey(interventions, mu_t_inter)
  return(interventions)
}

ProcessSheets <- function(sheets) {
  params <- sheets$`Parameters with Distributions`[, .(name = internal.name, mu = Mean, sigma = `Standard Deviation`)]
  params[, sigma := pmax(sigma, mu / 100)] #Stan crashes if sigma = 0
  frac_pui <- sheets$`PUI Details`[, .(name = internal.name, mu = Mean)]

  model.inputs <- TableToList(sheets$`Model Inputs`)
  model.inputs$total.population <- sheets$`Vaccine Distribution`[, sum(pop)]
  internal.args <- TableToList(sheets$Internal)
  interventions <- sheets$Interventions
  obs.data <- sheets$Data


  all.na <- rowAlls(is.na(as.matrix(obs.data[, -"date"])))
  obs.data <- obs.data[all.na == F]

  interventions[, sigma_t_inter := pmax(sigma_t_inter, 0.01)]
  interventions[, sigma_beta_inter := pmax(sigma_beta_inter, 0.001)]
  interventions[, sigma_len_inter := pmax(sigma_len_inter, 0.01)]
  if (internal.args$automatic.interventions) {
    if (nrow(interventions) == 0) {
      min.date <- internal.args$simulation.start.date
    } else {
      min.date <- min(interventions$mu_t_inter)
    }
    max.obs.date <- max(obs.data$date)
    future.interventions <- interventions[mu_t_inter > max.obs.date]
    past.interventions <- interventions[mu_t_inter <= max.obs.date]
    past.interventions <- AddInterventions(past.interventions, min.date = min.date, max.date = max.obs.date)
    interventions <- rbind(past.interventions, future.interventions)
    setkey(interventions, mu_t_inter)
  }

  if (is.na(internal.args$output.filestr)) {
    internal.args$output.filestr <- sub(".xlsx", " output", internal.args$input.file, fixed = T)
  }
  if (internal.args$add.timestamp.to.filestr) {
    internal.args$output.filestr <- paste0(internal.args$output.filestr, gsub(":", "-", as.character(date()), fixed = T))
  }
  if (is.null(internal.args$prior.infection.vaccine.scale)) {
    internal.args$prior.infection.vaccine.scale <- 1
  }
  if (is.null(internal.args$child.transmission.scale)) {
    internal.args$child.transmission.scale <- 1
  }
  if (is.null(internal.args$skip.second.dose.fraction)) {
    internal.args$skip.second.dose.fraction <- 0
  }
  internal.args$weights <- rep(1, length(DataTypes())) #TODO - make this an input?

  start_date <- internal.args$simulation.start.date + 1
  end_date <- model.inputs$end.date
  doses_future <- c(TableToList(sheets$`Vaccine Doses - Future`, "mrna"), TableToList(sheets$`Vaccine Doses - Future`, "jj"))

  if (sheets$`Vaccine Distribution`[, abs(sum(dose_proportion) - 1)] > 1e-4) {
    stop("The Percentage of Vaccinated column on Vaccine Distribution should sum to 100%")
  }

  vaccines.list <- GetVaccineParams(doses_actual = sheets$`Vaccine Doses - Observed`,
                                    doses_future = doses_future,
                    start_date = start_date,
                    end_date = end_date,
                    population = sheets$`Vaccine Distribution`,
                    variants = sheets$Variants,
                    child_transmission_scale = internal.args$child.transmission.scale,
                    skip_second_dose_fraction = internal.args$skip.second.dose.fraction)
  return(list(params = params, frac_pui = frac_pui, model.inputs = model.inputs, internal.args = internal.args, interventions = interventions, obs.data = obs.data, vaccines = vaccines.list$vaccines, vaccines_nonstan = vaccines.list$vaccines_nonstan))
}



