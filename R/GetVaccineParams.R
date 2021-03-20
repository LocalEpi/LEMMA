library(data.table)
library(matrixStats)


#currently uses total first and second doses per day with stable age distribution - could extend in future to day x age x dose_num
GetDoses <- function(doses_actual, doses_per_day_base, doses_per_day_increase, doses_per_day_maximum, start_increase_day, start_date, end_date, population, vax_uptake, max_second_dose_frac) {
  doses <- merge(doses_actual, data.table(date = seq(start_date, end_date, by ="day")), all.y = T)
  doses[is.na(dose1), dose1 := 0]
  doses[is.na(dose2), dose2 := 0]
  max_actual_vaccine_date <- doses_actual[, max(date)]
  doses[date > max_actual_vaccine_date, doses_available := doses_per_day_base]
  doses[date >= start_increase_day, doses_available := pmin(doses_per_day_maximum, doses_per_day_base + seq(from = 0, by = doses_per_day_increase, length.out = length(start_increase_day:end_date)))]

  total_max_vaccinated <- population[, sum(max_vax)]

  dose_spacing <- 24 #avg Moderna and Pfizer
  doses[, need_second_dose := NA]
  doses[1:dose_spacing, need_second_dose := 0]
  doses[, first_dose_done := F]
  nt <- nrow(doses)
  for (it in (dose_spacing + 1):nt) {
    doses$need_second_dose[it] <- pmax(0, sum(doses$dose1[1:(it - dose_spacing)]) - sum(doses$dose2[1:it]))
    if (doses[it, date] > max_actual_vaccine_date) {
      first_dose_done <- sum(doses$dose1[1:(it - 1)]) >= total_max_vaccinated
      if (first_dose_done) {
        max_second_dose_frac1 <- 1
        doses[it, first_dose_done := T]
      } else {
        max_second_dose_frac1 <- max_second_dose_frac[it]
      }
      doses[it, dose2 := pmin(doses_available * max_second_dose_frac1, need_second_dose)]
      if (first_dose_done) {
        doses[it, dose1 := 0]
      } else {
        doses[it, dose1 := doses_available - dose2]
      }
    }

  }
  return(doses[, .(date, dose1, dose2)])
}

GetVaccineParams <- function(doses_actual, doses_per_day_base, doses_per_day_increase, doses_per_day_maximum, start_increase_day, start_date, end_date, population, vax_uptake, max_second_dose_frac, variants, variant_day0) {
  population <- copy(population)
  population[, max_vax := vax_uptake * pop]
  population[dose_proportion == 0, max_vax := 0]

  doses <- GetDoses(doses_actual, doses_per_day_base, doses_per_day_increase, doses_per_day_maximum, start_increase_day, start_date, end_date, population, vax_uptake, max_second_dose_frac)

  lethality <- data.table(age = c(0, 5, 18, 30, 40, 50, 65, 75, 85),
                          hosp = c(1/4, 1/9, 1, 2, 3, 4, 5, 8, 13),
                          death = c(1/9, 1/16, 1, 4, 10, 30, 90, 220, 630))
  lethality[, icu_rate := sqrt(death/hosp)]
  lethality[, mort_rate := sqrt(death/hosp)]

  nt <- nrow(doses)
  num_variants <- nrow(variants)
  time_to_effect <- 10
  doses[, dose1_effective := shift(dose1, time_to_effect, fill = 0)]
  doses[, dose2_effective := shift(dose2, time_to_effect, fill = 0)]
  doses[, frac_fully := pmin(1, cumsum(dose2_effective) / cumsum(dose1_effective + 1e-10))]

  vaccinated_per_day <- doses[, dose1_effective]
  frac_fully <- doses[, frac_fully]

  it <- as.numeric(doses[, date] - variant_day0)
  growth <- rbind(matrix(variants$daily_growth_prior, nrow = sum(it <= 0), ncol = num_variants, byrow = T),
                  matrix(variants$daily_growth_future, nrow = sum(it > 0), ncol = num_variants, byrow = T))
  variant_frac <- matrix(variants$frac_on_day0, nrow = nt, ncol = num_variants, byrow = T)
  variant_frac <- variant_frac * growth ^ it #recycles it
  variant_frac <- variant_frac / rowSums(variant_frac) #recycles row sum

  vaccine_efficacy_against_progression <- vaccine_efficacy_for_susceptibility <- duration_vaccinated <- duration_natural  <- transmission_variant_multiplier <- rep(NA_real_, nt)
  for (it in 1:nt) {
    vaccine_efficacy_against_progression[it] <- sum((variants$vaccine_efficacy_against_progression_1 * (1 - frac_fully[it]) + variants$vaccine_efficacy_against_progression_2 * frac_fully[it]) * variant_frac[it, ])
    vaccine_efficacy_for_susceptibility[it] <- sum((variants$vaccine_efficacy_for_susceptibility_1 * (1 - frac_fully[it]) + variants$vaccine_efficacy_for_susceptibility_2 * frac_fully[it]) * variant_frac[it, ])

    duration_vaccinated[it] <- sum(variant_frac[it, ] * variants$duration_vaccinated_)
    duration_natural[it] <- sum(variant_frac[it, ] * variants$duration_natural_)
    transmission_variant_multiplier[it] <- sum(variant_frac[it, ] * variants$transmisson_mult)
  }

  #input is mort|infected, LEMMA uses icu|hosp and mort|icu
  var_hosp_mult <- variants$hosp_mult
  var_mort_mult <- var_icu_mult <- sqrt(variants$mort_mult / var_hosp_mult)

  #multiplier = rate_unvax / rate_pop
  frac_hosp_multiplier <- frac_icu_multiplier <- frac_mort_multiplier <- rep(NA_real_, nt)
  dt <- cbind(lethality, population)

  dt[, vax := 0]
  dt[, pop_frac := pop / sum(pop)]
  hosp_rate_pop <- dt[, sum(hosp * pop_frac)] #rate relative to 18-30
  icu_rate_pop <- dt[, sum(icu_rate * pop_frac)] #rate relative to 18-30
  mort_rate_pop <- dt[, sum(mort_rate * pop_frac)] #rate relative to 18-30

  for (it in 1:nt) {
    not_maxed_prop <- dt[vax < max_vax, sum(dose_proportion)]
    dt[, new_doses := (vax < max_vax) * dose_proportion / (1e-10 + not_maxed_prop) * vaccinated_per_day[it]]
    dt[, vax := vax + new_doses]
    stopifnot(!anyNA(dt))
    dt[, unvaccinated_pop := pmax(0, pop - vax)]
    dt[, unvaccinated_pop_frac := unvaccinated_pop / sum(unvaccinated_pop)] #num unvax/total unvax

    frac_hosp_multiplier[it] <- dt[, sum(hosp * unvaccinated_pop_frac)] / hosp_rate_pop * sum(var_hosp_mult * variant_frac[it, ])
    frac_icu_multiplier[it] <- dt[, sum(icu_rate * unvaccinated_pop_frac)] / icu_rate_pop * sum(var_icu_mult * variant_frac[it, ])
    frac_mort_multiplier[it] <- dt[, sum(mort_rate * unvaccinated_pop_frac)] / mort_rate_pop * sum(var_mort_mult * variant_frac[it, ])
  }

  date <- seq(start_date, end_date, by = "day")
  vaccines <- data.table(date, vaccinated_per_day, vaccine_efficacy_for_susceptibility, vaccine_efficacy_against_progression, duration_vaccinated, duration_natural, frac_hosp_multiplier, frac_icu_multiplier, frac_mort_multiplier, transmission_variant_multiplier)
  return(list(vaccines = vaccines, vaccines_nonstan = list(doses = doses, variant_frac = variant_frac, variants = variants)))
}

