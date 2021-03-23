library(data.table)
library(matrixStats)


#currently uses total first and second doses per day with stable age distribution - could extend in future to day x age x dose_num
GetDoses <- function(doses_actual, doses_future, population, start_date, end_date) {
  doses <- merge(doses_actual, data.table(date = seq(start_date, end_date, by ="day")), all.y = T)
  doses[is.na(dose1), dose1 := 0]
  doses[is.na(dose2), dose2 := 0]
  doses[is.na(doseJ), doseJ := 0]
  max_actual_vaccine_date <- doses_actual[, max(date)]
  doses[date > max_actual_vaccine_date, doses_available_mrna := doses_future$doses_per_day_base_mrna]
  doses[date >= doses_future$start_increase_day_mrna, doses_available_mrna := pmin(doses_future$doses_per_day_maximum_mrna, doses_future$doses_per_day_base_mrna + seq(from = 0, by = doses_future$doses_per_day_increase_mrna, length.out = length(doses_future$start_increase_day_mrna:end_date)))]

  doses[date > max_actual_vaccine_date, doses_available_jj := doses_future$doses_per_day_base_jj]
  doses[date >= doses_future$start_increase_day_jj, doses_available_jj := pmin(doses_future$doses_per_day_maximum_jj, doses_future$doses_per_day_base_jj + seq(from = 0, by = doses_future$doses_per_day_increase_jj, length.out = length(doses_future$start_increase_day_jj:end_date)))]

  dose_spacing <- 24 #avg Moderna and Pfizer
  nt <- nrow(doses)
  future_index <- which.max(doses$date > max_actual_vaccine_date)
  for (it in future_index:nt) {
    need_second_dose <- pmax(0, sum(doses$dose1[1:(it - dose_spacing)]) - sum(doses$dose2[1:(it - 1)]))
    doses[it, dose2 := pmin(doses_available_mrna, need_second_dose)]

    total_max_vax <- population[, sum(GetMaxVax(vax_uptake, date = start_date + it - 1, elligible_date, pop))]
    total_doses_given <- doses[1:it, sum(dose1) + sum(doseJ)]
    if (total_doses_given < total_max_vax) {
      #can overshoot total_max_vax for one day
      doses[it, dose1 := doses_available_mrna - dose2]
      doses[it, doseJ := doses_available_jj]
    }
  }
  return(doses[, .(date, dose1, dose2, doseJ)])
}

GetMaxVax <- function(vax_uptake, date, elligible_date, pop) {
  uptake <- vax_uptake
  uptake[date < elligible_date] <- 0
  return(uptake * pop)
}

GetVaccineParams <- function(doses_actual, doses_future, start_date, end_date, population, variants) {
  population <- copy(population)
  max_actual_vaccine_date <- doses_actual[, max(date)]
  doses <- GetDoses(doses_actual, doses_future, population, start_date, end_date)

  #CDC data is for 0 5 18 30 40, 50, 65, 75, 85
  #assume 5-11 and 12-15 are same as 5-17 and that 16-18 is same as 18-29 - to match vaccine elibibility brackets
  lethality <- data.table(age = c(0, 5, 12, 16, 30, 40, 50, 65, 75, 85),
                          hosp = c(1/4, 1/9, 1/9, 1, 2, 3, 4, 5, 8, 13),
                          death = c(1/9, 1/16, 1/16, 1, 4, 10, 30, 90, 220, 630))
  lethality[, icu_rate := sqrt(death/hosp)]
  lethality[, mort_rate := sqrt(death/hosp)]

  nt <- nrow(doses)
  num_variants <- nrow(variants)
  time_to_effect <- 10
  doses[, dose1_effective := shift(dose1, time_to_effect, fill = 0)]
  doses[, dose2_effective := shift(dose2, time_to_effect, fill = 0)]
  doses[, doseJ_effective := shift(doseJ, time_to_effect, fill = 0)]

  #frac1 = fraction with 1 mRNA dose
  #frac2 = fraction with 2 mRNA doses
  #fracJ = fraction with single J&J dose
  #total people with any dose: cumsum(dose1_effective) + cumsum(doseJ)  [because anyone with dose2_effective already had dose1_effective]
  doses[, any_dose := cumsum(dose1_effective) + cumsum(doseJ_effective)]
  doses[, frac1 := (cumsum(dose1_effective) - cumsum(dose2_effective)) / any_dose]
  doses[, frac2 := cumsum(dose2_effective) / any_dose]
  doses[, fracJ := cumsum(doseJ_effective) / any_dose]
  doses[any_dose == 0, frac1 := 1] #avoids errors
  doses[any_dose == 0, frac2 := 0]
  doses[any_dose == 0, fracJ := 0]
  doses$any_dose <- NULL

  # doses[, frac_fully := pmin(1, cumsum(dose2_effective) / cumsum(dose1_effective + 1e-10))]

  vaccinated_per_day <- doses[, dose1_effective + doseJ_effective]

  variant_frac <- matrix(nrow = nt, ncol = num_variants)
  for (i in 1:num_variants) {
    it_vec <- as.numeric(doses[, date] - variants$variant_day0[i])
    growth <- ifelse(it_vec <= 0, variants$daily_growth_prior[i], variants$daily_growth_future[i])
    variant_frac[, i] <- variants$frac_on_day0[i] * growth ^ it_vec
  }
  variant_frac <- variant_frac / rowSums(variant_frac) #recycles row sum

  vaccine_efficacy_against_progression <- vaccine_efficacy_for_susceptibility <- duration_vaccinated <- duration_natural  <- transmission_variant_multiplier <- rep(NA_real_, nt)
  for (it in 1:nt) {
    vaccine_efficacy_against_progression[it] <- sum(variant_frac[it, ] *
      (variants$vaccine_efficacy_against_progression_1 * doses$frac1[it] +
        variants$vaccine_efficacy_against_progression_2 * doses$frac2[it] +
        variants$vaccine_efficacy_against_progression_J * doses$fracJ[it]))
    vaccine_efficacy_for_susceptibility[it] <- sum(variant_frac[it, ] *
      (variants$vaccine_efficacy_for_susceptibility_1 * doses$frac1[it] +
        variants$vaccine_efficacy_for_susceptibility_2 * doses$frac2[it] +
        variants$vaccine_efficacy_for_susceptibility_J * doses$fracJ[it]))

    #input durations are in years, output in days
    duration_vaccinated[it] <- 365 * sum(variant_frac[it, ] * variants$duration_vaccinated_years)
    duration_natural[it] <- 365 * sum(variant_frac[it, ] * variants$duration_natural_years)
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

  date <- seq(start_date, end_date, by = "day")
  for (it in 1:nt) {
    dt[, max_vax := GetMaxVax(vax_uptake, date = start_date + it - 1, elligible_date, pop)]
    if (date[it] <= max_actual_vaccine_date) {
      dt[, new_dose_proportion := dose_proportion]
    } else {
      #after end of actual, average between previous proportion and population proportion
      dt[, new_dose_proportion := (dose_proportion + pop_frac) / 2]
    }
    not_maxed_prop <- dt[vax < max_vax, sum(new_dose_proportion)]
    dt[, new_doses := (vax < max_vax) * new_dose_proportion / (1e-10 + not_maxed_prop) * vaccinated_per_day[it]]
    dt[, vax := vax + new_doses]
    stopifnot(!anyNA(dt))
    dt[, unvaccinated_pop := pmax(0, pop - vax)]
    dt[, unvaccinated_pop_frac := unvaccinated_pop / sum(unvaccinated_pop)] #num unvax/total unvax

    frac_hosp_multiplier[it] <- dt[, sum(hosp * unvaccinated_pop_frac)] / hosp_rate_pop * sum(var_hosp_mult * variant_frac[it, ])
    frac_icu_multiplier[it] <- dt[, sum(icu_rate * unvaccinated_pop_frac)] / icu_rate_pop * sum(var_icu_mult * variant_frac[it, ])
    frac_mort_multiplier[it] <- dt[, sum(mort_rate * unvaccinated_pop_frac)] / mort_rate_pop * sum(var_mort_mult * variant_frac[it, ])
  }

  vaccines <- data.table(date, vaccinated_per_day, vaccine_efficacy_for_susceptibility, vaccine_efficacy_against_progression, duration_vaccinated, duration_natural, frac_hosp_multiplier, frac_icu_multiplier, frac_mort_multiplier, transmission_variant_multiplier)
  return(list(vaccines = vaccines, vaccines_nonstan = list(doses = doses, variant_frac = variant_frac, variants = variants)))
}

