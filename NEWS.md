# LEMMA 0.5.0.9001

## New Features
* Moved credibility interval calculation to Stan based on the Santa Cruz County COVID-19 Model (https://github.com/jpmattern/seir-covid19). Many thanks to Jann Paul Mattern and Mikala Caton!

# LEMMA 0.4.0.9000

## New Features
* Added fitting to ICU, deaths and cumulative hospital admissions

# LEMMA 0.3.0.9005

## Minor improvements and fixes
* CredibilityIntervalFromExcel now supports up to 9 interventions (see FAQ for details)


# LEMMA 0.3.0.9004 

## New Features
* “Best Guess” was sometimes misinterpreted so we changed “Best Guess” to “User’s Prior Projection”. User’s Prior Projection is only displayed if niter < 100.

* PDF now includes short-term (dates of observed data) and long-term (until Projection End Date) plots

* `plot.observed.data.long.term`, `plot.observed.data.short.term` added to internal.args (Internal sheet). These logicals control plotting the observed data bounds on the long-term and short-term plots.

* `lower.bound.label`, `upper.bound.label` added to internal.args (Internal sheet). For example, in San Francisco these are "Confirmed COVID19" and "Confirmed COVID19 + 30%PUI".

## Minor improvements and fixes
* Case Studies was updated to include plots
