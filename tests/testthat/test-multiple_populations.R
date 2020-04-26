test_that("two populations works", {
  #could be cleaned up and turned into a testthat 
  
  setwd("~/Dropbox/LEMMA_shared/JS code branch/LEMMA")
  library(devtools)
  library(data.table)
  
  RandInt <- function(n, min.val, max.val) {
    x <- floor(runif(n, min=min.val, max=(max.val + 1)))
  }
  if (F) {
    set.seed(NULL); sseed <- round(runif(1, 1, 10000)); set.seed(sseed); cat("sseed = ", sseed, "\n")
  } else {
    set.seed(2267)
  }
  
  start.date <- as.Date("2020/2/23") #short
  end.date <- as.Date("2020/3/25")
  
  num.params <- 3
  latent.period <- RandInt(num.params, 0, 6)
  hosp.length.of.stay1 <- RandInt(num.params, 5, 15)
  prop.hospitalized1 <- RandInt(num.params, 4, 6) / 100
  params <- data.table(
    latent.period = latent.period,
    illness.length.given.nonhosp = RandInt(num.params, 2, 10),
    exposed.to.hospital = RandInt(num.params, latent.period + 1, 12),
    hosp.length.of.stay1 = hosp.length.of.stay1,
    hosp.length.of.stay2 = hosp.length.of.stay1 + RandInt(num.params, 1, 5),
    prop.hospitalized1 = prop.hospitalized1,
    prop.hospitalized2 = prop.hospitalized1 + RandInt(num.params, 1, 3) / 100,
    prop.icu = 0,
    prop.vent = 0,
    patients.in.hosp.are.infectious = T)
  
  N1 <- 880000
  N2 <- 10000
  N <- N1 + N2
  params$r0.initial.11 <- 3
  
  #new infections in 1 caused by 2
  params$r0.initial.12 <- 0 #3 * N1/N 
  
  #new infections in 2 caused by 1
  params$r0.initial.21 <- 0 #3 * N2/N 
  
  params$r0.initial.22 <- 4
  
  num.init.exp <- 5
  initial.new.exposures <- matrix(runif(num.init.exp*num.params), nrow=num.params, ncol=num.init.exp)
  
  # load_all() 
  z.new <- SEIR(initial.new.exposures = initial.new.exposures, initial.conditions = c(875000, 8000), start.date, end.date, params)
  #compare vs running each init.exp separately (there was a bug, now fixed)
  for (i in 1:num.init.exp) {
    z.new1 <- SEIR(initial.new.exposures = initial.new.exposures[, i, drop = F], initial.conditions = c(875000, 8000), start.date, end.date, params)
    for (j in names(z.old)) {
      cat(i, j, "\n")
      stopifnot(all.equal(z.new[[j]][, , , i, drop = F], z.new1[[j]]))
    }
  }
  #compare vs running LEMMA v0.2.0.9001 with init.exp in each population 
  #modifications needed: see warning("temp!") in SEIR - disable truncation, set same init exp in both pops
  if (F) {
    unload()
    for (i in 1:2) {
      p <- copy(params)
      if (i == 1) {
        p$r0.initial <- params$r0.initial.11 
        p$hosp.length.of.stay <- params$hosp.length.of.stay1
        p$prop.hospitalized <- params$prop.hospitalized1
      } else {
        p$r0.initial <- params$r0.initial.22 
        p$hosp.length.of.stay <- params$hosp.length.of.stay2
        p$prop.hospitalized <- params$prop.hospitalized2
      }
      stopifnot(params$r0.initial.12 == 0, params$r0.initial.21 == 0)
      p$use.hosp.rate <- T
      z.old <- LEMMA:::SEIR(initial.new.exposures = initial.new.exposures, initial.conditions = c(875000, 8000)[i], start.date, end.date, p) 
      for (j in names(z.old)) {
        cat(i, j, "\n")
        stopifnot(all.equal(adrop(z.new[[j]][, , i, , drop = F], 3), z.old[[j]]))
      }
    }
  }
  
  
})
