#get mobility array from chris, get running on laptop with small iter, run on AWS with small iter, run on AWS with large iter (r6g.xlarge)

library(data.table)
output.filestr <- gsub(":| ", "-", paste0("Npop", Sys.time()))
# setwd("~/LEMMA/inst/extdata/")
setwd("~/Dropbox/jsLEMMA/inst/extdata/")
input.file <- "Npop input3.xlsx" #Npop3
sheets <- LEMMA:::ReadInputs(input.file)
inputs <- LEMMA:::ProcessSheets(sheets, input.file)

day0 <- inputs$internal.args$simulation.start.date
nt <- as.numeric(inputs$model.inputs$end.date - day0)

#mobility.array[i, j, t] is number of devices going to i from j on day t
mobility.array <- readRDS("CACountyMvmt2020-01-01to2020-07-27.rds") #npop x npop x nt
mobility.array.fips <- c("001", "003", "005", "007", "009", "011", "013", "015", "017",
                         "019", "021", "023", "025", "027", "029", "031", "033", "035",
                         "037", "039", "041", "043", "045", "047", "049", "051", "053",
                         "055", "057", "059", "061", "063", "065", "067", "069", "071",
                         "073", "075", "077", "079", "081", "083", "085", "087", "089",
                         "091", "093", "095", "097", "099", "101", "103", "105", "107",
                         "109", "111", "113", "115")
rownames(mobility.array) <- colnames(mobility.array) <- mobility.array.fips
dimnames(mobility.array)[[3]] <- as.Date("2020/1/1") + (1:dim(mobility.array)[3]) - 1
county.pop.dt <- readRDS("county population.RDS")
# dt <- unique(fread("~/Documents/GitHub/DL-COVID-19/DL-us-mobility-daterow.csv")[admin_level == 2 & admin1 == "California"])
# dt[, county := sub(" County", "", admin2)]
# dt[, ct_fips := substr(as.character(fips), 2, 4)]

# county.pop.dt <- county.pop.dt[population < 500000] #temp!

setkey(county.pop.dt, "ct_fips")
county.pop <- county.pop.dt$population
names(county.pop) <- county.pop.dt$county

index <- mobility.array.fips %in% county.pop.dt$ct_fips
mobility.array <- mobility.array[index, index, ]
stopifnot(all.equal(rownames(mobility.array), county.pop.dt$ct_fips))
rownames(mobility.array) <- colnames(mobility.array) <- county.pop.dt$county

county.set <- county.pop.dt$county

npop <- length(county.set)
mobility <- array(dim = c(nt, npop, npop))
for (it in 1:nt) {
  #mobility.array starts Jan 1; for dates beyond end of mobility.array, assume same as last mobility.array
  date.index <- pmin(as.integer(day0 - as.Date("2020/1/1")) + it, dim(mobility.array)[3])
  mobility.frac <- mobility.array[, , date.index]
  mobility.frac <- sweep(mobility.frac, 2, colSums(mobility.frac), "/")
for (ipop1 in 1:npop) {
  for (ipop2 in 1:ipop1) {
      if (ipop1 == ipop2) {
        mobility[it, ipop1, ipop2] <- mobility.frac[ipop1, ipop2]
      } else {
        mobility[it, ipop1, ipop2] <- mobility[it, ipop2, ipop1] <- (mobility.frac[ipop1, ipop2] * county.pop[ipop2] + mobility.frac[ipop2, ipop1] * county.pop[ipop1]) / (county.pop[ipop1] + county.pop[ipop2]) #mobility.frac[i, j] is fraction of devices going to i from j
      }
    }
  }
}


#mobility[it, ipop1, ipop2] = (# in pop1 moving to pop2 + # in pop2 moving to pop1) / (population1 + population2) [ipop1 != ipop2]
#mobility[it, ipop1, ipop1] = # in pop1 staying in ipop1 / population1


if (!exists("hosp.dt")) {
  hosp.dt <- fread("https://data.ca.gov/dataset/529ac907-6ba1-4cb7-9aae-8966fc96aeef/resource/42d33765-20fd-44b8-a978-b083b7542225/download/hospitals_by_county.csv")
  #hosp.dt <- hosp.dt[todays_date != "", .(county, date = as.Date(todays_date), hosp.conf = hospitalized_covid_confirmed_patients, hosp.pui = hospitalized_suspected_covid_patients)]
  hosp.dt <- hosp.dt[, .(county, date = as.Date(todays_date), hosp.conf = hospitalized_covid_confirmed_patients, hosp.pui = hospitalized_suspected_covid_patients)]

  hosp.dt <- hosp.dt[date > as.Date("2020/4/1"), .(county, date, hosp.conf, hosp.pui)]
}
counties <- hosp.dt[county %in% county.set]

obs.data.conf <- dcast(counties, date ~ county, value.var = "hosp.conf")
obs.data.pui <- dcast(counties, date ~ county, value.var = "hosp.pui")
setcolorder(obs.data.conf, c("date", county.set))
setcolorder(obs.data.pui, c("date", county.set))


mean.ini <- 1e-5 * county.pop
lambda_ini_exposed <- 1 / mean.ini

inputs$model.inputs$total.population <- county.pop
inputs$obs.data.conf <- obs.data.conf
inputs$obs.data.pui <- obs.data.pui
inputs$internal.args$lambda_ini_exposed <- lambda_ini_exposed

#update these to reflect changes in cell mobility and/or policy?
inputs$mu_beta_inter <- matrix(inputs$interventions$mu_beta_inter, nrow = length(inputs$interventions$mu_beta_inter), ncol = npop)
inputs$sigma_beta_inter <- matrix(inputs$interventions$sigma_beta_inter, nrow = length(inputs$interventions$sigma_beta_inter), ncol = npop)

inputs$mobility <- mobility
#reset obs.data.conf, obs.data.pui, lambda_ini_exposed, total.population, mu_beta_inter, sigma_beta_inter, mobility

inputs$internal.args$random.seed <- 10
inputs$internal.args$iter <- 1200
inputs$internal.args$adapt_delta <- 0.95
inputs$internal.args$max_treedepth <- 13 #last run was 10 - lots of exceeded warnings
inputs$internal.args$warmup <- round(inputs$internal.args$iter * 0.7)

saveRDS(inputs, file = paste0(output.filestr, "inputs.RDS"))
fit <- LEMMA:::CredibilityInterval(inputs)

library(matrixStats)
library(ggplot2)
library(data.table)

pdf(paste0(output.filestr, ".pdf"), width = 11, height = 8)
Rt <- rstan::extract(fit, pars = "Rt")[[1]]
rt <- as.data.table(Rt)
setnames(rt, c("V2", "V1", "value"), c("t", "county.num", "Rt"))
rt <- rt[, .(median.Rt = median(Rt)), by = c("t", "county.num")]

rt <- merge(rt, data.table(county.num = 1:length(county.set), county = county.set))
gg <- ggplot(rt, aes(x = as.Date("2020/2/15") + t, y = median.Rt, color = county)) + geom_line() + xlab("")
print(directlabels::direct.label(gg))

quantiles <- list()
hosp <- rstan::extract(fit, pars = "sim_data")[[1]]
for (i in seq_along(county.set)) {
  quant <- colQuantiles(hosp[, , i], probs = seq(0, 1, by = 0.05))
  rownames(quant) <- as.character(as.Date("2020/2/15") + 1:nrow(quant))
  quantiles[[county.set[i]]] <- quant
}

temp.data <- inputs$obs.data.conf
names(temp.data)[-1] <- paste0(names(temp.data)[-1], ".conf")
temp.data2 <- inputs$obs.data.pui
names(temp.data2)[-1] <- paste0(names(temp.data2)[-1], ".pui")
temp.data2$date <- NULL
inputs$obs.data <- cbind(temp.data, temp.data2)

expansion <- function() c(0, 0, 0, 0) #fixes a problem if AWS has old ggplot version
for (i in county.set) {

  LEMMA:::GetProjectionPlot(short.term = T, quantiles = quantiles, data.type = i, inputs = inputs)
}

for (p in c("duration_latent", "duration_rec_mild", "duration_pre_hosp",
            "duration_hosp_mod", "duration_hosp_icu", "frac_hosp", "ini_exposed",
            "sigma_obs", "r0", "t_inter", "len_inter",
            "mobility_coef_home_0", "mobility_coef_home_1",
            "mobility_coef_away_0", "mobility_coef_away_1")) {
  print(rstan::stan_hist(fit, p))
}
dev.off()

