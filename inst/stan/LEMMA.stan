// Modification of Santa Cruz County COVID-19 Model
// https://github.com/jpmattern/seir-covid19
// by Jann Paul Mattern and Mikala Caton

data {
  //////////////////////////////////////////
  // data required to run model

  int<lower=0> nobs_types;                           // number of observation data types
  int<lower=0> nobs[nobs_types];                     // number of timesteps with observations
  int <lower=0> nobs_max;                            // maximum of nobs
  int<lower=-1> tobs[nobs_types, nobs_max];           // obs times; -1 = NA
  matrix<lower=-1.0>[nobs_types, nobs_max] obs_data;  // observed confirmed (-1 = NA)
  int<lower=0> nt;                                    // number of time steps
  real<lower=0.0> npop;                               // total population

  //////////////////////////////////////////
  // prior parameter distributions
  real<lower=1.0> mu_duration_latent;     // mean duration in "exposed" stage
  real<lower=0.0> sigma_duration_latent;  // sd duration in "exposed" stage
  real<lower=1.0> mu_duration_rec_mild;   // mean duration in "infectious" stage for mild cases
  real<lower=0.0> sigma_duration_rec_mild;// sd duration in "infectious" stage for mild cases
  real<lower=1.0> mu_duration_pre_hosp;   // mean duration in "infectious" stage for hospitalized cases
  real<lower=0.0> sigma_duration_pre_hosp;// sd duration in "infectious" stage for hospitalized cases
  real<lower=1.0> mu_duration_hosp_mod;   // mean duration in hospital for non-ICU cases
  real<lower=0.0> sigma_duration_hosp_mod;// sd duration in hospital for non-ICU cases
  real<lower=1.0> mu_duration_hosp_icu;   // mean duration in hospital for ICU cases
  real<lower=0.0> sigma_duration_hosp_icu;// sd duration in hospital for ICU cases
  real<lower=1.0> mu_duration_mort_nonhosp;
  real<lower=0.0> sigma_duration_mort_nonhosp;


  real<lower=0.0> mu_r0;                  // mean initial beta estimate
  real<lower=0.0> sigma_r0;               // sd initial beta estimate

  real<lower=0.0> mu_frac_hosp;           // mean ICU + non-ICU
  real<lower=0.0> sigma_frac_hosp;        // sd ICU + non-ICU
  real<lower=0.0> mu_frac_icu;            // mean ICU as fraction of hosp
  real<lower=0.0> sigma_frac_icu;         // sd ICU as fraction of hosp
  real<lower=0.0> mu_frac_mort;           // mean mortality as fraction of ICU
  real<lower=0.0> sigma_frac_mort;        // sd mortality as fraction of ICU
  real<lower=0.0> mu_frac_mort_nonhosp;
  real<lower=0.0> sigma_frac_mort_nonhosp;

  real<lower=0.0> lambda_ini_exposed;     // parameter for initial conditions of "exposed"


  //////////////////////////////////////////
  // interventions

  int<lower=0> ninter;                      // number of interventions
  real<lower=1.0> mu_t_inter[ninter];       // mean start time of each interventions
  real<lower=0.0> sigma_t_inter[ninter];    // sd start time of each interventions
  real<lower=1.0> mu_len_inter[ninter];     // mean length of each intervention
  real<lower=0.0> sigma_len_inter[ninter];  // sd length of each intervention
  real<lower=0.0> mu_beta_inter[ninter];    // mean change in beta through intervention
  real<lower=0.0> sigma_beta_inter[ninter]; // sd change in beta through intervention

  real<lower=0.0> vaccinated_per_day[nt]; //vaccinated per day total
  real<lower=0.0, upper=1.0> vaccine_efficacy_for_susceptibility[nt];

  real<lower=1.0> duration_vaccinated[nt];
  real<lower=1.0> duration_natural[nt];
  real<lower=0.0> frac_hosp_multiplier_unvaccinated[nt];  //multiplier due to vaccines/variants
  real<lower=0.0> frac_hosp_multiplier_vaccinated[nt];  //multiplier due to vaccines/variants
  real<lower=0.0> frac_icu_multiplier_unvaccinated[nt];   //multiplier due to vaccines/variants
  real<lower=0.0> frac_icu_multiplier_vaccinated[nt];   //multiplier due to vaccines/variants
  real<lower=0.0> frac_mort_multiplier_unvaccinated[nt];   //multiplier due to vaccines/variants
  real<lower=0.0> frac_mort_multiplier_vaccinated[nt];   //multiplier due to vaccines/variants
  real<lower=0.0> frac_mort_nonhosp_multiplier_unvaccinated[nt];   //multiplier due to vaccines/variants
  real<lower=0.0> frac_mort_nonhosp_multiplier_vaccinated[nt];   //multiplier due to vaccines/variants

  real<lower=0.0> transmission_multiplier_unvaccinated[nt]; //multiplier due to variants and age (vaccine effect is in vaccine_efficacy_for_susceptibility)
  real<lower=0.0> transmission_multiplier_vaccinated[nt]; //multiplier due to variants and age (vaccine effect is in vaccine_efficacy_for_susceptibility)
   real<lower=0.0> prior_infection_vaccine_scale; //e.g. 0.8 = those with prior infection 20% less likely to get vaccinated than those without prior infection

  real<lower=0.0> sigma_obs_est_inv[nobs_types];

  real<lower=0.0> mu_frac_tested;         // mean of infected who test positive
  real<lower=0.0> sigma_frac_tested;      // sd of infected who test positiveof infected who test positive

}

transformed data {
  //assigning indices for state matrix x
  int Su = 1;
  int Sv_succ = 2;
  int Eu = 3;
  int Ev = 4;
  int Imildu = 5; //note: on compartment diagram this is labelled "Inonhospu"
  int Imildv = 6; //note: on compartment diagram this is labelled "Inonhospv"
  int Iprehu = 7;
  int Iprehv = 8;
  int Hmodu  = 9;
  int Hmodv  = 10;
  int Hicuu  = 11;
  int Hicuv  = 12;
  int Rliveu = 13;
  int Rlivev = 14;
  int Rpremort_nonhospu = 15;
  int Rpremort_nonhospv = 16;
  int Rmortu = 17;
  int Rmortv = 18;
  int Sv_fail = 19;

  int ncompartments = 19;

  int obs_hosp_census = 1;
  int obs_icu_census = 2;
  int obs_cum_deaths = 3;
  int obs_admits = 4;
  int obs_cases = 5;
  int obs_seroprev = 6;
}
parameters {
  real<lower=1.0> duration_latent; // duration is a minimum of 1 which is the stepsize of this model
  real<lower=1.0> duration_rec_mild;
  real<lower=1.0> duration_pre_hosp;
  real<lower=1.0> duration_hosp_mod;
  real<lower=1.0> duration_hosp_icu;
  real<lower=1.0> duration_mort_nonhosp;

  real<lower=0.005, upper=1.0> frac_hosp;
  real<lower=0.0, upper=1.0> frac_icu;
  real<lower=0.0, upper=1.0> frac_mort;
  real<lower=0.0, upper=1.0> frac_tested;
  real<lower=0.0, upper=1.0> frac_mort_nonhosp;

  real<lower=0.0> ini_exposed;
  real<lower=0.0> sigma_obs[nobs_types];
  real<lower=0.0> r0;
  real<lower=0.1> beta_multiplier[ninter];
  real<lower=1.0> t_inter[ninter];
  real<lower=1.0> len_inter[ninter];

}
transformed parameters {
  matrix<lower=0.0>[ncompartments,nt] x;
  matrix<lower=0.0>[nobs_types,nt] sim_data;
  real<lower=0.0> beta[nt];
  row_vector<lower=0.0>[nt] new_admits;
  row_vector<lower=0.0>[nt] new_admitsu;
  row_vector<lower=0.0>[nt] new_admitsv;
  row_vector<lower=0.0>[nt] new_cases;
  real<lower=0.0> total_cases[nt];
  real<lower=0.0> total_casesu[nt];
  real<lower=0.0> total_casesv[nt];

  {
    // variables in curly brackets will not have output, they are local variables
    real newEu;
    real newEv;
    real newIu;
    real newIv;
    real newrecu_mild;
    real newrecv_mild;
    real newrecu_mod;
    real newrecv_mod;
    real newhospu;
    real newhospv;
    real frac_hospu;
    real frac_hospv;
    real frac_icuu;
    real frac_icuv;
    real frac_mortu;
    real frac_mortv;
    real frac_mort_nonhospu;
    real frac_mort_nonhospv;
    real leave_icuu;
    real leave_icuv;
    real leave_premort_nonhosp;
    real leave_premort_nonhospu;
    real leave_premort_nonhospv;
    real beta_0;
    real vaccinated;
    real frac_vac_S;
    real newSv_succ;
    real newSv_fail;
    real newRlivev;
    real S_lostv;
    real R_lostv;
    real R_lostnatu;
    real R_lostnatv;
    real gained_protection;
    real infected_effective;

    //////////////////////////////////////////
    // Calculate beta for each time point
    beta_0 = r0 / (frac_hosp * duration_pre_hosp + (1 - frac_hosp) * duration_rec_mild);
    for (it in 1:nt) {
      beta[it] = beta_0;
      for (iinter in 1:ninter) {
        //k <- 2/s * qlogis(0.99) # = -2/s * qlogis(0.01) --> 2 * qlogis(0.99) = 9.19024
        //f <- m ^ plogis(k * (t - (d + s/2)))
        beta[it] = beta[it] * beta_multiplier[iinter] ^ inv_logit(9.19024 / len_inter[iinter] * (it - (t_inter[iinter] + len_inter[iinter] / 2)));
      }
    }

    // initial cond
    x[:,1] = rep_vector(0.0, ncompartments);
    total_cases[1] = 0.0;
    total_casesu[1] = 0.0;
    total_casesv[1] = 0.0;
    x[Eu,1] = ini_exposed;
    x[Su,1] = npop - ini_exposed;
    new_admits[1] = 0.0;
    new_admitsu[1] = 0.0;
    new_admitsv[1] = 0.0;
    new_cases[1] = 0.0;


    //////////////////////////////////////////
    // the SEIR model
    for (it in 1:nt-1){
      //////////////////////////////////////////
      // set transition variables
      infected_effective = transmission_multiplier_unvaccinated[it] * (x[Imildu,it] + x[Iprehu,it]) + transmission_multiplier_vaccinated[it] * (x[Imildv,it] + x[Iprehv,it]);

      newEu = fmin(x[Su,it],
      x[Su,it] * beta[it]/ npop * infected_effective);

      newEv = fmin(x[Sv_fail,it],
      x[Sv_fail,it] * beta[it] / npop * infected_effective);

      total_cases[it + 1] = newEu + newEv + total_cases[it];
      total_casesu[it + 1] = newEu + total_casesu[it];
      total_casesv[it + 1] = newEv + total_casesv[it];

      vaccinated = fmin(vaccinated_per_day[it], x[Su, it] + x[Eu, it] + x[Rliveu, it]);


      newIu = 1.0/duration_latent * x[Eu, it];
      newIv = 1.0/duration_latent * x[Ev, it];
      newhospu = 1.0/duration_pre_hosp * x[Iprehu, it];
      newhospv = 1.0/duration_pre_hosp * x[Iprehv, it];
      newrecu_mild = 1.0/duration_rec_mild * x[Imildu, it];
      newrecv_mild = 1.0/duration_rec_mild * x[Imildv, it];
      newrecu_mod = 1.0/duration_hosp_mod * x[Hmodu, it];
      newrecv_mod = 1.0/duration_hosp_mod * x[Hmodv, it];
      leave_icuu = 1.0/duration_hosp_icu * x[Hicuu, it];
      leave_icuv = 1.0/duration_hosp_icu * x[Hicuv, it];
      leave_premort_nonhospu = 1.0/duration_mort_nonhosp * x[Rpremort_nonhospu, it];
      leave_premort_nonhospv = 1.0/duration_mort_nonhosp * x[Rpremort_nonhospv, it];
      frac_vac_S = x[Su, it] / (x[Su, it] + x[Eu, it] + prior_infection_vaccine_scale * x[Rliveu, it]); //approx - ignores I and splits between S and Rlive - small magnitude
      newSv_succ = vaccine_efficacy_for_susceptibility[it] * vaccinated * frac_vac_S;
      newSv_fail = (1 - vaccine_efficacy_for_susceptibility[it]) * vaccinated * frac_vac_S;
      newRlivev = vaccine_efficacy_for_susceptibility[it] * vaccinated * (1 - frac_vac_S); //for now only failed vaccination of Rlive just stays in newRliveu [small problem for those who then lose natural immunity - will be counted as non-breakthrough but should be breakthrough]

      if (vaccine_efficacy_for_susceptibility[it + 1] > vaccine_efficacy_for_susceptibility[it]) {
         gained_protection = (vaccine_efficacy_for_susceptibility[it + 1] - vaccine_efficacy_for_susceptibility[it]) / (1 - vaccine_efficacy_for_susceptibility[it]) * x[Sv_fail, it];
      } else {
        gained_protection = - (vaccine_efficacy_for_susceptibility[it] - vaccine_efficacy_for_susceptibility[it + 1]) / vaccine_efficacy_for_susceptibility[it] * x[Sv_succ, it];
      }


      S_lostv = 1.0/duration_vaccinated[it] * x[Sv_succ, it];
      R_lostv = 1.0/duration_vaccinated[it] * x[Rlivev, it];
      R_lostnatu = 1.0/duration_natural[it] * x[Rliveu, it];
      R_lostnatv = 1.0/duration_natural[it] * x[Rlivev, it];

      frac_hospu = frac_hosp * frac_hosp_multiplier_unvaccinated[it];
      frac_hospv = frac_hosp * frac_hosp_multiplier_vaccinated[it];

      frac_icuu = frac_icu * frac_icu_multiplier_unvaccinated[it];
      frac_icuv = frac_icu * frac_icu_multiplier_vaccinated[it];

      frac_mortu = frac_mort * frac_mort_multiplier_unvaccinated[it];
      frac_mortv = frac_mort * frac_mort_multiplier_vaccinated[it];

      //these need to multiplier hosp*icu*mort because frac_mort is fraction died given icu, frac_icu is fraction icu given hosp
      frac_mort_nonhospu = frac_mort_nonhosp * frac_mort_nonhosp_multiplier_unvaccinated[it];
      frac_mort_nonhospv = frac_mort_nonhosp * frac_mort_nonhosp_multiplier_vaccinated[it];
      //////////////////////////////////////////
      // S -> E -> I

      x[Su, it+1] = x[Su, it] - newEu - newSv_succ - newSv_fail + R_lostnatu;
      x[Sv_succ, it+1] = x[Sv_succ, it] + newSv_succ - S_lostv + R_lostnatv + gained_protection;
      x[Sv_fail, it+1] = x[Sv_fail, it] - newEv + newSv_fail + S_lostv - gained_protection;
      x[Eu, it+1] = x[Eu, it] + newEu - newIu;
      x[Ev, it+1] = x[Ev, it] + newEv - newIv;
      x[Imildu, it+1] = x[Imildu, it] + newIu * (1 - frac_hospu) - newrecu_mild;
      x[Imildv, it+1] = x[Imildv, it] + newIv * (1 - frac_hospv) - newrecv_mild;
      x[Iprehu, it+1] = x[Iprehu, it] + newIu * frac_hospu - newhospu;
      x[Iprehv, it+1] = x[Iprehv, it] + newIv * frac_hospv - newhospv;
      x[Hmodu, it+1] = x[Hmodu, it] + newhospu * (1 - frac_icuu) - newrecu_mod + leave_icuu * (1 - frac_mortu);
      x[Hmodv, it+1] = x[Hmodv, it] + newhospv * (1 - frac_icuv) - newrecv_mod + leave_icuv * (1 - frac_mortv);
      x[Hicuu, it+1] = x[Hicuu, it] + newhospu * frac_icuu - leave_icuu;
      x[Hicuv, it+1] = x[Hicuv, it] + newhospv * frac_icuv - leave_icuv;
      x[Rliveu, it+1] = x[Rliveu, it] + newrecu_mild * (1 - frac_mort_nonhospu) + newrecu_mod - newRlivev + R_lostv - R_lostnatu;
      x[Rlivev, it+1] = x[Rlivev, it] + newrecv_mild * (1 - frac_mort_nonhospv) + newrecv_mod + newRlivev - R_lostv - R_lostnatv;
      x[Rpremort_nonhospu, it+1] = x[Rpremort_nonhospu, it] + newrecu_mild * frac_mort_nonhospu - leave_premort_nonhospu;
      x[Rpremort_nonhospv, it+1] = x[Rpremort_nonhospv, it] + newrecv_mild * frac_mort_nonhospv - leave_premort_nonhospv;
      x[Rmortu, it+1] = x[Rmortu, it] + leave_icuu * frac_mortu + leave_premort_nonhospu;
      x[Rmortv, it+1] = x[Rmortv, it] + leave_icuv * frac_mortv + leave_premort_nonhospv;

      //hospital admissions
      new_admits[it+1] = newhospu + newhospv;
      new_admitsu[it+1] = newhospu;
      new_admitsv[it+1] = newhospv;

      // new cases that are tested (assumes everyone tests positive on day transitions from E to I - not quite right but minor issue)
      new_cases[it+1] = (newIu + newIv) * frac_tested;

      //////////////////////////////////////////
      // test
      if (fabs(sum(x[:,it+1])-npop) > 1e-1){
        reject("Model is leaking, net gain: ", sum(x[:,it+1])-npop);
        // reject("Model is leaking, net gain: ", sum(x[:,it+1])-npop, "it= ", it, " ", x[:, it], x[:, it+1]);
      }
    }
  }

  // Data for fitting
  for (itype in 1:nobs_types) {
    if (itype == obs_hosp_census) {
      sim_data[itype] = x[Hmodu] + x[Hicuu] + x[Hmodv] + x[Hicuv];
    } else if (itype == obs_icu_census) {
      sim_data[itype] = x[Hicuu] + x[Hicuv];
    } else if (itype == obs_cum_deaths) {
      sim_data[itype] = x[Rmortu] + x[Rmortv];
    } else if (itype == obs_admits) {
      sim_data[itype] = new_admits;
    } else if (itype == obs_cases) {
      sim_data[itype] = new_cases;
    } else if (itype == obs_seroprev) {
      sim_data[itype] = (x[Sv_succ] + x[Rliveu] + x[Rlivev]) / npop;
    } else {
      reject("unexpected itype");
    }
  }
}
model {
  //////////////////////////////////////////
  // prior distributions
  r0 ~ normal(mu_r0, sigma_r0);

  for (iinter in 1:ninter) {
    beta_multiplier[iinter] ~ normal(mu_beta_inter[iinter], sigma_beta_inter[iinter]);
    t_inter[iinter] ~ normal(mu_t_inter[iinter], sigma_t_inter[iinter]);
    len_inter[iinter] ~ normal(mu_len_inter[iinter], sigma_len_inter[iinter]);
  }

  duration_latent ~ normal(mu_duration_latent, sigma_duration_latent);
  duration_rec_mild ~ normal(mu_duration_rec_mild, sigma_duration_rec_mild);
  duration_pre_hosp ~ normal(mu_duration_pre_hosp, sigma_duration_pre_hosp);
  duration_hosp_mod ~ normal(mu_duration_hosp_mod, sigma_duration_hosp_mod);
  duration_hosp_icu ~ normal(mu_duration_hosp_icu, sigma_duration_hosp_icu);
  duration_mort_nonhosp ~ normal(mu_duration_mort_nonhosp, sigma_duration_mort_nonhosp);

  frac_hosp ~ normal(mu_frac_hosp, sigma_frac_hosp);
  frac_icu ~ normal(mu_frac_icu, sigma_frac_icu);
  frac_mort ~ normal(mu_frac_mort, sigma_frac_mort);
  frac_tested ~ normal(mu_frac_tested, sigma_frac_tested);
  frac_mort_nonhosp ~ normal(mu_frac_mort_nonhosp, sigma_frac_mort_nonhosp);

  ini_exposed ~ exponential(lambda_ini_exposed);

  //////////////////////////////////////////
  // fitting observations
  sigma_obs ~ exponential(sigma_obs_est_inv);

  //deal with PUIs and NAs in R
  for (itype in 1:nobs_types) {
    if (nobs[itype] > 0) {
      obs_data[itype, 1:nobs[itype]] ~ normal(sim_data[itype, tobs[itype, 1:nobs[itype]]], sigma_obs[itype]);
    }
  }
}

generated quantities{
  real<lower=0.0> Rt[nt];
  {
    real frac_prehosp;
    real frac_vaccinated;
    real avg_duration;
    real avg_transmission_multiplier;
    for (it in 1:nt) {
      frac_vaccinated = (x[Imildv,it] + x[Iprehv,it]) / (1e-10 + x[Imildu,it] + x[Iprehu,it] + x[Imildv,it] + x[Iprehv,it]); //at t=1 converges to 0
      avg_transmission_multiplier = transmission_multiplier_unvaccinated[it] * (1 - frac_vaccinated) + transmission_multiplier_vaccinated[it] * frac_vaccinated;

      frac_prehosp = (1e-10 * frac_hosp + x[Iprehu, it] + x[Iprehv, it]) / (1e-10 + x[Iprehu, it] + x[Iprehv, it] + x[Imildu, it] + x[Imildv, it]); //at t=1 converges to frac_hosp
      avg_duration = frac_prehosp * duration_pre_hosp + (1 - frac_prehosp) * duration_rec_mild;
      Rt[it] = beta[it] * avg_transmission_multiplier * avg_duration * (x[Su, it] + x[Sv_fail, it]) / npop;
    }
  }
}
