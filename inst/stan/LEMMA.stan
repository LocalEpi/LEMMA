// Modification of Santa Cruz County COVID-19 Model
// https://github.com/jpmattern/seir-covid19
// by Jann Paul Mattern and Mikala Caton

data {
  //////////////////////////////////////////
  // data required to run model

  // obs_data_conf = total hosp census, ICU census, cumulative deaths, cumulative total hosp


  //want to have obs_data and tobs different sizes for each obs type

  int<lower=0> nobs_types;
  int<lower=0> nobs[nobs_types];                     // number of timesteps with observations
  int <lower=0> nobs_max;
  int<lower=-1> tobs[nobs_types, nobs_max];               // obs times; -1 = NA
  matrix<lower=-1.0>[nobs_types, nobs_max] obs_data;  // observed confirmed (-1 = NA)

  int<lower=0> nt;                       // number of time steps
  real<lower=0.0> npop;                             // total population
  int<lower=0, upper=1> extend;
  int<lower=0, upper=1> from_beginning;

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

  real<lower=0.0> mu_r0;                  // mean initial beta estimate
  real<lower=0.0> sigma_r0;               // sd initial beta estimate

  real<lower=0.0> mu_frac_hosp;           // mean ICU + non-ICU
  real<lower=0.0> sigma_frac_hosp;        // sd ICU + non-ICU
  real<lower=0.0> mu_frac_icu;            // mean ICU as fraction of hosp
  real<lower=0.0> sigma_frac_icu;         // sd ICU as fraction of hosp
  real<lower=0.0> mu_frac_mort;           // mean mortality as fraction of ICU
  real<lower=0.0> sigma_frac_mort;        // sd mortality as fraction of ICU

  //real<lower=0.0> lambda_ini_exposed;     // parameter for initial conditions of "exposed"
  real<lower=0.0> mu_iniE;
  real<lower=0.0> sigma_iniE;
  real<lower=0.0> mu_ini_Imild;
  real<lower=0.0> sigma_ini_Imild;
  real<lower=0.0> mu_ini_Ipreh;
  real<lower=0.0> sigma_ini_Ipreh;
  real<lower=0.0> mu_ini_Rlive;
  real<lower=0.0> sigma_ini_Rlive;

  //////////////////////////////////////////
  // interventions

  int<lower=0> ninter;                      // number of interventions
  real<lower=1.0> mu_t_inter[ninter];       // mean start time of each interventions
  real<lower=0.0> sigma_t_inter[ninter];    // sd start time of each interventions
  real<lower=1.0> len_inter[ninter];     // mean length of each intervention
  real<lower=0.0> mu_beta_inter[ninter];    // mean change in beta through intervention
  real<lower=0.0> sigma_beta_inter[ninter];    // sd change in beta through intervention

  real<lower=0.0> vaccinated_per_day[nt]; //vaccinated per day total
  real<lower=0.0, upper=0.0> vaccine_efficacy_for_susceptibility[nt];
  real<lower=0.0, upper=0.0> vaccine_efficacy_against_infectiousness[nt];
  real<lower=0.0, upper=0.0> vaccine_efficacy_against_progression[nt];

  real<lower=1.0> duration_vaccinated[nt];
  real<lower=1.0> duration_natural[nt];
  real<lower=0.0> frac_hosp_multiplier[nt];  //multiplier due to vaccines/variants
  real<lower=0.0> frac_icu_multiplier[nt];   //multiplier due to vaccines/variants
  real<lower=0.0> frac_mort_multiplier[nt];   //multiplier due to vaccines/variants
  real<lower=0.0> sigma_obs_est_inv[nobs_types];

  // VE_S vaccine efficacy for susceptibility
  // VE_I vaccine efficacy against infectiousness
  // VE_P vaccine efficacy against progression
  // "x" transmission reduction factor
  // x = 1 – (1 – VE_S)(1 – VE_I)

//VE_S | VE_P = VE_S & VE_P / VE_P = VE_S * VE_P | VE_S / VE_P = VE_S / VE_P - if not susceptible, can't progress VE_P | VE_S = 1
//VE_I | VE_P = VE_I & VE_P / VE_P = VE_I * VE_P | VE_I / VE_P = VE_I / VE_P - assume if can't transmit can't progress VE_P | VE_I = 1
// require VE_S < VE_P, VE_I < VE_P
}

transformed data {
  //assigning indices for state matrix x
  int Su = 1;
  int Sv = 2;
  int Eu = 3;
  int Ev = 4;
  int Imildu = 5;
  int Imildv = 6;
  int Iprehu = 7;
  int Iprehv = 8;
  int Hmodu  = 9;
  int Hmodv  = 10;
  int Hicuu  = 11;
  int Hicuv  = 12;
  int Rliveu = 13;
  int Rlivev = 14;
  int Rmort = 15;

  int ncompartments = 15;

  int obs_hosp_census = 1;
  int obs_icu_census = 2;
  int obs_cum_deaths = 3;
  int obs_cum_admits = 4;

  int nobs_notmissing = 0;

  real beta_limit;

  real<lower=0.0, upper=1.0> VE_IgP[nt]; //VE_I given VE_P
  real<lower=0.0, upper=1.0> VE_SgP[nt]; //VE_S given VE_P
  real<lower=0.0, upper=1.0> transmission_reduction_factor[nt];

  if (extend == 1) {
    beta_limit = 1e10; // no limit on beta if extending simulation
  } else {
    beta_limit = 2.0;
  }

  for (it in 1:nt) {
    VE_IgP[it] = vaccine_efficacy_against_infectiousness[it] / vaccine_efficacy_against_progression[it];
    VE_SgP[it] = vaccine_efficacy_for_susceptibility[it] / vaccine_efficacy_against_progression[it];
    transmission_reduction_factor[it] = 1 - (1 - vaccine_efficacy_against_infectiousness[it]) * (1 - vaccine_efficacy_for_susceptibility[it]);
  }


}
parameters {

  real<lower=1.0> duration_latent; // duration is a minimum of 1 which is the stepsize of this model
  real<lower=1.0> duration_rec_mild;
  real<lower=1.0> duration_pre_hosp;
  real<lower=1.0> duration_hosp_mod;
  real<lower=1.0> duration_hosp_icu;

  real<lower=0.005, upper=1.0> frac_hosp;
  real<lower=0.0, upper=1.0> frac_icu;
  real<lower=0.0, upper=1.0> frac_mort;

  real<lower=0.0> ini_E;
  real<lower=0.0> ini_Imild[1 - from_beginning];
  real<lower=0.0> ini_Ipreh[1 - from_beginning];
  real<lower=0.0> ini_Rlive[1 - from_beginning];
  real<lower=0.0> sigma_obs[nobs_types];


  real<lower=0.0> r0;
  real<lower=0.0> beta_multiplier[ninter];
  real<lower=1.0> t_inter[ninter];

}
transformed parameters {
  matrix<lower=0.0>[ncompartments,nt] x;
  matrix<lower=0.0>[nobs_types,nt] sim_data;
  real<lower=0.0, upper=beta_limit> beta[nt];
  row_vector<lower=0.0>[nt] Hadmits;
  real<lower=1e-10> newE_temp[nt-1];

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
    real leave_icuu;
    real leave_icuv;
    real beta_0;
    real zero;
    real vaccinated;
    real frac_vac_S;
    real newSv;
    real newRlivev;
    real S_lostv;
    real R_lostv;
    real R_lostnatu;
    real R_lostnatv;

    //////////////////////////////////////////
    // Calculate beta for each time point
    beta_0 = r0 / (frac_hosp * duration_pre_hosp + (1 - frac_hosp) * duration_rec_mild);
    for (it in 1:nt) {
      beta[it] = beta_0;
      for (iinter in 1:ninter) {
        //k <- 2/s * qlogis(0.99) # = -2/s * qlogis(0.01) --> 2 * qlogis(0.99) = 9.19024
        //f <- m ^ plogis(k * (t - (d + s/2)))
        beta[it] = beta[it] * beta_multiplier[iinter] ^ inv_logit(9.19024 / len_inter[iinter] * (it - (t_inter[iinter] + len_inter[iinter] / 2))); //TODO: document this; maybe could speed up too
      }
    }

    // initial cond
    zero = ini_E * 1e-15; //should be zero, hack for RStan (causes problems with RHat if constant)
    x[:,1] = rep_vector(zero, ncompartments); //: means all entries. puts a zero in x1-8 for initial entries

    x[Eu,1] = ini_E;
    if (from_beginning == 0) {
      x[Hmodu, 1] = obs_data[obs_hosp_census, 1]; //FIXME - works as special case only (needs non NA value)
      x[Hicuu, 1] = obs_data[obs_icu_census, 1]; //FIXME - works as special case only (needs non NA value)
      x[Rmort, 1] = obs_data[obs_cum_deaths, 1]; //FIXME - works as special case only (needs non NA value)
      x[Imildu, 1] = ini_Imild[1];
      x[Iprehu, 1] = ini_Ipreh[1];
      x[Rliveu, 1] = ini_Rlive[1];
    }

    x[Su,1] = npop - sum(x[2:ncompartments, 1]);

    Hadmits[1] = zero; //FIXME - this is wrong

    //////////////////////////////////////////
    // the SEIR model
    for (it in 1:nt-1){
      //////////////////////////////////////////
      // set transition variables
      newEu = fmin(x[Su,it],
        x[Su,it] * beta[it] / npop *
        (x[Imildu,it] + x[Iprehu,it] + x[Imildv,it] + x[Iprehv,it]));

      newEv = fmin(x[Sv,it],
        (1 - vaccine_efficacy_for_susceptibility[it]) * x[Sv,it] * beta[it] / npop *
        (x[Imildu,it] + x[Iprehu,it] + x[Imildv,it] + x[Iprehv,it]));


      if (it > 1 && it < 200 && extend == 0) {
        newE_temp[it] = newEu;
      } else {
        newE_temp[it] = 1 + zero; //ignore this
      }

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
      frac_vac_S = x[Su, it] / (x[Su, it] + x[Eu, it] + x[Rliveu, it]); //approx - ignores I and splits between S and Rlive - small magnitude
      newSv = vaccinated * frac_vac_S;
      newRlivev = vaccinated * (1 - frac_vac_S);

      S_lostv = 1.0/duration_vaccinated[it] * x[Sv, it];
      R_lostv = 1.0/duration_vaccinated[it] * x[Rlivev, it];
      R_lostnatu = 1.0/duration_natural[it] * x[Rliveu, it];
      R_lostnatv = 1.0/duration_natural[it] * x[Rlivev, it];

      //////////////////////////////////////////
      // S -> E -> I

      x[Su, it+1] = x[Su, it] - newEu - newSv + S_lostv + R_lostnatu;
      x[Sv, it+1] = x[Sv, it] - newEv + newSv - S_lostv + R_lostnatv;
      x[Eu, it+1] = x[Eu, it] + newEu - newIu;
      x[Ev, it+1] = x[Ev, it] + newEv - newIv;
      x[Imildu, it+1] = x[Imildu, it] + newIu * (1 - frac_hosp * frac_hosp_multiplier[it]) - newrecu_mild;
      x[Imildv, it+1] = x[Imildv, it] + newIv - newrecv_mild; //100% of vax are mild
      x[Ipreh, it+1] = x[Ipreh, it] + newIu * frac_hosp * frac_hosp_multiplier[it] - newhosp;
      x[Hmod, it+1] = x[Hmod, it] + newhosp * (1 - frac_icu * frac_icu_multiplier[it]) - newrec_mod;
      x[Hicu, it+1] = x[Hicu, it] + newhosp * frac_icu * frac_icu_multiplier[it] - leave_icu;
      x[Rliveu, it+1] = x[Rliveu, it] + newrecu_mild + newrec_mod + leave_icu * (1 - frac_mort * frac_mort_multiplier[it]) - newRlivev + R_lostv - R_lostnatu;
      x[Rlivev, it+1] = x[Rlivev, it] + newrecv_mild + newRlivev - R_lostv - R_lostnatv;
      x[Rmort, it+1] = x[Rmort, it] + leave_icu * frac_mort * frac_mort_multiplier[it];

      // cumulative hospital admissions
      Hadmits[it+1] = Hadmits[it] + newhosp;

      //////////////////////////////////////////
      // test
      if (fabs(sum(x[:,it+1])-npop) > 1e-1){
        reject("Model is leaking, net gain: ", sum(x[:,it+1])-npop)
      }
    }
  }

  // Data for fitting
  for (itype in 1:nobs_types) {
    if (itype == obs_hosp_census) {
      sim_data[itype] = x[Hmod] + x[Hicu];
    } else if (itype == obs_icu_census) {
      sim_data[itype] = x[Hicu];
    } else if (itype == obs_cum_deaths) {
      sim_data[itype] = x[Rmort];
    } else if (itype == obs_cum_admits) {
      sim_data[itype] = Hadmits;
    } else {
      reject("unexpected itype")
    }
  }
}
model {
  if (extend == 0) {
    //////////////////////////////////////////
    // prior distributions
    r0 ~ normal(mu_r0, sigma_r0);

    for (iinter in 1:ninter) {
      beta_multiplier[iinter] ~ normal(mu_beta_inter[iinter], sigma_beta_inter[iinter]);
      t_inter[iinter] ~ normal(mu_t_inter[iinter], sigma_t_inter[iinter]);
    }

    duration_latent ~ normal(mu_duration_latent, sigma_duration_latent);
    duration_rec_mild ~ normal(mu_duration_rec_mild, sigma_duration_rec_mild);
    duration_pre_hosp ~ normal(mu_duration_pre_hosp, sigma_duration_pre_hosp);
    duration_hosp_mod ~ normal(mu_duration_hosp_mod, sigma_duration_hosp_mod);
    duration_hosp_icu ~ normal(mu_duration_hosp_icu, sigma_duration_hosp_icu);

    frac_hosp ~ normal(mu_frac_hosp, sigma_frac_hosp);
    frac_icu ~ normal(mu_frac_icu, sigma_frac_icu);
    frac_mort ~ normal(mu_frac_mort, sigma_frac_mort);

    if (from_beginning == 1) {
      ini_E ~ exponential(1.0 / mu_iniE);
    } else {
      ini_E ~ normal(mu_iniE, sigma_iniE);
      ini_Imild ~ normal(mu_ini_Imild, sigma_ini_Imild);
      ini_Ipreh ~ normal(mu_ini_Ipreh, sigma_ini_Ipreh);
      ini_Rlive ~ normal(mu_ini_Rlive, sigma_ini_Rlive);
    }


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
}

generated quantities{
  real<lower=0.0> Rt[nt];
  real<lower=0.0> Rt_unvac[nt];
  real<lower=0.0> frac_vac;
  for (it in 1:nt) {
    frac_vac = (x[Sv, it] + x[Ev, it] + x[Imildv, it]) / (x[Su, it] + x[Eu, it] + x[Imildu, it] + x[Sv, it] + x[Ev, it] + x[Imildv, it] + x[Ipreh, it]);
    Rt_unvac[it] = beta[it] * (frac_hosp * frac_hosp_multiplier[it] * duration_pre_hosp + (1 - frac_hosp * frac_hosp_multiplier[it]) * duration_rec_mild) * (x[Su, it] + x[Sv, it]) / npop;
    Rt[it] = (1 - frac_vac * transmission_reduction_factor[it]) * Rt_unvac[it];
  }
}
