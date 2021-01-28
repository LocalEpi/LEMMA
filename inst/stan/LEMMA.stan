// Modification of Santa Cruz County COVID-19 Model
// https://github.com/jpmattern/seir-covid19
// by Jann Paul Mattern and Mikala Caton

data {
  //////////////////////////////////////////
  // data required to run model

  // obs_data_conf = total hosp census, ICU census, cumulative deaths, cumulative total hosp
  // obs_data_pui =  same

  int<lower=0> nobs_types;
  int<lower=0> nobs;                     // number of timesteps with observations
  int<lower=0> tobs[nobs];               // obs times; this is a vector
  matrix<lower=-1.0>[nobs_types, nobs] obs_data_conf;  // observed confirmed (-1 = NA)
  matrix<lower=-1.0>[nobs_types, nobs] obs_data_pui;   // observed PUI (-1 = NA)

  int<lower=0> nt;                       // number of time steps
  real npop;                             // total population
  int<lower=0, upper=1> extend;

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

  real<lower=0.0> mu_frac_pui[nobs_types];     // mean fraction of PUI that are true COVID+
  real<lower=0.0> sigma_frac_pui[nobs_types];  // sd fraction of PUI that are true COVID+

  real<lower=0.0> mu_frac_hosp;           // mean ICU + non-ICU
  real<lower=0.0> sigma_frac_hosp;        // sd ICU + non-ICU
  real<lower=0.0> mu_frac_icu;            // mean ICU as fraction of hosp
  real<lower=0.0> sigma_frac_icu;         // sd ICU as fraction of hosp
  real<lower=0.0> mu_frac_mort;           // mean mortality as fraction of ICU
  real<lower=0.0> sigma_frac_mort;        // sd mortality as fraction of ICU

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

  int<lower=0> nvac;
  real<lower=0.0> vaccinated_per_day[nvac];
  real<lower=0.0> vaccinated_t[nvac];
  real<lower=0.0> vaccine_efficacy_transmission;
  real<lower=0.0> vaccine_efficacy_susceptible;

}
transformed data {
  //assigning indices for state matrix x
  int Su = 1;
  int Sv = 2;
  int Eu = 3;
  int Ev = 4;
  int Imildu = 5;
  int Imildv = 6;
  int Ipreh = 7;
  int Hmod  = 8;
  int Hicu  = 9;
  int Rliveu = 10;
  int Rlivev = 11;
  int Rmort = 12;

  int ncompartments = 12;

  int obs_hosp_census = 1;
  int obs_icu_census = 2;
  int obs_cum_deaths = 3;
  int obs_cum_admits = 4;

  int nobs_notmissing = 0;

  real beta_limit;

  for (iobs in 1:nobs){
    for (itype in 1:nobs_types) {
      if (obs_data_conf[itype, iobs] > 0) {
        nobs_notmissing = nobs_notmissing + 1;
      }
    }
  }

  if (extend == 1) {
    beta_limit = 1e10; // no limit on beta if extending simulation
  } else {
    beta_limit = 2.0;
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

  real<lower=0> ini_exposed;

  real<lower=0> sigma_obs[nobs_types];


  real<lower=0.0> r0;
  real<lower=0.0> beta_multiplier[ninter];
  real<lower=1.0> t_inter[ninter];
 // real<lower=1.0> len_inter[ninter];

  // real<lower=0, upper=1> frac_PUI[nobs_types];
}
transformed parameters {
  matrix[ncompartments,nt] x;
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
    real newrec_mod;
    real newhosp;
    real leave_icu;
    real beta_0;
    real obs;
    real sim;
    real zero;
    real vaccinated;
    real frac_vac_S;
    real newSv;
    real newRlivev;

    //////////////////////////////////////////
    // Calculate beta for each time point
    beta_0 = r0 / (frac_hosp * duration_pre_hosp + (1 - frac_hosp) * duration_rec_mild);
    for (it in 1:nt) {
      beta[it] = beta_0;
      for (iinter in 1:ninter) {
        //k <- 2/s * qlogis(0.99) # = -2/s * qlogis(0.01) --> 2 * qlogis(0.99) = 9.19024
        //f <- m ^ plogis(k * (t - (d + s/2)))
        beta[it] = beta[it] * beta_multiplier[iinter] ^ inv_logit(9.19024 / mu_len_inter[iinter] * (it - (t_inter[iinter] + mu_len_inter[iinter] / 2))); //TODO: document this; maybe could speed up too
      }
    }

    // initial cond
    zero = ini_exposed * 1e-15; //should be zero, hack for RStan (causes problems with RHat if constant)
    x[:,1] = rep_vector(zero, ncompartments); //: means all entries. puts a zero in x1-8 for initial entries
    x[Su,1] = npop-ini_exposed;
    x[Eu,1] = ini_exposed;
    Hadmits[1] = zero;

    //////////////////////////////////////////
    // the SEIR model
    for (it in 1:nt-1){
      //////////////////////////////////////////
      // set transition variables
      newEu = fmin(x[Su,it],  x[Su,it] * beta[it]/npop * (x[Imildu,it] + x[Ipreh,it] + (1 - vaccine_efficacy_transmission) * x[Imildv,it]));
      newEv = fmin(x[Sv,it],  x[Sv,it] * beta[it]/npop * (x[Imildu,it] + x[Ipreh,it] + (1 - vaccine_efficacy_transmission) * x[Imildv,it]));


      if (it > 1 && it < 200 && extend == 0) {
        newE_temp[it] = newEu;
      } else {
        newE_temp[it] = 1 + zero; //ignore this
      }

      vaccinated = 0; //successfully vaccinated on this day
      for (ivac in 1:nvac) {
        if (it >= vaccinated_t[ivac]) {
          vaccinated = vaccinated + vaccine_efficacy_susceptible * vaccinated_per_day[ivac];
        }
      }
      vaccinated = fmin(vaccinated, x[Su, it] + x[Eu, it] + x[Rliveu, it]);


      newIu = 1.0/duration_latent * x[Eu, it];
      newIv = 1.0/duration_latent * x[Ev, it];
      newhosp = 1.0/duration_pre_hosp * x[Ipreh,it];
      newrecu_mild = 1.0/duration_rec_mild * x[Imildu,it];
      newrecv_mild = 1.0/duration_rec_mild * x[Imildv,it];
      newrec_mod = 1.0/duration_hosp_mod * x[Hmod,it];
      leave_icu = 1.0/duration_hosp_icu * x[Hicu, it];
      frac_vac_S = x[Su, it] / (x[Su, it] + x[Eu, it] + x[Rliveu, it]);
      newSv = vaccinated * frac_vac_S;
      newRlivev = vaccinated * (1 - frac_vac_S);

      //////////////////////////////////////////
      // S -> E -> I

      x[Su, it+1] = x[Su, it] - newEu - newSv;
      x[Sv, it+1] = x[Sv, it] - newEv + newSv;
      x[Eu, it+1] = x[Eu, it] + newEu - newIu;
      x[Ev, it+1] = x[Ev, it] + newEv - newIv;
      x[Imildu, it+1] = x[Imildu, it] + newIu * (1 - frac_hosp) - newrecu_mild;
      x[Imildv, it+1] = x[Imildv, it] + newIv - newrecv_mild; //100% of vax are mild
      x[Ipreh, it+1] = x[Ipreh, it] + newIu * frac_hosp - newhosp;
      x[Hmod, it+1] = x[Hmod, it] + newhosp * (1 - frac_icu) - newrec_mod;
      x[Hicu, it+1] = x[Hicu, it] + newhosp * frac_icu - leave_icu;
      x[Rliveu, it+1] = x[Rliveu, it] + newrecu_mild + newrec_mod + leave_icu * (1 - frac_mort) - newRlivev;
      x[Rlivev, it+1] = x[Rlivev, it] + newrecv_mild + newRlivev;
      x[Rmort, it+1] = x[Rmort, it] + leave_icu * frac_mort;

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
  //////////////////////////////////////////
  // prior distributions
  r0 ~ normal(mu_r0, sigma_r0);
  // frac_PUI ~ normal(mu_frac_pui, sigma_frac_pui);

  for (iinter in 1:ninter) {
    beta_multiplier[iinter] ~ normal(mu_beta_inter[iinter], sigma_beta_inter[iinter]);
    t_inter[iinter] ~ normal(mu_t_inter[iinter], sigma_t_inter[iinter]);
    // len_inter[iinter] ~ normal(mu_len_inter[iinter], sigma_len_inter[iinter]);
  }

  duration_latent ~ normal(mu_duration_latent, sigma_duration_latent);
  duration_rec_mild ~ normal(mu_duration_rec_mild, sigma_duration_rec_mild);
  duration_pre_hosp ~ normal(mu_duration_pre_hosp, sigma_duration_pre_hosp);
  duration_hosp_mod ~ normal(mu_duration_hosp_mod, sigma_duration_hosp_mod);
  duration_hosp_icu ~ normal(mu_duration_hosp_icu, sigma_duration_hosp_icu);

  frac_hosp ~ normal(mu_frac_hosp, sigma_frac_hosp);
  frac_icu ~ normal(mu_frac_icu, sigma_frac_icu);
  frac_mort ~ normal(mu_frac_mort, sigma_frac_mort);

  ini_exposed ~ exponential(lambda_ini_exposed);

  //////////////////////////////////////////
  // fitting observations
  sigma_obs ~ exponential(1.0);
  {
    vector[nobs_notmissing] error;
    real obs;
    real sim;
    int cnt;
    real scale = npop / 1000000;

    cnt = 0;
    for (iobs in 1:nobs){
      for (itype in 1:nobs_types) {
        if (obs_data_conf[itype, iobs] > 0) {
          cnt = cnt + 1;
          obs = obs_data_conf[itype, iobs];
          if (obs_data_pui[itype, iobs] > 0) {
            obs = obs + obs_data_pui[itype, iobs] * mu_frac_pui[itype];
          }
          sim = sim_data[itype, tobs[iobs]];
          error[cnt] = (obs - sim) / (sigma_obs[itype] * scale);
        }
      }
    }
    error ~ std_normal();
  }
}

generated quantities{
  real<lower=0.0> Rt[nt];
  real<lower=0.0> Rt_unvac[nt];
  real<lower=0.0> frac_inf_vac;
  for (it in 1:nt) {
    frac_inf_vac = x[Imildv, it] / (x[Imildv, it] + x[Imildu, it] + x[Ipreh, it]);
    Rt_unvac[it] = beta[it] * (frac_hosp * duration_pre_hosp + (1 - frac_hosp) * duration_rec_mild) * (x[Su, it] + x[Sv, it]) / npop; //not quite right because the fraction that goes to hospital changes as more are vaccinated (but duration_pre_hosp is ~7 and duration_rec_mild is ~5, not much difference)
    Rt[it] = (1 - frac_inf_vac * (1 - vaccine_efficacy_transmission)) * Rt_unvac[it];
  }
}
