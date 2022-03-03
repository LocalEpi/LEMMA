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

  real<lower=0.0> mu_frac_tested;         // mean of infected who test positive
  real<lower=0.0> sigma_frac_tested;      // sd of infected who test positiveof infected who test positive
  real<lower=0.0> mu_test_delay;
  real<lower=0.0> sigma_test_delay;
  real<lower=0.0> mu_severity;
  real<lower=0.0> sigma_severity;
  real<lower=0.0> mu_omicron_trans_multiplier;
  real<lower=0.0> sigma_omicron_trans_multiplier;
  real<lower=0.0> mu_hosp_delta;
  real<lower=0.0> sigma_hosp_delta;
  real<lower=0.0> mu_cases_delta;
  real<lower=0.0> sigma_cases_delta;
  real<lower=0.0> mu_duration_protection_infection;
  real<lower=0.0> sigma_duration_protection_infection;



  real<lower=0.0> lambda_initial_infected;     // parameter for initial conditions of "exposed"

  real<lower=0.0> sigma_obs_est_inv[nobs_types];

  real<lower=0.0> frac_hosp_lemma;
  real<lower=0.0> VE_infection;
  real<lower=0.0> VE_infection_delta;
  int<lower=0> holiday_start;
  int<lower=0> holiday_end;
  real<lower=0.0> holiday_multiplier;
  real<lower=0.0> omicron_recovered_booster_scale;
  real<lower=0.0> num_boosters[nt];
  real<lower=0.0> booster_VE_infection;
  real<lower=0.0> booster_VE_severe;
  real<lower=0.0> frac_incidental_omicron;
  real<lower=0.0> VE_severe_given_infection_0;
  real<lower=0.0> hosp_frac_delta_0;
  real<lower=0.0> case_frac_delta_0;
  real<lower=0.0> omicron_growth;

  int<lower=0> ninter;                      // number of interventions
  real<lower=1.0> t_inter[ninter];       // mean start time of each interventions
  real<lower=1.0> len_inter[ninter];     // mean length of each intervention
  real<lower=0.0> mu_beta_inter[ninter];    // mean change in beta through intervention
  real<lower=0.0> sigma_beta_inter[ninter]; // sd change in beta through intervention
}

transformed data {
  //assigning indices for state matrix x
  int S = 1;
  int E = 2;
  int Imild = 3;
  int Ipreh = 4;
  int Hmod  = 5;
  int P = 6;
  int Rlive = 7;

  int ncompartments = 7;

  int obs_hosp_census = 1;
  int obs_cases = 2;
}
parameters {
  real<lower=1.0> duration_latent; // duration is a minimum of 1 which is the stepsize of this model
  real<lower=1.0> duration_rec_mild;
  real<lower=1.0> duration_pre_hosp;
  real<lower=1.0> duration_hosp_mod;
  real<lower=1.0> duration_protection_infection;

  real<lower=0.0, upper=1.0> frac_tested;

  real<lower=0.0> initial_infected;
  real<lower=0.0> sigma_obs[nobs_types];

  real<lower=0.0> severity;
  real<lower=0.0> omicron_trans_multiplier;
  real<lower=1.0> test_delay; //days tested after exposed

  real<lower=0.0> hosp_delta;
  real<lower=0.0> cases_delta;

  // real<lower=0.0> beta_mult1;
  // real<lower=0.0> beta_mult2;
  // real<lower=0.0> beta_mult3;
  real<lower=0.0> beta_multiplier[ninter];
}
transformed parameters {
  matrix<lower=0.0>[ncompartments,nt] x;
  matrix<lower=0.0>[nobs_types,nt] sim_data;
  row_vector<lower=0.0>[nt] new_cases;
  row_vector<lower=0.0>[nt] soon_positive;
  real<lower=0.0> beta[nt];
  real<lower=0.0> frac_hosp_0;

  for (it in 1:nt) {
    for (itype in 1:nobs_types) {
      sim_data[itype, it] = 0.0;
    }
  }

  {
    // variables in curly brackets will not have output, they are local variables
    real newE;
    real newI;
    real frac_init_E;
    real beta_0;
    real avg_duration;
    real frac_hosp;
    real frac_boosters_to_susceptible;
    real new_protected;
    real increased_severity_protection;
    real frac_increased_severity_protection;
    real new_admits;
    real VE_severe_given_infection;
    real lost_protection_infection;
    row_vector[nt]  hosp_frac_delta;
    row_vector[nt]  case_frac_delta;

    for (it in 1:nt) {
      hosp_frac_delta[it] = 1 / (1 + hosp_frac_delta_0 * exp(it * omicron_growth));
      case_frac_delta[it] = 1 / (1 + case_frac_delta_0 * exp(it * omicron_growth));
    }

    VE_severe_given_infection = VE_severe_given_infection_0;

    frac_hosp_0 = frac_hosp_lemma * (1 - VE_severe_given_infection) * severity;
    avg_duration = frac_hosp_0 * duration_pre_hosp + (1 - frac_hosp_0) * duration_rec_mild;

    // print(frac_hosp_0, frac_hosp_lemma, VE_severe_given_infection, severity, avg_duration, duration_pre_hosp, duration_rec_mild)


    x[:,1] = rep_vector(0.0, ncompartments);
    frac_init_E = duration_latent / (duration_latent + avg_duration);
    x[E, 1] = frac_init_E * initial_infected;
    x[Imild, 1] = (1 - frac_init_E) * initial_infected * (1 - frac_hosp_0);
    x[Ipreh, 1] = (1 - frac_init_E) * initial_infected * frac_hosp_0;
    x[P, 1] = npop * VE_infection;
    x[S, 1] = npop - initial_infected - x[P, 1];

    soon_positive[1] = 0.0;
    new_cases[1] = 0.0; //not correct - pass NA obs cases at t = 1 so no fitting to this
    if (tobs[obs_cases, 1] == 1) {
      reject("Minimum tobs[obs_cases, :] is 2")
    }

    //beta_0 = omicron_trans_multiplier / avg_duration * rt_delta / (1 - VE_infection_delta);
    //for now assume rt_delta = 1
    beta_0 = omicron_trans_multiplier / avg_duration * 1.0 / (1 - VE_infection_delta);
    for (it in 1:nt) {
      beta[it] = beta_0;
      for (iinter in 1:ninter) {
        //k <- 2/s * qlogis(0.99) # = -2/s * qlogis(0.01) --> 2 * qlogis(0.99) = 9.19024
        //f <- m ^ plogis(k * (t - (d + s/2)))
        beta[it] = beta[it] * beta_multiplier[iinter] ^ inv_logit(9.19024 / len_inter[iinter] * (it - (t_inter[iinter] + len_inter[iinter] / 2)));
      }
    }

    for (it in 1:nt-1){
      newE = fmin(x[S, it], beta[it] * x[S, it] * (x[Imild, it] + x[Ipreh, it]) / npop);
      newI = x[E, it] / duration_latent;

      frac_boosters_to_susceptible = x[S, it] / (x[S, it] + omicron_recovered_booster_scale * x[Rlive, it]);
      new_protected = fmin(x[S, it] - newE, num_boosters[it] * frac_boosters_to_susceptible * booster_VE_infection);
      if (new_protected < 0) reject("new_protected < 0")


      increased_severity_protection = num_boosters[it] * frac_boosters_to_susceptible * (1 - booster_VE_infection); //boosted but not Protected from infection, still susceptible
      frac_increased_severity_protection = increased_severity_protection / fmax(x[S, it], 0.001); //avoid 0/0 problems

      //increase VE_severe_given_infection for boosters
      VE_severe_given_infection = booster_VE_severe * frac_increased_severity_protection +   VE_severe_given_infection * (1 - frac_increased_severity_protection);
      frac_hosp = frac_hosp_lemma * (1 - VE_severe_given_infection) * severity;
      new_admits = x[Ipreh, it] / duration_pre_hosp;

      lost_protection_infection = x[P, it] / duration_protection_infection;

      x[S, it + 1] = x[S, it] - new_protected - newE + lost_protection_infection;
      x[E, it + 1] = x[E, it] + newE - newI;
      x[Imild, it + 1] = x[Imild, it] + newI * (1 - frac_hosp) - x[Imild, it] / duration_rec_mild;
      x[Ipreh, it + 1] = x[Ipreh, it] + newI * frac_hosp - new_admits;
      x[Hmod, it + 1] = x[Hmod, it] + new_admits - x[Hmod, it] / duration_hosp_mod;
      x[Rlive, it + 1] = x[Rlive, it] + x[Hmod, it] / duration_hosp_mod + x[Imild, it] / duration_rec_mild;
      x[P, it + 1] = x[P, it] + new_protected - lost_protection_infection;

      soon_positive[it + 1] = soon_positive[it] + newE * frac_tested - soon_positive[it] / test_delay;
      new_cases[it + 1] = soon_positive[it + 1] / test_delay;
      if (is_nan(new_cases[it + 1])) {
        print("is nan cases:")
        print(it, " ", soon_positive[it], " ", soon_positive[it + 1], " ", newE, " ", frac_tested, " ", test_delay)
        reject("new_cases[it+1] is nan")
      }

      // if (is_nan(sum(x[:, it + 1]))) {
      //   print("is nan:")
      //   print(newE, " ", newI, " ", frac_init_E, " ", beta[it], " ", beta_0, " ", avg_duration, " ", frac_hosp_0, " ", frac_boosters_to_susceptible, " ", new_protected, " ", increased_severity_protection, " ", frac_increased_severity_protection, " ", frac_hosp, " ", new_admits, " ", VE_severe_given_infection, " ");
      //   reject("x is nan: it = ", it, " x[, it+1] = ", x[:, it + 1]);
      // }
      // test
      if (fabs(sum(x[:, it + 1])-npop) > 0.01) {
        reject("Model is leaking, net gain: ", sum(x[:,it+1])-npop);
        // reject("Model is leaking, net gain: ", sum(x[:,it+1])-npop, "it= ", it, " ", x[:, it], x[:, it+1]);
      }
    }

    // Data for fitting
    sim_data[obs_hosp_census] = x[Hmod] * (1 + frac_incidental_omicron) + hosp_delta * hosp_frac_delta;
    sim_data[obs_cases] = new_cases + cases_delta * case_frac_delta;
  }

  // print(x)
  // print(sim_data)
}

model {
  //////////////////////////////////////////
  // prior distributions
  //
  duration_latent ~ normal(mu_duration_latent, sigma_duration_latent);
  duration_rec_mild ~ normal(mu_duration_rec_mild, sigma_duration_rec_mild);
  duration_pre_hosp ~ normal(mu_duration_pre_hosp, sigma_duration_pre_hosp);
  duration_hosp_mod ~ normal(mu_duration_hosp_mod, sigma_duration_hosp_mod);
  duration_protection_infection ~ normal(mu_duration_protection_infection, sigma_duration_protection_infection);

  for (iinter in 1:ninter) {
    beta_multiplier[iinter] ~ normal(mu_beta_inter[iinter], sigma_beta_inter[iinter]);
  }

  frac_tested ~ normal(mu_frac_tested, sigma_frac_tested);
  severity ~ normal(mu_severity, sigma_severity);
  omicron_trans_multiplier ~ normal(mu_omicron_trans_multiplier, sigma_omicron_trans_multiplier);
  test_delay ~ normal(mu_test_delay, sigma_test_delay);
  hosp_delta ~ normal(mu_hosp_delta, sigma_hosp_delta);
  cases_delta ~ normal(mu_cases_delta, sigma_cases_delta);

  initial_infected ~ exponential(lambda_initial_infected);

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
  matrix<lower=0.0>[nobs_types,nt] sim_data_with_error;
  real<lower=0.0> Rt[nt];

  for (itype in 1:nobs_types) {
    if (nobs[itype] > 0) {
      for (it in 1:nt) {
        sim_data_with_error[itype, it] = fmax(0.0, normal_rng(sim_data[itype, it], sigma_obs[itype])); //can't vectorize - fmax not vectorized?
      }
    } else {
      sim_data_with_error[itype] = sim_data[itype];
    }
  }
  {
    real frac_prehosp;
    real avg_duration;
    for (it in 1:nt) {
      frac_prehosp = (1e-10 * frac_hosp_0 + x[Ipreh, it]) / (1e-10 + x[Ipreh, it] +  x[Imild, it]); //at t=1 converges to frac_hosp_0
      avg_duration = frac_prehosp * duration_pre_hosp + (1 - frac_prehosp) * duration_rec_mild;
      Rt[it] = beta[it] * avg_duration * x[S, it] / npop;
    }
  }
}
