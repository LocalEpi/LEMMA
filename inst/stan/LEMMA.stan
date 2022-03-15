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
  real<lower=1.0> mu_duration_latent1;     // mean duration in "exposed" stage
  real<lower=0.0> sigma_duration_latent1;  // sd duration in "exposed" stage
  real<lower=1.0> mu_duration_latent2;     // mean duration in "exposed" stage
  real<lower=0.0> sigma_duration_latent2;  // sd duration in "exposed" stage
  real<lower=1.0> mu_duration_rec_mild;   // mean duration in "infectious" stage for mild cases
  real<lower=0.0> sigma_duration_rec_mild;// sd duration in "infectious" stage for mild cases
  real<lower=1.0> mu_duration_pre_hosp;   // mean duration in "infectious" stage for hospitalized cases
  real<lower=0.0> sigma_duration_pre_hosp;// sd duration in "infectious" stage for hospitalized cases
  real<lower=1.0> mu_duration_hosp_mod1;   // mean duration in hospital for non-ICU cases
  real<lower=0.0> sigma_duration_hosp_mod1;// sd duration in hospital for non-ICU cases
  real<lower=1.0> mu_duration_hosp_mod2;   // mean duration in hospital for non-ICU cases
  real<lower=0.0> sigma_duration_hosp_mod2;// sd duration in hospital for non-ICU cases

  real<lower=0.0> mu_frac_tested;         // mean of infected who test positive
  real<lower=0.0> sigma_frac_tested;      // sd of infected who test positiveof infected who test positive
  real<lower=0.0> mu_test_delay;
  real<lower=0.0> sigma_test_delay;
  real<lower=0.0> mu_beta0;
  real<lower=0.0> sigma_beta0;
  real<lower=0.0> mu_frac_hosp1_naive;
  real<lower=0.0> sigma_frac_hosp1_naive;
  real<lower=0.0> mu_trans_multiplier;
  real<lower=0.0> sigma_trans_multiplier;
  real<lower=0.0> mu_frac_hosp_multiplier;
  real<lower=0.0> sigma_frac_hosp_multiplier;

  real<lower=0.0> mu_duration_protection_infection;
  real<lower=0.0> sigma_duration_protection_infection;

  real<lower=0.0> lambda_initial_infected1;     // parameter for initial conditions of "exposed"
  real<lower=0.0> lambda_initial_infected2;     // parameter for initial conditions of "exposed"
  // real<lower=0.0> initial_infected2_fraction;

  real<lower=0.0> sigma_obs_est_inv[nobs_types];

  real<lower=0.0> VE_infection1_init;
  real<lower=0.0> VE_infection2_init;
  real<lower=0.0> init_hosp1;

  // real<lower=0.0> omicron_recovered_booster_scale;
  // real<lower=0.0> num_boosters[nt];
  // real<lower=0.0> booster_VE_infection;
  // real<lower=0.0> booster_VE_severe;
  real<lower=0.0> frac_incidental1;
  real<lower=0.0> frac_incidental2;
  real<lower=0.0> VE_severe_given_infection1_init;
  real<lower=0.0> VE_severe_given_infection2_init;

  int<lower=0> variant2_introduction;

  int<lower=0> ninter;                      // number of interventions
  real<lower=1.0> t_inter[ninter];       // mean start time of each interventions
  real<lower=1.0> len_inter[ninter];     // mean length of each intervention
  real<lower=0.0> mu_beta_inter[ninter];    // mean change in beta through intervention
  real<lower=0.0> sigma_beta_inter[ninter]; // sd change in beta through intervention

  int fit_to_data;
}

transformed data {
  // general assumption - protection from variant 2 implies protection from variant 1

  //assigning indices for state matrix x
  int S = 1; //susceptible to both variants
  int E1 = 2; //exposed to variant 1
  int E2 = 3; //exposed to variant 2
  int Imild1 = 4; //infected with variant 1 - will not be admitted
  int Imild2 = 5; //infected with variant 2 - will not be admitted
  int Ipreh1 = 6; //infected with variant 1 - will be admitted
  int Ipreh2 = 7; //infected with variant 2 - will be admitted
  int Hmod1  = 8; //hospitalized with variant 1
  int Hmod2  = 9; //hospitalized with variant 2
  int P1 = 10; //protected from infection by variant 1 by vaccines or old infection (short-term)
  int P12 = 11; //protected from infection by variants 1 and 2 by vaccines or old infection (short-term)
  int Rlive1 = 12; //protected from infection by variant 1 by recent infection with variant 1 (long-term); some variant 1 infections protect only from 1
  int Rlive12 = 13; //protected from infection by variant 1 and 2 by recent infection with either variant (long-term) - all variant 2 infections protect from both; some variant 1 infections protect from both

  int ncompartments = 13;

  int obs_hosp_census = 1; //combined hospitalized
  int obs_cases1 = 2;
  int obs_cases2 = 3;
}
parameters {
  real<lower=1.0> duration_latent1; // duration is a minimum of 1 which is the stepsize of this model
  real<lower=1.0> duration_latent2; // duration is a minimum of 1 which is the stepsize of this model
  real<lower=1.0, upper = 999.0> duration_rec_mild; //was getting weird errors with duration_rec_mild=Inf
  real<lower=1.0> duration_pre_hosp;
  real<lower=1.0> duration_hosp_mod1;
  real<lower=1.0> duration_hosp_mod2;
  real<lower=1.0> duration_protection_infection;
  real<lower=0.001, upper = 1.0> frac_hosp1_naive;

  real<lower=0.001> frac_hosp_multiplier;

  real<lower=0.001, upper=1.0> frac_tested;

  real<lower=0.001, upper=0.1*npop> initial_infected1;
  real<lower=0.001, upper=0.1*npop> initial_infected2;
  real<lower=0.001> sigma_obs[nobs_types];

  real<lower=1.0> test_delay; //days tested after exposed

  real<lower=0.001, upper=9.9> beta_0;
  real<lower=0.001, upper=3.0> beta_multiplier[ninter];
  real<lower=0.001, upper=3.09> trans_multiplier;
}
transformed parameters {
  matrix<lower=0.0>[ncompartments,nt] x;
  matrix<lower=0.0>[nobs_types,nt] sim_data;
  row_vector<lower=0.0>[nt] new_cases1;
  row_vector<lower=0.0>[nt] new_cases2;
  row_vector<lower=0.0>[nt] soon_positive1;
  row_vector<lower=0.0>[nt] soon_positive2;
  real<lower=0.0> beta[nt];
  real Rt1[nt];
  real Rt2[nt];
  // real<lower=0.0> frac_hosp_0;
  real<lower=0.0, upper=1.0> frac_hosp2_naive;

  for (it in 1:nt) {
    for (itype in 1:nobs_types) {
      sim_data[itype, it] = 0.0;
    }
  }

  {
    // variables in curly brackets will not have output, they are local variables
    real newE1;
    real newE2;
    real newI1;
    real newI2;
    real frac_init_E;
    real avg_duration;
    real frac_boosters_to_susceptible;
    real new_protected1;
    real new_protected12;
    real increased_severity_protection;
    real frac_increased_severity_protection;
    real new_admits1;
    real new_admits2;
    real VE_severe_given_infection1;
    real VE_severe_given_infection2;
    real lost_protection_infection1;
    real lost_protection_infection12;
    real frac_hosp_init;
    real S2;
    real frac_E2_from_S;
    real frac_E2_from_P1;
    real frac_E2_from_Rlive1;
    real frac_hosp1;
    real frac_hosp2;

    frac_hosp2_naive = frac_hosp1_naive * frac_hosp_multiplier;

    // VE_severe_given_infection = VE_severe_given_infection_1_init;
    VE_severe_given_infection1 = VE_severe_given_infection1_init;
    VE_severe_given_infection2 = VE_severe_given_infection2_init;

    frac_hosp_init = frac_hosp1_naive * (1 - VE_severe_given_infection1);
    avg_duration = frac_hosp_init * duration_pre_hosp + (1 - frac_hosp_init) * duration_rec_mild;

    x[:,1] = rep_vector(0.0, ncompartments);
    frac_init_E = duration_latent1 / (duration_latent1 + avg_duration);
    x[E1, 1] = frac_init_E * initial_infected1;
    if (is_nan(x[E1, 1])) {
      print("x[E1, 1] is nan: ", frac_hosp1_naive, " ", VE_severe_given_infection1, " ", frac_hosp_init, " ", duration_pre_hosp, " ", duration_rec_mild, " ", avg_duration, " ", duration_latent1, " ", frac_init_E)
    }
    x[Imild1, 1] = (1 - frac_init_E) * initial_infected1 * (1 - frac_hosp_init);
    if (x[Imild1, 1] < 0) {
      print("x[Imild1, 1]=",x[Imild1, 1]," frac_init_E=",frac_init_E," initial_infected1=", initial_infected1, " frac_hosp_init=", frac_hosp_init, " duration_latent1=",duration_latent1, " duration_rec_mild=",duration_rec_mild)
    }
    x[Ipreh1, 1] = (1 - frac_init_E) * initial_infected1 * frac_hosp_init;
    x[Hmod1, 1] = init_hosp1;
    x[P1, 1] = npop * (VE_infection1_init - VE_infection2_init);
    x[P12, 1] = npop * VE_infection2_init;

    x[S, 1] = npop - initial_infected1 - x[P1, 1] - x[P12, 1] - x[Hmod1, 1];

if (x[S, 1] < -10) {
  print("x[S, 1] < -10 x[S, 1]=", x[S, 1], " initial_infected1=", initial_infected1, " x= ", x[:, 1])
}

// Rt1[1] = 1; //to avoid nan error
// Rt2[1] = 1; //to avoid nan error

    soon_positive1[1] = 0.0;
    new_cases1[1] = 0.0; //not correct - pass NA obs cases at t = 1 so no fitting to this
    if (tobs[obs_cases1, 1] == 1) {
      reject("Minimum tobs[obs_cases, :] is 2")
    }
    soon_positive2[1] = 0.0;
    new_cases2[1] = 0.0; //not correct - pass NA obs cases at t = 1 so no fitting to this
    if (tobs[obs_cases2, 1] == 1) {
      reject("Minimum tobs[obs_cases, :] is 2")
    }

    for (it in 1:nt) {
      beta[it] = beta_0;
      for (iinter in 1:ninter) {
        //k <- 2/s * qlogis(0.99) # = -2/s * qlogis(0.01) --> 2 * qlogis(0.99) = 9.19024
        //f <- m ^ plogis(k * (t - (d + s/2)))
        beta[it] = beta[it] * beta_multiplier[iinter] ^ inv_logit(9.19024 / len_inter[iinter] * (it - (t_inter[iinter] + len_inter[iinter] / 2)));
      }
    }

    for (it in 1:nt-1){
      newE1 = fmin(x[S, it], beta[it] * x[S, it] * (x[Imild1, it] + x[Ipreh1, it]) / npop);
      S2 = x[S, it] + x[P1, it] + x[Rlive1, it];

      frac_E2_from_S = x[S, it] / S2;
      frac_E2_from_P1 = x[P1, it] / S2;
      frac_E2_from_Rlive1 = 1 - (frac_E2_from_S + frac_E2_from_P1);

      if (it == variant2_introduction) {
        newE2 = initial_infected2;
      } else {
        newE2 = fmin(S2, beta[it] * trans_multiplier * S2 * (x[Imild2, it] + x[Ipreh2, it]) / npop);

      }
      Rt1[it] = beta[it] * x[S, it] / npop * duration_rec_mild;
      Rt2[it] = beta[it] * trans_multiplier * S2 / npop * duration_rec_mild;

      newI1 = x[E1, it] / duration_latent1;
      newI2 = fmin(x[E2, it], x[E2, it] / duration_latent2);

      //this all needs work
      // frac_boosters_to_susceptible = x[S, it] / (x[S, it] + omicron_recovered_booster_scale * x[Rlive, it]);
      // new_protected = fmin(x[S, it] - newE, num_boosters[it] * frac_boosters_to_susceptible * booster_VE_infection);
      // if (new_protected < 0) reject("new_protected < 0")
      //
      //
      // increased_severity_protection = num_boosters[it] * frac_boosters_to_susceptible * (1 - booster_VE_infection); //boosted but not Protected from infection, still susceptible
      // frac_increased_severity_protection = increased_severity_protection / fmax(x[S, it], 0.001); //avoid 0/0 problems
      //
      // //increase VE_severe_given_infection for boosters
      // VE_severe_given_infection = booster_VE_severe * frac_increased_severity_protection +   VE_severe_given_infection * (1 - frac_increased_severity_protection);

      frac_hosp1 = frac_hosp1_naive * (1 - VE_severe_given_infection1);
      frac_hosp2 = frac_hosp2_naive * (1 - VE_severe_given_infection2);

      //new_protected1 was subtracted from S but this seems more complicated now
      new_protected1 = 0;
      new_protected12 = 0;

      new_admits1 = x[Ipreh1, it] / duration_pre_hosp;
      new_admits2 = x[Ipreh2, it] / duration_pre_hosp;

      lost_protection_infection1 = x[P1, it] / duration_protection_infection;
      lost_protection_infection12 = x[P12, it] / duration_protection_infection;

      x[S, it + 1] = x[S, it] - newE1 + lost_protection_infection1 - newE2 * frac_E2_from_S;
      x[E1, it + 1] = x[E1, it] + newE1 - newI1;
      x[E2, it + 1] = x[E2, it] + newE2 - newI2;

      // if (it == variant2_introduction || it == (variant2_introduction+1)) {
        //   print("it = ", it, " newE2 = ", newE2, " newI2 = ", newI2, " x[E2, it] = ", x[E2, it], " x[E2, it + 1] = ", x[E2, it + 1], " duration_latent2 = ", duration_latent2)
        // }

        x[Imild1, it + 1] = x[Imild1, it] + newI1 * (1 - frac_hosp1) - x[Imild1, it] / duration_rec_mild;
        x[Imild2, it + 1] = x[Imild2, it] + newI2 * (1 - frac_hosp2) - x[Imild2, it] / duration_rec_mild;

        x[Ipreh1, it + 1] = x[Ipreh1, it] + newI1 * frac_hosp1 - new_admits1;
        x[Ipreh2, it + 1] = x[Ipreh2, it] + newI2 * frac_hosp2 - new_admits2;

        x[Hmod1, it + 1] = x[Hmod1, it] + new_admits1 - x[Hmod1, it] / duration_hosp_mod1;
        x[Hmod2, it + 1] = x[Hmod2, it] + new_admits2 - x[Hmod2, it] / duration_hosp_mod2;

        x[P1, it + 1] = x[P1, it] + new_protected1 - lost_protection_infection1 +  lost_protection_infection12 - newE2 * frac_E2_from_P1;
        x[P12, it + 1] = x[P12, it] + new_protected12 - lost_protection_infection12;

        x[Rlive1, it + 1] = x[Rlive1, it] + x[Hmod1, it] / duration_hosp_mod1 + x[Imild1, it] / duration_rec_mild - newE2 * frac_E2_from_Rlive1;
        x[Rlive12, it + 1] = x[Rlive12, it] + x[Hmod2, it] / duration_hosp_mod2 + x[Imild2, it] / duration_rec_mild;



        soon_positive1[it + 1] = soon_positive1[it] + newE1 * frac_tested - soon_positive1[it] / test_delay;
        soon_positive2[it + 1] = soon_positive2[it] + newE2 * frac_tested - soon_positive2[it] / test_delay;
        new_cases1[it + 1] = soon_positive1[it + 1] / test_delay;
        new_cases2[it + 1] = soon_positive2[it + 1] / test_delay;

        // if (is_nan(new_cases1[it + 1])) {
          //   print("is nan cases:")
          //   print(it, " ", soon_positive[it], " ", soon_positive[it + 1], " ", newE, " ", frac_tested, " ", test_delay)
          //   reject("new_cases[it+1] is nan")
          // }

          // if (is_nan(sum(x[:, it + 1]))) {
            //   print("is nan:")
            //   print(newE, " ", newI, " ", frac_init_E, " ", beta[it], " ", beta_0, " ", avg_duration, " ", frac_hosp_0, " ", frac_boosters_to_susceptible, " ", new_protected, " ", increased_severity_protection, " ", frac_increased_severity_protection, " ", frac_hosp, " ", new_admits, " ", VE_severe_given_infection, " ");
            //   reject("x is nan: it = ", it, " x[, it+1] = ", x[:, it + 1]);
            // }
            for (ii in 1:ncompartments) {
              if (x[ii, it + 1] < -100) {
                print("beta[it] * trans_multiplier * S2 * (x[Imild2, it] + x[Ipreh2, it]) / npop=", beta[it] * trans_multiplier * S2 * (x[Imild2, it] + x[Ipreh2, it]) / npop)
                print("beta[it]=", beta[it], " trans_multiplier=", trans_multiplier, " S2=", S2, " x[Imild2, it] + x[Ipreh2, it]=", x[Imild2, it] + x[Ipreh2, it], " npop=", npop)

              reject("x[ii, it + 1] < -100: it= ", it, " ii=", ii, " variant2_introduction=",variant2_introduction,  " newE2=", newE2, " x[Imild2, it]=", x[Imild2, it]," x[Ipreh2, it]=", x[Ipreh2, it], " x[:, it]=", x[:, it], " ==== x[:, it+1]=  ", x[:, it+1]);
            }
            }

            // test
            if (fabs(sum(x[:, it + 1])-npop) > 0.01) {
              // reject("Model is leaking, net gain: ", sum(x[:,it+1])-npop);
              reject("Model is leaking, net gain: ", sum(x[:,it+1])-npop, "it= ", it, " ", x[:, it], " ====  ", x[:, it+1]);
            }
    }

    // Data for fitting
    sim_data[obs_hosp_census] = x[Hmod1] * (1 + frac_incidental1) + x[Hmod2] * (1 + frac_incidental2);
    sim_data[obs_cases1] = new_cases1;
    sim_data[obs_cases2] = new_cases2;
  }

  // print(x)
  // print(sim_data)
}

model {
  //////////////////////////////////////////
  // prior distributions
  //
  duration_latent1 ~ normal(mu_duration_latent1, sigma_duration_latent1);
  duration_latent2 ~ normal(mu_duration_latent2, sigma_duration_latent2);
  duration_rec_mild ~ normal(mu_duration_rec_mild, sigma_duration_rec_mild);
  duration_pre_hosp ~ normal(mu_duration_pre_hosp, sigma_duration_pre_hosp);
  duration_hosp_mod1 ~ normal(mu_duration_hosp_mod1, sigma_duration_hosp_mod1);
  duration_hosp_mod2 ~ normal(mu_duration_hosp_mod2, sigma_duration_hosp_mod2);
  duration_protection_infection ~ normal(mu_duration_protection_infection, sigma_duration_protection_infection);

  beta_0 ~ normal(mu_beta0, sigma_beta0);
  for (iinter in 1:ninter) {
    beta_multiplier[iinter] ~ normal(mu_beta_inter[iinter], sigma_beta_inter[iinter]);
  }

  frac_tested ~ normal(mu_frac_tested, sigma_frac_tested);
  frac_hosp1_naive ~ normal(mu_frac_hosp1_naive, sigma_frac_hosp1_naive);
  frac_hosp_multiplier ~ normal(mu_frac_hosp_multiplier, sigma_frac_hosp_multiplier);
  trans_multiplier ~ normal(mu_trans_multiplier, sigma_trans_multiplier);
  test_delay ~ normal(mu_test_delay, sigma_test_delay);

  // initial_infected1 ~ exponential(lambda_initial_infected1);
  initial_infected1 ~ normal(1/lambda_initial_infected1, 0.1 * 1/lambda_initial_infected1);
  initial_infected2 ~ normal(1/lambda_initial_infected2, 0.1 * 1/lambda_initial_infected2);

  //////////////////////////////////////////
  // fitting observations

  if (fit_to_data) {
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
  matrix[nobs_types,nt] sim_data_with_error;
  // real<lower=0.0> Rt[nt];

  if (fit_to_data) {
    for (itype in 1:nobs_types) {
      if (nobs[itype] > 0) {
        for (it in 1:nt) {
          sim_data_with_error[itype, it] = fmax(0.0, normal_rng(sim_data[itype, it], sigma_obs[itype])); //can't vectorize - fmax not vectorized?
        }
      } else {
        sim_data_with_error[itype] = sim_data[itype];
      }
    }
  } else {
    for (itype in 1:nobs_types) {
      sim_data_with_error[itype] = sim_data[itype];
    }
  }
  {
    real frac_prehosp;
    real avg_duration;
    // for (it in 1:nt) {
      //   frac_prehosp = (1e-10 * frac_hosp_0 + x[Ipreh, it]) / (1e-10 + x[Ipreh, it] +  x[Imild, it]); //at t=1 converges to frac_hosp_0
      //   avg_duration = frac_prehosp * duration_pre_hosp + (1 - frac_prehosp) * duration_rec_mild;
      //   Rt[it] = beta[it] * avg_duration * x[S, it] / npop;
      // }
  }
}
