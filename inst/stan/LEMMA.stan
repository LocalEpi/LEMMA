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
  real<lower=0.0> mu_Rt1_init;
  real<lower=0.0> sigma_Rt1_init;
  real<lower=0.0> mu_frac_hosp1_naive;
  real<lower=0.0> sigma_frac_hosp1_naive;
  real<lower=0.0> mu_trans_multiplier;
  real<lower=0.0> sigma_trans_multiplier;
  real<lower=0.0> mu_frac_hosp_multiplier;
  real<lower=0.0> sigma_frac_hosp_multiplier;

  real<lower=0.0> mu_duration_protection_infection;
  real<lower=0.0> sigma_duration_protection_infection;
  real<lower=0.0> mu_VE_infection1;
  real<lower=0.0> sigma_VE_infection1;
  real<lower=0.0> mu_immune_evasion;
  real<lower=0.0> sigma_immune_evasion;

  // real <lower=0.0> mu_initial_exposed1;
  // real <lower=0.0> sigma_initial_exposed1;
  real <lower=0.0> mu_initial_exposed2_frac;
  real <lower=0.0> sigma_initial_exposed2_frac;
  // real <lower=0.0> mu_initial_infected1;
  // real <lower=0.0> sigma_initial_infected1;
  // real <lower=0.0> mu_initial_preh1;
  // real <lower=0.0> sigma_initial_preh1;

  real<lower=0.0> sigma_obs_est_inv[nobs_types];

  real<lower=0.0> frac_case2_growth_obs;
  real<lower=0.0> sigma_frac_case2_growth_obs;
  real<lower=0.0> variant2_crossover_days_obs;
  real<lower=0.0> sigma_variant2_crossover_days;

  real<lower=0.0> init_hosp1;

  real<lower=0.0> num_boosters[nt];
  real<lower=0.0> booster_VE_infection1;
  real<lower=0.0> booster_VE_infection2;
  real<lower=0.0> booster_VE_severe_given_infection1;
  real<lower=0.0> booster_VE_severe_given_infection2;
  real<lower=0.0> frac_incidental1;
  real<lower=0.0> frac_incidental2;
  real<lower=0.0> VE_severe_given_infection1_init;
  real<lower=0.0> VE_severe_given_infection2_init;

  int<lower=0> variant2_introduction;

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
  int P1 = 10; //protected from infection by variant 1 by vaccines or infection
  int P12 = 11; //protected from infection by variants 1 and 2 by vaccines or infection
  //some variant 1 infections protect only from 1, some variant 1 infections protect only from both

  int ncompartments = 11;

  int obs_hosp_census = 1; //combined hospitalized
  int obs_cases = 2; //combined cases

  int days_delay_measure_growth = 5; //number of days after variant2_introduction to start calculating frac2 growth
  int nfrac_growth = nt - (variant2_introduction + days_delay_measure_growth) + 1;
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

  // real<lower=0.001, upper=0.1*npop> initial_exposed1;
  real<lower=1e-6, upper=1> initial_exposed2_frac;
  // real<lower=0.001, upper=0.1*npop> initial_infected1;
  // real<lower=0.001, upper=0.1*npop> initial_preh1;
  real<lower=0.001> sigma_obs[nobs_types];

  real<lower=1.0> test_delay; //days tested after exposed

  real<lower=0.001, upper=9.9> Rt1_init;
  real<lower=0.001, upper=3.0> trans_multiplier;
  real<lower=0.001, upper=1.0> VE_infection1;
  real<lower=0.001, upper=1.0> immune_evasion;
}
transformed parameters {
  matrix<lower=0.0>[ncompartments,nt] x;
  matrix<lower=0.0>[nobs_types,nt] sim_data;
  row_vector<lower=0.0>[nt] new_cases1;
  row_vector<lower=0.0>[nt] new_cases2;
  row_vector<lower=0.0>[nt] soon_positive1;
  row_vector<lower=0.0>[nt] soon_positive2;

  real<lower=0.0, upper=1.0> frac_hosp2_naive;
  real frac_case2_growth;
  real variant2_crossover_days;
  vector<lower=0.0, upper=1.0>[nfrac_growth] frac_case2;
  vector[2] frac_case2_growth_coef;
  real<lower=0.0> beta_0;

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
    real increased_severity_protection;
    real frac_increased_severity_protection;
    real new_admits1;
    real new_admits2;
    real VE_severe_given_infection1;
    real VE_severe_given_infection2;
    real lost_protection_infection1;
    real lost_protection_infection12;
    real frac_hosp_init;
    real frac_preh_init;
    real S2;
    real frac_E2_from_S;
    real frac_E2_from_P1;
    real frac_hosp1;
    real frac_hosp2;
    real frac_boosters_to_S;
    real frac_boosters_to_P1;
    real new_P12_from_S;
    real new_P1_from_S;
    real new_P12_from_P1;
    real S_increased_severe_protection;
    real P1_increased_severe_protection;
    real frac_increased_severity_protection1;
    real frac_increased_severity_protection2;
    int max_compartment;
    int index;
    real recovered1;
    matrix[nfrac_growth, 2] X;

    frac_hosp2_naive = frac_hosp1_naive * frac_hosp_multiplier;
    if (frac_hosp2_naive > 1) {
      reject("frac_hosp2_naive > 1")
    }

    VE_severe_given_infection1 = VE_severe_given_infection1_init;
    VE_severe_given_infection2 = VE_severe_given_infection2_init;

    frac_hosp_init = frac_hosp1_naive * (1 - VE_severe_given_infection1);
    frac_preh_init = frac_hosp_init * duration_pre_hosp / (duration_pre_hosp + duration_rec_mild);
    avg_duration = frac_hosp_init * duration_pre_hosp + (1 - frac_hosp_init) * duration_rec_mild;

    x[:,1] = rep_vector(0.0, ncompartments);
    frac_init_E = duration_latent1 / (duration_latent1 + avg_duration);
    x[E1, 1] = duration_latent1 / duration_hosp_mod1 * init_hosp1 / frac_hosp_init; //frac_init_E * initial_infected1;

    x[Ipreh1, 1] = duration_pre_hosp / duration_hosp_mod1 * init_hosp1; //(1 - frac_init_E) * initial_infected1 * frac_preh_init;
    x[Imild1, 1] = duration_rec_mild / (Rt1_init * duration_latent1) * x[E1, 1] - x[Ipreh1, 1]; //initial_infected1; //(1 - frac_init_E) * initial_infected1 * (1 - frac_preh_init);
    x[Hmod1, 1] = init_hosp1;
    x[P1, 1] =  npop * VE_infection1 * immune_evasion; //npop * (VE_infection1 - VE_infection2);
    x[P12, 1] = npop * VE_infection1 * (1 - immune_evasion); //npop * VE_infection2;

    x[S, 1] = npop - (x[E1, 1] + x[Imild1, 1] + x[Ipreh1, 1] + x[Hmod1, 1] + x[P1, 1] + x[P12, 1]);

    //Rt1[it] = beta_0 * x[S, it] / npop * duration_rec_mild;
    beta_0 = Rt1_init / (x[S, 1] / npop * duration_rec_mild);
    if (beta_0 < 0) {
      reject("beta_0 < 0") //beta_0 has lower=0.0 but this avoids going through the 1:nt-1 loop first
    }

    if (x[S, 1] < -10) {
      // print("x[S, 1] < -10 x[S, 1]=", x[S, 1], " initial_infected1=", initial_infected1, " x= ", x[:, 1])
    }

    soon_positive1[1] = 0.0;
    new_cases1[1] = 0.0; //not correct - pass NA obs cases at t = 1 so no fitting to this
    if (tobs[obs_cases, 1] == 1) {
      reject("Minimum tobs[obs_cases, :] is 2")
    }
    soon_positive2[1] = 0.0;
    new_cases2[1] = 0.0; //not correct - pass NA obs cases at t = 1 so no fitting to this



    if (booster_VE_infection1 - booster_VE_infection2 < 0) {
      reject("booster_VE_infection1 - booster_VE_infection2 < 0")
    }

    for (it in 1:nt-1){
      newE1 = fmin(x[S, it], beta_0 * x[S, it] * (x[Imild1, it] + x[Ipreh1, it]) / npop);
      S2 = x[S, it] + x[P1, it];

      newI1 = fmin(x[E1, it], x[E1, it] / duration_latent1); //fmin shouldn't be needed since duration_latent1 >= 1 but can't hurt
      if (it == variant2_introduction) {
        newE2 = (x[E1, it] + newE1 - newI1) * initial_exposed2_frac;
      } else {
        newE2 = fmin(S2, beta_0 * trans_multiplier * S2 * (x[Imild2, it] + x[Ipreh2, it]) / npop);
      }
      newE2 = fmax(newE2, 0.0); //fmax shouldn't be needed but can't hurt
      newI2 = fmin(x[E2, it], x[E2, it] / duration_latent2);

      //boosters:
      //S -> S or P1 or P12
      //P1 -> P1 or P12
      //P12 -> P12

      //needs double checking and code for 0/0, making sure new_X_from_Y is not more than Y (accounting for newE)
      frac_boosters_to_S = x[S, it] / (x[S, it] + x[P1, it] + x[P12, it]);
      frac_boosters_to_P1 = x[P1, it] / (x[S, it] + x[P1, it] + x[P12, it]);

      new_P12_from_S = num_boosters[it] * frac_boosters_to_S * booster_VE_infection2;
      new_P1_from_S = num_boosters[it] * frac_boosters_to_S * (booster_VE_infection1 - booster_VE_infection2);
      S_increased_severe_protection = num_boosters[it] * frac_boosters_to_S * (1 - booster_VE_infection1); //boosted but still susceptible to both, increased protection from severe disease
      new_P12_from_P1 = num_boosters[it] * frac_boosters_to_P1 * booster_VE_infection2;
      P1_increased_severe_protection = num_boosters[it] * frac_boosters_to_P1 * (1 - booster_VE_infection2) + new_P1_from_S; //boosted but still susceptible to variant2, increased protection from severe disease

      frac_increased_severity_protection1 = S_increased_severe_protection / (x[S, it] + 1e-8); //avoid divide by 0
      frac_increased_severity_protection2 = (S_increased_severe_protection + P1_increased_severe_protection) / S2;

      VE_severe_given_infection1 = booster_VE_severe_given_infection1 * frac_increased_severity_protection1 +   VE_severe_given_infection1 * (1 - frac_increased_severity_protection1);
      VE_severe_given_infection2 = booster_VE_severe_given_infection2 * frac_increased_severity_protection2 +   VE_severe_given_infection2 * (1 - frac_increased_severity_protection2);

      frac_hosp1 = frac_hosp1_naive * (1 - VE_severe_given_infection1);
      frac_hosp2 = frac_hosp2_naive * (1 - VE_severe_given_infection2);

      new_admits1 = x[Ipreh1, it] / duration_pre_hosp;
      new_admits2 = x[Ipreh2, it] / duration_pre_hosp;

      lost_protection_infection1 = x[P1, it] / duration_protection_infection;
      lost_protection_infection12 = x[P12, it] / duration_protection_infection;

      frac_E2_from_S = x[S, it] / S2;
      frac_E2_from_P1 = 1 - frac_E2_from_S;

      x[S, it + 1] = x[S, it] - newE1 + lost_protection_infection1 - newE2 * frac_E2_from_S - new_P12_from_S - new_P1_from_S;

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

        recovered1 = x[Hmod1, it] / duration_hosp_mod1 + x[Imild1, it] / duration_rec_mild;

        x[P1, it + 1] = x[P1, it] - lost_protection_infection1 +  lost_protection_infection12 - newE2 * frac_E2_from_P1 + new_P1_from_S - new_P12_from_P1 + recovered1 * immune_evasion;
        x[P12, it + 1] = x[P12, it] - lost_protection_infection12 + new_P12_from_S + new_P12_from_P1 +  x[Hmod2, it] / duration_hosp_mod2 + x[Imild2, it] / duration_rec_mild + recovered1 * (1 - immune_evasion);

        if (is_nan(x[P12, it + 1])) {
          print(x[P12, it], " ", lost_protection_infection12, " ", new_P12_from_S," ",  new_P12_from_P1, " ",  x[Hmod2, it], " ", duration_hosp_mod2," ",  x[Imild2, it], " ", duration_rec_mild, " ", recovered1, " ", immune_evasion)
        }

        soon_positive1[it + 1] = soon_positive1[it] + newE1 * frac_tested - soon_positive1[it] / test_delay;
        soon_positive2[it + 1] = soon_positive2[it] + newE2 * frac_tested - soon_positive2[it] / test_delay;
        new_cases1[it + 1] = soon_positive1[it + 1] / test_delay;
        new_cases2[it + 1] = soon_positive2[it + 1] / test_delay;

        // if (is_nan(new_cases1[it + 1])) {
          //   print("is nan cases:")
          //   print(it, " ", soon_positive[it], " ", soon_positive[it + 1], " ", newE, " ", frac_tested, " ", test_delay)
          //   reject("new_cases[it+1] is nan")
          // }

          if (is_nan(sum(x[:, it + 1]))) {
            print("is nan:")
            // print(newE, " ", newI, " ", frac_init_E, " ", beta_0, " ", avg_duration, " ", frac_hosp_0, " ", frac_boosters_to_susceptible, " ", new_protected, " ", increased_severity_protection, " ", frac_increased_severity_protection, " ", frac_hosp, " ", new_admits, " ", VE_severe_given_infection, " ");
            reject("x is nan: it = ", it, "x[, it] = ", x[:, it], " x[, it+1] = ", x[:, it + 1]);
          }
          if (min(x[:, it]) < 0) {
            for (ii in 1:ncompartments) {
              if (x[ii, it + 1] < 0) {
                max_compartment = 1;
                for (ii2 in 1:ncompartments) {
                  if (x[ii2, it + 1] > x[max_compartment, it + 1]) {
                    max_compartment = ii2;
                  }
                }
                if (x[ii, it + 1] < -100) {
                  print("negative compartment: x[", ii, ", ", it+1, "] = ", x[ii, it + 1], " moving to compartment ", max_compartment, "x[max_compartment, it + 1] = ", x[max_compartment, it + 1])
                  if (ii == 2) {
                    print("duration_latent1 = ", duration_latent1, " x[E1, it+1] = ", x[E1, it+1], "x[E1, it] = ", x[E1, it], " newE1 = ", newE1, " newI1 = ", newI1)
                    print("beta_0 = ", beta_0)
                    print(x[, it])

                  }
                  // print("negative compartment: x[", ii, ", ", it+1, "] = ", x[ii, it + 1], " moving to compartment ", max_compartment, "x[max_compartment, it + 1] = ", x[max_compartment, it + 1], " x[:, it]=", x[:, it], " ==== x[:, it+1]=  ", x[:, it+1])
                }
                x[max_compartment, it + 1] = x[max_compartment, it + 1] + x[ii, it + 1];
                x[ii, it + 1] = 0;
              }
            }
          }
          // for (ii in 1:ncompartments) {
            //   if (x[ii, it + 1] < -100) {
              //     print("beta[it] * trans_multiplier * S2 * (x[Imild2, it] + x[Ipreh2, it]) / npop=", beta[it] * trans_multiplier * S2 * (x[Imild2, it] + x[Ipreh2, it]) / npop)
              //     print("beta[it]=", beta[it], " trans_multiplier=", trans_multiplier, " S2=", S2, " x[Imild2, it] + x[Ipreh2, it]=", x[Imild2, it] + x[Ipreh2, it], " npop=", npop)
              //
              //     reject("x[ii, it + 1] < -100: it= ", it, " ii=", ii, " variant2_introduction=",variant2_introduction,  " newE2=", newE2, " x[Imild2, it]=", x[Imild2, it]," x[Ipreh2, it]=", x[Ipreh2, it], " x[:, it]=", x[:, it], " ==== x[:, it+1]=  ", x[:, it+1]);
              //   }
              // }

              // test
              if (fabs(sum(x[:, it + 1])-npop) > 0.01) {
                // reject("Model is leaking, net gain: ", sum(x[:,it+1])-npop);
                reject("Model is leaking, net gain: ", sum(x[:,it+1])-npop, "it= ", it, " ", x[:, it], " ====  ", x[:, it+1]);
              }
    }

    for (it in (variant2_introduction + days_delay_measure_growth):nt) {
      //soon after variant2_introduction+1, new_cases2 = 0 (because not enough time to have any infected, only exposed)
      index = it - (variant2_introduction + days_delay_measure_growth) + 1;
      frac_case2[index] = fmin(1-1e-6, fmax(1e-6, new_cases2[it] / (new_cases1[it] + new_cases2[it] + 1e-8)));
      X[index, 1] = 1.0;
      X[index, 2] = it;
    }

    //     nt <- new_variant_frac[, .N]
    // X <- cbind(rep(1, nt), 1:nt)
    // y <- qlogis(new_variant_frac$variant2)
    // # b <- solve(t(X) %*% X) %*% t(X) %*% y
    // b <- solve(t(X) %*% X, t(X) %*% y)


    frac_case2_growth_coef = mdivide_left_spd(X' * X, X' * logit(frac_case2));
    frac_case2_growth = frac_case2_growth_coef[2];
    variant2_crossover_days = -frac_case2_growth_coef[1] / frac_case2_growth_coef[2];

    // Data for fitting
    sim_data[obs_hosp_census] = x[Hmod1] * (1 + frac_incidental1) + x[Hmod2] * (1 + frac_incidental2);
    sim_data[obs_cases] = new_cases1 + new_cases2;
  }
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

  Rt1_init ~ normal(mu_Rt1_init, sigma_Rt1_init);

  frac_tested ~ normal(mu_frac_tested, sigma_frac_tested);
  frac_hosp1_naive ~ normal(mu_frac_hosp1_naive, sigma_frac_hosp1_naive);
  frac_hosp_multiplier ~ normal(mu_frac_hosp_multiplier, sigma_frac_hosp_multiplier);
  trans_multiplier ~ normal(mu_trans_multiplier, sigma_trans_multiplier);
  test_delay ~ normal(mu_test_delay, sigma_test_delay);
  VE_infection1 ~ normal(mu_VE_infection1, sigma_VE_infection1);
  immune_evasion ~ normal(mu_immune_evasion, sigma_immune_evasion);

  // initial_exposed1 ~ normal(mu_initial_exposed1, sigma_initial_exposed1);
  initial_exposed2_frac ~ normal(mu_initial_exposed2_frac, sigma_initial_exposed2_frac);
  // initial_infected1 ~ normal(mu_initial_infected1, sigma_initial_infected1);
  // initial_preh1 ~ normal(mu_initial_preh1, sigma_initial_preh1);

  //////////////////////////////////////////
  // fitting observations

  sigma_obs ~ exponential(sigma_obs_est_inv);
  if (fit_to_data == 1) {
    //deal with PUIs and NAs in R
    for (itype in 1:nobs_types) {
      if (nobs[itype] > 0) {
        obs_data[itype, 1:nobs[itype]] ~ normal(sim_data[itype, tobs[itype, 1:nobs[itype]]], sigma_obs[itype]);
      }
    }
    frac_case2_growth_obs ~ normal(frac_case2_growth, sigma_frac_case2_growth_obs);
    variant2_crossover_days_obs ~ normal(variant2_crossover_days, sigma_variant2_crossover_days);
  } else if (fit_to_data == 0) {
    frac_case2_growth_obs ~ normal(0, 1);
    variant2_crossover_days_obs ~ normal(0, 1);
  } else if (fit_to_data == 2) {
    for (itype in 1:nobs_types) {
      if (nobs[itype] > 0) {
        obs_data[itype, 1:nobs[itype]] ~ normal(sim_data[itype, tobs[itype, 1:nobs[itype]]], sigma_obs[itype]);
      }
      frac_case2_growth_obs ~ normal(0, 1);
      variant2_crossover_days_obs ~ normal(0, 1);
    }
  } else {
    reject("bad fit_to_data")
  }
}

generated quantities{
  matrix[nobs_types,nt] sim_data_with_error;
  real<lower=0.0> Rt1[nt];
  real<lower=0.0> Rt2[nt];

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
    real S2;
    for (it in 1:nt) {
      S2 = x[S, it] + x[P1, it];
      Rt1[it] = beta_0 * x[S, it] / npop * duration_rec_mild;
      Rt2[it] = beta_0 * trans_multiplier * S2 / npop * duration_rec_mild;
    }
  }
}
