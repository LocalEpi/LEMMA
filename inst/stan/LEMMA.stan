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
  int<lower=0> npops;
  real<lower=-1.0> obs_data_conf[nobs_types, nobs, npops];  // observed confirmed (-1 = NA)
  real<lower=-1.0> obs_data_pui[nobs_types, nobs, npops];   // observed PUI (-1 = NA)
  int<lower=0> nt;                       // number of time steps
  real<lower=0> population[npops];                             // total population
  int<lower=0, upper=1> extend;

  //////////////////////////////////////////
  // prior parameter distributions
  real<lower=1.0> mu_duration_latent;     // mean duration in "exposed" stage
  real<lower=0.0> sigma_duration_latent;  // sd duration in "exposed" stage
  real<lower=1.0> mu_duration_rec_mild;   // mean duration in "infectious" stage for mild cases
  real<lower=0.0> sigma_duration_rec_mild;// sd duration in "infectious" stage for mild cases
  real<lower=1.0> mu_duration_pre_hosp[npops];   // mean duration in "infectious" stage for hospitalized cases
  real<lower=0.0> sigma_duration_pre_hosp[npops];// sd duration in "infectious" stage for hospitalized cases
  real<lower=1.0> mu_duration_hosp_mod[npops];   // mean duration in hospital for non-ICU cases
  real<lower=0.0> sigma_duration_hosp_mod[npops];// sd duration in hospital for non-ICU cases
  real<lower=1.0> mu_duration_hosp_icu[npops];   // mean duration in hospital for ICU cases
  real<lower=0.0> sigma_duration_hosp_icu[npops];// sd duration in hospital for ICU cases

  real<lower=0.0> mu_r0[npops];                  // mean initial beta estimate
  real<lower=0.0> sigma_r0[npops];               // sd initial beta estimate

  real<lower=0.0> mu_frac_pui[nobs_types];     // mean fraction of PUI that are true COVID+
  real<lower=0.0> sigma_frac_pui[nobs_types];  // sd fraction of PUI that are true COVID+

  real<lower=0.0> mu_frac_hosp[npops];           // mean ICU + non-ICU
  real<lower=0.0> sigma_frac_hosp[npops];        // sd ICU + non-ICU
  real<lower=0.0> mu_frac_icu[npops];            // mean ICU as fraction of hosp
  real<lower=0.0> sigma_frac_icu[npops];         // sd ICU as fraction of hosp
  real<lower=0.0> mu_frac_mort[npops];           // mean mortality as fraction of ICU
  real<lower=0.0> sigma_frac_mort[npops];        // sd mortality as fraction of ICU

  real<lower=0.0> lambda_ini_exposed[npops];     // parameter for initial conditions of "exposed"

  //////////////////////////////////////////
  // interventions

  int<lower=0> ninter;                      // number of interventions
  real<lower=1.0> mu_t_inter[ninter];       // mean start time of each interventions
  real<lower=0.0> sigma_t_inter[ninter];    // sd start time of each interventions
  real<lower=1.0> mu_len_inter[ninter];     // mean length of each intervention
  real<lower=1.0> sigma_len_inter[ninter];  // sd length of each intervention
  real<lower=0.0> mu_beta_inter[ninter, npops];    // mean change in beta through intervention
  real<lower=0.0> sigma_beta_inter[ninter, npops]; // sd change in beta through intervention

  real<lower=0.0> len_inter_age;
  real<lower=0.0> t_inter_age;
  real<lower=0.0> mu_frac_hosp_multiplier;
  real<lower=0.0> sigma_frac_hosp_multiplier;
  real<lower=0.0> mu_frac_icu_multiplier;
  real<lower=0.0> sigma_frac_icu_multiplier;
  real<lower=0.0> mu_frac_mort_multiplier;
  real<lower=0.0> sigma_frac_mort_multiplier;
}
transformed data {
  //assigning indices for state matrix x
  int S = 1;
  int E = 2;
  int Imild = 3;
  int Ipreh = 4;
  int Hmod  = 5;
  int Hicu  = 6;
  int Rlive = 7;
  int Rmort = 8;

  int ncompartments = 8;

  int obs_hosp_census = 1;
  int obs_icu_census = 2;
  int obs_cum_deaths = 3;
  int obs_cum_admits = 4;

  int nobs_notmissing = 0;
  real beta_limit;

  for (ipop in 1:npops) {
    for (iobs in 1:nobs){
      for (itype in 1:nobs_types) {
        if (obs_data_conf[itype, iobs, ipop] > 0) {
          nobs_notmissing = nobs_notmissing + 1;
        }
      }
    }
  }

  if (extend == 1) {
    beta_limit = 9999; // no limit on beta if extending simulation
  } else {
    beta_limit = 2.0;
  }
}
parameters {

  real<lower=1.0> duration_latent; // duration is a minimum of 1 which is the stepsize of this model
  real<lower=1.0> duration_rec_mild;

  real<lower=1.0> duration_pre_hosp[npops];
  real<lower=1.0> duration_hosp_mod[npops];
  real<lower=1.0> duration_hosp_icu[npops];

  real<lower=0.005, upper=1.0> frac_hosp_0[npops];
  real<lower=0.0, upper=1.0> frac_icu_0[npops];
  real<lower=0.0, upper=1.0> frac_mort_0[npops];


  real<lower=0> ini_exposed[npops];

  real<lower=0> sigma_obs[nobs_types,npops];

  real<lower=0.0> r0[npops];
  real<lower=0.0> beta_multiplier[ninter,npops];
  real<lower=1.0> t_inter[ninter];
  real<lower=1.0> len_inter[ninter];

  real<lower=0.0> frac_hosp_multiplier;
  real<lower=0.0> frac_icu_multiplier;
  real<lower=0.0> frac_mort_multiplier;

  real<lower=0, upper=1> frac_PUI[nobs_types];

  real<lower=0, upper=1> m[npops];
  simplex[4] s1;
  simplex[4] s2;

}
transformed parameters {

  real<lower=0.0> x[ncompartments,nt,npops];
  real<lower=0.0> sim_data[nobs_types,nt,npops];
  real<lower=0.0, upper=2.0> beta[nt,npops];
  real<lower=0.0> beta_mat[nt-1,npops,npops]; //could make local, drop nt
  real<lower=0.0> Hadmits[nt,npops];
  real<lower=0.0, upper=1.0> frac_hosp[nt, npops];
  real<lower=0.0, upper=1.0> frac_icu[nt, npops];
  real<lower=0.0, upper=1.0> frac_mort[nt, npops];

  {
    // variables in curly brackets will not have output, they are local variables

    real newE[npops];
    real newI[npops];
    real newrec_mild[npops];
    real newrec_mod[npops];
    real newhosp[npops];
    real leave_icu[npops];
    real beta_0[npops];
    real zero;
    real N1 = population[1];
    real N2 = population[2];
    // real beta_mat[npops,npops];

    //////////////////////////////////////////
    for (it in 1:nt) {

    }
    for (ipop in 1:npops) {
      for (it in 1:nt) {
        frac_hosp[it, ipop] = frac_hosp_0[ipop] * frac_hosp_multiplier ^ inv_logit(9.19024 / len_inter_age * (it - (t_inter_age + len_inter_age / 2)));
        frac_icu[it, ipop] = frac_icu_0[ipop] * frac_icu_multiplier ^ inv_logit(9.19024 / len_inter_age * (it - (t_inter_age + len_inter_age / 2)));
        frac_mort[it, ipop] = frac_icu_0[ipop] * frac_mort_multiplier ^ inv_logit(9.19024 / len_inter_age * (it - (t_inter_age + len_inter_age / 2)));
      }
    }


    // Calculate beta for each time point
    for (ipop in 1:npops) {
      beta_0[ipop] = r0[ipop] / (frac_hosp[1, ipop] * duration_pre_hosp[ipop] + (1 - frac_hosp[1, ipop]) * duration_rec_mild);
      for (it in 1:nt) {
        beta[it, ipop] = beta_0[ipop];
        for (iinter in 1:ninter) {
          //k <- 2/s * qlogis(0.99) # = -2/s * qlogis(0.01) --> 2 * qlogis(0.99) = 9.19024
          //f <- m ^ plogis(k * (t - (d + s/2)))
          beta[it, ipop] = beta[it, ipop] * beta_multiplier[iinter, ipop] ^ inv_logit(9.19024 / len_inter[iinter] * (it - (t_inter[iinter] + len_inter[iinter] / 2))); //TODO: document this; maybe could speed up too
        }
      }
    }



    // initial cond
    zero = ini_exposed[1] * 1e-15; //should be zero, hack for RStan (causes problems with RHat if constant)
    for (ipop in 1:npops) {
      for (icompartment in 1:ncompartments) {
        x[icompartment,1,ipop] = zero;
      }
      x[S,1,ipop] = population[ipop] - ini_exposed[ipop];
      x[E,1,ipop] = ini_exposed[ipop];
      Hadmits[1,ipop] = zero;
    }

    //////////////////////////////////////////
    // the SEIR model
    //print(frac_icu, change_frac_icu, frac_icu_0)
    for (it in 1:nt-1){
      //////////////////////////////////////////
      // set transition variables

      beta_mat[it, 1, 1] = beta[it, 1] * ((N1 + N2) / N1 - (1 - m[1]) * N2 / N1);
      beta_mat[it, 2, 2] = beta[it, 2] * ((N1 + N2) / N2 - (1 - m[2]) * N1 / N2);
      beta_mat[it, 1, 2] = beta[it, 1] * (s1[1] * (1 - m[1]) + s1[2] * (1 - m[2])) + beta[it, 2] * (s1[3] * (1 - m[1]) + s1[4] * (1 - m[2]));
      beta_mat[it, 2, 1] = beta[it, 1] * (s2[1] * (1 - m[1]) + s2[2] * (1 - m[2])) + beta[it, 2] * (s2[3] * (1 - m[1]) + s2[4] * (1 - m[2]));

      for (ipop1 in 1:npops) {
        newE[ipop1] = 0;
      }
      for (ipop1 in 1:npops) {
        for (ipop2 in 1:npops) {
          newE[ipop1] = newE[ipop1] + beta_mat[it, ipop1, ipop2] * x[S,it,ipop1] * (x[Imild,it,ipop2] + x[Ipreh,it,ipop2]) / (N1 + N2);
        }
      }
      for (ipop in 1:npops) {
        newE[ipop] = fmin(newE[ipop], x[S,it,ipop]);
      }

      for (ipop in 1:npops) {
        newI[ipop] = 1.0/duration_latent * x[E, it, ipop];
        newhosp[ipop] = 1.0/duration_pre_hosp[ipop] * x[Ipreh,it,ipop];
        newrec_mild[ipop] = 1.0/duration_rec_mild * x[Imild,it,ipop];
        newrec_mod[ipop] = 1.0/duration_hosp_mod[ipop] * x[Hmod,it,ipop];
        leave_icu[ipop] = 1.0/duration_hosp_icu[ipop] * x[Hicu, it,ipop];

        //////////////////////////////////////////
        // S -> E -> I

        x[S, it+1, ipop] = x[S, it, ipop] - newE[ipop];
        x[E, it+1, ipop] = x[E, it, ipop] + newE[ipop] - newI[ipop];
        x[Imild, it+1, ipop] = x[Imild, it, ipop] + newI[ipop] * (1 - frac_hosp[it, ipop]) - newrec_mild[ipop];
        x[Ipreh, it+1, ipop] = x[Ipreh, it, ipop] + newI[ipop] * frac_hosp[it, ipop] - newhosp[ipop];
        x[Hmod, it+1, ipop] = x[Hmod, it, ipop] + newhosp[ipop] * (1 - frac_icu[it, ipop]) - newrec_mod[ipop];
        x[Hicu, it+1, ipop] = x[Hicu, it, ipop] + newhosp[ipop] * frac_icu[it, ipop] - leave_icu[ipop];
        x[Rlive, it+1, ipop] = x[Rlive, it, ipop] + newrec_mild[ipop] + newrec_mod[ipop] + leave_icu[ipop] * (1 - frac_mort[it, ipop]);
        x[Rmort, it+1, ipop] = x[Rmort, it, ipop] + leave_icu[ipop] * frac_mort[it, ipop];

        // cumulative hospital admissions
        Hadmits[it+1, ipop] = Hadmits[it, ipop] + newhosp[ipop];

        //////////////////////////////////////////
        // test
        if (fabs(sum(x[:,it+1,ipop])-population[ipop]) > 1e-1){
          reject("Model is leaking, net gain: ", sum(x[:,it+1,ipop])-population[ipop])
        }
      }
    }
  }

  // // Data for fitting
  for (ipop in 1:npops) {
    for (it in 1:nt) {
      for (itype in 1:nobs_types) {
        if (itype == obs_hosp_census) {
          sim_data[itype, it, ipop] = x[Hmod, it, ipop] + x[Hicu, it, ipop];
        } else if (itype == obs_icu_census) {
          sim_data[itype, it, ipop] = x[Hicu, it, ipop];
        } else if (itype == obs_cum_deaths) {
          sim_data[itype, it, ipop] = x[Rmort, it, ipop];
        } else if (itype == obs_cum_admits) {
          sim_data[itype, it, ipop] = Hadmits[it, ipop];
        } else {
          reject("unexpected itype")
        }
      }
    }
  }

}
model {
  ////////////////////////////////////////
  // prior distributions
  for (ipop in 1:npops) {
    r0[ipop] ~ normal(mu_r0[ipop], sigma_r0[ipop]);

    duration_pre_hosp[ipop] ~ normal(mu_duration_pre_hosp[ipop], sigma_duration_pre_hosp[ipop]);
    duration_hosp_mod[ipop] ~ normal(mu_duration_hosp_mod[ipop], sigma_duration_hosp_mod[ipop]);
    duration_hosp_icu[ipop] ~ normal(mu_duration_hosp_icu[ipop], sigma_duration_hosp_icu[ipop]);

    frac_hosp_0[ipop] ~ normal(mu_frac_hosp[ipop], sigma_frac_hosp[ipop]);
    frac_icu_0[ipop] ~ normal(mu_frac_icu[ipop], sigma_frac_icu[ipop]);
    frac_mort_0[ipop] ~ normal(mu_frac_mort[ipop], sigma_frac_mort[ipop]);

    ini_exposed[ipop] ~ exponential(lambda_ini_exposed[ipop]);
  }

  duration_latent ~ normal(mu_duration_latent, sigma_duration_latent);
  duration_rec_mild ~ normal(mu_duration_rec_mild, sigma_duration_rec_mild);

  frac_PUI ~ normal(mu_frac_pui, sigma_frac_pui);

  for (iinter in 1:ninter) {
    for (ipop in 1:npops) {
      beta_multiplier[iinter,ipop] ~ normal(mu_beta_inter[iinter,ipop], sigma_beta_inter[iinter,ipop]);
    }
    t_inter[iinter] ~ normal(mu_t_inter[iinter], sigma_t_inter[iinter]);
    len_inter[iinter] ~ normal(mu_len_inter[iinter], sigma_len_inter[iinter]);
  }



  frac_hosp_multiplier ~ normal(mu_frac_hosp_multiplier, sigma_frac_hosp_multiplier);
  frac_icu_multiplier ~ normal(mu_frac_icu_multiplier, sigma_frac_icu_multiplier);
  frac_mort_multiplier ~ normal(mu_frac_mort_multiplier, sigma_frac_mort_multiplier);


  // //////////////////////////////////////////
  // // fitting observations
  for (ipop in 1:npops) {
    for (itype in 1:nobs_types) {
      sigma_obs[itype, ipop] ~ exponential(1.0);
    }
  }

  {
    vector[nobs_notmissing] error; //could try to vectorize this again
    real obs;
    real sim;
    int cnt = 0;

    for (ipop in 1:npops) {
      for (iobs in 1:nobs){
        for (itype in 1:nobs_types) {
          if (obs_data_conf[itype, iobs, ipop] > 0) {
            cnt = cnt + 1;
            obs = obs_data_conf[itype, iobs, ipop];
            if (obs_data_pui[itype, iobs, ipop] > 0) {
              obs = obs + obs_data_pui[itype, iobs, ipop] * frac_PUI[itype];
            }
            sim = sim_data[itype, tobs[iobs], ipop];
            error[cnt] = (obs - sim) / sigma_obs[itype, ipop];
          }
        }
      }
    }
    error ~ std_normal();
  }
}
generated quantities{
  // real<lower=0.0> Rt[nt];
  // for (it in 1:nt) {
    //   Rt[it] = beta[it] * (frac_hosp * duration_pre_hosp + (1 - frac_hosp) * duration_rec_mild) * x[S, it] / population;
    // }
}
