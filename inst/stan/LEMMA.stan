// Modification of Santa Cruz County COVID-19 Model
// https://github.com/jpmattern/seir-covid19
// by Jann Paul Mattern and Mikala Caton


data {
  //////////////////////////////////////////
  // data required to run model

  // obs_data_conf = total hosp census, ICU census, cumulative deaths, cumulative total hosp
  // obs_data_pui =  same

  int<lower=0> nobs;                     // number of timesteps with observations
  int<lower=0> tobs[nobs];               // obs times; this is a vector
  int<lower=0> npops;
  real<lower=-1.0> obs_data_conf[nobs, npops];  // observed confirmed (-1 = NA)
  real<lower=-1.0> obs_data_pui[nobs, npops];   // observed PUI (-1 = NA)
  int<lower=0> nt;                       // number of time steps
  real<lower=0, upper = 1> mobility[nt, npops, npops];
  real<lower=0> population[npops];                             // total population
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

  real<lower=0.0> frac_pui;     //fraction of PUI that are true COVID+

  real<lower=0.0> mu_frac_hosp;           // mean ICU + non-ICU
  real<lower=0.0> sigma_frac_hosp;        // sd ICU + non-ICU
  real<lower=0.0> frac_icu;            // mean ICU as fraction of hosp
  real<lower=0.0> frac_mort;           // mean mortality as fraction of ICU

  real<lower=0.0> lambda_ini_exposed[npops];     // parameter for initial conditions of "exposed"

  //////////////////////////////////////////
  // interventions

  int<lower=0> ninter;                      // number of interventions
  real<lower=1.0> mu_t_inter[ninter];       // mean start time of each interventions
  real<lower=0.0> sigma_t_inter[ninter];    // sd start time of each interventions
  real<lower=1.0> mu_len_inter[ninter];     // mean length of each intervention
  real<lower=0.0> sigma_len_inter[ninter];  // sd length of each intervention
  real<lower=0.0> mu_beta_inter[ninter, npops];    // mean change in beta through intervention
  real<lower=0.0> sigma_beta_inter[ninter, npops]; // sd change in beta through intervention

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
      if (obs_data_conf[iobs, ipop] > 0) {
        nobs_notmissing = nobs_notmissing + 1;
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

  real<lower=1.0> duration_pre_hosp;
  real<lower=1.0> duration_hosp_mod;
  real<lower=1.0> duration_hosp_icu;

  real<lower=0.005, upper=1.0> frac_hosp;


  real<lower=0> ini_exposed[npops];

  real<lower=0> sigma_obs;

  real<lower=0.0> r0;
  real<lower=0.0> beta_multiplier[ninter,npops];
  real<lower=1.0> t_inter[ninter];
  real<lower=1.0> len_inter[ninter];

  real mobility_coef0;
  real mobility_coef1;
}
transformed parameters {

  real<lower=0.0> x[ncompartments,nt,npops];
  real<lower=0.0> sim_data[nt,npops];
  real<lower=0.0, upper=2.0> beta[nt,npops];
  real<lower=0.0> beta_mat[nt-1,npops,npops]; //could make local, drop nt

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
    real total_population = sum(population);
    real m;
    // real beta_mat[npops,npops];

    //////////////////////////////////////////

    // Calculate beta for each time point
    for (ipop in 1:npops) {
      beta_0[ipop] = r0 / (frac_hosp * duration_pre_hosp + (1 - frac_hosp) * duration_rec_mild);
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
    }

    //////////////////////////////////////////
    // the SEIR model
    //print(frac_icu, change_frac_icu, frac_icu_0)
    for (it in 1:nt-1){
      //////////////////////////////////////////
      // set transition variables
      for (ipop1 in 1:npops) {
        newE[ipop1] = 0;

        for (ipop2 in 1:npops) {
          m = inv_logit(mobility_coef0 + mobility_coef1 * mobility[it, ipop1, ipop2]);
          if (ipop1 == ipop2) {
            beta_mat[it, ipop1, ipop2] = beta[it, ipop1] * (1 + m * (total_population / population[ipop1] - 1));
          } else {
            beta_mat[it, ipop1, ipop2] = 0.5 * (beta[it, ipop1] + beta[it, ipop2]) * (1 - m);
          }
        }
      }

      for (ipop1 in 1:npops) {
        for (ipop2 in 1:npops) {
          newE[ipop1] = newE[ipop1] + beta_mat[it, ipop1, ipop2] * x[S,it,ipop1] * (x[Imild,it,ipop2] + x[Ipreh,it,ipop2]) / total_population;
        }
      }

      for (ipop in 1:npops) {
        newE[ipop] = fmin(newE[ipop], x[S,it,ipop]);
        newI[ipop] = 1.0/duration_latent * x[E, it, ipop];
        newhosp[ipop] = 1.0/duration_pre_hosp * x[Ipreh,it,ipop];
        newrec_mild[ipop] = 1.0/duration_rec_mild * x[Imild,it,ipop];
        newrec_mod[ipop] = 1.0/duration_hosp_mod * x[Hmod,it,ipop];
        leave_icu[ipop] = 1.0/duration_hosp_icu * x[Hicu, it,ipop];

        //////////////////////////////////////////
        // S -> E -> I

        x[S, it+1, ipop] = x[S, it, ipop] - newE[ipop];
        x[E, it+1, ipop] = x[E, it, ipop] + newE[ipop] - newI[ipop];
        x[Imild, it+1, ipop] = x[Imild, it, ipop] + newI[ipop] * (1 - frac_hosp) - newrec_mild[ipop];
        x[Ipreh, it+1, ipop] = x[Ipreh, it, ipop] + newI[ipop] * frac_hosp - newhosp[ipop];
        x[Hmod, it+1, ipop] = x[Hmod, it, ipop] + newhosp[ipop] * (1 - frac_icu) - newrec_mod[ipop];
        x[Hicu, it+1, ipop] = x[Hicu, it, ipop] + newhosp[ipop] * frac_icu - leave_icu[ipop];
        x[Rlive, it+1, ipop] = x[Rlive, it, ipop] + newrec_mild[ipop] + newrec_mod[ipop] + leave_icu[ipop] * (1 - frac_mort);
        x[Rmort, it+1, ipop] = x[Rmort, it, ipop] + leave_icu[ipop] * frac_mort;

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
      sim_data[it, ipop] = x[Hmod, it, ipop] + x[Hicu, it, ipop];
    }
  }
}
model {
  ////////////////////////////////////////
  // prior distributions

  r0 ~ normal(mu_r0, sigma_r0);
  mobility_coef0 ~ normal(0, 5);
  mobility_coef1 ~ normal(0, 5);

  duration_pre_hosp ~ normal(mu_duration_pre_hosp, sigma_duration_pre_hosp);
  duration_hosp_mod ~ normal(mu_duration_hosp_mod, sigma_duration_hosp_mod);
  duration_hosp_icu ~ normal(mu_duration_hosp_icu, sigma_duration_hosp_icu);

  frac_hosp ~ normal(mu_frac_hosp, sigma_frac_hosp);

  for (ipop in 1:npops) {
    ini_exposed[ipop] ~ exponential(lambda_ini_exposed[ipop]);
  }

  duration_latent ~ normal(mu_duration_latent, sigma_duration_latent);
  duration_rec_mild ~ normal(mu_duration_rec_mild, sigma_duration_rec_mild);

  for (iinter in 1:ninter) {
    for (ipop in 1:npops) {
      beta_multiplier[iinter,ipop] ~ normal(mu_beta_inter[iinter,ipop], sigma_beta_inter[iinter,ipop]);
    }
    t_inter[iinter] ~ normal(mu_t_inter[iinter], sigma_t_inter[iinter]);
    len_inter[iinter] ~ normal(mu_len_inter[iinter], sigma_len_inter[iinter]);
  }

  // //////////////////////////////////////////
  // // fitting observations
  sigma_obs ~ exponential(1.0);

  {
    vector[nobs_notmissing] error; //could try to vectorize this again
    real obs;
    real sim;
    int cnt = 0;

    for (ipop in 1:npops) {
      for (iobs in 1:nobs){
        if (obs_data_conf[iobs, ipop] > 0) {
          cnt = cnt + 1;
          obs = obs_data_conf[iobs, ipop];
          if (obs_data_pui[iobs, ipop] > 0) {
            obs = obs + obs_data_pui[iobs, ipop] * frac_pui;
          }
          sim = sim_data[tobs[iobs], ipop];
          error[cnt] = (obs - sim) / sigma_obs;
        }
      }
    }
    error ~ std_normal();
  }
}
generated quantities{
}
