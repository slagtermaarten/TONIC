data {
  int prior_only;  // should the likelihood be ignored?
  int use_baseline; // should the baseline expression be used?
  int use_clinical_response; // should clinical response be used?
  int use_TL; // should the number of prior treatment lines be used?
  int use_LN; // should lymph node only be used?
  int<lower=1> N_arms; // number of induction arms
  int<lower=1> N_obs;  // total number of observations
  int arm[N_obs];
  int<lower=0> df_sigma_arm;
  int<lower=0> df_sigma_FC;
  int<lower=0> df_FC_mu;
  int<lower=0> df_b;
  int<lower=0> df_response_contrib;
  vector[N_obs] FC;  // response variable
  vector[N_obs] baseline;  // baseline (confounding) variable
  vector[N_obs] clinical_response;  // Did a response occur?
  vector[N_obs] TL;  // Disease free interval between primary and first met
  vector[N_obs] LN;  // Patients with lymph node only metastases
}
parameters {
  real FC_mu;  
  vector[N_arms] FC_arm;  // group-level effects
  real<lower=0> sigma_arm;  // group-level standard deviations
  real<lower=0> sigma_FC;  // residual SD at the individual-level
  real bc[use_baseline]; // contribution of baseline value to observed FC
  real rc[use_clinical_response]; // contribution of baseline value to observed FC
  real tc[use_TL]; // contribution of TL to observed FC
  real lc[use_LN]; // contribution of lymph-node only metastasis to observed FC
}
model {
  // priors including all constants
  if (use_baseline) {
    bc ~ student_t(df_b, 0, 10);
  }
  if (use_clinical_response) {
    rc ~ student_t(df_response_contrib, 0, 10);
  }
  if (use_TL) {
    // Use the same amount of degrees of freedom on purpose
    tc ~ student_t(df_response_contrib, 0, 10);
  }
  if (use_LN) {
    // Use the same amount of degrees of freedom on purpose
    lc ~ student_t(df_response_contrib, 0, 10);
  }
  sigma_arm ~ student_t(df_sigma_arm, 0, 10);
  sigma_FC ~ student_t(df_sigma_FC, 0, 10);
  FC_mu ~ student_t(df_FC_mu, 0, 10);
  FC_arm ~ student_t(sigma_arm, FC_mu, 10);

  // likelihood including all constants
  if (!prior_only) {
    vector[N_obs] mu = FC_arm[arm];
    if (use_baseline) {
      mu += bc[1] * baseline;
    }
    if (use_clinical_response) {
      mu += rc[1] * clinical_response;
    }
    if (use_TL) {
      mu += tc[1] * TL;
    }
    if (use_LN) {
      mu += lc[1] * LN;
    }
    target += student_t_lpdf(FC | sigma_FC, mu, 10);
  }
}
generated quantities {
  // 2019-02-01 09:18 Compilation problems unless all vars are declared first
  vector[N_arms - 1] FC_arm_normalized;
  vector[N_obs] log_lik;
  vector[N_obs] mu = FC_arm[arm];

  FC_arm_normalized = FC_arm[2:5] - FC_arm[1];

  if (!prior_only) {
    if (use_baseline) {
      mu += bc[1] * baseline;
    }
    if (use_clinical_response) {
      mu += rc[1] * clinical_response;
    }
    if (use_TL) {
      mu += tc[1] * TL;
    }
    if (use_LN) {
      mu += lc[1] * LN;
    }
    for (n in 1:N_obs) {
      log_lik[n] = student_t_lpdf(FC[n] | sigma_FC, mu[n], 10);
    }
  }
}

// vim: set ft=c:
