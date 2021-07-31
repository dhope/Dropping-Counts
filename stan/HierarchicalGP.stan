data {
  int n_doy;
  int n_year;
  int n_d2s;
  int n_groups; // Number of groups (4)
  int<lower=0> n; // Number of rows
  int d2s[n]; // Distance to shore
  int doy[n]; // Day of Year
  int year[n]; // Year
  real exposure[n]; // Minutes of exposure
  int DC[n]; // Counts summed by 15 quadrats
  matrix[n_doy, n_doy] D_doy; // Distance between points in day of year
  matrix[n_d2s, n_d2s] D_d2s;
  int gr[n]; // Age-Year group, 1 =Adult2007;Juvenile2007;Adult2008;Juvenile2008

}

transformed data {
  real d2s_i[n_d2s];
  vector[3] counts;
  int n_comps;
  for (t in 1:n_d2s)
    d2s_i[t] = t;
  n_comps = rows(counts);
  for (i in 1:n_comps)
    counts[i] = 2;
}
parameters {
  matrix[n_d2s,n_groups] GP_d2s_group_std;
  vector[n_d2s] d2s_std;
  vector[n_groups] group_std;
  real<lower=0> tot_var;
  real alpha;
  simplex[n_comps] prop_var;
  real<lower=0> length_GP_group_long;
  // real<lower=0> length_GP_group_short;
}
transformed parameters {
    vector[n] lambda;
  matrix[n_d2s,n_groups] GP_d2s_group;
  vector[n_d2s] d2s_re;
  vector[n_groups] groups_re;
  vector[n_comps] vars;
  real<lower=0> sigma_d2s;
  real<lower=0> sigma_group;
  real<lower=0> sigma_GP_group_long;
  // real<lower=0> sigma_GP_group_short;
  vars = n_comps * prop_var * tot_var;
  sigma_d2s = sqrt(vars[1]);
  sigma_group = sqrt(vars[2]);
  sigma_GP_group_long = sqrt(vars[3]);
  // sigma_GP_group_short = sqrt(vars[4]);
  d2s_re = sigma_d2s * d2s_std;
  groups_re = sigma_group * group_std;
  {
    matrix[n_d2s, n_d2s] cov_group_d2s;
    matrix[n_d2s, n_d2s] L_cov_group_d2s;
    cov_group_d2s = gp_exp_quad_cov(d2s_i, sigma_GP_group_long,
    length_GP_group_long);
    // + gp_exp_quad_cov(d2s_i, sigma_GP_group_short,
    // length_GP_group_short);
    // for (g in 1:n_d2s) {
    //   cov_group_d2s[g, g] = cov_group_d2s[g, g] + 1e-12;
    // }
    L_cov_group_d2s = cholesky_decompose(cov_group_d2s);
    GP_d2s_group = L_cov_group_d2s * GP_d2s_group_std;
  }

  for (i in 1:n){
      lambda[n] = alpha
      + d2s_re[d2s[i]]
      + exposure[i]
      + groups_re[gr[i]]
      // + GP_doy[doy[i]]
      + GP_d2s_group[d2s[i], gr[i]];
  }
  }
model {

  alpha ~ normal(0,0.1);
  tot_var ~ gamma(3, 3);
  prop_var ~ dirichlet(counts);
  to_vector(GP_d2s_group_std) ~ normal(0,1);
  d2s_std ~ normal(0,1);
  group_std ~ normal(0,1);
  length_GP_group_long ~ weibull(30,3);
  // length_GP_group_short ~ weibull(30,2);

  DC ~ poisson_log(lambda);

}

// generated quantities {
//   matrix[N_years,N_states] y_new;
//   matrix[N_years,N_states] y_new_pred;
//   {
//     real level;
//     level = normal_rng(0.0, sigma_year);
//     for (state in 1:N_states) {
//       for (t in 1:N_years) {
//         if (t < 12) {
//           y_new[t,state] = state_re[state]
//           + region_re[state_region_ind[state]]
//           + GP_state[t,state]
//           + GP_region[t,state_region_ind[state]]
//           + mu + year_re[t];
//           } else {
//             y_new[t,state] = state_re[state]
//             + region_re[state_region_ind[state]]
//             + GP_state[t,state]
//             + GP_region[t,state_region_ind[state]]
//             + level;
//           }
//           y_new_pred[t,state] =
//           beta_rng(inv_logit(y_new[t,state]) * nu,
//           nu * (1 - inv_logit(y_new[t,state])));
//           }
//         }
//       }
//     }
