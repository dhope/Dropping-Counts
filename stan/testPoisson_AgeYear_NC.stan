// Model of Counts of Shorebird Droppings at Boundary Bay
// This is wrong as we want the pattern to differ between migratory periods and right now
// it is a single point.
// cov_GPL2 macro extracted from ulam object with get_stancode
functions{
    matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + delta;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + delta;
        return K;
    }
}


data {
  int n_doy;
  int n_year;
  int n_d2s;
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



// Parameters
parameters {
  real alpha;
  vector[n_year] X;
  vector<lower=0>[5] etasq;
  vector<lower=0>[5] rhosq;
  vector[n_d2s] z_d2s[4];
  vector[n_doy] z_doy;


}

transformed parameters{
  vector[n_doy] X1;
  vector[n_d2s] X2[4];
  matrix[n_doy, n_doy] SIGMA_doy;
  matrix[n_doy, n_doy] L_SIGMA_doy;
  matrix[n_d2s, n_d2s] SIGMA_d2s[4];
  matrix[n_d2s, n_d2s] L_SIGMA_d2s[4];
  vector[n] lambda;

  SIGMA_doy = cov_GPL2(D_doy, etasq[5], rhosq[5], 0.01);
  L_SIGMA_doy = cholesky_decompose(SIGMA_doy);
  X1 = L_SIGMA_doy * z_doy;
  for (k in 1:4) {
    SIGMA_d2s[k] = cov_GPL2(D_d2s, etasq[k], rhosq[k], 0.01);
    L_SIGMA_d2s[k] = cholesky_decompose(SIGMA_d2s[k]);
    X2[k] = L_SIGMA_d2s[k] * z_d2s[k];
    }
  for(i in 1:n) lambda[i] = alpha +  X[year[i]] +
                X1[doy[i]] + X2[gr[i]][d2s[i]] + exposure[i];
}

// Models
model {
    alpha ~ normal(0, .1);
    X ~ normal(0, 0.1);
    rhosq ~ exponential( 0.5 );
    etasq ~ exponential( 2 );
    for(i in 1:4) z_d2s[i] ~ normal( 0 , 1 );
    z_doy ~ normal( 0 , 1 );

   target +=  poisson_log_lpmf(DC | lambda);
}

