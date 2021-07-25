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
  vector<lower=0>[7] etasq;
  vector<lower=0>[7] rhosq;
  vector[n_doy] X1[2];
  vector[n_d2s] X2[4];

}

transformed parameters{
  matrix[n_doy, n_doy] SIGMA_doy[2];
  matrix[n_d2s, n_d2s] SIGMA_d2s[4];
  vector[n] lambda;
  for(y in 1:n_year) SIGMA_doy[y] = cov_GPL2(D_doy, etasq[4+y], rhosq[4+y], 0.01);
  for (k in 1:4) SIGMA_d2s[k] = cov_GPL2(D_d2s, etasq[k], rhosq[k], 0.01);
  for(i in 1:n) lambda[i] = alpha +  X[year[i]] +
      X1[year[i]][doy[i]] + X2[gr[i]][d2s[i]] + exposure[i];
}

// Models
model {
  alpha ~ normal(0, .1);
  X ~ normal(0, 0.1);
  rhosq ~ exponential( 0.5 );
  etasq ~ exponential( 2 );
  for(y in 1:2) X1[y] ~ multi_normal(rep_vector(0,n_doy), SIGMA_doy[y]);
  for (z in 1:4) X2[z] ~ multi_normal(rep_vector(0,n_d2s), SIGMA_d2s[z]);

  target +=  poisson_log_lpmf(DC | lambda);
}

