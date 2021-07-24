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

   int num_zeros(int[] y) {
    int s = 0;
    for (n in 1:size(y))
      s += (y[n] == 0);
    return s;
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
  vector[n] exposure; // Minutes of exposure
  int DC[n]; // Counts summed by 15 quadrats
  matrix[n_doy, n_doy] D_doy; // Distance between points in day of year
  matrix[n_d2s, n_d2s] D_d2s;

}
transformed data {
  vector[n] log_exposure;
  int<lower = 0> N_zero = num_zeros(DC);
  int<lower = 1> y_nonzero[n - N_zero];
  int N_nonzero = 0;
  for (i in 1:n) {
    if (DC[i] == 0) continue;
    N_nonzero += 1;
    y_nonzero[N_nonzero] = DC[i];
  }
  log_exposure = log(exposure);
}


// Parameters
parameters {
  real<lower=0, upper=1> theta;
  real alpha;
  vector[n_year] X;
  vector<lower=0>[2] etasq;
  vector<lower=0>[2] rhosq;
  vector[n_doy] X1;
  vector[n_d2s] X2;

}

transformed parameters{
  matrix[n_doy, n_doy] SIGMA_doy;
  matrix[n_d2s, n_d2s] SIGMA_d2s;
  vector[n] lambda;
  SIGMA_doy = cov_GPL2(D_doy, etasq[1], rhosq[1], 0.01);
  SIGMA_d2s = cov_GPL2(D_d2s, etasq[2], rhosq[2], 0.01);
  lambda = exp(alpha +  X[year] + X1[doy] + X2[d2s] + log_exposure);
}

// Models
model {
    theta ~ beta(1,1); // prior on zero inflated param
    alpha ~ normal(0, .1);
    X ~ normal(0, 0.1);
    rhosq ~ exponential( 0.5 );
    etasq ~ exponential( 2 );
    X1 ~ multi_normal(rep_vector(0,n_doy), SIGMA_doy);
    X2 ~ multi_normal(rep_vector(0,n_d2s), SIGMA_d2s);

   target
     += N_zero
          * log_sum_exp(bernoulli_lpmf(1 | theta),
                        bernoulli_lpmf(0 | theta)
                          + poisson_lpmf(0 | lambda));
   target += N_nonzero * bernoulli_lpmf(0 | theta);
   target += poisson_lpmf(y_nonzero | lambda);
}

