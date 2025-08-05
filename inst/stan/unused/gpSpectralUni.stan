functions {
  // Returns N-length row. FIXED; use in td{}.
  vector basis_phi(real L, int m, vector x) {
    int N = rows(x);
    vector[N] basis_phi = 1 / sqrt(L) * sin((pi() * m) / (2 * L) * (x + L));
    return basis_phi;
  }
  
  // Returns NxM. FIXED; use in td{}.
  // This is the eigenfunctions on X.
  matrix basis_phis(real L, int M, vector x) {
    int N = rows(x);
    matrix[N, M] phis;
    for (m in 1 : M) {
      phis[ : , m] = basis_phi(L, m, x);
    }
    return phis;
  }
  
  // Returns spectral density.
  real spd(real alpha, real rho, real lambda) {
    real spd = (alpha ^ 2) * sqrt(2 * pi()) * rho
               * exp(-.5 * (rho ^ 2) * lambda ^ 2);
    return spd;
  }
  
  // Returns M-length vector of spectral densities.
  vector spds(real alpha, real rho, vector lambdas) {
    int M = rows(lambdas);
    vector[M] spds;
    for (m in 1 : M) {
      spds[m] = spd(alpha, rho, lambdas[m]);
    }
    return spds;
  }
  
  // Computes lambda (eigenfunction's eigenvalue). This is FIXED; use td{}
  vector lambdas(real L, int M) {
    vector[M] lambdas;
    for (m in 1 : M) {
      lambdas[m] = ((pi() * m) / (2 * L)) ^ 2;
    }
    return lambdas;
  }
  
  // Returns NxM matrix. Needs to be multiplied by "beta" (gp_z).
  // gp_z should *not* be N-length, but M-length.
  matrix spd_gp(matrix x_phi, real alpha, real rho, vector lambdas) {
    int M = cols(x_phi);
    vector[M] spdsOut = spds(alpha, rho, lambdas);
    return diag_post_multiply(x_phi, sqrt(spdsOut));
  }
}
data {
  int N;
  vector[N] x;
  vector[N] y;
  int M; // Number of basis functions
}
transformed data {
  real L = max(sqrt(x .* x)) * 5 / 2; // -L < x < L
  matrix[N, M] gp_x_phi = basis_phis(L, M, x);
  vector[M] gp_lambdas = lambdas(L, M);
  print(dims(gp_x_phi));
}
parameters {
  // Reg
  real<lower=0> sigma;
  
  // GP
  real<lower=0> alpha;
  real<lower=0> rho;
  vector[M] gp_z;
}
transformed parameters {
  vector[N] yhat = spd_gp(gp_x_phi, alpha, rho, gp_lambdas) * gp_z;
}
model {
  // Prior
  alpha ~ student_t(3, 0, 1);
  rho ~ normal(0, 1);
  gp_z ~ std_normal();
  
  y ~ normal(yhat, sigma);
}
generated quantities {
  vector[N] yrep;
  for (n in 1 : N) {
    yrep[n] = yhat[n] + normal_rng(0, sigma);
  }
}

