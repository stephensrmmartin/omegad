functions {
  // Computes lambda (eigenfunction's eigenvalue). This is FIXED; use td{}
  vector lambdas(real L, int M){
    vector[M] lambdas;
    for(m in 1:M){
      lambdas[m] = ((pi()*m)/(2*L))^2;
    }
    return(lambdas);
  }
  // Returns N-length row. FIXED; use in td{}.
  vector basis_phi(real L, int m, vector x){
    int N = rows(x);
    vector[N] basis_phi = 1/sqrt(L) * sin((pi()*m)/(2*L) * (x + L));
    return(basis_phi);
  }

  // Returns NxM. FIXED; use in td{}.
  // This is the eigenfunctions on X.
  matrix basis_phis(real L, int M, vector x){
    int N = rows(x);
    matrix[N,M] phis;
    for(m in 1:M){
      phis[,m] = basis_phi(L,m,x);
    }
    return(phis);
  }

  // Returns spectral density.
  real spd(real alpha, real rho, real lambda){
    real spd = (alpha^2)*sqrt(2*pi())*rho*exp(-.5 * (rho^2) * lambda);
    return(spd);
  }

  // Returns M-length vector of spectral densities.
  vector spds(real alpha, real rho, vector lambdas){
    int M = rows(lambdas);
    vector[M] spds;
    for(m in 1:M){
      spds[m] = spd(alpha, rho, lambdas[m]);
    }
    return(spds);
  }

  // Returns NxM matrix. Needs to be multiplied by "beta" (gp_z).
  // gp_z should *not* be N-length, but M-length.
  matrix spd_gp(matrix x_phi, real alpha, real rho, vector lambdas){
    int M = cols(x_phi);
    vector[M] spdsOut = spds(alpha, rho, lambdas);
    return(diag_post_multiply(x_phi,sqrt(spdsOut)));
  }

  // Returns N-length vector of gp output. may be slightyl faster.
  vector spd_gp_fast(matrix x_phi, real alpha, real rho, vector lambdas, vector gp_z){
    int M = cols(x_phi);
    vector[M] spdsOut = spds(alpha, rho, lambdas);
    return(x_phi * (sqrt(spdsOut) .* gp_z));
  }

  // Returns NxF matrix of omegas
  matrix omega_one(matrix lambda_loc_mat,int[,] F_inds, int[] F_inds_num, matrix shat){
    int N = rows(shat);
    int F = rows(lambda_loc_mat);
    int J = cols(shat);
    matrix[N,J] vhat = shat .* shat;
    matrix[N,F] vhat_sum;
    vector[F] numerator;
    matrix[N, F] omega;
    for(f in 1:F){
      vhat_sum[,f] = vhat[,F_inds[f,1:F_inds_num[f]]] * rep_vector(1,F_inds_num[f]);
      numerator[f] = sum(lambda_loc_mat[f,])^2;
      omega[,f] = numerator[f] ./ (numerator[f] + vhat_sum[,f]);
    }
    return(omega);

  }

  // Creates FxJ matrix of binary codes; 1 where F has a loading.
  matrix loadings_to_ones(int[,] F_inds,int[] F_inds_num){
    int F = size(F_inds);
    int J = size(F_inds[1]);
    matrix[F,J] lambda_ones;
    for(f in 1:F){
      for(j in 1:J){
        lambda_ones[f,j] = 0; // Init to zero
      }
      lambda_ones[f,F_inds[f,1:F_inds_num[f]]] = rep_row_vector(1,F_inds_num[f]);
    }
    return(lambda_ones);
  }

  vector ones(int num){
    vector[num] ones;
    for(n in 1:num){
      ones[n] = 1;
    }
    return(ones);
  }

  // Returns NxF matrix of omegas
  matrix omega_two(matrix lambda_loc_mat, int[,] F_inds, int[] F_inds_num, matrix theta_cor_L, matrix shat){
    int N = rows(shat);
    int F = rows(lambda_loc_mat);
    int J = cols(shat);
    matrix[N,J] vhat = shat .* shat;
    matrix[F,F] theta_loc_cov = multiply_lower_tri_self_transpose(theta_cor_L);
    // matrix[F*2,F*2] theta_cov = multiply_lower_tri_self_transpose(theta_cor_L);
    // matrix[F,F] theta_loc_cov = theta_cov[1:F,1:F];
    matrix[J,J] implied_cov_fixed = lambda_loc_mat' * theta_loc_cov * lambda_loc_mat;
    matrix[F,J] lambda_ones = loadings_to_ones(F_inds,F_inds_num);
    vector[F] numerator;
    matrix[N,F] omega;
    for(f in 1:F){
      numerator[f] = sum(lambda_loc_mat[f,])^2;
      omega[,f] = numerator[f] ./ (lambda_ones[f]*implied_cov_fixed*lambda_ones[f]' + vhat[,F_inds[f,1:F_inds_num[f]]]*ones(F_inds_num[f]));
    }
    return(omega);
  }

  matrix omega_total(matrix lambda_loc_mat, matrix theta_cor_L, matrix shat) {
    int N = rows(shat);
    int F = rows(lambda_loc_mat);
    int J = cols(shat);
    matrix[N,J] vhat = shat .* shat;
    matrix[F,F] theta_loc_cov = multiply_lower_tri_self_transpose(theta_cor_L);
    vector[J] one = ones(J);
    real numerator = one' * lambda_loc_mat' * theta_loc_cov * lambda_loc_mat * one;
    matrix[N, 1] omega_total;
    omega_total[,1]= numerator ./ (numerator + vhat*one);
    return(omega_total);
  }

}
data {
  int N;
  int J;
  int F;
  int F_inds[F,J];
  int M; // GP basis functions

  matrix[N,J] x;

  int P; // Number of exogenous covariates for theta_sca
  matrix[N,P] exo_x; // Design matrix for exogenous covariates.

}

transformed data {
  real L = 3.0*5.0/2.0; // gp Constant. May change later to vary based on current estimated latent values.
  vector[M] gp_lambdas = lambdas(L,M);
  matrix[N,M] gp_exo_phi[P-1];
  int F_inds_num[F];
  int N_loadings;
  for(f in 1:F){
    F_inds_num[f] = 0;
    for(j in 1:J){
      if(F_inds[f,j] != 0){
        F_inds_num[f] += 1;
      }
    }
  }
  N_loadings = sum(F_inds_num);

  for(p in 1:(P-1)){
    gp_exo_phi[p] = basis_phis(L, M, exo_x[,(p+1)]);
  }

}

parameters {
  // Loadings
  row_vector<lower=0>[N_loadings] lambda_loc; // Pos. constraint!
  row_vector<lower=0>[N_loadings] lambda_sca; // Pos. constraint!
  // Intercepts
  row_vector[J] nu_loc;
  row_vector[J] nu_sca;
  //Latents
  matrix[N,F] theta_loc_z;
  matrix[N,F] theta_sca_z;
  cholesky_factor_corr[F] theta_cor_L;

  // GP
  row_vector[F] gp_linear;
  vector<lower=0>[F] gp_rho;
  vector<lower=0>[F] gp_alpha;
  matrix[M,F] gp_z;


  //Exogenous model for scale factors
  //Using independent GPs (so kernel is based on distance on each variable, not both).
  matrix[P-1, F] exo_gp_linear;
  matrix<lower=0>[P-1,F] exo_gp_rho;
  matrix<lower=0>[P-1,F] exo_gp_alpha;
  matrix[M,F] exo_gp_z[P-1];

}

transformed parameters {
  matrix[N,F] theta_loc = theta_loc_z*theta_cor_L';
  matrix[N,M] gp_theta_phi[F];
  matrix[N,F] theta_sca = theta_sca_z;
  matrix[F,J] lambda_loc_mat;
  matrix[F,J] lambda_sca_mat;
  matrix[N,J] yhat;
  matrix[N,J] shat;
  matrix[N,F*2] theta;
  /* Latent GP */
  for(f in 1:F){
    gp_theta_phi[f] = basis_phis(L, M, theta_loc[,f]);
    theta_sca[,f] += theta_loc[,f] * gp_linear[f]; // Linear
    theta_sca[,f] += spd_gp_fast(gp_theta_phi[f], gp_alpha[f], gp_rho[f], gp_lambdas, gp_z[,f]); // SPD-BF-RBF
  }
  /* Exogenous GP */
  for(p in 1:(P-1)){
    for(f in 1:F){
      theta_sca[,f] += exo_x[,(p+1)]*exo_gp_linear[p,f]; // Linear
      theta_sca[,f] += spd_gp_fast(gp_exo_phi[p], exo_gp_alpha[p,f], exo_gp_rho[p,f],gp_lambdas, exo_gp_z[p,1:M,f]);
    }
  }

  /* Repackage */
  theta[,1:F] = theta_loc;
  theta[,(F+1):(2*F)] = theta_sca;

  // Init Loadings to zero
  for(f in 1:F){
    for(j in 1:J){
      lambda_loc_mat[f,j] = 0;
      lambda_sca_mat[f,j] = 0;
    }
  }
  // Unroll lambda_loc and lambda_sca
  {
    int count = 1;
    for(f in 1:F){
      lambda_loc_mat[f,F_inds[f,1:F_inds_num[f]]] = lambda_loc[count:(count - 1 + F_inds_num[f])];
      lambda_sca_mat[f,F_inds[f,1:F_inds_num[f]]] = lambda_sca[count:(count - 1 + F_inds_num[f])];
      count += F_inds_num[f];
    }

  }
  //Predictions
  yhat = rep_matrix(nu_loc,N) + theta_loc*lambda_loc_mat;
  shat = exp(rep_matrix(nu_sca,N) + theta_sca*lambda_sca_mat);
}

model {
  // Priors
  /* Measurement */
  // lambda_loc ~ std_normal();
  // lambda_sca ~ std_normal();
  lambda_loc ~ gamma(6,10);
  lambda_sca ~ lognormal(0,1);
  nu_loc ~ std_normal();
  // nu_sca ~ std_normal();
  nu_sca ~ normal(-.5, 1);
  to_vector(theta_loc_z) ~ std_normal();
  to_vector(theta_sca_z) ~ std_normal();
  theta_cor_L ~ lkj_corr_cholesky(1);

  /* GP */
  gp_linear ~ std_normal();
  to_vector(gp_z) ~ std_normal();
  gp_alpha ~ student_t(3,0,2);
  gp_rho ~ std_normal();

  /* Exogenous GPs */
  to_vector(exo_gp_linear) ~ std_normal();
  to_vector(exo_gp_alpha) ~ student_t(3,0,2);
  to_vector(exo_gp_rho) ~ std_normal();
  for(p in 1:(P-1)){
    to_vector(exo_gp_z[p]) ~ std_normal();
  }
  // rho ~ normal(0,2);
  // rho ~ gamma(5,10)
  // rho ~ inv_gamma(5,5);


  // Likelihood
  to_vector(x) ~ normal(to_vector(yhat),to_vector(shat));
}

generated quantities {
  //matrix omega_one(matrix lambda_loc_mat,int[,] F_inds, int[] F_inds_num, matrix shat){
  // matrix[N,F] omega1 = omega_one(lambda_loc_mat, F_inds, F_inds_num, shat);
  // matrix omega_two(matrix lambda_loc_mat, int[,] F_inds, int[] F_inds_num, matrix theta_cor_L, matrix shat){
  // matrix[N,F] omega2 = omega_two(lambda_loc_mat, F_inds, F_inds_num, theta_cor_L, shat);
  matrix[1,F] omega1_expected = omega_one(lambda_loc_mat,F_inds,F_inds_num,exp(rep_matrix(nu_sca,1)));
  matrix[1,F] omega2_expected = omega_two(lambda_loc_mat,F_inds,F_inds_num,theta_cor_L,exp(rep_matrix(nu_sca,1)));
  matrix[1,1] omega_total_expected = omega_total(lambda_loc_mat, theta_cor_L, exp(rep_matrix(nu_sca,1)));
  matrix[F,F] theta_cor = multiply_lower_tri_self_transpose(theta_cor_L);
}
