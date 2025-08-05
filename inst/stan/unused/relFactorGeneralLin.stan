functions {
  matrix rbf_kernel(vector f, real alpha, real rho) {
    int N = rows(f);
    array[N] real f_array = to_array_1d(f);
    matrix[N, N] out = gp_exp_quad_cov(f_array, alpha, rho);
    for (n in 1 : N) {
      out[n, n] += 1e-9;
    }
    return out;
  }
  
  // Returns NxF matrix of omegas
  matrix omega_one(matrix lambda_loc_mat, array[,] int F_inds,
                   array[] int F_inds_num, matrix shat) {
    int N = rows(shat);
    int F = rows(lambda_loc_mat);
    int J = cols(shat);
    matrix[N, J] vhat = shat .* shat;
    matrix[N, F] vhat_sum;
    vector[F] numerator;
    matrix[N, F] omega;
    for (f in 1 : F) {
      vhat_sum[ : , f] = vhat[ : , F_inds[f, 1 : F_inds_num[f]]]
                         * rep_vector(1, F_inds_num[f]);
      numerator[f] = sum(lambda_loc_mat[f,  : ]) ^ 2;
      omega[ : , f] = numerator[f] ./ (numerator[f] + vhat_sum[ : , f]);
    }
    return omega;
  }
  
  // Creates FxJ matrix of binary codes; 1 where F has a loading.
  matrix loadings_to_ones(array[,] int F_inds, array[] int F_inds_num) {
    int F = size(F_inds);
    int J = size(F_inds[1]);
    matrix[F, J] lambda_ones;
    for (f in 1 : F) {
      for (j in 1 : J) {
        lambda_ones[f, j] = 0; // Init to zero
      }
      lambda_ones[f, F_inds[f, 1 : F_inds_num[f]]] = rep_row_vector(1,
                                                                    F_inds_num[f]);
    }
    return lambda_ones;
  }
  
  vector ones(int num) {
    vector[num] ones;
    for (n in 1 : num) {
      ones[n] = 1;
    }
    return ones;
  }
  
  // Returns NxF matrix of omegas
  matrix omega_two(matrix lambda_loc_mat, array[,] int F_inds,
                   array[] int F_inds_num, matrix theta_cor_L, matrix shat) {
    int N = rows(shat);
    int F = rows(lambda_loc_mat);
    int J = cols(shat);
    matrix[N, J] vhat = shat .* shat;
    matrix[F, F] theta_loc_cov = multiply_lower_tri_self_transpose(theta_cor_L);
    // matrix[F*2,F*2] theta_cov = multiply_lower_tri_self_transpose(theta_cor_L);
    // matrix[F,F] theta_loc_cov = theta_cov[1:F,1:F];
    matrix[J, J] implied_cov_fixed = lambda_loc_mat' * theta_loc_cov
                                     * lambda_loc_mat;
    matrix[F, J] lambda_ones = loadings_to_ones(F_inds, F_inds_num);
    vector[F] numerator;
    matrix[N, F] omega;
    for (f in 1 : F) {
      numerator[f] = sum(lambda_loc_mat[f,  : ]) ^ 2;
      omega[ : , f] = numerator[f]
                      ./ (lambda_ones[f] * implied_cov_fixed
                          * lambda_ones[f]'
                          + vhat[ : , F_inds[f, 1 : F_inds_num[f]]]
                            * ones(F_inds_num[f]));
    }
    return omega;
  }
}
data {
  int N;
  int J;
  int F;
  array[F, J] int F_inds;
  
  matrix[N, J] x;
}
transformed data {
  array[F] int F_inds_num;
  int N_loadings;
  for (f in 1 : F) {
    F_inds_num[f] = 0;
    for (j in 1 : J) {
      if (F_inds[f, j] != 0) {
        F_inds_num[f] += 1;
      }
    }
  }
  N_loadings = sum(F_inds_num);
}
parameters {
  // Loadings
  row_vector<lower=0>[N_loadings] lambda_loc; // Pos. constraint!
  row_vector<lower=0>[N_loadings] lambda_sca; // Pos. constraint!
  // Intercepts
  row_vector[J] nu_loc;
  row_vector[J] nu_sca;
  //Latents
  matrix[N, F] theta_loc_z;
  matrix[N, F] theta_sca_z;
  cholesky_factor_corr[F] theta_cor_L;
  
  // GP
  row_vector[F] linear_beta;
  // vector<lower=0>[F] rho;
  // vector<lower=0>[F] alpha;
  // matrix[N,F] gp_z;
}
transformed parameters {
  matrix[N, F] theta_loc = theta_loc_z * theta_cor_L';
  matrix[N, F] theta_sca = theta_sca_z;
  matrix[F, J] lambda_loc_mat;
  matrix[F, J] lambda_sca_mat;
  matrix[N, J] yhat;
  matrix[N, J] shat;
  matrix[N, F * 2] theta;
  // theta_sca += theta_loc .* rep_matrix(linear_beta,N);
  // theta[,(F+1):(F*2)] += theta[,1:F] .* rep_matrix(linear_beta,N);
  for (f in 1 : F) {
    // theta_sca[,f] += cholesky_decompose(rbf_kernel(theta_loc[,f],alpha[f],rho[f]))*gp_z[,f];
    theta_sca[ : , f] += theta_loc[ : , f] * linear_beta[f];
    theta[ : , f] = theta_loc[ : , f];
    theta[ : , F + f] = theta_sca[ : , f];
    // theta[,(f+F)] += cholesky_decompose(rbf_kernel(theta[,f],alpha[f],rho[f]))*gp_z[,f];
  }
  
  // Init to zero
  for (f in 1 : F) {
    for (j in 1 : J) {
      lambda_loc_mat[f, j] = 0;
      lambda_sca_mat[f, j] = 0;
    }
  }
  // Unroll lambda_loc and lambda_sca
  {
    int count = 1;
    for (f in 1 : F) {
      lambda_loc_mat[f, F_inds[f, 1 : F_inds_num[f]]] = lambda_loc[count : (count
                                                                    - 1
                                                                    + F_inds_num[f])];
      lambda_sca_mat[f, F_inds[f, 1 : F_inds_num[f]]] = lambda_sca[count : (count
                                                                    - 1
                                                                    + F_inds_num[f])];
      count += F_inds_num[f];
    }
  }
  //Predictions
  yhat = rep_matrix(nu_loc, N) + theta_loc * lambda_loc_mat;
  shat = exp(rep_matrix(nu_sca, N) + theta_sca * lambda_sca_mat);
}
model {
  // Priors
  // lambda_loc ~ std_normal();
  // lambda_sca ~ std_normal();
  lambda_loc ~ gamma(6, 10);
  lambda_sca ~ lognormal(0, 1);
  nu_loc ~ std_normal();
  nu_sca ~ std_normal();
  to_vector(theta_loc_z) ~ std_normal();
  to_vector(theta_sca_z) ~ std_normal();
  theta_cor_L ~ lkj_corr_cholesky(1);
  
  linear_beta ~ std_normal();
  
  // to_vector(gp_z) ~ std_normal();
  // alpha ~ student_t(3,0,2);
  // rho ~ normal(0,2);
  
  // Likelihood
  to_vector(x) ~ normal(to_vector(yhat), to_vector(shat));
}
generated quantities {
  //matrix omega_one(matrix lambda_loc_mat,int[,] F_inds, int[] F_inds_num, matrix shat){
  matrix[N, F] omega1 = omega_one(lambda_loc_mat, F_inds, F_inds_num, shat);
  // matrix omega_two(matrix lambda_loc_mat, int[,] F_inds, int[] F_inds_num, matrix theta_cor_L, matrix shat){
  matrix[N, F] omega2 = omega_two(lambda_loc_mat, F_inds, F_inds_num,
                                  theta_cor_L, shat);
  matrix[1, F] omega1_expected = omega_one(lambda_loc_mat, F_inds,
                                           F_inds_num,
                                           exp(rep_matrix(nu_sca, 1)));
  matrix[1, F] omega2_expected = omega_two(lambda_loc_mat, F_inds,
                                           F_inds_num, theta_cor_L,
                                           exp(rep_matrix(nu_sca, 1)));
  matrix[F, F] theta_cor = multiply_lower_tri_self_transpose(theta_cor_L);
}

