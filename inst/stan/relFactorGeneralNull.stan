functions {
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
    matrix[F*2,F*2] theta_cov = multiply_lower_tri_self_transpose(theta_cor_L);
    matrix[F,F] theta_loc_cov = theta_cov[1:F,1:F];
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
    matrix[F,F] theta_loc_cov = multiply_lower_tri_self_transpose(theta_cor_L)[1:F, 1:F];
    vector[J] one = ones(J);
    real numerator = one' * lambda_loc_mat' * theta_loc_cov * lambda_loc_mat * one;
    matrix[N, 1] omega_total;
    omega_total[,1] = numerator ./ (numerator + vhat*one);
    return(omega_total);
  }

}
data {
  int N;
  int J;
  int F;
  int F_inds[F,J];

  matrix[N,J] x;

  int P; // Number of exogenous covariates for theta_sca
  matrix[N,P] exo_x; // Design matrix for exogenous covariates.

}

transformed data {
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


}

parameters {
  // Loadings
  row_vector<lower=0>[N_loadings] lambda_loc; // Pos. constraint!
  // Intercepts
  row_vector[J] nu_loc;
  row_vector[J] nu_sca;
  //Latents
  matrix[N,F] theta_z;
  cholesky_factor_corr[F] theta_cor_L;

}

transformed parameters {
  matrix[N,F] theta = theta_z*theta_cor_L';
  matrix[F,J] lambda_loc_mat;
  matrix[N,J] yhat;
  matrix[N,J] shat;

  // Init to zero
  for(f in 1:F){
    for(j in 1:J){
      lambda_loc_mat[f,j] = 0;
    }
  }
  // Unroll lambda_loc and lambda_sca
  {
    int count = 1;
    for(f in 1:F){
      lambda_loc_mat[f,F_inds[f,1:F_inds_num[f]]] = lambda_loc[count:(count - 1 + F_inds_num[f])];
      count += F_inds_num[f];
    }

  }
  //Predictions
  yhat = rep_matrix(nu_loc,N) + theta[,1:F]*lambda_loc_mat;
  shat = exp(rep_matrix(nu_sca, N));
}

model {
  // Priors
  /* Measurement */
  lambda_loc ~ std_normal();
  nu_loc ~ std_normal();
  //nu_sca ~ std_normal();
  nu_sca ~ normal(-.5, 1);
  to_vector(theta_z) ~ std_normal();
  theta_cor_L ~ lkj_corr_cholesky(1);

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
  matrix[F*2,F*2] theta_cor = multiply_lower_tri_self_transpose(theta_cor_L);
}
