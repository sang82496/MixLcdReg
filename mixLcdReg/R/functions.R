# Generated from _main.Rmd: do not edit by hand

#' Generate data from mixture of log-concave regressions
#' 
#' The mixture of log-concave regressions model is defined as follows:
#' 
#' At time t there are n_t points generated as follows:
#' 
#' Given Z_i^{(t)} = k, Y_i^{(t)}  = theta_{k0} + theta_k^T X^{(t)} + epsilon_i^{(t)}
#' 
#' where epsilon_i^{(t)} ~ exp(g_k)
#' 
#' and P(Z_i^{(t)} =k |X^{(t)} ;alpha) = \frac{exp(alpha_{k0} + alpha_k^T X^{(t)} )}{\sum_{l=1}^K exp(alpha_{l0} + alpha_l^T X^{(t)} )}
#' 
#' This function generates Y and Z
#' 
#' For simplicity we here assume alpha = 0, K = 2, p = 2 and d = 1
#'
#' @param X a T-by-p matrix of covariates, where `X[[t]]` being the p-vector of independent variable at time t
#' @param theta0 length K list, with `theta0[[k]]` being the estimate for the intercept coefficient
#' of the regression for kth group
#' @param theta length K list, with `theta[[k]]` being the p-by-1 vector. Initial estimate for coefficients
#' of the regression for kth group
#' @param alpha (p+1)-by-K matrix. The coefficients for the cluster probabilities.
#' @param num_points T-vector of integers giving the number of points n_t to
#' generate at each time point t.
#' @param err_type a K-vector of strings, each of which is the type of error distribution
#' for the each group. Now, it is one of 'Gaussian', 'exp'
#' @return Y length T list with `Y[[t]]` being a n_t-by-d matrix
#' @return Z length T list with `Z[[t]]` being a n_t vector of cluster memberships
#' @export
generate_mix_lcd_reg <- function(X,
                                 theta0,
                                 theta,
                                 alpha,
                                 num_points,
                                 err_type = c('Gaussian', 'exp')) {
  TT = dim(X)[1]
  p = dim(X)[2]
  Y = list()
  Z = list()
  for (t in 1:TT){
    tmp = exp(alpha[1,] + X[t,] %*% alpha[2:(p+1),])
    prob = tmp / sum(tmp)
    Z[[t]] = rbinom(num_points[t], 1, prob[2])
    n1 = sum(Z[[t]] == 0)
    n2 = sum(Z[[t]] == 1)
    Xt1 = theta0[[1]] + X[t,] %*% as.matrix(theta[[1]]) 
    Xt2 = theta0[[2]] + X[t,] %*% as.matrix(theta[[2]])     
    
    if (err_type[1] == 'Gaussian') {
      e1 = rnorm(n1)
    } else {
      e1 = rexp(n1, 1) - 1
    }
    if (err_type[2] == 'Gaussian') {
      e2 = rnorm(n2)
    } else {
      e2 = rexp(n2, 1) - 1
    }
    
    Y[[t]] = rep(NA, num_points[t])
    Y[[t]][Z[[t]] == 0] = rep(Xt1, n1) + e1
    Y[[t]][Z[[t]] == 1] = rep(Xt2, n2) + e2    
  }
  return(list(Y = Y,
              Z = Z))
}

#' Plot raw data when `d = 1`
#' 
#' @param Y length T list with `Y[[t]]` being a n_t-by-d matrix
#' @export
plot_data <- function(Y) {
  par(mfrow = c(1, 2), pty = 's', cex = 0.7)
  plot(c(0, 0), xlim = c(1, TT), ylim = range(unlist(Y)), col = 'white')
  for (t in 1:TT){
    points(rep(t, length(Y[[t]])), Y[[t]])
  }
}

#' Plot data colored by cluster assignment with cluster means when `d=1`
#' 
#' @param X a T-by-p matrix of covariates, where `X[[t]]` being the p-vector of independent variable at time t
#' @param Y length T list with `Y[[t]]` being a n_t-by-d matrix
#' @param Z a length T list with `Z[[t]]` being a n_t vector of cluster assignments
#' @param theta0 length K list, with `theta0[[k]]` being the estimate for the intercept coefficient
#' of the regression for kth group
#' @param theta length K list, with `theta[[k]]` being the p-by-1 vector. Initial estimate for coefficients
#' of the regression for kth group
#' @export
plot_data_and_model <- function(X, Y, Z, theta0, theta) {
  plot(c(0, 0), xlim = c(1, TT), ylim = c(-2, 10), col = 'white')
  for (t in 1:TT){
    points(rep(t, sum(Z[[t]]==0)), Y[[t]][Z[[t]]==0], col = 'red')
    points(rep(t, sum(Z[[t]]==1)), Y[[t]][Z[[t]]==1], col = 'blue')
  }
  lines(1:TT, rep(theta0[[1]], TT) + X %*% theta[[1]], col = 'black', lwd = 2)
  lines(1:TT, rep(theta0[[2]], TT) + X %*% theta[[2]], col = 'black', lwd = 2)
}


#' Binning on Y
#' 
#' @param X a T-by-p matrix of covariates, where `X[[t]]` being the p-vector of independent variable at time t
#' @param Y length T list with `Y[[t]]` being a n_t-by-d matrix
#' @param number of binning
#' @param min_count_ratio min count ratio
#' @return Y_bin N_bin-by-d matrix, indicating the center of the bins
#' @return Y_count N_bin-vector indicating how many points belong to each bin
#' @return X_bin N_bin-by-p matrix, which is expanded version of X so that its dimension agrees with Y_bin's 
#' @return time_indicator length N_bin factor indicating the time of the data points 
#' @return N the total number of observations of Y_i^(t) in Y 
#' @export
binningY_d1 = function(X,
                       Y, 
                       B = 0,
                       min_count_ratio = 0){ # now, it only works for d = 1
  TT = length(Y)
  p = dim(X)[2]
  
  if (B == 0){
    Y_bin = matrix(unlist(Y), ncol = 1)
    nt = unlist(lapply(Y, length))
    N = 0
    time_indicator = c()
    X_bin = c()
    for (t in 1:TT){
      N = N + nt[t]
      time_indicator = c(time_indicator, rep(t, nt[t]))
      new_row = matrix(rep(X[t,], nt[t]), ncol = p, byrow = T)
      X_bin = rbind(X_bin, new_row)
    }
    Y_count = rep(1L, N)
    
  } else {
    
    Y_range = unlist(lapply(Y, range))
    minY = min(Y_range)
    maxY = max(Y_range)
    bin = seq(from=minY, to=maxY, length= B + 1)
    binnedY = lapply(Y, findInterval, bin, rightmost.closed = T)
    binnedY = lapply(binnedY, factor, levels = 1:B)
    count = lapply(binnedY, table)
    N = do.call(sum, lapply(count, sum))
    mid = rep(0, B)
    for (i in 1:B){
      mid[i] = (bin[i+1] + bin[i])/2
    }

    # resizing
    Y_bin = c()
    Y_count = c()
    X_bin = c()
    time_indicator = c()
    for (t in 1:TT){
      tbl = count[[t]]
      count_nonzero = 0
      for (j in 1:B){
        if (tbl[j] > N*min_count_ratio/TT){
          Y_bin = c(Y_bin, c(mid[j]))
          Y_count = c(Y_count, tbl[j])
          count_nonzero = count_nonzero + 1
        }
      }
      time_indicator = c(time_indicator, rep(t, count_nonzero))
      new_row = t(matrix(rep(X[t,], count_nonzero), nrow = p))
      X_bin = rbind(X_bin, new_row)
    }
    Y_bin = matrix(Y_bin, ncol = 1)
    time_indicator = factor(time_indicator)
  }
  
  return(list(Y_bin = Y_bin, 
         Y_count = Y_count,
         X_bin = X_bin,
         time_indicator = time_indicator,
         N = N
         ))
}

#' Mstep_shift (d = 1)
#' 
#' @param X_bin N_bin-by-p matrix, which is expanded version of X so that its dimension agrees with Y_bin's 
#' @param Y_bin N_bin-by-d matrix, indicating the center of the bin
#' @param theta length K list, with `theta[[k]]` being the p-by-1 vector. Initial estimate for coefficients
#' of the regression for kth group
#' @param weight N_bin by K matrix. weight[i, k] represents weight of ith residual for kth group
#' @param idx N_bin by K logical matrix. idx[i, k] represents whether the corresponding weight is above r_bar
#' @param K number of clusters
#' @return theta0 length K list, with `theta0[[k]]` being the estimate for the intercept coefficient
#' of the regression for kth group
#' @export
Mstep_shift = function(X_bin,
                       Y_bin,
                       theta, 
                       weight,
                       idx,
                       K){
  theta0 = list()
  for (k in 1:K){
    idx_k = idx[,k] 
    # except for the intercept term
    theta0[[k]] = t(Y_bin[idx_k,] - X_bin[idx_k,] %*% as.matrix(theta[[k]])) %*% weight[idx_k,k] / sum(weight[idx_k,k])
  }
  return(theta0)
}

#' updating alpha
#' 
#' @param X a T-by-p matrix of covariates, where `X[[t]]` being the p-vector of independent variable at time t
#' @param weight N_bin by K matrix. weight[i, k] represents weight of ith residual for kth group
#' @param lambda_alpha penalty parameter for alpha
#' @param time_indicator length N_bin factor indicating the time of the data points 
#' @return alpha (p+1)-by-K matrix. The coefficients for the cluster probabilities.
#' @export
Mstep_alpha = function(X,
                       weight,
                       lambda_alpha,
                       time_indicator){
  lambda_max = lambda_alpha * 100
  lambdas = exp(seq(from = log(lambda_max), to = log(lambda_alpha), length = 30))
  weight.sum = apply(weight, 2, tapply, time_indicator, sum)
  fit = glmnet::glmnet(x = X,
                        y = weight.sum,
                        lambda = lambdas,
                        family = "multinomial",
                        intercept = TRUE)
  coefs = glmnet::coef.glmnet(fit, s = lambda_alpha)
  alpha = as.matrix(do.call(cbind, coefs))
  return(alpha)  # (p+1) by K matrix
}

#' calculating residuals for each k
#' 
#' @param X_bin N_bin-by-p matrix, which is expanded version of X so that its dimension agrees with Y_bin's 
#' @param Y_bin N_bin-by-d matrix, indicating the center of the bin
#' @param theta0 length K list, with `theta0[[k]]` being the estimate for the intercept coefficient
#' of the regression for kth group
#' @param theta length K list, with `theta[[k]]` being the p-by-1 vector. Initial estimate for coefficients
#' of the regression for kth group
#' @param K number of clusters
#' @return length K list, with `resi0[[k]]` N_bin-by-d residual matrix
#' @export
calc_resi = function(X_bin,
                     Y_bin,
                     theta0, 
                     theta,
                     K){
  resi = list()
  for (k in 1:K){
    resi[[k]] = Y_bin - rep(theta0[[k]], dim(X_bin)[1]) - X_bin %*% matrix(theta[[k]])
  }
  return(resi) # length K list, with `theta0[[k]]` N_bin-by-d residual matrix
}

#' updating `mlelcd` g's
#' 
#' @param resi length K list, with `resi[[k]]` N_bin-by-d residual matrix
#' @param weight N_bin by K matrix. weight[i, k] represents weight of ith residual for kth group
#' @param idx N_bin by K logical matrix. idx[i, k] represents whether the corresponding weight is above r_bar
#' @param K number of clusters
#' @return g length K list, with `g[[k]]` being `mlelcd` for the kth group
#' @export 
Mstep_g = function(resi, 
                   weight, 
                   idx,
                   K){
  g = list()
  for (k in 1:K){
    idx_k = idx[,k]
    suppressWarnings({g[[k]] = LogConcDEAD::mlelcd(resi[[k]][idx_k,], w = weight[idx_k,k])})
  }
  return(g)
}

#' calculating \pi_{tk}(\alpha) = P(Z_i^{(t)} = k| X^(t))
#' 
#' @param X_bin N_bin-by-p matrix, which is expanded version of X so that its dimension agrees with Y_bin's 
#' @param alpha (p+1)-by-K matrix. The coefficients for the cluster probabilities.
#' @return pi_k N_bin-by-K matrix with `pi_k[i,k]` indicates P(Z_i = k| X^(t))
#' @export
pi_k = function(X_bin, alpha){
  p = dim(X_bin)[2]
  tmp = exp(alpha[1,] + X_bin %*% alpha[2:(p+1),])
  pi_k = tmp / rowSums(tmp) # N_bin by K matrix
  return(pi_k) 
}

#' calculating the surrogate loglikelihood Q
#' 
#' @param X_bin N_bin-by-p matrix, which is expanded version of X so that its dimension agrees with Y_bin's 
#' @param g length K list, with `g[[k]]` being `mlelcd` for the kth group
#' @param resi length K list, with `resi[[k]]` N_bin-by-d residual matrix
#' @param theta length K list, with `theta[[k]]` being the p-by-1 vector. Initial estimate for coefficients
#' of the regression for kth group
#' @param alpha (p+1)-by-K matrix. The coefficients for the cluster probabilities.
#' @param weight N_bin by K matrix. weight[i, k] represents weight of ith residual for kth group
#' @param idx N_bin by K logical matrix. idx[i, k] represents whether the corresponding weight is above r_bar
#' @param lambda_alpha penalty parameter for alpha
#' @param lambda_theta penalty parameter for theta
#' @param K number of clusters
#' @param N the total number of observations of Y_i^(t) in Y 
#' @return Q the surrogate loglikelihood of the current parameters
#' @export
# calculating the surrogate
calc_surr = function(X_bin,
                     g, 
                     resi, 
                     theta, 
                     alpha, 
                     weight, 
                     idx,
                     lambda_alpha,
                     lambda_theta,
                     K,
                     N){
  P = pi_k(X_bin, alpha) # N_bin by K matrix
  p = dim(X_bin)[2]
  sum = c()
  for (k in 1:K){
    tmp = weight[,k] * (LogConcDEAD::dlcd(resi[[k]], g[[k]], uselog = T) + log(P[,k]))
    tmp[!idx[,k]] = 0
    sum = cbind(sum, tmp) # N_bin by K matrix
  }
  Q = mean(rowSums(sum)) - lambda_alpha * sum(abs(alpha[2:(p+1),])) / N - lambda_theta * sum(abs(unlist(theta))) / N
  return(Q)
}

#' Updating responsibility
#' 
#' @param X_bin N_bin-by-p matrix, which is expanded version of X so that its dimension agrees with Y_bin's 
#' @param Y_count N_bin-vector indicating how many points belong to each bin
#' @param resi length K list, with `resi[[k]]` N_bin-by-d residual matrix
#' @param alpha (p+1)-by-K matrix. The coefficients for the cluster probabilities.
#' @param g length K list, with `g[[k]]` being `mlelcd` for the kth group
#' @param K number of clusters
#' @param r_bar threshold for responsibility
#' @return resp N_bin by K matrix. `resp[i, k]` being the current estimate for P(Z_i = k|X, Y) 
#' @return weight N_bin by K matrix. weight[i, k] represents weight of ith residual for kth group
#' @return idx N_bin by K logical matrix. idx[i, k] represents whether the corresponding weight is above r_bar
#' @export
E_step = function(X_bin,
                  Y_count,
                  resi, 
                  alpha, 
                  g,
                  K,
                  r_bar){
  P = pi_k(X_bin, alpha)
  sum = c()
  for (k in 1:K){
    tmp = LogConcDEAD::dlcd(resi[[k]], g[[k]]) * P[,k]
    sum = cbind(sum, tmp)
  }
  resp = sum/rowSums(sum)
  
  if (sum(is.nan(resp)) > 0) {
  print("There are some NaN's in the new responsibilties")
  }
  
  weight = resp * Y_count
  idx = weight > r_bar
  
  return(list(resp = resp,
              weight = weight,
              idx = idx
              ))
}

#' Updating theta
#' 
#' @param X_bin N_bin-by-p matrix, which is expanded version of X so that its dimension agrees with Y_bin's 
#' @param X a T-by-p matrix of covariates, where `X[[t]]` being the p-vector of independent variable at time t
#' @param Y_bin N_bin-by-d matrix, indicating the center of the bin
#' @param g length K list, with `g[[k]]` being `mlelcd` for the kth group
#' @param weight N_bin by K matrix. weight[i, k] represents weight of ith residual for kth group
#' @param idx_old N_bin by K logical matrix, but used in the previous iteration. 
#' idx[i, k] represents whether the corresponding previous weight is above r_bar
#' @param theta0 length K list, with `theta0[[k]]` being the estimate for the intercept coefficient
#' of the regression for kth group
#' @param theta length K list, with `theta[[k]]` being the p-by-1 vector. Initial estimate for coefficients
#' of the regression for kth group
#' @param lambda_theta penalty parameter for theta
#' @param K number of clusters
#' @param time_indicator length N_bin factor indicating the time of the data points 
#' @return theta0 length K list, with `theta0[[k]]` being the estimate for the intercept coefficient
#' of the regression for kth group
#' @return theta length K list, with `theta[[k]]` being the p-by-1 vector of estimate for coefficients
#' of the regression for kth group
#' @export
Mstep_theta = function(X_bin,
                       X,
                       Y_bin,
                       g, 
                       weight, 
                       idx_old, 
                       theta0, 
                       theta,
                       lambda_theta,
                       K,
                       time_indicator){
  theta0_new = list()
  theta_new = list()
  for (k in 1:K){
    tmp = LP_d1(X_bin, X, Y_bin, g[[k]], weight[,k], idx_old[,k], theta0[[k]], theta[[k]], lambda_theta, time_indicator)
    theta0_new[[k]] = tmp$theta0_k
    theta_new[[k]] = tmp$theta_k
  }
  return(list('theta0' = theta0_new , 'theta' = theta_new ))
}

#' Updating theta for each k
#' 
#' @param X_bin N_bin-by-p matrix, which is expanded version of X so that its dimension agrees with Y_bin's 
#' @param X a T-by-p matrix of covariates, where `X[[t]]` being the p-vector of independent variable at time t
#' @param Y_bin N_bin-by-d matrix, indicating the center of the bin
#' @param g_k `mlelcd` for the kth group
#' @param weight_k N_bin vector of weight for kth group
#' @param idx_old_k N_bin logical vector, but used in the previous iteration. 
#' idx_old_k[i] represents whether the corresponding previous weight for kth group is above r_bar
#' @param theta0_k an intercept coefficient for the regression for the kth group
#' @param theta_k a p-vector of coefficients for the regression for the kth group
#' @param lambda_theta penalty parameter for theta
#' @param time_indicator length N_bin factor indicating the time of the data points 
#' @return theta0_k estimate for the intercept coefficient of the regression for kth group
#' @return theta_k p-vector of estimates for coefficients of the regression for kth group
#' @export
LP_d1 = function(X_bin,
                 X,
                 Y_bin,
                 g_k, 
                 weight_k, 
                 idx_old_k, 
                 theta0_k, 
                 theta_k,
                 lambda_theta,
                 time_indicator){
  weight = weight_k[idx_old_k]
  Y_idx = Y_bin[idx_old_k]
  X_idx = X_bin[idx_old_k,]
  p = dim(X_idx)[2]
  n = dim(X_idx)[1]   # the number of points in C_n
  J = length(g_k$beta) # the number of affine functions
  resi = Y_idx - rep(theta0_k, n) - X_idx %*% as.matrix(theta_k)
  
  L = min(resi)
  U = max(resi) 
  const_mat = matrix(0, nrow = J*n + 2*TT, ncol = 2*(n+p+1)) 
  const_vec = rep(0, J*n + 2*TT)
  
  # epigraph part
  for (i in 1:n) {
    for (j in 1:J) {
      const_mat[(i-1)*J + j, i] = 1
      const_mat[(i-1)*J + j, n + i] = -1
    }
    const_mat[((i-1)*J + 1):((i-1)*J + J), (2*n+1)] = g_k$b
    const_mat[((i-1)*J + 1):((i-1)*J + J), (2*n+2)] = -g_k$b
    const_mat[((i-1)*J + 1):((i-1)*J + J), (2*n+3):(2*n+p+2)] = g_k$b %*% X_idx[i,]
    const_mat[((i-1)*J + 1):((i-1)*J + J), (2*n+p+3):(2*(n+p+1))] = 
      - const_mat[((i-1)*J + 1):((i-1)*J + J), (2*n+3):(2*n+p+2)]
    const_vec[((i-1)*J + 1):((i-1)*J + J)] = Y_idx[i] * g_k$b - g_k$beta
  }
  
    # feasibility part
  const_mat[(J*n+1):(J*n+TT), 2*n+1] = rep(1, TT)
  const_mat[(J*n+1):(J*n+TT), 2*n+2] = -rep(1, TT)
  const_mat[(J*n+1):(J*n+TT), 2*(n+1)+ (1:p)] = X 
  const_mat[(J*n+1):(J*n+TT), 2*(n+1)+p+ (1:p)] = -X
  const_mat[(J*n+TT+1):(J*n+TT+TT), 2*n+(1:(2*(p+1)))] = -const_mat[(J*n+1):(J*n+TT), 2*n+(1:(2*(p+1)))]
  
  for (t in 1:TT){  
    const_vec[J*n + t] = min(Y_bin[time_indicator == t & idx_old_k])- L
    const_vec[J*n+TT + t] = U - max(Y_bin[time_indicator == t & idx_old_k])
  }
  
  
  obj_coef = c(weight, -weight, 0, 0, rep(-lambda_theta, 2*p))
  const_dir = rep("<=", J*n + 2*TT)
  
  # solving LP
  lp_res = Rsymphony::Rsymphony_solve_LP(obj = obj_coef, mat = const_mat, dir = const_dir, rhs = const_vec,  max = T)
  theta0_k = lp_res$solution[2*n+1] - lp_res$solution[2*n+2]
  theta_k = as.matrix(lp_res$solution[(2*n+3):((2*n+p+2))] - lp_res$solution[(2*n+p+3):(2*(n+p+1))])
  return(list(theta0_k = theta0_k, theta_k = theta_k)) #theta
}

#' Mixture of log-concave regression
#' 
#' @param X a T-by-p matrix of covariates, where `X[[t]]` being the p-vector of independent variable at time t
#' @param Y length T list with `Y[[t]]` being a n_t-by-d matrix
#' @param K number of clusters
#' @param number of binning
#' @param min_count_ratio min count ratio
#' @param r_bar threshold for responsibility
#' @param lambda_alpha penalty parameter for alpha
#' @param lambda_theta penalty parameter for theta
#' @param max_iter number of maximum iterations of EM to perform
#' @param iter_eta threshold for the iteration. If the increment of loglikelihood is smaller than iter_eta,
#' we terminate the iterations.
#' @return X a T-by-p matrix of covariates, where `X[[t]]` being the p-vector of independent variable at time t
#' @return Y length T list with `Y[[t]]` being a n_t-by-d matrix
#' @return Y_count N_bin-vector indicating how many points belong to each bin
#' @return N the total number of observations of Y_i^(t) in Y 
#' @return resp_init N_bin by K matrix. `resp_init[i, k]` being the initial estimate for P(Z_i = k|X, Y) 
#' @return weight_init N_bin by K matrix. weight_init[i, k] represents initial weight of ith residual for kth group
#' @return idx_init N_bin by K logical matrix. idx_init[i, k] represents whether the corresponding initial weight is above r_bar
#' @return theta0_init length K list, with `theta0_init[[k]]` being the initial estimate for the intercept coefficient of the regression for kth group
#' @return theta_init length K list, with `theta_init[[k]]` being the p-by-1 vector. 
#' Initial estimates for coefficients of the regression for kth group
#' @return alpha_init (p+1)-by-K matrix. The coefficients for the initial cluster probabilities.
#' @return resi_init length K list with `resi_init[[k]]` N_bin-by-d initial residual matrix
#' @return g_init length K list, with `g_init[[k]]` being the initial `mlelcd` for the kth group
#' @return resp N_bin by K matrix. `resp[i, k]` being the final estimate for P(Z_i = k|X, Y) 
#' @return weight N_bin by K matrix. weight[i, k] represents initial weight of ith residual for kth group
#' @return idx N_bin by K logical matrix. idx[i, k] represents whether the corresponding final weight is above r_bar
#' @return theta0 length K list, with `theta0[[k]]` being the final estimate for the intercept coefficient of the regression for kth group
#' @return theta length K list, with `theta[[k]]` being the p-by-1 vector. 
#' Final estimates for coefficients of the regression for kth group
#' @return alpha (p+1)-by-K matrix. The coefficients for the cluster probabilities.
#' @return resi length K list with `resi[[k]]` N_bin-by-d final residual matrix
#' @return g length K list, with `g[[k]]` being the final `mlelcd` for the kth group
#' @return Q a vector of the surrogate loglikelihoods of the parameters, stored at each iteration
#' @return Q_every a vector of the surrogate loglikelihoods of the parameters, stored at each step
#' @export
mixLcdReg <- function(X, 
                      Y, 
                      K, 
                      B = 0, 
                      min_count_ratio = 0,
                      r_bar, 
                      lambda_alpha, 
                      lambda_theta, 
                      max_iter = 100, 
                      iter_eta = 1e-6) {
  
  # preprocessing
  binnedY = binningY_d1(X, Y, B, min_count_ratio)
  Y_bin = binnedY$Y_bin
  Y_count = binnedY$Y_count
  X_bin = binnedY$X_bin
  time_indicator = binnedY$time_indicator
  N = binnedY$N
  p = dim(X_bin)[2]
  
  # initialization
  # note that it is for d = 1

  # use flexmix to get initial resp, theta:
  flex = flexmix::flexmix(Y_bin ~ X_bin, k = K, weights = as.vector(Y_count))
  
  # initial responsibility
  resp_init = flexmix::posterior(flex) 
  
  # thresholding after resp * Y_count
  weight_init = resp_init * Y_count
  idx_init = weight_init > r_bar
  
  # theta initialization
  theta0_init = list()
  theta_init = list()
  for (k in 1:K){
    theta0_init[[k]] = flexmix::parameters(flex, comp = k)[1]
    theta_init[[k]] = matrix(flexmix::parameters(flex, comp = k)[2:(p+1)])
  }

    ## initial alpha
    alpha_init = Mstep_alpha(X, weight_init, lambda_alpha, time_indicator)
    
    ## initial theta shift
    theta0_init = Mstep_shift(X_bin, Y_bin, theta_init, weight_init, idx_init, K)
    resi_init = calc_resi(X_bin, Y_bin, theta0_init, theta_init, K)
    
    ## initial g
    g_init = Mstep_g(resi_init, weight_init, idx_init, K)
  
  ## initial Q
  Q = calc_surr(X_bin, g_init, resi_init, theta_init, alpha_init, weight_init, idx_init, lambda_alpha, lambda_theta, K, N)
  Q_every = Q
  
  idx_old = idx_init
  resi_old = resi_init
  alpha_old = alpha_init
  theta0_old = theta0_init
  theta_old = theta_init
  g_old = g_init
  
  # iteration
  for (i in seq(max_iter)) {
    
    ## E-step
  ## E-step: update responsibilities
  Estep = E_step(X_bin, Y_count, resi_old, alpha_old, g_old, K, r_bar)
  resp_new = Estep$resp
  weight_new = Estep$weight
  idx_new = Estep$idx
    Q_every = append(Q_every, calc_surr(X_bin, g_old, resi_old, theta_old, alpha_old, weight_new, idx_new, lambda_alpha, lambda_theta, K, N))
    
    ## M-step
  ## M-step: update estimates of (alpha,theta,g)

    ### Mstep_alpha
    alpha_new = Mstep_alpha(X, weight_new, lambda_alpha, time_indicator)
    Q_every = append(Q_every, calc_surr(X_bin, g_old, resi_old, theta_old, alpha_new, weight_new, idx_new, lambda_alpha, lambda_theta, K, N))
    
    ### Mstep_theta
    M_theta = Mstep_theta(X_bin, X, Y_bin, g_old, weight_new, idx_old, theta0_old, theta_old, lambda_theta, K, time_indicator)
    theta0_new = M_theta$theta0
    theta_new = M_theta$theta
    
    ### Mstep_shift
    theta0_new = Mstep_shift(X_bin, Y_bin, theta_new, weight_new, idx_new, K)
    resi_new = calc_resi(X_bin, Y_bin, theta0_new, theta_new, K)
    Q_every = append(Q_every, calc_surr(X_bin, g_old, resi_new, theta_new, alpha_new, weight_new, idx_new, lambda_alpha, lambda_theta, K, N))
    
    ### Mstep_g
    g_new = Mstep_g(resi_new, weight_new, idx_new, K)
    Q_every = append(Q_every, calc_surr(X_bin, g_new, resi_new, theta_new, alpha_new, weight_new, idx_new, lambda_alpha, lambda_theta, K, N))
    
    ### loglikelihood
    Q = append(Q, calc_surr(X_bin, g_new, resi_new, theta_new, alpha_new, weight_new, idx_new, lambda_alpha, lambda_theta, K, N))
    
    ## termination criteria
    print(i)
    if (abs((Q[i+1]-Q[i])/(Q[i])) <= iter_eta  | i==max_iter){
      break;
    } else {
      idx_old = idx_new
      resi_old = resi_new
      alpha_old = alpha_new
      theta0_old = theta0_new
      theta_old = theta_new
      g_old = g_new
    }
  }
  print(i)
  return(list(X = X,
              X_bin = X_bin,
              Y = Y,
              Y_bin = Y_bin,
              Y_count = Y_count,
              N = N,
              resp_init = resp_init,
              weight_init = weight_init,
              idx_init = idx_init,
              theta0_init = theta0_init,
              theta_init = theta_init,
              alpha_init = alpha_init,
              resi_init = resi_init,
              g_init = g_init, 
              resp = resp_new,
              weight = weight_new,
              idx = idx_new,
              theta0 = theta0_new,
              theta = theta_new,
              alpha = alpha_new,
              resi = resi_new,
              g = g_new, 
              Q = Q,
              Q_every = Q_every))
}
