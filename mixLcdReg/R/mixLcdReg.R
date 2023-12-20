# Generated from _main.Rmd: do not edit by hand

#' Mixture of log-concave regression
#' 
#' @param X a T-by-p matrix of covariates, where `X[[t]]` being the p-vector of independent variable at time t
#' @param Y length T list with `Y[[t]]` being a n_t-by-d matrix
#' @param K number of clusters
#' @param number of binning
#' @param min_count_ratio 
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
                      B = 40, 
                      min_count_ratio = 0,
                      r_bar, 
                      lambda_alpha, 
                      lambda_theta, 
                      max_iter = 100, 
                      iter_eta = 1e-6) {
  
  # preprocessing
  binnedY = binningY_d1(X, Y, B)
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
    if (Q[i+1]-Q[i] <= iter_eta  | i==max_iter){
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
