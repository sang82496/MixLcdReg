# Generated from _main.Rmd: do not edit by hand

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
