# Generated from _main.Rmd: do not edit by hand

#' Plot raw data when `d = 1`
#' 
#' @param Y length T list with `Y[[t]]` being a n_t-by-d matrix
#' 
#' @export
#' 
plot_data <- function(Y) {
  par(mfrow = c(1, 2), pty = 's', cex = 0.7)
  plot(c(0, 0), xlim = c(1, TT), ylim = range(unlist(Y)), col = 'white')
  for (t in 1:TT){
    points(rep(t, length(Y[[t]])), Y[[t]])
  }
}
