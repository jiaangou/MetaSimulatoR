#' @title Weighted frequencies

#' @description Weights neighborhood densities according to an exponential function

#' @param N densities of neighborhood cells
#' @param delta coefficient determining strength of frequency-dependence (higher delta results in higher weights in high densities)

#' @return returns the weights of each cell
#' @export



weighted_freq <- function(N, delta){

  w <- exp(delta*(N+1))

  w_freq <- w/sum(w)

  return(w_freq)
}

