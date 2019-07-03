#' @title Generates an instance of the SHO function as well as a math function funcSHO. 
#' @param f
#' @param f0
#' @param Q
#' @param k
#' @return ...
#' @export
getSHO <- function(f, f0, Q, k) {
  C = funcSHO(c(f0, Q, k), f); 
  abs_C = abs(C); 
  theta_C = unwrap(Arg(C))*180/pi;
  return(list(C = C, abs_C = abs_C, theta_C = theta_C))
}

# helper functions
#' @export
funcSHO <- function(c, x) {
  1/c[3] * 1/ ( (1-(x/c[1])^2) + 1i*x/c[1]/c[2])
} 
#' @export
funcSHO_abs <- function(c,x) {
  1/c[3] * 1/ sqrt((abs((1-(x/c[1])^2))^2 + abs((x/c[1]/c[2]))^2))
}
