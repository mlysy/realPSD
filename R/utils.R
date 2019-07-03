#' @title Unwrap radian phases by adding multiples of 2*pi as appropriate to
#' remove jumps greater than tol. tol defaults to pi.
#' @param a Radian angle.
#' @param tol Tolerance.
#' @param dim Dimension.
#' @export
unwrap <- function(a, tol = pi, dim = 1) {
  sz = dim(a)
  nd = length(sz)
  if (nd == 0) {
    sz = length(a)
    nd = 1
  }
  if (! (length(dim) == 1 && dim == round(dim)) && dim > 0 && dim < (nd + 1))
    stop("unwrap: dim must be an integer and valid dimension")
  # Find the first non-singleton dimension
  while (dim < (nd + 1) && sz[dim] == 1)
    dim = dim + 1
  if (dim > nd)
    dim = 1
  # Don't let anyone use a negative value for TOL.
  tol = abs(tol)
  rng = 2*pi
  m = sz[dim]
  # Handle case where we are trying to unwrap a scalar, or only have
  # one sample in the specified dimension.
  if (m == 1)       
    return(a)
  # Take first order difference to see so that wraps will show up
  # as large values, and the sign will show direction.
  idx = list()
  for (i  in  1:nd) 
    idx[[i]] = 1:sz[i]
  idx[[dim]] = c(1,1:(m-1))
  d = a[unlist(idx)] - a
  # Find only the peaks, and multiply them by the range so that there
  # are kronecker deltas at each wrap point multiplied by the range
  # value.
  p =  rng * (((d > tol) > 0) - ((d < -tol) > 0))
  # Now need to "integrate" this so that the deltas become steps.
  if (nd == 1)
    r = cumsum(p)
  else  
    r = apply(p, MARGIN = dim, FUN = cumsum)
  # Now add the "steps" to the original data and put output in the
  # same shape as originally.
  a + r
} 

# test P1 == P2 == P3
# Z = c(1 - 1i, 2 + 1i, 3 - 1i, 4 + 1i, 
#   1 + 2i, 2 - 2i, 3 + 2i, 4 - 2i, 
#   1 - 3i, 2 + 3i, 3 - 3i, 4 + 3i, 
#   1 + 4i, 2 - 4i, 3 + 4i, 4 - 4i)
# P1 = Arg(Z)
# P2 = atan2(Im(Z), Re(Z))
# P3 = Im(log(Z))
