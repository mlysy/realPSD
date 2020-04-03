# carfima polynomials

Rcpp::sourceCpp("ComplexPoly.cpp")

# square norm of complex polynomial sum( alpha * (1i * x)^(0:d) )
# vectorize in argument x
sqr_poly <- function(alpha, x) {
  pol <- sapply(x, function(xx) sum( alpha * (1i * xx)^(1:length(alpha)-1) ))
  pol
  ## abs(pol)^2
}

# now the for-loop way
sqr_poly2 <- function(alpha, x) {
  deg <- length(alpha)-1 # degree
  polr <- poli <- 0 # real and imaginary parts
  sigr <- sigi <- 1 # real and imaginary signs
  for(j in 0:deg) {
    if(j %% 2 == 0) {
      polr <- polr + sigr * alpha[j+1] * x^j
      sigr <- sigr * -1
    } else {
      poli <- poli + sigi * alpha[j+1] * x^j
      sigi <- sigi * -1
    }
  }
  polr + 1i * poli
}

# now separating even and odd polynomials
sqr_poly3 <- function(alpha, x) {
  deg <- length(alpha)-1 # degree
  polr <- poli <- 0 # real and imaginary parts
  sig <- 1 # sign
  no <- floor((deg+1)/2) # number of odd terms
  me <- (deg+1) > 2*no # whether there is one more even term
  x2 <- x*x # square of argument
  j <- 0
  if(no > 0) {
    for(j in 1:no) {
      # update real
      polr <- polr + sig * alpha[2*(j-1)+1] * x2^(j-1)
      # update imaginary
      poli <- poli + sig * alpha[2*(j-1)+2] * x2^(j-1)
      # update sign
      sig <- sig * -1
    }
  }
  if(me) {
    # extra real term
    j <- j+1
    polr <- polr + sig * alpha[2*(j-1)+1] * x2^(j-1)
  }
  polr + 1i * x * poli
}

deg <- sample(c(0, 1, 2, 5, 10), 1)
n <- 10
alpha <- rnorm(deg+1)
x <- rnorm(n)
p1 <- sqr_poly(alpha, x)
p2 <- sqr_poly3(alpha, x)
## max(abs(Re(p1 - p2)))
## max(abs(Im(p1 - p2)))
max(abs(p1 - p2))

deg <- 5
n <- 1
alpha <- rnorm(deg+1)
x <- rnorm(n)
p1 <- sqr_poly(alpha, x)
p2 <- cpoly(x[1], alpha)
p1
p2
abs(p1-p2)
