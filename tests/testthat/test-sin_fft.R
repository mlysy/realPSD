context("sin_fft")

test_that("sin_fft via geometric series is same as direct summation.", {
  ntest <- 10
  for(ii in 1:ntest) {
    fs <- runif(1)
    N <- sample(10:200, 1)
    A <- rexp(1)
    xi <- rexp(1)
    phi <- rexp(1)
    full_freq <- sample(c(TRUE, FALSE), 1)
    tseq <- (1:N-1)/fs
    Xt <- A * sin(2*pi * xi * tseq + phi)
    if(full_freq) {
      freq <- (1:N-1)/N*fs
      Xf1 <- fft(Xt) # calculate via FFT
    } else {
      freq <- rnorm(sample(10:200, 1))
      # calculate via direct summation
      Xf1 <- sapply(freq, function(f) {
        sum(exp(-2i*pi * (1:N-1)/fs * f) * A * sin(2*pi * xi * (1:N-1)/fs + phi))
      })
    }
    # calculate via geometric series
    Xf2 <- sin_fft(freq, fs, N, A, xi, phi)
    expect_equal(Xf1, Xf2)
  }
})
