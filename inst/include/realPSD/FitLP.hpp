/// @file FitLP.hpp

/// Optimal value of `zeta = log(sigma^2)` given `phi`.
///
/// @param[in] Zbar Vector of log of binned periodograms, potentially including the bias term `C_B`.
/// @param[in] logUbar Vector of log of normalized PSD at bin-average frequencies.
///
/// @return Scalar estimate of `zeta`.
template <class Type>
Type zeta_nls(const array<Type>& Zbar, const array<Type>& logUbar) {
  return (Zbar - logUbar).sum() / Zbar.size();
}

/// Objective function for the LP method.
///
/// This corresponds to the negative loglikelihood of logs of binned periodograms, times `B/2`, i.e., half the bin size.
///
/// @param[in] Zbar Vector of log of binned periodograms, potentially including the bias term `C_B`.
/// @param[in] logUbar Vector of log of normalized PSD at bin-average frequencies.
///
/// @return Scalar value of the objective function.
template <class Type>
Type nll_lp(const array<Type>& Zbar, const array<Type>& logUbar) {
  Type zeta_hat = zeta_lp(Zbar, logUbar);
  return (Zbar - zeta_hat - logUbar).matrix().squaredNorm();
}
