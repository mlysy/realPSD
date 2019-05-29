/// @file FitNLS.hpp

/// Optimal value of `tau = sigma^2` given `phi`.
///
/// @param[in] Ybar Vector of binned periodograms.
/// @param[in] Ubar Vector of normalized PSD at bin-average frequencies.
///
/// @return Scalar estimate of `tau`.
template <class Type>
Type tau_nls(const array<Type>& Ybar, const array<Type>& Ubar) {
  Type tau = (Ybar * Ubar).sum();
  return tau/(Ubar * Ubar).sum();
}

/// Objective function for the NLS method.
///
/// This corresponds to unweighted sum-of-squares between the binned periodogram and the PSD at the average value.
///
/// @param[in] Ybar Vector of binned periodograms.
/// @param[in] Ubar Vector of normalized PSD at bin-average frequencies.
///
/// @return Scalar value of the objective function.
template <class Type>
Type nll_nls(const array<Type>& Ybar, const array<Type>& Ubar) {
  Type tau_hat = tau_nls(Ybar, Ubar);
  return (Ybar - tau_hat * Ubar).matrix().squaredNorm();
}
