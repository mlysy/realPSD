/// @file FitMLE.hpp

/// Optimal value of `tau = sigma^2` given `phi`.
///
/// @param[in] YU Vector corresponding to `Y/U(phi)`.
///
/// @return Scalar estimate of `tau`.
template <class Type>
Type tau_mle(const array<Type>& YU) {
  return YU.mean();
}

/// Objective function for the MLE method.
///
/// This corresponds to the negative Whittle loglikelihood.
///
/// @param[in] YU Vector corresponding to `Y/U(phi)`.
/// @param[in] U Corresponding vector of normalized PSD values `U(phi)`.
///
/// @return Scalar value of the objective function.
template <class Type>
Type nll_mle(const array<Type>& YU, const array<Type>& U) {
  Type tau_hat = tau_mle(YU);
  return (YU/tau_hat + U.log()).sum() + U.size() * log(tau_hat);
}
