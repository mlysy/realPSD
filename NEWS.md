# **realPSD** 0.0.0.9000 (development version)

## Breaking changes

- `periodogram()` has different inputs, drops the Nyquist frequency, and no longer folds by default.  See documentation for details.
- `fs` argument to all fit methods is dropped.
- All fit methods have a new `scale` argument.
- All fit methods are parametrized wrt `zeta = log(sigma^2)`.
- The `nll()` methods now include the correct scale factor such that post-TMB scaling is no longer required.
