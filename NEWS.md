# **realPSD** 0.0.0.9000 (development version)

## Breaking changes

- `periodogram()` has different inputs, drops the Nyquist frequency, and no longer folds by default.  See documentation for details.
- `fs` argument to all fit methods is dropped.
- In `show_psd()`, `showf_psd()`, etc., replace arguments `fseq` by `freq`.
- All TMB fit methods have a new `scale` argument.
- All fit methods are parametrized wrt `zeta = log(sigma^2)`.
- The `nll()` methods now include the correct scale factor such that post-TMB scaling is no longer required.
- Simplified arguments to `make_psd_model()`.  Have not updated vignette correspondingly...

## Major updates

- Added a base class `psd_model` which performs all the `TMB::MakeADFun()` calls under the hood.
