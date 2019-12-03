## TODO List

- [ ] Settle on some naming conventions.

	- `tseq`: Vector of time points.
	- `dt`: Interobservation time.
	- `Xt`: Vector of observations in the time domain.
	- `fseq`: Vector of frequencies.
	- `fs`: Sampling frequency.
	- `Xf`: Vector of observations in the frequency domain.
	- `Yf`: Vector of periodogram ordinates in the frequency domain, i.e., `Yf = abs(Xf)^2/N`.
	- `psd_c`: PSD of continuous-time signal.
	- `psd_fs`: PSD of discrete-time signal.

- [x] Create the R functions `tsSim`, `periodogram`, and `fisherGstat` mentioned in the vignette.  To do this:

    - [x] Port the corresponding MATLAB code from [**realSHO**](https://github.com/mlysy/realSHO) into R (`periodogram` is called `get_periodogram` there).
	
	- [x] Use the R package [**fftw**](https://CRAN.R-project.org/package=fftw) for all FFT calculations.
	
	- [x] Please document everything properly with **roxygen**.
	
	- [x] Add unit tests where applicable.
	
	The usual MO for github contributions is:
	
    1. Clone the repo locally.
    2. Create your own branch, e.g., `ferris-devel`.
    3. Make modifications on your branch.
    4. Push your branch to remote and create PR of `ferris-devel` to `master`.
	
- [x] Unit test the new `res` methods for `LP_nlp` and `NLS_nlp`.

- [x] Refactor the `res` methods.  Currently, they are part of a `SIMULATE` branch which also calculates the objective function itself.  If only residuals are desired, then we shouldn't compute the objective function as well.  The `res` method of `LP_nlp` is particularly inefficient, as to compute `zeta` currently the residuals get calculated twice (logs and all).

    The relevant methods are now called `LP_res` and `NLS_res`.

- [x] Add new models to the package:

	- [x] `SHOW_nat`: The "natural" parametrization, i.e., that which Bryan original used.
	- [ ] `SHOW_comp`: The "computational" parametrizatoin, i.e., the one we originally used for **realPSD** (probably still the one documented in vignette).
	- [x] `SHOW_log`: The natural parametrization, but with each parameter input on the log scale.  This is probably the one which will yield the best results.

	Now that `SHOWFit.hpp` has been properly refactored, adding new models can be done with minimal copy-pasting.  We'll probably need to refactor the unit tests though, lest we want to copy-paste like crazy.

	A more permanent solution that will allow users to readily create their own models can be done via templating, e.g., R package  [**whisker**](https://CRAN.R-project.org/package=whisker).  This is exactly the approach used under the hood by [**usethis**](https://CRAN.R-project.org/package=usethis).
	
- [ ] Fix `LP_nll` method.  Right now, we don't pass in `C_B` constant, which means that if the Hessian of `LP_nll` is calculated it won't give the right standard errors.  One way around this is to include `C_B` in the `fs` input.  However, this means we have to wastefully calculate `log(exp(C_B)`.  A more effiicent alternative is perhaps for the `LP` methods to accept argument `logUBar`, rather than calculate it internally.

	On second thought, adding to the `fs` argument is probably easiest.  Besides, we're only talking about a single extra `exp`, and roundoff error is negligible.  To avoid confusion, should rename the argment, or add an extra argument?
	
- [x] `_nlp` methods should optionally `REPORT` `tau`/`zeta`.  It's a little extra work if we don't actually want `_nlp`, but creating a whole new computational graph just to get `tau`/`zeta` alone seems like a lot of extra work.

- [ ] Add unit tests for gradients.  Currently only `fn()` method is checked, but with `TMB_OBJECTIVE_PTR` getting passed around so many times it would be a nice sanity check.

- [ ] Prune git history of large objects -- repo is almost 30Mb and there's no way that's all source code!  Somewhere along the way we must have committed e.g., object files, PDFs, etc.

	This is a really annoying task, somewhat simplified using [this method](https://rtyley.github.io/bfg-repo-cleaner/).

