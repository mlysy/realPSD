## TODO List

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
	
- [ ] Unit test the new `res` methods for `LP_nlp`, `LP_nll`, and `NLS_nll`.

- [ ] Refactor the `res` methods.  Currently, they are part of a `SIMULATE` branch which also calculates the objective function itself.  If only residuals are desired, then we shouldn't compute the objective function as well.  

- [ ] The `res` method of `LP_nlp` is particularly inefficient, as to compute `zeta` currently the residuals get calculated twice (logs and all).

- [ ] Add new models to the package:

	- [ ] `SHOW_nat`: The "natural" parametrization, i.e., that which Bryan original used.
	- [ ] `SHOW_comp`: The "computational" parametrizatoin, i.e., the one we originally used for **realPSD** (probably still the one documented in vignette).
	- [ ] `SHOW_log`: The natural parametrization, but with each parameter input on the log scale.  This is probably the one which will yield the best results.

	Now that `SHOWFit.hpp` has been properly refactored, adding new models can be done with minimal copy-pasting.  We'll probably need to refactor the unit tests though, lest we want to copy-paste like crazy.

	A more permanent solution that will allow users to readily create their own models can be done via templating, e.g., R package  [**whisker**](https://CRAN.R-project.org/package=whisker).  This is exactly the approach used under the hood by [**usethis**](https://CRAN.R-project.org/package=usethis).
	

- [ ] Prune git history of large objects -- repo is almost 30Mb and there's no way that's all source code!  Somewhere along the way we must have committed e.g., object files, PDFs, etc.

	This is a really annoying task, somewhat simplified using [this method](https://rtyley.github.io/bfg-repo-cleaner/).

