## TODO List

- [ ] Create the R functions `tsSim`, `periodogram`, and `fisherGstat` mentioned in the vignette.  To do this:

    - [ ] Port the corresponding MATLAB code from [**realSHO**](https://github.com/mlysy/realSHO) into R (`periodogram` is called `get_periodogram` there).
	
	- [ ] Use the R package [**fftw**](https://CRAN.R-project.org/package=fftw) for all FFT calculations.
	
	- [ ] Please document everything properly with **roxygen**.
	
	- [ ] Add unit tests where applicable.
	
	The usual MO for github contributions is:
	
	    1. Clone the repo locally.
		2. Create your own branch, e.g., `ferris-devel`.
		3. Make modifications on your branch.
		4. Push your branch to remote and create PR of `ferris-devel` to `master`.
