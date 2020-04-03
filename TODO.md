## TODO List

- [ ] Settle on some naming conventions.

	- `tseq`: Vector of time points.
	- `dt`: Interobservation time.
	- `Xt`: Vector of observations in the time domain.
	- `fseq`: Vector of frequencies.
	- `fs`: Sampling frequency.
	- `Xf`: Vector of observations in the frequency domain.
	- `Yf`: Vector of periodogram ordinates in the frequency domain, i.e., `Yf = abs(Xf)^2/N`.
	- `psd`: PSD of continuous-time signal.
	- `psd_fs`: PSD of discrete-time signal.
	- Let's always use `snake_case` instead of `camelCase` or `PascalCase`.  The only exception is C++ class names, which are `PascalCase`.  Most importantly, never use `cran.case` for R objects.

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
	
	- [ ] `CARFIMA`: The general CARFIMA model defined in the vignette.  
	
	Now that `SHOWFit.hpp` has been properly refactored, adding new models can be done with minimal copy-pasting.  We'll probably need to refactor the unit tests though, lest we want to copy-paste like crazy.

	A more permanent solution that will allow users to readily create their own models can be done via templating, e.g., R package  [**whisker**](https://CRAN.R-project.org/package=whisker).  This is exactly the approach used under the hood by [**usethis**](https://CRAN.R-project.org/package=usethis).  See below for a potential API.
		
- [ ] Fix `LP_nll` method.  Right now, we don't pass in `C_B` constant, which means that if the Hessian of `LP_nll` is calculated it won't give the right standard errors.  One way around this is to include `C_B` in the `fs` input.  However, this means we have to wastefully calculate `log(exp(C_B)`.  A more effiicent alternative is perhaps for the `LP` methods to accept argument `logUBar`, rather than calculate it internally.

	On second thought, adding to the `fs` argument is probably easiest.  Besides, we're only talking about a single extra `exp`, and roundoff error is negligible.  To avoid confusion, should rename the argment, or add an extra argument?
	
- [x] `_nlp` methods should optionally `REPORT` `tau`/`zeta`.  It's a little extra work if we don't actually want `_nlp`, but creating a whole new computational graph just to get `tau`/`zeta` alone seems like a lot of extra work.

- [ ] Add unit tests for gradients.  Currently only `fn()` method is checked, but with `TMB_OBJECTIVE_PTR` getting passed around so many times it would be a nice sanity check.

- [ ] Prune git history of large objects -- repo is almost 30Mb and there's no way that's all source code!  This involves two things:

	1.  Remove everything from the repo that shouldn't be part of the R package, e.g., PDF vignettes, internal R Markdown reports on the Technometrics simulation study, etc.
	
	2.  Even if we remove these objects, they remain in the Git history and therefore add a lot of memory to the repo.  Therefore, it is necessary to rewrite the repo history without these objects committed.  This is a difficult task, somewhat simplified using [this method](https://rtyley.github.io/bfg-repo-cleaner/).

## API for Adding New Models

Here's the C++ code specifying the model:
```c
namespace MyModel {

template <class Type>
class UFun {
  private:
  // Convenience typedefs for TMB/Eigen matrix inputs.
  // Might eventually make these typedefs part of realPSD namespace,
  // but this works for now.
  /// Typedef equivalent to `matrix<Type>`.
  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXd_t;
  /// Typedef equivalent to `Ref <matrix<Type> >`.
  typedef Eigen::Ref <MatrixXd_t> RefMatrix_t;
  /// Typedef equivalent to `const Ref <const matrix<Type> >`.
  typedef const Eigen::Ref <const MatrixXd_t> cRefMatrix_t;
  // more internal variables as needed.
  public:
  /// Constructor.
  UFun(int N, type1& arg1, type2& arg2, ...);
  /// TMB-specific constructor.
  UFun(int N, objective_function<Type>* obj);
  /// Set frequency vector.
  void set_f(cRefMatrix_t& f);
  /// Evaluate the normalized PSD.
  void eval(RefMatrix_t U, cRefMatrix_t& phi);  
};

}
```

The challenge with this approach is that the arguments of `UFun()` aside from `N` are model-dependent.  So TMB needs to know how to construct an object of this class.  Here's an idea:

```c
namespace MyModel {

class UFun {...}; // class definition

// must enclose constructor definition in these #statements
// for R -> C++ argument passing to work
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
inline UFun<Type>::UFun(int N, objective_function<Type>* obj) {
  DATA_VECTOR(arg1);
  DATA_SCALAR(arg2); // etc
  UFun(N, arg1, arg2);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
}
```

Then, we would modify the corresponding line of `FitMethods` to e.g.,

```c
template<class Type, class UFun>
Type FitMethods(objective_function<Type>* obj) {
  // pick method
  DATA_STRING(method);
  if(method == "UFun") {
    // data
    DATA_MATRIX(f);
    // parameters
    PARAMETER_MATRIX(phi);
    // calculate U
    int N = f.size();
    UFun Ufun(N, obj);
    Ufun.set_f(f);
    ...
  }
  ...
}
```

To check if this works:

1.  Rewrite the current TMB models `SHOW_nat` and `SHOW_log` in this form, and test to make sure it all works.  Also, please add gradient unit tests for these as suggested above.

2.  Let's try a test model  `SHOW_test::UFun()`.  The constructor will have an extra argument `double mult_factor` which simply scales `UFun.eval()` by a constant factor.  Of course this is useless, but it allows us to test the new `obj` constructor approach.

3.  Implement the CARFIMA model.  In this case, the additional constructor arguments are `p` and `q`.  Please see `vignettes/realPSD.Rmd`, `tests/dontrun/ComplexPoly.cpp` and accompanying `tests/dontrun/carfima-poly.R` for how to evaluate "PSD" polynomials efficiently using only real operations.
