# TODO List

# Current

- [ ] `tmb_model_` does not need to be a mandatory member of `psd_model`.  It is really just for compatibility with models provided by packages built with **TMBtools**.  Could achieve the same outcome (telling `TMB::MakeADFun()` which **TMBtools** `model` to look for) using the `ctor_args` of the derived class constructor.

- [ ] The `set_f()` member of `{PSD}_Model.hpp` currently seems to be redundant.  That is, at no point in the current implementation do we need to "reset" the frequency basis, so users are currently duplicating code to set frequencies in `set_f()` and the constructor.

	Should probably wait for some time before making a breaking change to the interface of the package.
	
	But one thing we can do now is call `set_f()` from within the constructor in all provided models and examples.  At least this way there's no code duplication, i.e., we're showing users best programming practices.

- [ ] `sim_time()` should support thinning and truncating of the generated time series.  That is, if you want `N` observations at sampling frequency `fs`, generate `N_fac * N` observations at sampling frequency `fs_fac * fs`, then thin/truncate appropriately.  The larger `N_fac` and `fs_fac`, the more accurate the approximation invoked by `sim_time()`.

- [ ] Add examples for R functions.  Do this with `@example`, not `@examples`.  See [here](https://r-pkgs.org/man.html) for details.

- [ ] Add CARFIMA model to package.

# Depreciated

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

**Update:** Turns out the method above doesn't work because the **TMB** macros `DATA_MATRIX()`, etc. don't work inside class member functions.  Instead we define an external constructor as follows:

```c
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
UFun<Type> MyModel_ctor(objective_function<Type>* obj) {
  DATA_VECTOR(arg1);
  DATA_SCALAR(arg2); // etc
  return UFun<Type>(N, arg1, arg2);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
```

This works fine, except there's a potential copy of the `UFun<Type>` object constructed inside the function when it is returned.

### General Paradigm

- User writes `MyModel` class which defines two types of public members:

	1.  Generic public members which have the same name and signature for each user defined class `MyModelA`, `MyModelB`, etc.
	
	2.  A class constructor which has model-specific input parameters.

- Package `GenPkg` contains generic code which can operate on an instantiated object `my_model` of type `MyModel<Type>`.

- Package automatically creates a **TMB** file out of user + package code, i.e., which can be compiled on-the-fly or added to the user's R/**TMB** package created with **TMBtools**.

Here's what different pieces might look like:

```c
/// @file {{{Model}}}_Generics.hpp
/// @brief Generic code which supplies everything needed to the **TMB** compiler for `{{{Model}}}`.

#ifndef {{{Model}}}_Generics_hpp
#define {{{Model}}}_Generics_hpp 1

#include "{{{Header}}}.hpp" // model class definition
#include "GenPkg.hpp" // the generic package code

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type {{{Model}}}_Generics(objective_function<Type>* obj) {
  return GenericMethods<Type, {{{Class}}}<Type> >(obj, {{{Ctor}}}<Type>);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
```

The replacements here are:

- `Model`: The name of the model to use.
- `Header`: The name of the header file which declares it.
- `Class`: The name of the class and possibly the namespace containing it.
- `Ctor`: The name of the external constructor and possibly the namespace containing it.

Upon replacing these terms, the file `MyModel_Generics.hpp` can be placed in the `src/TMB` folder of the user's package and work as expected.  Alternatively, if we change the extension to `cpp` and append the following lines to the bottom:

```
template<class Type>
Type objective_function<Type>::operator() () {
  return MyModel_Generics<Type>(this);
}
```

Then we can compile `MyModel_Generics.cpp` on-the-fly within an R session.  (Note: We also need to add `#include <TMB.hpp>` at the top of the file and remove include guards to switch to on-the-fly mode.)

```c
/// @file GenericMethods.hpp
/// @brief Metafile to include all generic TMB functionality that gets attached to `MyModel`.

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type, class MyModel>
Type GenericMethods(objective_function<Type>* obj,
                    MyModel (*make_MyModel)(objective_function<Type>*)) {
  // construct MyModel object
  MyModel my_model = make_MyModel(obj);
  // do other things
  return ...;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
```

Note that `MyModel` and `make_MyModel` are called without the template `Type` argument, but it must be supplied when calling the function, i.e.,

```c
val = GenericMethods<Type, MyNS::MyModel<Type> >(obj, MyNS::make_MyModel<Type>);
```
It might be possible to use `template<class> class MyModel` as the second template parameter, in which case the `<Type>` should be added to `MyModel` inside the function and removed from its call above.

So, please write the following function:

```r
#' Create a generic TMB file from user's model definition.
#'
#' @param model Name of model.
#' @param header Name of header file including extension.  If missing defaults to `{model}.hpp`.
#' @param class Name of class definition, with enclosing namespace if it exists.  If missing defaults to `{model}`.
#' @param ctor Name of external constructor, with enclosing namespace if it exists.  If missing defaults to `make_{model}`.
#' @param standalone If `TRUE` creates a standalone `cpp` file to pass to [TMB::compile()].  Otherwise, creates an `hpp` header file to be placed in a package created with **TMBtools**.
#' @param path Name of output file.  If missing defaults to `{model}_Generics.{cpp/hpp}` depending on the value of `standalone`.  If `NULL` doesn't create file but prints contents to console.  Generates an error if `path` file already exists.
#'
#' @return The name of the file (invisibly), or nothing but prints the output to the console.
make_psd_model <- funtion(name, include, class, ctor, 
                          standalone = TRUE, path) {}
```

As for the template file `Model_Generics.hpp`, it should be placed in `inst/include/realPSD`.

