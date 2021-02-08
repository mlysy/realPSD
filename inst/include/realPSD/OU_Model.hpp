/// @file OU_Model.hpp

#include "realPSD/utils.hpp"

namespace ou {
  using namespace realPSD; // import functions from utils.hpp

  /// The Ornstein-Uhlenbeck model class.
  ///
  /// The PSD of the Ornstein-Uhlenbeck (OU) model is defined as
  ///
  /// \f[
  /// U(f \mid \alpha) = \frac{1}{(2\pi f)^2 + \alpha^2}. 
  /// \f]
  /// The parameter is supplied as \f$\phi = \log \alpha\f$.
  template <class Type>
  class UFun {
  private:
    // internal variables
    int N_; ///> problem dimensions
    matrix<Type> f_; ///>column vector of frequencies.
    matrix<Type> f2_; ///> column vector of squared frequencies.
    Type scale_; ///> Scaling factor for frequencies.
    matrix<Type> extra_arg_; // vector of additional arguments
    /// Set frequency vector.
    void set_f(cRefMatrix<Type>& f);
  public:
    /// Constructor.
    UFun(cRefMatrix<Type>& f, cRefMatrix<Type>& extra_arg);
    /// Evaluate the normalized PSD.
    void eval(RefMatrix<Type> U, cRefMatrix<Type>& phi);
  };

  /// @param[in] f Frequency vector of size `N`.
  /// @param[in] extra_arg Additional argument vector.
  template<class Type>
  inline UFun<Type>::UFun(cRefMatrix<Type>& f, cRefMatrix<Type>& extra_arg) {
    // N_ = N;
    set_f(f); // note that this will also set `N_`
    scale_ = Type(4.0 * EIGEN_PI * EIGEN_PI);
    extra_arg_ = extra_arg; // additional constructor arguments
  }

  /// @param[in] f Frequency vector of size `N`.
  template<class Type>
  inline void UFun<Type>::set_f(cRefMatrix<Type>& f) {
    N_ = f.size();
    f_ = zero_matrix<Type>(N_,1);
    f_ = f;
    f2_ = zero_matrix<Type>(N_,1);
    f2_ = f.cwiseProduct(f);
    return;
  }

  /// @param[out] U Normalized PSD vector of size `N`.
  /// @param[in] phi Computational basis parameters.  Here we have `phi = log(alpha)`.
  template <class Type>
  inline void UFun<Type>::eval(RefMatrix<Type> U, cRefMatrix<Type>& phi) {
    U = (scale_ * f2_).array() + exp(Type(2.0)*phi(0,0));
    // U = exp(Type(2.0)*phi(1,0)) / U.array();
    U = U.array().inverse();
    return;
  }

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
  /// External constructor for `ou::UFun` objects.
  ///
  /// The arguments to this function are always `f` and `obj`.  Inside the function, we can specify additional TMB macros (`DATA_VECTOR`, etc.), to obtain inputs to the `UFun` constructor.
  ///
  /// @param[in] N Number of frequency/psd observations.
  /// @param[in] obj Pointer to the TMB object.
  /// @return An `ou::UFun` object.
  template<class Type>
  UFun<Type> make_Ufun(cRefMatrix<Type>& f, objective_function<Type>* obj) {
    DATA_MATRIX(extra_arg); // get extra_arg from R
    return UFun<Type>(f, extra_arg);
  }
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
  
} // end namespace ou
