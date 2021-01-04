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
  public:
    /// Constructor.
    UFun(int N);
    /// Set frequency vector.
    void set_f(cRefMatrix<Type>& f);
    /// Evaluate the normalized PSD.
    void eval(RefMatrix<Type> U, cRefMatrix<Type>& phi);
  };

  /// @param[in] N Number of frequency/PSD observations.
  template<class Type>
  inline UFun<Type>::UFun(int N) {
    N_ = N;
    f_ = zero_matrix<Type>(N_, 1);
    f2_ = zero_matrix<Type>(N_,1);
    scale_ = Type(4.0 * EIGEN_PI * EIGEN_PI);
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
  /// The arguments to this function are always `N` and `obj`.  Inside the function, we can specify additional TMB macros (`DATA_VECTOR`, etc.), to obtain inputs to the `UFun` constructor.
  ///
  /// @param[in] N Number of frequency/psd observations.
  /// @param[in] obj Pointer to the TMB object.
  /// @return An `ou::UFun` object.
  template<class Type>
  UFun<Type> make_Ufun(int N, objective_function<Type>* obj) {
    return UFun<Type>(N);
  }
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
  
} // end namespace ou
