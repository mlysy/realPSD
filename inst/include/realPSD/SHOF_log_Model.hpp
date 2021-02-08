/// @file SHOF_log_Model.hpp

#ifndef realPSD_SHOF_log_Model_hpp
#define realPSD_SHOF_log_Model_hpp 1

#include "utils.hpp"

namespace realPSD {

  /// Scale-free PSD for the Simple Harmonic Oscillator + 1/f Noise Model Without White Noise.
  ///
  /// The Scale-free PSD `U(phi)` for the SHOF model is
  ///
  /// \f[
  /// U(f \mid f_0, Q, R_f, alpha) = R_f/f^alpha + \frac{1}{[(f/f_0)^2-1]^2 + [f/(f_0 Q)]^2}.
  /// \f]
  ///
  /// The computational basis for the normalized PSD parameters is \f$\boldsymbol{\phi} = (\log f_0, \log Q, \log R_f, alpha)$\f.
  namespace SHOF_log {
    template <class Type>
    class UFun {
    private:
      // internal variables
      int N_; ///> problem dimensions
      matrix<Type> f2_; ///> Vector of squared frequencies.
      matrix<Type> logf_; ///> Vector of log frequencies.
      matrix<Type> extra_arg_; ///> Vector of extra arguments.
      /// Set frequency vector.
      void set_f(cRefMatrix<Type>& f);
    public:
      /// Constructor.
      UFun(cRefMatrix<Type>& f, cRefMatrix<Type>& extra_arg);
      // /// TMB-specific constructor.
      // UFun(int N, objective_function<Type>* obj);
      /// Evaluate the normalized PSD.
      void eval(RefMatrix<Type> U, cRefMatrix<Type>& phi);
    };

    template<class Type>
    inline UFun<Type>::UFun(cRefMatrix<Type>& f, cRefMatrix<Type>& extra_arg) {
      // N_ = N;
      // f2_ = zero_matrix<Type>(N_,1);
      // logf_ = zero_matrix<Type>(N_,1);
      set_f(f);
      extra_arg_ = extra_arg;
    }

    template<class Type>
    inline void UFun<Type>::set_f(cRefMatrix<Type>& f) {
      N_ = f.size();
      f2_ = zero_matrix<Type>(N_,1);
      logf_ = zero_matrix<Type>(N_,1);
      f2_ = f.cwiseProduct(f);
      logf_.array() = f.array().log();
      return;
    }

    template <class Type>
    inline void UFun<Type>::eval(RefMatrix<Type> U, cRefMatrix<Type>& phi) {
      U = (f2_/exp(Type(2.0) * phi(0,0))).array() - Type(1.0);
      U = U.cwiseProduct(U);
      U += f2_/exp(Type(2.0) * (phi(0,0) + phi(1,0)));
      U = 1.0/U.array() + (phi(2,0) - phi(3,0) * logf_.array()).exp();
      return;
    }

    #undef TMB_OBJECTIVE_PTR
    #define TMB_OBJECTIVE_PTR obj
    /// TMB-style constructor.
    ///
    /// The `DATA_*` macros don't not work properly inside class methods, 
    /// only regular functions.  Therefore, the following "external" constructor is used.
    template<class Type>
    UFun<Type> make_Ufun(cRefMatrix<Type>& f, objective_function<Type>* obj) {
      DATA_MATRIX(extra_arg);
      return UFun<Type>(f, extra_arg);
    }

    #undef TMB_OBJECTIVE_PTR
    #define TMB_OBJECTIVE_PTR this

  } // end namespace SHOF_log

} // end namespace realPSD

#endif
