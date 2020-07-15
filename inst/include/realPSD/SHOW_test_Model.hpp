/// @file SHOW_test_Model.hpp

#ifndef realPSD_SHOW_test_Model_hpp
#define realPSD_SHOW_test_Model_hpp 1

#include "utils.hpp"

namespace realPSD {

  namespace SHOW_test {

    /// Scale-free PSD for the Simple Harmonic Oscillator + White Noise (SHOW) model.
    ///
    /// The Scale-free PSD (`UFun`) for the SHOW model is
    ///
    /// \f[
    /// U(f \mid f_0, Q, R_w) = R_w + \frac{1}{[(f/f_0)^2-1]^2 + [f/(f_0 Q)]^2}.
    /// \f]
    ///
    /// The parameters are supplied as the vector \f$\boldsymbol{\phi} = (\log f_0, \log Q, \log R_w)\f$.
    template <class Type>
    class UFun {
    private:
      //typedefs
      /// Typedef equivalent to `matrix<Type>`.
      typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXd_t;
      /// Typedef equivalent to `Ref <matrix<Type> >`.
      typedef Eigen::Ref <MatrixXd_t> RefMatrix_t;
      /// Typedef equivalent to `const Ref <const matrix<Type> >`.
      typedef const Eigen::Ref <const MatrixXd_t> cRefMatrix_t;
      // internal variables
      int N_; ///> problem dimensions
      matrix<Type> f2_; ///> Vector of squared frequencies.
      Type alpha_; ///> Multiplicative factor.
    public:
      /// Constructor.
      UFun(int N, Type alpha);
      /// TMB-specific constructor.
      // UFun(int N, objective_function<Type>* obj);
      /// Set frequency vector.
      void set_f(cRefMatrix_t& f);
      /// Evaluate the normalized PSD.
      void eval(RefMatrix_t U, cRefMatrix_t& phi);
    };

    template<class Type>
    inline UFun<Type>::UFun(int N, Type alpha) {
      N_ = N;
      f2_ = zero_matrix<Type>(N_,1);
      alpha_ = alpha;
    }

    template<class Type>
    inline void UFun<Type>::set_f(cRefMatrix_t& f) {
      N_ = f.size();
      f2_ = zero_matrix<Type>(N_,1);
      f2_ = f.cwiseProduct(f);
      return;
    }

    /// Parameters are: `phi = (f0, Q, Rw = Aw/sigma^2)`.
    // psd = 1/((f^2 / f0^2 - 1)^2 + f^2/(f0*Q)^2) + Rw
    template <class Type>
    inline void UFun<Type>::eval(RefMatrix_t U, cRefMatrix_t& phi) {
      U = (f2_/exp(Type(2.0) * phi(0,0))).array() - Type(1.0);
      U = U.cwiseProduct(U);
      U += f2_/exp(Type(2.0) * (phi(0,0) + phi(1,0)));
      U = 1.0/U.array() + exp(phi(2,0));
      U *= alpha_;
      return;
    }

    #undef TMB_OBJECTIVE_PTR
    #define TMB_OBJECTIVE_PTR obj
    /// TMB-style constructor.
    ///
    /// The `PARAMETER` macro does not work properly inside class methods, 
    /// only regular functions.  Therefore, the following "external" constructor is used.
    template<class Type>
    UFun<Type> make_Ufun(int N, objective_function<Type>* obj) {
      PARAMETER(alpha);
      UFun<Type> Ufun(N, alpha);
      return Ufun;
    }

    #undef TMB_OBJECTIVE_PTR
    #define TMB_OBJECTIVE_PTR this

  } // end namespace SHOW_test

  #undef TMB_OBJECTIVE_PTR
  #define TMB_OBJECTIVE_PTR obj
  /// TMB method dispatch.
  ///
  /// The external constructor is passed to `FitMethods` via a C++ callback, i.e., a function-pointer argument to a function.
  template<class Type, class UFun>
  Type FitMethods(objective_function<Type>* obj,
      UFun (*make_Ufun)(int, objective_function<Type>*)) {
    // data
    DATA_MATRIX(f);
    // parameters
    PARAMETER_MATRIX(phi);
    // PARAMETER_MATRIX(alpha);
    // calculate U
    int N = f.size();
    // UFun<Type> Ufun(N, alpha);
    // UFun<Type> Ufun(N, this);
    UFun Ufun = make_Ufun(N, obj);
    Ufun.set_f(f);
    SIMULATE {
      matrix<Type> U(N,1);
      Ufun.eval(U, phi);
      REPORT(U);
    }
    return Type(0.0);
  }
  #undef TMB_OBJECTIVE_PTR
  #define TMB_OBJECTIVE_PTR this

} // end namespace realPSD

#endif
