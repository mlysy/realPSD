/// @file SHOWF_log_Model.hpp
///
/// @brief Scale-free PSD for the Simple Harmonic Oscillator + White Noise (SHOW) + 1/f Noise Model.
///
/// The Scale-free PSD (`UFun`) for the SHOWF model is
///
/// \f[
/// U(f \mid f_0, Q, R_w) = R_w + R_f/f^alpha + \frac{1}{[(f/f_0)^2-1]^2 + [f/(f_0 Q)]^2}.
/// \f]
///
/// The log parameterization is given by the vector \f$\boldsymbol{\phi} = (\log f_0, \log Q, \log R_w, \log R_f, \log alpha)$\f
/// where R_w = A_w/tau, R_f = A_f/tau, tau = sigma^2.

#ifndef realPSD_SHOWF_log_Model_hpp
#define realPSD_SHOWF_log_Model_hpp 1

#include "utils.hpp"

namespace realPSD {

  namespace SHOWF_log {
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
    public:
      /// Constructor.
      UFun(int N);
      // /// TMB-specific constructor.
      // UFun(int N, objective_function<Type>* obj);
      /// Set frequency vector.
      void set_f(cRefMatrix_t& f);
      /// Evaluate the normalized PSD.
      void eval(RefMatrix_t U, cRefMatrix_t& phi);
    };

    template<class Type>
    inline UFun<Type>::UFun(int N) {
      N_ = N;
      f2_ = zero_matrix<Type>(N_,1);
    }

    template<class Type>
    inline void UFun<Type>::set_f(cRefMatrix_t& f) {
      N_ = f.size();
      f2_ = zero_matrix<Type>(N_,1);
      f2_ = f.cwiseProduct(f);
      return;
    }

    /// Parameters are: `phi = (log f0, log Q, log Rw, log Rf, log alpha)`.
    // psd = 1/((f^2 / f0^2 - 1)^2 + f^2/(f0*Q)^2) + Rw + Rf/f^alpha
    template <class Type>
    inline void UFun<Type>::eval(RefMatrix_t U, cRefMatrix_t& phi) {
      U = (f2_/exp(Type(2.0) * phi(0,0))).array() - Type(1.0);
      U = U.cwiseProduct(U);
      U += f2_/exp(Type(2.0) * (phi(0,0) + phi(1,0)));
      U = 1.0/U.array() + exp(phi(2,0)) + exp(phi(3,0))/f2_.array().pow(exp(phi(4,0)));
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
      UFun<Type> Ufun(N);
      return Ufun;
    }

    #undef TMB_OBJECTIVE_PTR
    #define TMB_OBJECTIVE_PTR this

  } // end namespace SHOWF_log

} // end namespace realPSD

#endif
