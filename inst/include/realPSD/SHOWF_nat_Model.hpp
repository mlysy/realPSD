/// @file SHOWF_nat_Model.hpp
///
/// @brief Scale-free PSD for the Simple Harmonic Oscillator + White Noise (SHOW) + 1/f Noise Model.
///
/// The Scale-free PSD (`UFun`) for the SHOWF model is
///
/// \f[
/// U(f \mid f_0, Q, R_w) = R_w + R_f/f^alpha + \frac{1}{[(f/f_0)^2-1]^2 + [f/(f_0 Q)]^2}.
/// \f]
///
/// The natural parameterization is given by the vector \f$\boldsymbol{\phi} = (f_0, Q, R_w, R_f, alpha)$\f
/// where R_w = A_w/tau, R_f = A_f/tau, tau = sigma^2.

#ifndef realPSD_SHOWF_nat_Model_hpp
#define realPSD_SHOWF_nat_Model_hpp 1

#include "utils.hpp"

namespace realPSD {

  namespace SHOWF_nat {
    template <class Type>
    class UFun {
    private:
      // //typedefs
      // /// Typedef equivalent to `matrix<Type>`.
      // typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXd_t;
      // /// Typedef equivalent to `Ref <matrix<Type> >`.
      // typedef Eigen::Ref <MatrixXd_t> RefMatrix_t;
      // /// Typedef equivalent to `const Ref <const matrix<Type> >`.
      // typedef const Eigen::Ref <const MatrixXd_t> cRefMatrix_t;
      // internal variables
      int N_; ///> problem dimensions
      matrix<Type> f2_; ///> Vector of squared frequencies.
      matrix<Type> extra_arg_; ///> Vector of extra arguments.
      /// Set frequency vector.
      void set_f(cRefMatrix<Type>& f);
    public:
      /// Constructor.
      UFun(cRefMatrix<Type>& f);
      UFun(cRefMatrix<Type>& f, cRefMatrix<Type>& extra_arg);
      // /// TMB-specific constructor.
      // UFun(int N, objective_function<Type>* obj);
      /// Evaluate the normalized PSD.
      void eval(RefMatrix<Type> U, cRefMatrix<Type>& phi);
    };

    template<class Type>
    inline UFun<Type>::UFun(cRefMatrix<Type>& f) {
      set_f(f);
    }

    template<class Type>
    inline UFun<Type>::UFun(cRefMatrix<Type>& f, cRefMatrix<Type>& extra_arg) {
      set_f(f);
      extra_arg_ = extra_arg;
    }

    template<class Type>
    inline void UFun<Type>::set_f(cRefMatrix<Type>& f) {
      N_ = f.size();
      f2_ = zero_matrix<Type>(N_,1);
      f2_ = f.cwiseProduct(f);
      return;
    }

    /// Parameters are: `phi = (f0, Q, Rw, Rf, alpha)`.
    // psd = 1/((f^2 / f0^2 - 1)^2 + f^2/(f0*Q)^2) + Rw + Rf/f^alpha
    template <class Type>
    inline void UFun<Type>::eval(RefMatrix<Type> U, cRefMatrix<Type>& phi) {
      U = (f2_/(phi(0,0) * phi(0,0))).array() - Type(1.0);
      U = U.cwiseProduct(U);
      U += f2_/(phi(0,0) * phi(1,0) * phi(0,0) * phi(1,0));
      U = 1.0/U.array() + phi(2,0) + phi(3,0)/f2_.array().pow(phi(4,0)/2.0);
      return;
    }

    #undef TMB_OBJECTIVE_PTR
    #define TMB_OBJECTIVE_PTR obj
    /// TMB-style constructor.
    ///
    /// The `PARAMETER` macro does not work properly inside class methods, 
    /// only regular functions.  Therefore, the following "external" constructor is used.
    template<class Type>
    UFun<Type> make_Ufun(cRefMatrix<Type>& f, objective_function<Type>* obj) {
      // DATA_MATRIX(extra_arg);
      // return UFun<Type>(f, extra_arg);
      return UFun<Type>(f);
    }

    #undef TMB_OBJECTIVE_PTR
    #define TMB_OBJECTIVE_PTR this

  } // end namespace SHOWF_nat

} // end namespace realPSD

#endif
