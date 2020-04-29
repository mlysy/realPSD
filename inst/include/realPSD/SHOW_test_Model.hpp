/// @file SHOW_test_Model.hpp
///
/// Test a new way of passing R obj to C++

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
    /// The parameters are supplied as the vector \f$\boldsymbol{\phi} = (f_0, Q, R_w)\f$.
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
      Type mult_factor_; /// a uesless obj for testing purposes only
    public:
      /// Constructor.
      UFun(int N, Type mult_factor);
      /// TMB-specific constructor.
      UFun(int N, objective_function<Type>* obj);
      /// Set frequency vector.
      void set_f(cRefMatrix_t& f);
      /// Evaluate the normalized PSD.
      void eval(RefMatrix_t U, cRefMatrix_t& phi);
    };

    template<class Type>
    inline UFun<Type>::UFun(int N, Type mult_factor) {
      N_ = N;
      f2_ = zero_matrix<Type>(N_,1);
      mult_factor_ = mult_factor;
    }

    #undef TMB_OBJECTIVE_PTR
    #define TMB_OBJECTIVE_PTR obj

    template<class Type>
    inline UFun<Type>::UFun(int N, objective_function<Type>* obj) {
      DATA_SCALAR(mult_factor);
      (UFun(N, mult_factor));
    }

    #undef TMB_OBJECTIVE_PTR
    #define TMB_OBJECTIVE_PTR this

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
      U = (f2_/(phi(0,0) * phi(0,0))).array() - Type(1.0);
      U = U.cwiseProduct(U);
      U += f2_/(phi(0,0) * phi(1,0) * phi(0,0) * phi(1,0));
      U = 1.0/U.array() + phi(2,0);
      // U = U * mult_factor_; // would lead to errors in the unit test
      U = U.array() * mult_factor_;
      return;
    }

  } // end namespace SHOW_test

} // end namespace realPSD

#endif
