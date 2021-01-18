/// @file SHOW_comp_Model.hpp
///

#ifndef realPSD_SHOW_comp_Model_hpp
#define realPSD_SHOW_comp_Model_hpp 1

#include "utils.hpp"

namespace realPSD {

  namespace SHOW_comp {

    /// Scale-free PSD for the Simple Harmonic Oscillator + White Noise (SHOW) model.
    ///
    /// The Scale-free PSD (`UFun`) for the SHOW model is
    ///
    /// \f[
    /// U(f \mid f_0, Q, R_w) = R_w + \frac{1}{[(f/f_0)^2-1]^2 + [f/(gamma)]^2}.
    /// \f]
    ///
    /// The parameters are supplied as the vector \f$\boldsymbol{\phi} = (f_0, gamma, R_w)\f$ where gamma = f0*Q.
    template <class Type>
    class UFun {
    private:
      //typedefs
      /// Typedef equivalent to `matrix<Type>`.
      // typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXd_t;
      /// Typedef equivalent to `Ref <matrix<Type> >`.
      // typedef Eigen::Ref <MatrixXd_t> RefMatrix_t;
      /// Typedef equivalent to `const Ref <const matrix<Type> >`.
      // typedef const Eigen::Ref <const MatrixXd_t> cRefMatrix_t;
      // internal variables
      int N_; ///> problem dimensions
      matrix<Type> f2_; ///> Vector of squared frequencies.
      /// Set frequency vector.
      void set_f(cRefMatrix<Type>& f);
    public:
      /// Constructor.
      UFun(int N, cRefMatrix<Type>& f);
      // /// TMB-specific constructor.
      // UFun(int N, objective_function<Type>* obj);
      /// Evaluate the normalized PSD.
      void eval(RefMatrix<Type> U, RefMatrix<Type>& phi);
    };

    template<class Type>
    inline UFun<Type>::UFun(int N, cRefMatrix<Type>& f) {
      N_ = N;
      // f2_ = zero_matrix<Type>(N_,1);
      set_f(f);
    }

    // #undef TMB_OBJECTIVE_PTR
    // #define TMB_OBJECTIVE_PTR obj

    // template<class Type>
    // inline UFun<Type>::UFun(int N, objective_function<Type>* obj) {
    //   (UFun(N));
    // }

    // #undef TMB_OBJECTIVE_PTR
    // #define TMB_OBJECTIVE_PTR this

    template<class Type>
    inline void UFun<Type>::set_f(cRefMatrix<Type>& f) {
      N_ = f.size();
      f2_ = zero_matrix<Type>(N_,1);
      f2_ = f.cwiseProduct(f);
      return;
    }

    /// Parameters are: `phi = (f0, gamma = f0*Q, Rw = Aw/sigma^2)`.
    // psd = 1/((f^2 / f0^2 - 1)^2 + f^2/(f0*Q)^2) + Rw
    template <class Type>
    inline void UFun<Type>::eval(RefMatrix<Type> U, cRefMatrix<Type>& phi) {
      U = (f2_/(phi(0,0) * phi(0,0))).array() - Type(1.0);
      U = U.cwiseProduct(U);
      U += f2_/(phi(1,0) * phi(1,0));
      U = 1.0/U.array() + phi(2,0);
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
      // pick method
      DATA_STRING(method);
      // assign f or fbar based on different methods
      if(method == "UFun" || method.find("MLE_") == 0) {
        DATA_MATRIX(f);
        return UFun<Type>(N, f);
      } else if(method.find("LP_") == 0 || method.find("NLS_") == 0) {
        DATA_MATRIX(fbar);
        return UFun<Type>(N, fbar);
      } else {
        error("Unknown method."); 
      }
    }

    #undef TMB_OBJECTIVE_PTR
    #define TMB_OBJECTIVE_PTR this

  } // end namespace SHOW_comp

} // end namespace realPSD

#endif
