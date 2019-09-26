/// @file SHOWModel.hpp
///
/// **TODO:**
/// - `Ref` class for inputs.
/// - Parameter transformations.

#ifndef realPSD_SHOWModel_hpp
#define realPSD_SHOWModel_hpp 1

#include "utils.hpp"

namespace realPSD {

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
    int N_; // problem dimensions
    matrix<Type> f2_; ///> Vector of squared frequencies.
  public:
    /// Constructor.
    UFun(int N);
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

  /// Parameters are: `phi = (f0, Q, Rw = Aw/sigma^2)`.
  // psd = 1/((f^2 / f0^2 - 1)^2 + f^2/(f0*Q)^2) + Rw
  template <class Type>
  inline void UFun<Type>::eval(RefMatrix_t U, cRefMatrix_t& phi) {
    U = (f2_/(phi(0,0) * phi(0,0))).array() - Type(1.0);
    U = U.cwiseProduct(U);
    U += f2_/(phi(0,0) * phi(1,0) * phi(0,0) * phi(1,0));
    U = 1.0/U.array() + phi(2,0);
    return;
  }

}

#endif
