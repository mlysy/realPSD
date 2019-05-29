/// @file SHOWModel.hpp
///
/// **TODO:**
/// - `Ref` class for inputs.
/// - Parameter transformations.

#include "utils.hpp"

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
  int n_; // problem dimensions
  matrix<Type> f2_; ///> Vector of squared frequencies.
public:
  /// Constructor.
  UFun(int n);
  /// Set frequency vector.
  void set_f(cRefMatrix_t& f);
  /// Evaluate the normalized PSD.
  void eval(RefMatrix_t U, cRefMatrix_t& phi);
};

template<class Type>
inline UFun<Type>::UFun(int n) {
  n_ = n;
  f2_ = zero_matrix<Type>(n_,1);
}

template<class Type>
inline void UFun<Type>::set_f(cRefMatrix_t& f) {
  n_ = f.size();
  f2_ = zero_matrix<Type>(n_,1);
  f2_ = f.cwiseProduct(f);
  return;
}

/// Parameters are: `phi = (f0, gamma = f0/Q, Rw = Aw/sigma^2)`.
template <class Type>
inline void UFun<Type>::eval(RefMatrix_t U, cRefMatrix_t& phi) {
  U = (f2_/(phi(0,0) * phi(0,0))).array() - Type(1.0);
  U = U.cwiseProduct(U);
  U += f2_/(phi(1,0) * phi(1,0));
  U = 1.0/U.array() + phi(2,0);
  return;
}
