/// @file utils.hpp

#ifndef realPSD_utils_hpp
#define realPSD_utils_hpp 1

namespace realPSD {

  /// Typedef equivalent to `Ref <MatrixXd<Type> >`.
  template <class Type>
  using RefMatrix = Eigen::Ref <Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> >;
  /// Typedef equivalent to `const Ref <const MatrixXd<Type> >`
  template <class Type>
  using cRefMatrix = const Eigen::Ref <const Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> >;
  /// Typedef equivalent to `Ref <VectorXd<Type> >`.
  template <class Type>
  using RefVector = Eigen::Ref <Eigen::Matrix<Type, Eigen::Dynamic, 1> >;
  /// Typedef equivalent to `const Ref <const VectorXd<Type> >`
  template <class Type>
  using cRefVector = const Eigen::Ref <const Eigen::Matrix<Type, Eigen::Dynamic, 1> >;

  /// Create `matrix<Type>` of zeros.
  ///
  /// @param[in] n Integer number of rows.
  /// @param[in] p Integer number of columns.
  /// @return A `matrix<Type>` of size `n x p` initialized with zeros.
  template<class Type>
  matrix<Type> zero_matrix(int n, int p) {
    matrix<Type> out(n,p);
    out.setZero();
    return out;
  }

}

#endif
