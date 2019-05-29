/// @file utils.hpp

#ifndef realPSD_utils_hpp
#define realPSD_utils_hpp 1

namespace realPSD {

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
