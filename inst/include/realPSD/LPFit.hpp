/// @file LPFit.hpp

#ifndef realPSD_LPFit_hpp
#define realPSD_LPFit_hpp 1

#include "utils.hpp"

namespace realPSD {

  /// Class for the LP fitting method.
  template <class Type>
  class LP {
  private:
    // Typedefs
    /// Typedef equivalent to `MatrixXd<Type>`.
    typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXd_t;
    /// Typedef equivalent to `Ref <MatrixXd<Type> >`
    typedef Eigen::Ref <MatrixXd_t> RefMatrix_t;
    /// Typedef equivalent to `const Ref <const MatrixXd<Type> >`
    typedef const Eigen::Ref <const MatrixXd_t> cRefMatrix_t;
    // internal variables
    int N_;
    matrix<Type> Zbar_;
    matrix<Type> ZLU_;
    Type zeta_;
    // allocate internal memory
    void init(int N);
  public:
    /// Constructor.
    LP(int N);
    /// Setter for Zbar.
    void set_Zbar(cRefMatrix_t& Zbar);
    /// Optimal value of `zeta = log(sigma^2)` given `phi`.
    Type zeta(cRefMatrix_t& logUbar);
    /// Objective function for the LP method.
    Type nll(cRefMatrix_t& logUbar);
  };

  template <class Type>
  inline LP<Type>::LP(int N) {
    init(N);
  }

  template <class Type>
  inline void LP<Type>::init(int N) {
    N_ = N;
    Zbar_ = zero_matrix<Type>(N_, 1);
    ZLU_ = zero_matrix<Type>(N_, 1);
    return;
  }

  template <class Type>
  inline void LP<Type>::set_Zbar(cRefMatrix_t& Zbar) {
    if(Zbar.size() != N_) init(Zbar.size());
    Zbar_ = Zbar;
    return;
  }
  
  /// @param[in] Zbar Vector of log of binned periodograms, potentially including the bias term `C_B`.
  /// @param[in] logUbar Vector of log of normalized PSD at bin-average frequencies.
  ///
  /// @return Scalar estimate of `zeta`.
  template <class Type>
  inline Type LP<Type>::zeta(cRefMatrix_t& logUbar) {
    ZLU_ = Zbar_ - logUbar;
    return ZLU_.sum() / N_;
  }

  /// This corresponds to the negative loglikelihood of logs of binned periodograms, times `B/2`, i.e., half the bin size.
  ///
  /// @param[in] Zbar Vector of log of binned periodograms, potentially including the bias term `C_B`.
  /// @param[in] logUbar Vector of log of normalized PSD at bin-average frequencies.
  ///
  /// @return Scalar value of the objective function.
  template <class Type>
  inline Type LP<Type>::nll(cRefMatrix_t& logUbar) {
    zeta_ = zeta(logUbar);
    ZLU_.array() -= zeta_;
    return ZLU_.squaredNorm();
  }

}
  
#endif

