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
    matrix<Type> logUbar_;
    Type zeta_;
    /// allocate internal memory
    void init(int N);
    /// Objective function with internal `ZLU`.
    Type nll(const Type zeta);
  public:
    /// Constructor.
    LP(int N);
    /// Setter for `Zbar`.
    void set_Zbar(cRefMatrix_t& Zbar);
    /// Optimal value of `zeta = log(sigma^2)` given `phi`.
    Type zeta(cRefMatrix_t& Ubar);
    /// Objective function for the LP method.
    Type nll(cRefMatrix_t& Ubar, const Type zeta);
    /// Profiled objective function for the LP method.
    Type nlp(cRefMatrix_t& Ubar);
  };

  /// @param[in] N Length of `Zbar`.
  template <class Type>
  inline LP<Type>::LP(int N) {
    init(N);
  }

  /// Initializes `Zbar_` and `ZLU_` as one-column matrix of size `N x 1`.
  ///
  /// @param[in] N Length of `Zbar`.
  template <class Type>
  inline void LP<Type>::init(int N) {
    N_ = N;
    Zbar_ = zero_matrix<Type>(N_, 1);
    ZLU_ = zero_matrix<Type>(N_, 1);
    return;
  }
  
  /// Resets the internal value of `Zbar`.  Optionally reallocates memory if `Zbar.size() != Zbar_size()`.
  ///
  /// @param[in] Zbar Vector of log of periodogram bin averages, potentially including the bias term `C_B`.
  template <class Type>
  inline void LP<Type>::set_Zbar(cRefMatrix_t& Zbar) {
    if(Zbar.size() != N_) init(Zbar.size());
    Zbar_ = Zbar;
    return;
  }
  
  /// @param[in] Ubar Vector of bin-averaged normalized PSDs.
  ///
  /// @return Scalar estimate of `zeta`.
  template <class Type>
  inline Type LP<Type>::zeta(cRefMatrix_t& Ubar) {
    logUbar_ = Ubar.array().log();
    ZLU_ = Zbar_ - logUbar_;
    return ZLU_.sum() / N_;
  }

  /// @param[in] Ubar Vector of bin-averaged normalized PSDs.
  ///
  /// @return Value of the objective function (scalar).
  template <class Type>
  inline Type LP<Type>::nlp(cRefMatrix_t& Ubar) {
    // logUbar_ = Ubar.array().log();
    zeta_ = zeta(Ubar);
    return nll(zeta_);
    // ZLU_.array() -= zeta_;
    // return ZLU_.squaredNorm();
  }

  /// The LP objective function is the negative loglikelihood of logs of binned periodograms, times `B/2`, i.e., half the bin size.
  ///
  /// @param[in] Ubar Vector of bin-averaged normalized PSDs.
  /// @param[in] zeta Log of the PSD scale factor, `zeta = log(sigma^2)`.
  /// @return Value of the objective function (scalar).
  template <class Type>
  inline Type LP<Type>::nll(cRefMatrix_t& Ubar, const Type zeta) {
    // return logUbar(0,0);
    logUbar_ = Ubar.array().log();
    ZLU_ = Zbar_ - logUbar_;
    // ZLU_.array() -= zeta_;
    // return ZLU_.squaredNorm();
    // return ((Zbar_ - logUbar).array() - zeta).matrix().squaredNorm();
    // return logUbar(0,0);
    // ZLU_ = logUbar;
    // return logUbar(0,0);
    // return ZLU_.sum();
    return nll(zeta);
  }

  /// @param[in] zeta Log of the PSD scale factor, `zeta = log(sigma^2)`.
  /// @return Value of the objective function (scalar).
  template <class Type>
  inline Type LP<Type>::nll(const Type zeta) {
    ZLU_.array() -= zeta;
    return ZLU_.squaredNorm();
  }


}
  
#endif

