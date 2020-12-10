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
    /// Setter for `ZLU`.
    void set_ZLU(cRefMatrix_t& Ubar);
  public:
    /// Constructor.
    LP(int N);
    /// Setter for internal `Zbar`.
    void set_Zbar(cRefMatrix_t& Zbar);
    /// Getter for internal `zeta`.
    Type get_zeta();
    /// Optimal value of `zeta = log(sigma^2)` given `Ubar`.
    Type zeta(cRefMatrix_t& Ubar);
    /// Objective function for the LP method.
    Type nll(cRefMatrix_t& Ubar, const Type zeta);
    /// Profiled objective function for the LP method.
    Type nlp(cRefMatrix_t& Ubar);
    /// Residual vector for the LP method.
    void res(RefMatrix_t R, cRefMatrix_t& Ubar, const Type zeta);
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
    // logUbar_ = Ubar.array().log();
    // ZLU_ = Zbar_ - logUbar_;
    set_ZLU(Ubar);
    return ZLU_.sum() / N_;
  }

  /// @return Scalar value of `zeta_`.
  /// @warning Must be called after a call to `LP.nlp`.
  template <class Type>
  inline Type LP<Type>::get_zeta() {
    return zeta_;
  }

  /// The LP objective function is given by
  ///
  /// \f[
  /// Q(\bar{\boldsymbol{U}}, \zeta) = \sum_{m=1}^{N_B} (\bar Z_m - \zeta - \log \bar U_m)^2.
  /// \f]
  ///
  /// It is the negative loglikelihood of a delta-method approximation of the distribution of logs of binned periodogram values, up to a factor of `B/2`, where `B` is the bin size.
  ///
  /// @param[in] Ubar Vector of bin-averaged normalized PSDs.
  /// @param[in] zeta Log of the PSD scale factor, `zeta = log(sigma^2)`.
  ///
  /// @return Value of the objective function (scalar).
  template <class Type>
  inline Type LP<Type>::nll(cRefMatrix_t& Ubar, const Type zeta) {
    // return logUbar(0,0);
    // logUbar_ = Ubar.array().log();
    // ZLU_ = Zbar_ - logUbar_;
    set_ZLU(Ubar);
    // ZLU_.array() -= zeta_;
    // return ZLU_.squaredNorm();
    // return ((Zbar_ - logUbar).array() - zeta).matrix().squaredNorm();
    // return logUbar(0,0);
    // ZLU_ = logUbar;
    // return logUbar(0,0);
    // return ZLU_.sum();
    return nll(zeta);
  }

  /// Sets the internal value `ZLU_ = Zbar - log(Ubar)`.
  ///
  /// @param[in] Ubar Vector of bin-averaged normalized PSDs.
  template <class Type>
  inline void LP<Type>::set_ZLU(cRefMatrix_t& Ubar) {
    logUbar_ = Ubar.array().log();
    ZLU_ = Zbar_ - logUbar_;
    return;
  }

  /// Calculates the LP residual vector
  ///
  /// \f[
  /// R = \bar Z - \zeta - \log \bar U.
  /// \f]
  ///
  /// @param[out] R Vector of residuals.
  /// @param[in] Ubar Vector of bin-averaged normalized PSDs.
  /// @param[in] zeta Log of the PSD scale factor, `zeta = log(sigma^2)`.
  template <class Type>
  inline void LP<Type>::res(RefMatrix_t R, cRefMatrix_t& Ubar,
			    const Type zeta) {
    // zeta_ = zeta(Ubar);
    set_ZLU(Ubar);
    R = ZLU_.array() - zeta;
    // logUbar_ = Ubar.array().log();
    // R = Zbar_ - logUbar_;
    return;
  }

  /// @param[in] zeta Log of the PSD scale factor, `zeta = log(sigma^2)`.
  /// 
  /// @return Value of the objective function (scalar).
  template <class Type>
  inline Type LP<Type>::nll(const Type zeta) {
    ZLU_.array() -= zeta;
    return ZLU_.squaredNorm();
  }

  /// Calculates the objective function given `Ubar` at the optimal value of `zeta(Ubar)`.
  ///
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


}
  
#endif

