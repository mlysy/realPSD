/// @file NLSFit.hpp

#ifndef realPSD_NLSFit_hpp
#define realPSD_NLSFit_hpp

#include "utils.hpp"

namespace realPSD {

  /// Class for the NLS fitting method.
  template <class Type>
  class NLS {
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
    matrix<Type> Ybar_;
    matrix<Type> YU_;
    Type tau_;
    Type wgt_;
    // allocate internal memory
    void init(int N);
  public:
    /// Constructor
    NLS(int N);
    /// setter for internal `Ybar`.
    void set_Ybar(cRefMatrix_t& Ybar);
    /// Optimal value of `tau = sigma^2` given `Ubar`.
    Type tau(cRefMatrix_t& Ubar);
    /// Objective function for the NLS method.
    Type nll(cRefMatrix_t& Ubar, const Type tau);
    /// Profiled objective function for the NLS method.
    Type nlp(cRefMatrix_t& Ubar);
    /// Residual vector for the NLS method.
    void res(RefMatrix_t R, cRefMatrix_t& Ubar);
  };

  /// @param[in] N Length of `Ybar`.
  template <class Type>
  inline NLS<Type>::NLS(int N) {
    init(N);
  }

  /// Initializes `Ybar_` and `YU_` as a one-column matrix of size `N x 1`.
  /// @param[in] N Length of `Ybar`.
  template <class Type>
  inline void NLS<Type>::init(int N) {
    N_ = N;
    Ybar_ = zero_matrix<Type>(N_, 1);
    YU_ = zero_matrix<Type>(N_, 1);
    return;
  }


  /// Resets the internal value of `Ybar`.  Optionally reallocates memory if `Ybar.size() != Ybar_size()`.
  ///
  /// @param[in] Ybar Vector of bin-averaged periodogram values.
  template <class Type>
  inline void NLS<Type>::set_Ybar(cRefMatrix_t& Ybar) {
    if(Ybar.size() != N_) init(Ybar.size());
    Ybar_ = Ybar;
    return;
  }
  
  
  /// Optimal value of `tau = sigma^2` given `Ubar`.
  ///
  /// @param[in] Ubar Vector of normalized PSD values at the bin-average frequencies.
  ///
  /// @return Scalar estimate of `tau`.
  template <class Type>
  inline Type NLS<Type>::tau(cRefMatrix_t& Ubar) {
    wgt_ = Ybar_.cwiseProduct(Ubar).sum();
    return wgt_ / Ubar.cwiseProduct(Ubar).sum();
  }

  /// The NLS objective function is given by
  ///
  /// \f[
  /// Q(\bar{\boldsymbol{U}}, \tau) = \sum_{m=1}^{N_B} (\bar Y_m - \tau \cdot \bar U_m)^2.
  /// \f]
  ///
  ///
  /// @param[in] Ubar Vector of normalized PSD values at the bin-average frequencies.
  /// @param[in] tau PSD scale factor `tau = sigma^2`.
  ///
  /// @return Value of the objective function (scalar).
  template <class Type>
  inline Type NLS<Type>::nll(cRefMatrix_t& Ubar, const Type tau) {
    // res(YU_, Ubar, tau);
    YU_ = Ybar_ - tau * Ubar;
    return YU_.squaredNorm();
  }

  /// Calculates the NLS residual vector
  ///
  /// \f[
  /// R = \bar Y - \hat \tau \cdot \bar U,
  /// \f]
  ///
  /// where \f$\hat \tau\f$ is the conditional optimizer of \f$\tau\f$ calculated by `tau(Ubar)`.
  ///
  /// @param[out] R Vector of residuals.
  /// @param[in] Ubar Vector of normalized PSD values at the bin-average frequencies.
  /// @param[in] tau PSD scale factor `tau = sigma^2`.
  template <class Type>
  inline void NLS<Type>::res(RefMatrix_t R, cRefMatrix_t& Ubar) {
    tau_ = tau(Ubar);
    R = Ybar_ - tau_ * Ubar;
    return;
    // tau_ = tau(Ubar);
    // YU_ = Ybar_ - tau_ * Ubar;
    // return YU_;
  }


  /// Calculates the objective function given `Ubar` at the optimal value of `tau(Ubar)`.
  ///
  /// @param[in] Ubar Vector of normalized PSD values at the bin-average frequencies.
  ///
  /// @return Value of the objective function (scalar).
  template <class Type>
  inline Type NLS<Type>::nlp(cRefMatrix_t& Ubar) {
    tau_ = tau(Ubar);
    return nll(Ubar, tau_);
  }

}

#endif

