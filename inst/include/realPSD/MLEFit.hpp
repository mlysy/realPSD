/// @file MLEFit.hpp

#ifndef realPSD_MLEFit_hpp
#define realPSD_MLEFit_hpp 1

#include "utils.hpp"

namespace realPSD {

  /// Class for the MLE fitting method.
  template <class Type>
  class MLE {
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
    matrix<Type> Y_;
    matrix<Type> YU_;
    matrix<Type> logU_;
    Type tau_;
    // allocate internal memory
    void init(int N);
    // Objective function with internal `YU` and `logU`.
    Type nll(const Type tau);
  public:
    /// Constructor
    MLE(int N);
    /// setter for internal `Y`.
    void set_Y(cRefMatrix_t& Y);
    /// Optimal value of `tau = sigma^2` given `U`.
    Type tau(cRefMatrix_t& U);
    /// Objective function for the MLE method.
    Type nll(cRefMatrix_t& U, const Type tau);
    /// Profiled objective function for the MLE method.
    Type nlp(cRefMatrix_t& U);
  };

  /// @param[in] N Length of `Y`.
  template <class Type>
  inline MLE<Type>::MLE(int N) {
    init(N);
  }

  /// Initializes `Y_`, `YU_` and `logU_` as one-column matrices of size `N x 1`.
  /// @param[in] N Length of `Y`.
  template <class Type>
  inline void MLE<Type>::init(int N) {
    N_ = N;
    Y_ = zero_matrix<Type>(N_, 1);
    YU_ = zero_matrix<Type>(N_, 1);
    logU_ = zero_matrix<Type>(N_, 1);
    return;
  }


  /// Resets the internal value of `Y`.  Optionally reallocates memory if `Y.size() != Y_size()`.
  ///
  /// @param[in] Y Vector of periodogram values.
  template <class Type>
  inline void MLE<Type>::set_Y(cRefMatrix_t& Y) {
    if(Y.size() != N_) init(Y.size());
    Y_ = Y;
    return;
  }

  /// @param[in] U Vector of normalized PSD values.
  ///
  /// @return Scalar estimate of `tau`.
  template <class Type>
  inline Type MLE<Type>::tau(cRefMatrix_t& U) {
    YU_ = Y_.array() / U.array();
    return YU_.sum() / N_;
  }

  /// The MLE objective function is given by
  ///
  /// \f[
  /// Q(\boldsymbol{U}, \tau) = K \log \tau + \sum_{k=1}^{K} \frac{Y_k}{\tau U_k} + \log U_k.
  /// \f]
  ///
  /// It is the negative of the Whittle loglikelihood.
  ///
  /// @param[in] U Vector of normalized PSD values.
  /// @param[in] tau PSD scale factor `tau = sigma^2`.
  ///
  /// @return Value of the objective function (scalar).
  template <class Type>
  inline Type MLE<Type>::nll(cRefMatrix_t& U, const Type tau) {
    YU_ = Y_.array() / U.array();
    logU_ = U.array().log();
    return nll(tau);
  }

  /// @param[in] tau PSD scale factor `tau = sigma^2`.
  ///
  /// @return Value of the objective function (scalar).
  template <class Type>
  inline Type MLE<Type>::nll(const Type tau) {
    YU_ = YU_/tau + logU_;
    return YU_.sum() + N_ * log(tau);
  }

  /// Calculates the objective function given `U` at the optimal value of `tau(U)`.
  /// @param[in] U Vector of normalized PSD values.
  ///
  /// @return Value of the objective function (scalar).
  template <class Type>
  inline Type MLE<Type>::nlp(cRefMatrix_t& U) {
    tau_ = tau(U);
    logU_ = U.array().log();
    return N_ * (1.0 + log(tau_)) + logU_.sum();
  }  

}

#endif
