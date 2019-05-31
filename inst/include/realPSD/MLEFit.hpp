/// @file MLEFit.hpp

#ifndef realPSD_MLEFit_hpp
#define realPSD_MLEFit_hpp 1

#include "utils.hpp"

namespace realPSD {

  /// Class for the MLE fitting method
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
    // Objective function with internal `YU`.
    Type nll(const Type tau);
  public:
    /// Constructor
    MLE(int N);
    /// setter for internal `Y`.
    void set_Y(cRefMatrix_t& Y);
    /// Optimal value of `tau = sigma^2` given `phi`.
    Type tau(cRefMatrix_t& U);
    /// Objective function for the MLE method.
    Type nll(cRefMatrix_t& U, const Type tau);
    /// Profiled objective function for the MLE method.
    Type nlp(cRefMatrix_t& U);
  };

  template <class Type>
  inline MLE<Type>::MLE(int N) {
    init(N);
  }

  template <class Type>
  inline void MLE<Type>::init(int N) {
    N_ = N;
    Y_ = zero_matrix<Type>(N_, 1);
    YU_ = zero_matrix<Type>(N_, 1);
    logU_ = zero_matrix<Type>(N_, 1);
    return;
  }


  template <class Type>
  inline void MLE<Type>::set_Y(cRefMatrix_t& Y) {
    if(Y.size() != N_) init(Y.size());
    Y_ = Y;
    return;
  }

  /// @param[in] Vector of normalized PSD values.
  ///
  /// @return Scalar estimate of `tau`.
  template <class Type>
  inline Type MLE<Type>::tau(cRefMatrix_t& U) {
    YU_ = Y_.array() / U.array();
    return YU_.sum() / N_;
  }

  /// Objective function for the MLE method.
  ///
  /// This corresponds to the negative Whittle loglikelihood.
  ///
  /// @param[in] Vector of normalized PSD values.
  ///
  /// @return Scalar value of the objective function.
  template <class Type>
  inline Type MLE<Type>::nll(cRefMatrix_t& U, const Type tau) {
    YU_ = Y_.array() / U.array();
    logU_ = U.array().log();
    return nll(tau);
  }

  template <class Type>
  inline Type MLE<Type>::nll(const Type tau) {
    YU_ = YU_/tau + logU_;
    return YU_.sum() + N_ * log(tau);
  }

  template <class Type>
  inline Type MLE<Type>::nlp(cRefMatrix_t& U) {
    tau_ = tau(U);
    logU_ = U.array().log();
    return N_ * (1.0 + log(tau_)) + logU_.sum();
  }  

}

#endif
