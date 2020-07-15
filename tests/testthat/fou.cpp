/// @file fou.cpp
///
/// @brief fractional OU model, passing data to TMB

#include <TMB.hpp>

namespace fou {
  /// create `matrix<Type>` of zeros.
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

  /// A template class for fractional OU model
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
    int N_; ///> problem dimensions
    matrix<Type> f2_; ///> column vector of squared frequencies.
  public:
    /// Constructor.
    UFun(int N);
    /// Set frequency vector.
    void set_f(cRefMatrix_t& f);
    /// Evaluate the normalized PSD.
    void eval(RefMatrix_t U, cRefMatrix_t& phi, cRefMatrix_t& f);
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

  /// Parameters are: `phi = (H, gamma)`
  // psd = |f|^(1-2H) / (f^2 + gamma^2)
  template <class Type>
  inline void UFun<Type>::eval(RefMatrix_t U, cRefMatrix_t& phi, cRefMatrix_t& f) {
    U = f2_.array() + phi(1,0) * phi(1,0);
    U = (f.array().abs()).pow(1-Type(2.0)*phi(0,0)) / U.array(); // use ArrayBase::pow() and ::abs() functions
    return;
  }

  #undef TMB_OBJECTIVE_PTR
  #define TMB_OBJECTIVE_PTR obj
    template<class Type>
    UFun<Type> make_Ufun(int N, objective_function<Type>* obj) {
      UFun<Type> Ufun(N);
      return Ufun;
    }
  #undef TMB_OBJECTIVE_PTR
  #define TMB_OBJECTIVE_PTR this
} // end namespace fou


#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type, class UFun>
Type FitMethods(objective_function<Type>* obj,
    UFun (*make_Ufun)(int, objective_function<Type>*)) {
  // data
  DATA_MATRIX(f);
  // parameters
  PARAMETER_MATRIX(phi);
  // calculate U
  int N = f.size();
  UFun Ufun = make_Ufun(N, obj);
  Ufun.set_f(f);
  SIMULATE {
    matrix<Type> U(N,1);
    Ufun.eval(U, phi, f);
    REPORT(U);
  }
  return Type(0.0);
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

template<class Type>
Type objective_function<Type>::operator() () {
  using namespace fou;
  return FitMethods<Type, UFun<Type> >(this, make_Ufun<Type>);
}

