/// @file MLEMethods.hpp

#ifndef realPSD_MLEMethods_hpp
#define realPSD_MLEMethods_hpp 1

#include "MLEFit.hpp"

namespace realPSD {
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

  template<class Type>
  Type MLE_tau(objective_function<Type>* obj) {
    // data
    DATA_MATRIX(f);
    DATA_MATRIX(Y);
    DATA_VECTOR(fs);
    // parameters
    PARAMETER_MATRIX(phi);
    // intermediate variables
    int N = f.size();
    UFun<Type> Ufun(N);
    MLE<Type> mle(N);
    matrix<Type> U(N,1);
    // calculate U
    Ufun.set_f(f);
    mle.set_Y(Y);
    Ufun.eval(U, phi);
    U = U * fs;
    // calculate tau_hat
    return mle.tau(U);
  }

  template<class Type>
  Type MLE_nll(objective_function<Type>* obj) {
    // data
    DATA_MATRIX(f);
    DATA_MATRIX(Y);
    DATA_VECTOR(fs);
    // parameters
    PARAMETER_MATRIX(phi);
    PARAMETER(tau);
    // intermediate variables
    int N = f.size();
    UFun<Type> Ufun(N);
    MLE<Type> mle(N);
    matrix<Type> U(N,1);
    // calculate U
    Ufun.set_f(f);
    mle.set_Y(Y);
    Ufun.eval(U, phi);
    U = U * fs;
    // calculate nll
    return mle.nll(U, tau);
  }

  template<class Type>
  Type MLE_nlp(objective_function<Type>* obj) {
    // data
    DATA_MATRIX(f);
    DATA_MATRIX(Y);
    DATA_VECTOR(fs);
    // parameters
    PARAMETER_MATRIX(phi);
    // intermediate variables
    int N = f.size();
    UFun<Type> Ufun(N);
    MLE<Type> mle(N);
    matrix<Type> U(N,1);
    // calculate U
    Ufun.set_f(f);
    mle.set_Y(Y);
    Ufun.eval(U, phi);
    U = U * fs;
    // calculate nlp
    return mle.nlp(U);
  }

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
} // end namespace realPSD

#endif
