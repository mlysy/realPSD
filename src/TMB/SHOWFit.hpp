#include "realPSD/SHOWModel.hpp"
#include "realPSD/LPMethods.hpp"
#include "realPSD/MLEFit.hpp"
#include "realPSD/NLSMethods.hpp"
// #include "realPSD/LPFit.hpp"
// #include "realPSD/NLSFit.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type SHOWFit(objective_function<Type>* obj) {
  using namespace realPSD;
  // pick method
  DATA_STRING(method);
  if(method == "UFun") {
    // data
    DATA_MATRIX(f);
    // parameters
    PARAMETER_MATRIX(phi);
    // evaluate normalized PSD
    int N = f.size();
    UFun<Type> Ufun(N);
    Ufun.set_f(f);
    SIMULATE {
      // calculate U
      matrix<Type> U(N,1);
      Ufun.eval(U, phi);
      REPORT(U);
    }
    return Type(0);
  } else if(method == "LP_zeta") {
    return LP_zeta(obj);
  } else if(method == "LP_nlp") {
    return LP_nlp(obj);
  } else if(method == "LP_nll") {
    return LP_nll(obj);
  } else if(method == "MLE_tau") {
    // FIXME: factor out these methods into MLEMethods.hpp
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
    Ufun.set_f(f);
    mle.set_Y(Y);
    // calculate U
    Ufun.eval(U, phi);
    U = U * fs;
    // calculate tau_hat
    return mle.tau(U);
  } else if(method == "MLE_nll") {
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
    Ufun.set_f(f);
    mle.set_Y(Y);
    // calculate U
    Ufun.eval(U, phi);
    U = U * fs;
    // calculate nll
    return mle.nll(U, tau);
  } else if(method == "MLE_nlp") {
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
    Ufun.set_f(f);
    mle.set_Y(Y);
    // calculate U
    Ufun.eval(U, phi);
    U = U * fs;
    // calculate nlp
    return mle.nlp(U);
  } else if(method == "NLS_tau") {
    return NLS_tau(obj);
  } else if(method == "NLS_nll") {
    return NLS_nll(obj);
  } else if(method == "NLS_nlp") {
    return NLS_nlp(obj);
  } else {
    error("Unknown method.");
  }
  return Type(0.0);
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
