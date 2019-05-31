#include "realPSD/SHOWModel.hpp"
#include "realPSD/LPFit.hpp"
#include "realPSD/MLEFit.hpp"
#include "realPSD/NLSFit.hpp"

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
    // data
    DATA_MATRIX(fbar);
    DATA_MATRIX(Zbar);
    // parameters
    PARAMETER_MATRIX(phi);
    // intermediate variables
    int N = fbar.size();
    UFun<Type> Ufun(N);
    LP<Type> lp(N);
    matrix<Type> Ubar(N,1);
    Ufun.set_f(fbar);
    lp.set_Zbar(Zbar);
    // calculate Ubar
    Ufun.eval(Ubar, phi);
    // calculate zeta_hat
    return lp.zeta(Ubar);
  } else if(method == "LP_nlp") {
    // data
    DATA_MATRIX(fbar);
    DATA_MATRIX(Zbar);
    // parameters
    PARAMETER_MATRIX(phi);
    // intermediate variables
    int N = fbar.size();
    UFun<Type> Ufun(N);
    LP<Type> lp(N);
    matrix<Type> Ubar(N,1);
    Ufun.set_f(fbar);
    lp.set_Zbar(Zbar);
    // calculate Ubar
    Ufun.eval(Ubar, phi);
    // calculate nlp
    return lp.nlp(Ubar);
  } else if(method == "LP_nll") {
    // data
    DATA_MATRIX(fbar);
    DATA_MATRIX(Zbar);
    // parameters
    PARAMETER_MATRIX(phi);
    PARAMETER(zeta);
    // intermediate variables
    int N = fbar.size();
    UFun<Type> Ufun(N);
    LP<Type> lp(N);
    matrix<Type> Ubar(N,1);
    Ufun.set_f(fbar);
    lp.set_Zbar(Zbar);
    // calculate Ubar
    Ufun.eval(Ubar, phi);
    return lp.nll(Ubar, zeta);
  } else if(method == "MLE_tau") {
    // data
    DATA_MATRIX(f);
    DATA_MATRIX(Y);
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
    // calculate tau_hat
    return mle.tau(U);
  } else if(method == "MLE_nll") {
    // data
    DATA_MATRIX(f);
    DATA_MATRIX(Y);
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
    // calculate nll
    return mle.nll(U, tau);
  } else if(method == "MLE_nlp") {
    // data
    DATA_MATRIX(f);
    DATA_MATRIX(Y);
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
    // calculate nlp
    return mle.nlp(U);
  } else if(method == "NLS_tau") {
    // data
    DATA_MATRIX(fbar);
    DATA_MATRIX(Ybar);
    // parameters
    PARAMETER_MATRIX(phi);
    // intermediate variables
    int N = fbar.size();
    UFun<Type> Ufun(N);
    NLS<Type> nls(N);
    matrix<Type> Ubar(N,1);
    Ufun.set_f(fbar);
    nls.set_Ybar(Ybar);
    // calculate Ubar
    Ufun.eval(Ubar, phi);
    // calculate tau_hat
    return nls.tau(Ubar);
  } else if(method == "NLS_nll") {
    // data
    DATA_MATRIX(fbar);
    DATA_MATRIX(Ybar);
    // parameters
    PARAMETER_MATRIX(phi);
    PARAMETER(tau);
    // intermediate variables
    int N = fbar.size();
    UFun<Type> Ufun(N);
    NLS<Type> nls(N);
    matrix<Type> Ubar(N,1);
    Ufun.set_f(fbar);
    nls.set_Ybar(Ybar);
    // calculate Ubar
    Ufun.eval(Ubar, phi);
    // calculate nll
    return nls.nll(Ubar, tau);
  } else if(method == "NLS_nlp") {
    // data
    DATA_MATRIX(fbar);
    DATA_MATRIX(Ybar);
    // parameters
    PARAMETER_MATRIX(phi);
    // intermediate variables
    int N = fbar.size();
    UFun<Type> Ufun(N);
    NLS<Type> nls(N);
    matrix<Type> Ubar(N,1);
    Ufun.set_f(fbar);
    nls.set_Ybar(Ybar);
    // calculate Ubar
    Ufun.eval(Ubar, phi);
    // calculate nlp
    return nls.nlp(Ubar);
  } else {
    error("Unknown method.");
  }
  return Type(0.0);
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
