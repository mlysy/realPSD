#include "realPSD/SHOWModel.hpp"
#include "realPSD/LPFit.hpp"

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
  } else if(method == "zeta") {
    // data
    DATA_MATRIX(fbar);
    DATA_MATRIX(Zbar);
    // parameters
    PARAMETER_MATRIX(phi);
    // intermediate variables
    int N = fbar.size();
    UFun<Type> Ufun(N);
    LP<Type> lp(N);
    matrix<Type> logUbar(N,1);
    Ufun.set_f(fbar);
    lp.set_Zbar(Zbar);
    // calculate log(Ubar)
    Ufun.eval(logUbar, phi);
    logUbar = logUbar.array().log();
    // calculate zeta_hat
    return lp.zeta(logUbar);
  } else if(method == "nlp") {
    // data
    DATA_MATRIX(fbar);
    DATA_MATRIX(Zbar);
    // parameters
    PARAMETER_MATRIX(phi);
    // intermediate variables
    int N = fbar.size();
    UFun<Type> Ufun(N);
    LP<Type> lp(N);
    matrix<Type> logUbar(N,1);
    Ufun.set_f(fbar);
    lp.set_Zbar(Zbar);
    // calculate log(Ubar)
    Ufun.eval(logUbar, phi);
    logUbar = logUbar.array().log();
    // calculate nlp
    return lp.nlp(logUbar);
  } else if(method == "nll") {
    // data
    DATA_MATRIX(fbar);
    DATA_MATRIX(Zbar);
    // parameters
    PARAMETER_MATRIX(phi);
    PARAMETER(zeta);
    // PARAMETER_MATRIX(logUbar);
    // intermediate variables
    int N = fbar.size();
    UFun<Type> Ufun(N);
    LP<Type> lp(N);
    matrix<Type> logUbar(N,1);
    // matrix<Type> ZLU(N,1);
    // Type ZLU_sum;
    // logUbar.setZero();
    Ufun.set_f(fbar);
    lp.set_Zbar(Zbar);
    // calculate log(Ubar)
    Ufun.eval(logUbar, phi);
    logUbar = logUbar.array().log();
    // ZLU = logUbar;
    // ZLU_sum = ZLU.sum();
    // REPORT(logUbar);
    // REPORT(ZLU);
    // REPORT(ZLU_sum);
    // calculate nll
    // ZLU_sum = ((Zbar - logUbar).array() - zeta).matrix().squaredNorm();
    // ZLU_sum = logUbar(0,0);
    // REPORT(ZLU_sum);
    return lp.nll(logUbar, zeta);
  } else {
    error("Unknown method.");
  }
  return Type(0.0);
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
