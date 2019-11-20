#include "realPSD/SHOWModel.hpp"
#include "realPSD/UFunMethods.hpp"
#include "realPSD/MLEMethods.hpp"
#include "realPSD/LPMethods.hpp"
#include "realPSD/NLSMethods.hpp"
// #include "realPSD/MLEFit.hpp"
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
    return UFun_eval(obj);
  } else if(method == "LP_zeta") {
    return LP_zeta(obj);
  } else if(method == "LP_nlp") {
    return LP_nlp(obj);
  } else if(method == "LP_nll") {
    return LP_nll(obj);
  } else if(method == "MLE_tau") {
    return MLE_tau(obj);
  } else if(method == "MLE_nll") {
    return MLE_nll(obj);
  } else if(method == "MLE_nlp") {
    return MLE_nlp(obj);
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
