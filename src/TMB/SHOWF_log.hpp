#include "realPSD/SHOWF_log_Model.hpp"
#include "realPSD/FitMethods.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type SHOWF_log(objective_function<Type>* obj) {
  using namespace realPSD;
  return FitMethods<Type, SHOWF_log::UFun<Type> >(obj, SHOWF_log::make_Ufun<Type>);
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
