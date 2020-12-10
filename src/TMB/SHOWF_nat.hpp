#include "realPSD/SHOWF_nat_Model.hpp"
#include "realPSD/FitMethods.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type SHOWF_nat(objective_function<Type>* obj) {
  using namespace realPSD;
  return FitMethods<Type, SHOWF_nat::UFun<Type> >(obj, SHOWF_nat::make_Ufun<Type>);
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
