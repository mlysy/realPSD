#include "realPSD/SHOW_nat_Model.hpp"
#include "realPSD/FitMethods.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type SHOW_nat(objective_function<Type>* obj) {
  using namespace realPSD;
  return FitMethods<Type, SHOW_nat::UFun<Type> >(obj);
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
