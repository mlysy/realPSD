#include "realPSD/SHOW_comp_Model.hpp"
#include "realPSD/FitMethods.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type SHOW_comp(objective_function<Type>* obj) {
  using namespace realPSD;
  return FitMethods<Type, SHOW_comp::UFun<Type> >(obj);
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
