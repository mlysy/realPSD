#include "realPSD/SHOW_test_Model.hpp"
#include "realPSD/FitMethods.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type SHOW_test(objective_function<Type>* obj) {
  using namespace realPSD;
  return FitMethods<Type, SHOW_test::UFun<Type> >(obj);
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
