#include "realPSD/SHOW_test_Model.hpp"
#include "realPSD/FitMethods.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type SHOW_test(objective_function<Type>* obj) {
  using namespace realPSD;
  return FitMethods<Type, SHOW_test::UFun<Type> >(obj, SHOW_test::make_Ufun<Type>);
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

// template<class Type>
// Type objective_function<Type>::operator() () {
//   using namespace realPSD;
//   return FitMethods<Type, SHOW_test::UFun<Type> >(this, SHOW_test::make_Ufun<Type>);
// }
