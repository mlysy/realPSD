#include <TMB.hpp>

#include "{{PSDTemplate}}.hpp"
#include "realPSD/FitMethods.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type {{PSDTemplate}}(objective_function<Type>* obj) {
  using namespace realPSD;
  return FitMethods<Type, {{PSDTemplate}}::UFun<Type> >(obj, {{PSDTemplate}}::make_Ufun<Type>);
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "{{PSDTemplate}}") {
    return {{PSDTemplate}}(this);
  } else {
    error("Unknown model.");
  }
  return 0;
}
