#include <TMB.hpp>

#include "fou.hpp"
#include "realPSD/FitMethods.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type fou(objective_function<Type>* obj) {
  using namespace realPSD;
  return FitMethods<Type, fou::UFun<Type> >(obj, fou::make_Ufun<Type>);
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "fou") {
    return fou(this);
  } else {
    error("Unknown model.");
  }
  return 0;
}
