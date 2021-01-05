/// @file OU_Generics.hpp
/// @brief Generic code which supplies everything needed to the **TMB** compiler for `OU`.


#include <TMB.hpp>
#include "realPSD/FitMethods.hpp"
#include "OU_Model.hpp" // model class definition
// #include "" // the generic package code

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type OU_Generics(objective_function<Type>* obj) {
  return realPSD::FitMethods<Type, ou::UFun<Type> >(obj, ou::make_Ufun<Type>);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "OU") {
    return OU_Generics<Type>(this);
  } else {
    error("Unknown model.");
  }
  return Type(0.0);
}


