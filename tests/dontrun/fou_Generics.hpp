/// @file fou_Generics.hpp
/// @brief Generic code which supplies everything needed to the **TMB** compiler for `fou`.

#ifndef fou_Generics_hpp
#define fou_Generics_hpp 1

#include "fou.hpp" // model class definition
#include "realPSD/FitMethods.hpp" // the generic package code

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type fou_Generics(objective_function<Type>* obj) {
  return realPSD::FitMethods<Type, fou::UFun<Type> >(obj, fou::make_Ufun<Type>);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this


#endif

