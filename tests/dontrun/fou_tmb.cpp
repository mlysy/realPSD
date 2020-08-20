/// @file fou_tmb.hpp
/// @brief header file to pass user defined model to TMB interface in R

// #ifndef fou_tmb_hpp
// #define fou_tmb_hpp 1

#include <TMB.hpp>
#include "fou.hpp"
#include "realPSD/FitMethods.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type> 
Type fou_tmb(objective_function<Type>* obj) {
  return realPSD::FitMethods<Type, fou::UFun<Type> >(obj, fou::make_Ufun<Type>);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

template<class Type>
Type objective_function<Type>::operator() () {
  return fou_tmb<Type>(this);
}

// #endif
