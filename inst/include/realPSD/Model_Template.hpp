/// @file {{{Model}}}_Generics.hpp
/// @brief Generic code which supplies everything needed to the **TMB** compiler for `{{{Model}}}`.

{{^Standalone}}
#ifndef {{{Model}}}_Generics_hpp
#define {{{Model}}}_Generics_hpp 1
{{/Standalone}}

{{#Standalone}}
#include <TMB.hpp>
{{/Standalone}}
#include "realPSD/FitMethods.hpp"
#include "{{{Header}}}" // model class definition
// #include "{{{Include}}}" // the generic package code

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type {{{Model}}}_Generics(objective_function<Type>* obj) {
  return {{{GenericMethods}}}<Type, {{{Class}}}<Type> >(obj, {{{Ctor}}}<Type>);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

{{#Standalone}}
template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "{{{Model}}}") {
    return {{{Model}}}_Generics<Type>(this);
  } else {
    error("Unknown model.");
  }
  return Type(0.0);
}
{{/Standalone}}

{{^Standalone}}
#endif
{{/Standalone}}
