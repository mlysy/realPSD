/// @file {{{Model}}}_Generics.hpp
/// @brief Generic code which supplies everything needed to the **TMB** compiler for `{{{Model}}}`.

{{#Hide}}
#ifndef {{{Model}}}_Generics_hpp
#define {{{Model}}}_Generics_hpp 1
{{/Hide}}

{{#Show}}
#include <TMB.hpp>
{{/Show}}
#include "{{{Header}}}.hpp" // model class definition
#include "{{{Include}}}.hpp" // the generic package code

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type {{{Model}}}_Generics(objective_function<Type>* obj) {
  return {{{GenericMethods}}}<Type, {{{Class}}}<Type> >(obj, {{{Ctor}}}<Type>);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

{{#Show}}
template<class Type>
Type objective_function<Type>::operator() () {
  return {{{Model}}}_Generics<Type>(this);
}
{{/Show}}

{{#Hide}}
#endif
{{/Hide}}
