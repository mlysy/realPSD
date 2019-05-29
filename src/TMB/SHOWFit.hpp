#include "realPSD/SHOWModel.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type SHOWFit(objective_function<Type>* obj) {
  // using namespace realPSD;
  // data
  DATA_MATRIX(f);
  // parameters
  PARAMETER_MATRIX(phi);
  // evaluate normalized PSD
  int N = f.size();
  UFun<Type> Ufun(N);
  Ufun.set_f(f);
  SIMULATE {
    matrix<Type> U(N,1);
    Ufun.eval(U, phi);
    REPORT(U);
  }
  return Type(0);
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
