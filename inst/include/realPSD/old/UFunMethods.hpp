/// @file UFunMethods.hpp

#ifndef realPSD_UFunMethods_hpp
#define realPSD_UFunMethods_hpp 1

namespace realPSD {
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

  template<class Type>
  Type UFun_eval(objective_function<Type>* obj,
		 const matrix<Type>& f, matrix<Type>& phi) {
    // evaluate normalized PSD
    int N = f.size();
    UFun<Type> Ufun(N, obj);
    Ufun.set_f(f);
    SIMULATE {
      // calculate U
      matrix<Type> U(N,1);
      Ufun.eval(U, phi);
      REPORT(U);
    }
    return Type(0);
  }

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
} // end namespace realPSD

#endif
