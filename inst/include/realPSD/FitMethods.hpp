/// @file FitMethods.hpp

#ifndef realPSD_FitMethods_hpp
#define realPSD_FitMethods_hpp 1

// #include "LPMethods.hpp"
// #include "NLSMethods.hpp"
// #include "MLEMethods.hpp"
#include "LPFit.hpp"
#include "NLSFit.hpp"
#include "MLEFit.hpp"

namespace realPSD {
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

  template<class Type, class UFun>
  Type FitMethods(objective_function<Type>* obj,
		  UFun (*make_Ufun)(int, objective_function<Type>*)) {
    // pick method
    DATA_STRING(method);
    if(method == "UFun") {
      // data
      DATA_MATRIX(f);
      // parameters
      PARAMETER_MATRIX(phi);
      // calculate U
      int N = f.size();
      // UFun Ufun(N, obj);
      UFun Ufun = make_Ufun(N, obj);
      Ufun.set_f(f);
      SIMULATE {
        matrix<Type> U(N,1);
        Ufun.eval(U, phi);
        REPORT(U);
      }
      return Type(0);
    } else if(method.find("LP_") == 0) {
      // data
      DATA_MATRIX(fbar);
      DATA_MATRIX(Zbar);
      DATA_VECTOR(fs);
      // parameters
      PARAMETER_MATRIX(phi);
      // calculate Ubar
      int N = fbar.size();
      matrix<Type> Ubar(N,1);
      // UFun Ufun(N, obj);
      UFun Ufun = make_Ufun(N, obj);
      Ufun.set_f(fbar);
      Ufun.eval(Ubar, phi);
      Ubar = Ubar * fs;
      // return LP_methods(obj, method, fbar, Zbar, Ubar);
      LP<Type> lp(N);
      lp.set_Zbar(Zbar);
      if(method == "LP_zeta") {
      	return lp.zeta(Ubar);
      } else if(method == "LP_nlp") {
      	Type nlp = lp.nlp(Ubar);
      	SIMULATE {
      	  Type zeta = lp.get_zeta();
      	  REPORT(zeta);
      	}
      	return nlp;
      } else if(method == "LP_res") {
      	PARAMETER(zeta);
      	matrix<Type> res(N,1); // output variable
      	lp.res(res, Ubar, zeta); // calculate residuals
      	ADREPORT(res); // set them to be differentiable
      	return Type(0.0);
      } else if(method == "LP_nll") {
      	PARAMETER(zeta);
      	return lp.nll(Ubar, zeta);
      } else {
      	// std::string msg = "Unknown method " + method + ".";
      	error("Unknown LP method.");
      }
    } else if(method.find("MLE_") == 0) {
      // data
      DATA_MATRIX(f);
      DATA_MATRIX(Y);
      DATA_VECTOR(fs);
      // parameters
      PARAMETER_MATRIX(phi);
      // calculate U
      int N = f.size();
      matrix<Type> U(N,1);
      // UFun Ufun(N, obj);
      UFun Ufun = make_Ufun(N, obj);
      Ufun.set_f(f);
      Ufun.eval(U, phi);
      U = U * fs;
      // return MLE_methods(obj, method, f, Y, U);
      MLE<Type> mle(N);
      mle.set_Y(Y);
      if(method == "MLE_tau") {
      	return mle.tau(U);
      } else if(method == "MLE_nlp") {
      	Type nlp = mle.nlp(U);
      	SIMULATE {
      	  Type tau = mle.get_tau();
      	  REPORT(tau);
      	}
      	return nlp;
      } else if(method == "MLE_nll") {
      	PARAMETER(tau);
      	return mle.nll(U, tau);
      } else {
      	error("Unknown MLE method.");
      }
    } else if(method.find("NLS_") == 0) {
      // data
      DATA_MATRIX(fbar);
      DATA_MATRIX(Ybar);
      DATA_VECTOR(fs);
      // parameters
      PARAMETER_MATRIX(phi);
      // calculate Ubar
      int N = fbar.size();
      matrix<Type> Ubar(N,1);
      // UFun Ufun(N, obj);
      UFun Ufun = make_Ufun(N, obj);
      Ufun.set_f(fbar);
      Ufun.eval(Ubar, phi);
      Ubar = Ubar * fs;
      // return NLS_methods(obj, method, fbar, Ybar, Ubar);
      NLS<Type> nls(N);
      nls.set_Ybar(Ybar);
      if(method == "NLS_tau") {
      	return nls.tau(Ubar);
      } else if(method == "NLS_nlp") {
      	Type nlp = nls.nlp(Ubar);
      	SIMULATE {
      	  Type tau = nls.get_tau();
      	  REPORT(tau);
      	}
      	return nlp;
      } else if(method == "NLS_res") {
      	PARAMETER(tau);
      	// calculate residuals
      	matrix<Type> res(N,1);
      	nls.res(res, Ubar, tau);
      	ADREPORT(res); // make residuals differentiable
      	return Type(0.0);
      } else if(method == "NLS_nll") {
      	PARAMETER(tau);
      	return nls.nll(Ubar, tau);
      } else {
      	error("Unknown NLS method.");
      }
    } else {
      error("Unknown method.");
    }
    return Type(0.0);
  }

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
} // end namespace realPSD

#endif
