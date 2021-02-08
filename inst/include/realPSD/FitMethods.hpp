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
      UFun Ufun = make_Ufun(f, obj);
      // Ufun.set_f(f);
      matrix<Type> U(N,1);
      Ufun.eval(U, phi);
      ADREPORT(U);
      return Type(0.0);
    } else if(method.find("LP_") == 0) {
      // data
      DATA_MATRIX(fbar);
      DATA_MATRIX(Zbar);
      // DATA_VECTOR(fs);
      DATA_VECTOR(scale);
      // parameters
      PARAMETER_MATRIX(phi);
      // calculate Ubar
      int N = fbar.size();
      matrix<Type> Ubar(N,1);
      // UFun Ufun(N, obj);
      UFun Ufun = make_Ufun(fbar, obj);
      // Ufun.set_f(fbar);
      Ufun.eval(Ubar, phi);
      // Ubar = Ubar * fs;
      // return LP_methods(obj, method, fbar, Zbar, Ubar);
      LP<Type> lp(N);
      lp.set_Zbar(Zbar);
      if(method == "LP_ufun") {
	ADREPORT(Ubar);
	return Type(0.0);
      } else if(method == "LP_zeta") {
      	return lp.zeta(Ubar) + scale(0);
      } else if(method == "LP_nlp") {
      	Type nlp = scale(1) * lp.nlp(Ubar);
      	SIMULATE {
      	  Type zeta = lp.get_zeta() + scale(0);
      	  REPORT(zeta);
      	}
      	return nlp;
      } else if(method == "LP_res") {
      	PARAMETER(zeta);
      	matrix<Type> res(N,1); // output variable
      	lp.res(res, Ubar, zeta - scale(0)); // calculate residuals
      	ADREPORT(res); // set them to be differentiable
      	return Type(0.0);
      } else if(method == "LP_nll") {
      	PARAMETER(zeta);
      	return scale(1) * lp.nll(Ubar, zeta - scale(0));
      } else {
      	// std::string msg = "Unknown method " + method + ".";
      	error("Unknown LP method.");
      }
    } else if(method.find("MLE_") == 0) {
      // data
      DATA_MATRIX(f);
      DATA_MATRIX(Y);
      // DATA_VECTOR(fs);
      DATA_SCALAR(scale);
      // parameters
      PARAMETER_MATRIX(phi);
      // calculate U
      int N = f.size();
      matrix<Type> U(N,1);
      // UFun Ufun(N, obj);
      UFun Ufun = make_Ufun(f, obj);
      // Ufun.set_f(f);
      Ufun.eval(U, phi);
      // U = U * fs;
      // return MLE_methods(obj, method, f, Y, U);
      MLE<Type> mle(N);
      mle.set_Y(Y);
      if(method == "MLE_ufun") {
	ADREPORT(U);
	return Type(0.0);
      }
      else if(method == "MLE_zeta") {
      	return log(mle.tau(U)) + scale;
      } else if(method == "MLE_nlp") {
      	Type nlp = mle.nlp(U) + N * scale;
      	SIMULATE {
      	  Type tau = mle.get_tau();
	  Type zeta = log(tau) + scale;
      	  REPORT(zeta);
      	}
      	return nlp;
      } else if(method == "MLE_nll") {
      	PARAMETER(zeta);
	Type tau = exp(zeta - scale);
      	return mle.nll(U, tau) + N * scale;
      } else if(method == "MLE_res") {
      	PARAMETER(zeta);
      	// calculate residuals
	Type tau = exp(zeta - scale);
      	matrix<Type> res(N,1);
      	mle.res(res, U, tau);
      	ADREPORT(res); // make residuals differentiable
      	return Type(0.0);
      } else {
      	error("Unknown MLE method.");
      }
    } else if(method.find("NLS_") == 0) {
      // data
      DATA_MATRIX(fbar);
      DATA_MATRIX(Ybar);
      // DATA_VECTOR(fs);
      DATA_SCALAR(scale);
      // parameters
      PARAMETER_MATRIX(phi);
      // calculate Ubar
      int N = fbar.size();
      matrix<Type> Ubar(N,1);
      // UFun Ufun(N, obj);
      UFun Ufun = make_Ufun(fbar, obj);
      // Ufun.set_f(fbar);
      Ufun.eval(Ubar, phi);
      // Ubar = Ubar * fs;
      // return NLS_methods(obj, method, fbar, Ybar, Ubar);
      NLS<Type> nls(N);
      nls.set_Ybar(Ybar);
      if(method == "NLS_ufun") {
	ADREPORT(Ubar);
	return Type(0.0);
      } else if(method == "NLS_zeta") {
      	return log(nls.tau(Ubar)) + scale;
      } else if(method == "NLS_nlp") {
      	Type nlp = nls.nlp(Ubar);
      	SIMULATE {
      	  Type tau = nls.get_tau();
	  Type zeta = log(tau) + scale;
      	  REPORT(zeta);
	  // REPORT(tau);
      	}
      	return nlp;
      } else if(method == "NLS_res") {
      	PARAMETER(zeta);
      	// calculate residuals
	Type tau = exp(zeta - scale);
      	matrix<Type> res(N,1);
      	nls.res(res, Ubar, tau);
      	ADREPORT(res); // make residuals differentiable
      	return Type(0.0);
      } else if(method == "NLS_nll") {
      	PARAMETER(zeta);
	Type tau = exp(zeta - scale);
	SIMULATE {
	  REPORT(tau);
	  REPORT(Ubar);
	  REPORT(Ybar);
	}
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
