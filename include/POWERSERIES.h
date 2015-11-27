/*---------------------------------------------------------
 * Class template for powerseries (infinite multivariate polynomials)
 * A powerseries in d dimensions is given by a function f:|N^d -> T
 ----------------------------------------------------------*/
#ifndef POWERSERIES_H
#define POWERSERIES_H
#include "POLYNOMIAL.h"
#include <functional>
#include <cstdarg>
#include "combinatorics.h"

namespace iRRAM{
  namespace PWRSERIES_IMPL{ // implementation details for having special type for dimension 0
    template <unsigned int n, class T>
    class POWERSERIES;
    // define coeff type 
    template <unsigned int n, class T>
    struct COEFF_TYPE{
      typedef POWERSERIES<n-1,T> type;   
    };

    template <class T>
    struct COEFF_TYPE<1,T>{
      typedef T type;   
    };

    template <unsigned int n, class T>
    class POWERSERIES{
      using coeff_type = typename COEFF_TYPE<n,T>::type;
      using seq = std::function<coeff_type(const unsigned long)>;
      using sq_ptr = std::shared_ptr<seq>;
    private:
      sq_ptr f;
    public:
      POWERSERIES(const T& x): POWERSERIES([x] (unsigned int i) {
	  if(i==0){
	    return coeff_type(x);
	  } 
	  else{
	    return coeff_type(T());
	  }
	}) {}; 
      POWERSERIES() : POWERSERIES(T()) {};
      POWERSERIES(sq_ptr f): f(f) {};
      POWERSERIES(seq f) : f(sq_ptr(new seq(f))) {}; 
      // constructor from function N^n -> T
      template <typename... ARGS>
      POWERSERIES(std::function<T(const unsigned long, ARGS...)>);
      // get the coefficient with given index
      //POWERSERIES(seq  f): POWERSERIES(sq_ptr(new seq(f))) {};
      //T get(const std::vector<unsigned long>& v) const {
      //  return (*f)(v);
      //}
      coeff_type  operator[](const unsigned long) const;

      int known_coeffs() const{
	return 0;
      }
    };

    template <unsigned int n, class T>
    template <typename... ARGS>
    POWERSERIES<n,T>::POWERSERIES(std::function<T(const unsigned long, ARGS...)> f){
      using seq = std::function<coeff_type(const unsigned long)>;
      using sq_ptr = std::shared_ptr<seq>;
      auto sq_fun = [f] (const unsigned long i) {
	std::function<T(ARGS...)> pwr_fun = [f,i] (ARGS... rest) {return f(i, rest...);};
	return POWERSERIES<n-1,T>(pwr_fun);
      };
      this->f = sq_ptr(new seq(sq_fun));
    }

// fix the first parameter and return a POWERSERIES with decreased dimension
    template <unsigned int n, class T>
    typename COEFF_TYPE<n,T>::type POWERSERIES<n,T>::operator[](const unsigned long i) const{
      return (*f)(i);
    }
    // add coefficients
    template <unsigned int n, class T>
    POWERSERIES<n,T> operator+(const POWERSERIES<n,T>& lhs, const POWERSERIES<n,T>& rhs){
      return POWERSERIES<n,T>([lhs, rhs] (const unsigned int i) {return lhs[i]+rhs[i];});
    }

    // subtract coefficients
    template <unsigned int n, class T>
    POWERSERIES<n,T> operator-(const POWERSERIES<n,T>& lhs, const POWERSERIES<n,T>& rhs){
      return POWERSERIES<n,T>([lhs, rhs] (const unsigned int i) {return lhs[i]-rhs[i];});
    }


// multiply coefficients
    template <unsigned int n, class T>
    POWERSERIES<n,T> operator*(const POWERSERIES<n,T>& lhs, const POWERSERIES<n,T>& rhs){
      return POWERSERIES<n,T>([lhs, rhs] (const unsigned long i) {
	  typename COEFF_TYPE<n,T>::type c;
	  for(int l=0; l<=i;l++){
	    auto coeff = lhs[l]*rhs[i-l];
	    c = c+coeff;
	  }
	  return c;
	});
    }
  } // PWRSERIES_IMPL namespace
  template <unsigned int n, class T>
  using POWERSERIES = PWRSERIES_IMPL::POWERSERIES<n,T>;
} // iRRAM namespace

#endif
