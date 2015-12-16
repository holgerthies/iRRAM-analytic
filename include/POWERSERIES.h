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
    // forward declarations
    template <unsigned int n, class T>
    class POWERSERIES;

    template <unsigned int n, class T, typename... ARGS>
    T get_coeff(const POWERSERIES<n,T>& pwr, const unsigned long i, ARGS... rest);

    template<class T>
    T get_coeff(const T& x);

    template <unsigned int n, class T>
    T constant_coefficient(const POWERSERIES<n,T>& pwr);

    template<class T>
    T constant_coefficient(const T& x);

    template<unsigned int n, class T>
    POWERSERIES<n,T> inverse(const POWERSERIES<n,T>&);

    template<unsigned int n, class T>
    POWERSERIES<n,T> derivative(const POWERSERIES<n,T>&, const unsigned int, const unsigned int);

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

      // helper class to invert a powerseries
      // and cache coefficients to speed up recursion
      class INVERSION{
      private:
	mutable std::vector<coeff_type> cache;
	std::shared_ptr<POWERSERIES> pwr;
      public:
	INVERSION(const POWERSERIES& pwr):
	  pwr(std::make_shared<POWERSERIES>(pwr))
	  {};
	coeff_type get_coeff(const unsigned long) const;
      };

      friend POWERSERIES inverse<>(const POWERSERIES&);
      sq_ptr f;
      mutable std::vector<coeff_type> cache; // precomputed coefficients
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
      // copy constructor
      POWERSERIES(const POWERSERIES& pwr) {
	// do not copy cache
	this->f = pwr.f;
      }
      // constructor from function N^n -> T
      template <typename... ARGS>
      POWERSERIES(std::function<T(const unsigned long, ARGS...)>);
      // get the coefficient with given index
      template <typename... ARGS>
      T get(ARGS...);
      coeff_type  operator[](const unsigned long) const;
      // unary -
      POWERSERIES operator-();

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

    template <unsigned int n, class T, typename... ARGS>
    T get_coeff(const POWERSERIES<n,T>& pwr, const unsigned long i, ARGS... rest){
      return get_coeff(pwr[i], rest...);
    }

    template<class T>
    T get_coeff(const T& x){
      return x;
    }

    template<unsigned int n, class T>
    T constant_coefficient(const POWERSERIES<n,T>& pwr){
      return constant_coefficient(pwr[0]);
    }

    template<class T>
    T constant_coefficient(const T& x){
      return x;
    }

    template <unsigned int n, class T>
    template <typename... ARGS>
    T POWERSERIES<n,T>::get(ARGS... args){
      return get_coeff(*this, args...);
    }


// fix the first parameter and return a POWERSERIES with decreased dimension
    template <unsigned int n, class T>
    typename COEFF_TYPE<n,T>::type POWERSERIES<n,T>::operator[](const unsigned long i) const{
      if(cache.size() <= i){
	auto sz=cache.size();
	cache.resize(i+1);
	// fill cache
	for(int j=sz; j<=i;j++ )
	{
	  cache[j] = (*f)(j);
	}
      }
      return cache[i];
    }

    // composition of formal powerseries
    // insert for each variable one powerseries
    // recursively applying the faa di bruno formula
    // for each of the inserted powerseries all coefficients
    // with first parameter 0 have to be 0
    template<unsigned int m, class T>
    POWERSERIES<m,T> compose(const POWERSERIES<1,T>& p, const POWERSERIES<m,T>& q){
      using coeff_type = typename COEFF_TYPE<m,T>::type;
      auto series = [p,q] (const unsigned int i) -> coeff_type {
	if(i==0) return p[0];
	std::vector<coeff_type> c(i+1);
	for(int k=1; k<=i; k++){
	  for(auto& P : partitions(i,k)){
	    coeff_type cp = p[k];
	    for(auto j : P){
	      cp = cp*q[j];
	    }
	    c[k] = c[k-1]+cp;
	  }
	}
	return c[i];
      };
      return POWERSERIES<m,T>(series);
    }


    
    // unary -
    template <unsigned int n, class T>
    POWERSERIES<n,T> POWERSERIES<n,T>::operator-(){
      auto thisptr=std::make_shared<POWERSERIES<n,T>>(*this);
      return POWERSERIES<n,T>([thisptr] (const unsigned long i) {return -(*thisptr)[i];});
    }
    // add coefficients
    template <unsigned int n, class T>
    POWERSERIES<n,T> operator+(const POWERSERIES<n,T>& lhs, const POWERSERIES<n,T>& rhs){
      auto lhsptr = std::make_shared<POWERSERIES<n,T>>(lhs);
      auto rhsptr = std::make_shared<POWERSERIES<n,T>>(rhs);
      return POWERSERIES<n,T>([lhsptr, rhsptr] (const unsigned int i) {return (*lhsptr)[i]+(*rhsptr)[i];});
    }

    // subtract coefficients
    template <unsigned int n, class T>
    POWERSERIES<n,T> operator-(const POWERSERIES<n,T>& lhs, const POWERSERIES<n,T>& rhs){
      auto lhsptr = std::make_shared<POWERSERIES<n,T>>(lhs);
      auto rhsptr = std::make_shared<POWERSERIES<n,T>>(rhs);
      return POWERSERIES<n,T>([lhsptr, rhsptr] (const unsigned int i) {return (*lhsptr)[i]-(*rhsptr)[i];});
    }

    // division of power series
    template <unsigned int n, class T>
    POWERSERIES<n,T> operator/(const POWERSERIES<n,T>& lhs, const POWERSERIES<n,T>& rhs){
      return lhs*inverse(rhs);
    }


    // multiply coefficients
    template <unsigned int n, class T>
    POWERSERIES<n,T> operator*(const POWERSERIES<n,T>& lhs, const POWERSERIES<n,T>& rhs){
      auto lhsptr = std::make_shared<POWERSERIES<n,T>>(lhs);
      auto rhsptr = std::make_shared<POWERSERIES<n,T>>(rhs);
      return POWERSERIES<n,T>([lhsptr, rhsptr] (const unsigned long i) {
	  typename COEFF_TYPE<n,T>::type c;
	  for(int l=0; l<=i;l++){
	    auto coeff = (*lhsptr)[l]*(*rhsptr)[i-l];
	    c = c+coeff;
	  }
	  return c;
	});
    }

    // scalar addition
    template <unsigned int n, class T>
    POWERSERIES<n,T> operator+(const T& lhs, const POWERSERIES<n,T>& rhs){
      auto lhsptr = std::make_shared<POWERSERIES<n,T>>(lhs);
      auto rhsptr = std::make_shared<POWERSERIES<n,T>>(rhs);
      return POWERSERIES<n,T>([lhsptr, rhsptr] (const unsigned int i) {
	  if(i==0)
	    return (*lhsptr)+(*rhsptr)[0];
	  return (*rhsptr)[i];
	});
    }

    template <unsigned int n, class T>
    POWERSERIES<n,T> operator+(const POWERSERIES<n,T>& lhs, const T& rhs){
      return rhs+lhs;  
    }
    // scalar multiplication
    template <unsigned int n, class T>
    POWERSERIES<n,T> operator*(const T& lhs, const POWERSERIES<n,T>& rhs){
      auto lhsptr = std::make_shared<T>(lhs);
      auto rhsptr = std::make_shared<POWERSERIES<n,T>>(rhs);
      return POWERSERIES<n,T>([lhsptr, rhsptr] (const unsigned int i) {return (*lhsptr)*(*rhsptr)[i];});
    }
    template <unsigned int n, class T>
    POWERSERIES<n,T> operator*(const POWERSERIES<n,T>& lhs, const T& rhs){
      return rhs*lhs;  
    }

    template <unsigned int n, class T>
    POWERSERIES<n,T> inverse(const POWERSERIES<n,T>& pwr){
      using inv_class = typename POWERSERIES<n,T>::INVERSION;
      auto coeff_getter = inv_class(pwr);
      return POWERSERIES<n,T>([coeff_getter] (const unsigned int i) {
	  return coeff_getter.get_coeff(i);
	});
    }


    template<class T>
    T derivative(const T& x, const unsigned int i, const unsigned int d){
      return x;
    }
    // d-th formal derivative w.r.t. the i-th variable
    template<unsigned int n, class T>
    POWERSERIES<n,T> derivative(const POWERSERIES<n,T>& pwr, const unsigned int i, const unsigned int d){
      auto pwrptr = std::make_shared<POWERSERIES<n,T>>(pwr);
      if(i==0){
	return POWERSERIES<n,T>([pwrptr, d] (unsigned int j){
	    T fact=1; 
	    for(int k=j+1; k<=j+d; k++)
	    {
	      fact *= T(k);
	    }
	    return fact*(*pwrptr)[d+j];  
	  } );
      }
      else{
	return POWERSERIES<n,T>([pwrptr, i, d] (unsigned int j){
	    return derivative((*pwrptr)[j], i-1, d);
	  } );
	
      }
    }
    

    template<class T>
    T inverse(const T& x){
      return 1/x;
    }

    // get coefficient for inverse power series
    template <unsigned int n, class T> 
    typename COEFF_TYPE<n,T>::type POWERSERIES<n,T>::INVERSION::get_coeff(const unsigned long i) const{
      auto sz=cache.size();
      if(i >= sz)
      {
	cache.resize(i+1);
      }
      for(unsigned long j=sz; j<=i; j++){
	if(j==0){
	  //std::cout << "printing " << std::endl;
	  //print((*pwr)[0]);
	  cache[j] = inverse((*pwr)[0]);
	} else{
	  for(int k=0; k<=j; k++){
	    cache[j] = cache[j] + (*pwr)[k]*cache[j-k];
	  }
	  cache[j] = -cache[0]*cache[j];
	}
      }
      return cache[i];
    }
    
  } // PWRSERIES_IMPL namespace
  template <unsigned int n, class T>
  using POWERSERIES = typename PWRSERIES_IMPL::COEFF_TYPE<n+1,T>::type;
} // iRRAM namespace

#endif
