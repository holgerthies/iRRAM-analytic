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
  // forward declarations

  namespace PWRSERIES_IMPL{ // implementation details for having special type for dimension 0
    // forward declarations
    template <unsigned int n, class T>
    class POWERSERIES;

    template <unsigned int n, class T>
    std::shared_ptr<POWERSERIES<n,T>> add(std::shared_ptr<POWERSERIES<n,T>>, std::shared_ptr<POWERSERIES<n,T>>);

    template <unsigned int n, class T>
    std::shared_ptr<POWERSERIES<n,T>> subtract(std::shared_ptr<POWERSERIES<n,T>>, std::shared_ptr<POWERSERIES<n,T>>);

  template <unsigned int n, class T>
    std::shared_ptr<POWERSERIES<n,T>> multiply(std::shared_ptr<POWERSERIES<n,T>>, std::shared_ptr<POWERSERIES<n,T>>);

  template <unsigned int n, class T>
    std::shared_ptr<POWERSERIES<n,T>> unary_minus(std::shared_ptr<POWERSERIES<n,T>>);

    std::shared_ptr<REAL> add(std::shared_ptr<REAL> lhs, std::shared_ptr<REAL> rhs);

    std::shared_ptr<REAL> multiply(std::shared_ptr<REAL> lhs, std::shared_ptr<REAL> rhs);

    std::shared_ptr<REAL> add(std::shared_ptr<REAL> lhs, std::shared_ptr<REAL> rhs);

    std::shared_ptr<REAL> subtract(std::shared_ptr<REAL> lhs, std::shared_ptr<REAL> rhs);

    std::shared_ptr<REAL> unary_minus(std::shared_ptr<REAL> x);

    template <unsigned int n, class T>
      std::shared_ptr<POWERSERIES<n,T>> scalar_multiply(const T& lhs, const std::shared_ptr<POWERSERIES<n,T>>& rhs);

    std::shared_ptr<REAL> scalar_multiply(const REAL& lhs, const std::shared_ptr<REAL>& rhs);

    template <unsigned int n, class T, typename... ARGS>
    T get_coeff(const POWERSERIES<n,T>& pwr, const unsigned long i, ARGS... rest);

    template<class T>
    T get_coeff(const T& x);

    template <unsigned int n, class T>
      T constant_coefficient(const std::shared_ptr<POWERSERIES<n,T>>& pwr);

    template<class T>
      T constant_coefficient(const std::shared_ptr<T>& x);

    template<unsigned int n, class T>
    std::shared_ptr<POWERSERIES<n,T>> inverse(const std::shared_ptr<POWERSERIES<n,T>>&);

    template<unsigned int n, class T>
      std::shared_ptr<POWERSERIES<n,T>> derivative(const std::shared_ptr<POWERSERIES<n,T>>&, const unsigned int, const unsigned int);

    // define coeff type 
    template <unsigned int n, class T>
    struct COEFF_TYPE{
      typedef POWERSERIES<n-1,T> type;   
    };

    template <class T>
    struct COEFF_TYPE<1,T>{
      typedef T type;   
    };

    
    using std::make_shared;

    template <unsigned int n, class T>
    class POWERSERIES{
      using coeff_type = typename COEFF_TYPE<n,T>::type;
      using coeff_ptr = std::shared_ptr<typename COEFF_TYPE<n,T>::type>;
      using seq = std::function<coeff_ptr(const unsigned long)>;
      using sq_ptr = std::shared_ptr<seq>;
    private:

      // helper class to invert a powerseries
      // and cache coefficients to speed up recursion
      class INVERSION{
      private:
	mutable std::vector<coeff_ptr> cache;
	std::shared_ptr<POWERSERIES> pwr;
      public:
	INVERSION(const std::shared_ptr<POWERSERIES>& pwr):
	  pwr(pwr)
	  {};
	INVERSION(const POWERSERIES& pwr):
	  pwr(make_shared<POWERSERIES>(pwr))
	  {};
	coeff_ptr get_coeff(const unsigned long) const;
      };

      friend std::shared_ptr<POWERSERIES> inverse<>(const std::shared_ptr<POWERSERIES>&);
      sq_ptr f;
      mutable std::vector<coeff_ptr> cache; // precomputed coefficients
    public:
      POWERSERIES(const T& x): POWERSERIES([x] (unsigned int i) {
	  if(i==0){
	    return make_shared<coeff_type>(x);
	  } 
	  else{
	    return make_shared<coeff_type>(T());
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
      coeff_ptr  operator[](const unsigned long) const;
      // unary -
      POWERSERIES operator-();

      int known_coeffs() const{
	return 0;
      }
    };

    template<unsigned int n,class T>
    std::shared_ptr<T> make_recursive(std::function<T()> f){
      return make_shared<T>(f());
    }

    template<unsigned int n, class T, typename... ARGS>
    std::shared_ptr<typename COEFF_TYPE<n,T>::type> make_recursive(std::function<T(const unsigned long, ARGS...)> f){
      auto sq_fun = [f] (const unsigned long i)  {
	std::function<T(ARGS...)> pwr_fun = [f,i] (ARGS... rest) {return f(i, rest...);};
	return make_recursive<n-1,T>(pwr_fun);
      };
      return make_shared<typename COEFF_TYPE<n,T>::type>(sq_fun);
    } 



    template <unsigned int n, class T>
    template <typename... ARGS>
    POWERSERIES<n,T>::POWERSERIES(std::function<T(const unsigned long, ARGS...)> f){
      auto sq_fun = [f] (const unsigned long i) -> coeff_ptr {
	std::function<T(ARGS...)> pwr_fun = [f,i] (ARGS... rest) {return f(i, rest...);};
	return make_recursive<n,T>(pwr_fun);
      };
      this->f = sq_ptr(new seq(sq_fun));
    }


    template <unsigned int n, class T, typename... ARGS>
    T get_coeff(const POWERSERIES<n,T>& pwr, const unsigned long i, ARGS... rest){
      return get_coeff(*(pwr[i]), rest...);
    }

    template<class T>
    T get_coeff(const std::shared_ptr<T>& x){
      return *x;
    }

    template<unsigned int n, class T>
      T constant_coefficient(const std::shared_ptr<POWERSERIES<n,T>>& pwr){
      return constant_coefficient((*pwr)[0]);
    }

    template<class T>
    T constant_coefficient(const std::shared_ptr<T>& x){
      return *x;
    }

    template <unsigned int n, class T>
    template <typename... ARGS>
    T POWERSERIES<n,T>::get(ARGS... args){
      return get_coeff(*this, args...);
    }


// fix the first parameter and return a POWERSERIES with decreased dimension
    template <unsigned int n, class T>
      std::shared_ptr<typename COEFF_TYPE<n,T>::type> POWERSERIES<n,T>::operator[](const unsigned long i) const{
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
//      using coeff_type = typename COEFF_TYPE<m,T>::type;
//      auto series = [p,q] (const unsigned int i) -> coeff_type {
//	if(i==0) return p[0];
//	std::vector<coeff_type> c(i+1);
//	for(int k=1; k<=i; k++){
//	  for(auto& P : partitions(i,k)){
//	    coeff_type cp = p[k];
//	    for(auto j : P){
//	      cp = cp*q[j];
//	    }
//	    c[k] = c[k-1]+cp;
//	  }
//	}
//	return c[i];
//      };
//      return POWERSERIES<m,T>(series);
    }


    
    // unary -
    template <unsigned int n, class T>
    POWERSERIES<n,T> POWERSERIES<n,T>::operator-(){
      auto thisptr=std::make_shared<POWERSERIES<n,T>>(*this);
      return *(unary_minus(thisptr));
    }
    // add coefficients
    template <unsigned int n, class T>
    POWERSERIES<n,T> operator+(const POWERSERIES<n,T>& lhs, const POWERSERIES<n,T>& rhs){
      auto lhsptr = std::make_shared<POWERSERIES<n,T>>(lhs);
      auto rhsptr = std::make_shared<POWERSERIES<n,T>>(rhs);
      return *(add<n,T>(lhsptr, rhsptr));
    }

    // subtract coefficients
    template <unsigned int n, class T>
    POWERSERIES<n,T> operator-(const POWERSERIES<n,T>& lhs, const POWERSERIES<n,T>& rhs){
      auto lhsptr = std::make_shared<POWERSERIES<n,T>>(lhs);
      auto rhsptr = std::make_shared<POWERSERIES<n,T>>(rhs);
      return *(subtract<n,T>(lhsptr, rhsptr));
    }

    // division of power series
    template <unsigned int n, class T>
    POWERSERIES<n,T> operator/(const POWERSERIES<n,T>& lhs, const POWERSERIES<n,T>& rhs){
      //return lhs*inverse(rhs);
    }


    // multiply coefficients
    template <unsigned int n, class T>
    POWERSERIES<n,T> operator*(const POWERSERIES<n,T>& lhs, const POWERSERIES<n,T>& rhs){
      auto lhsptr = std::make_shared<POWERSERIES<n,T>>(lhs);
      auto rhsptr = std::make_shared<POWERSERIES<n,T>>(rhs);
      return *(multiply(lhsptr,rhsptr));
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
    std::shared_ptr<POWERSERIES<n,T>> inverse(const std::shared_ptr<POWERSERIES<n,T>>& pwr){
      using inv_class = typename POWERSERIES<n,T>::INVERSION;
      auto coeff_getter = inv_class(pwr);
      return std::make_shared<POWERSERIES<n,T>>([coeff_getter] (const unsigned int i) {
	  return coeff_getter.get_coeff(i);
	});
    }


    std::shared_ptr<REAL> derivative(const std::shared_ptr<REAL>& x, const unsigned int i, const unsigned int d){
      return x;
    }
    // d-th formal derivative w.r.t. the i-th variable
    template<unsigned int n, class T>
      std::shared_ptr<POWERSERIES<n,T>> derivative(const std::shared_ptr<POWERSERIES<n,T>>& pwr, const unsigned int i, const unsigned int d){
      if(i==0){
	return make_shared<POWERSERIES<n,T>>([pwr, d] (unsigned int j){
	    T fact=1; 
	    for(int k=j+1; k<=j+d; k++)
	    {
	      fact *= T(k);
	    }
	    return scalar_multiply(fact,(*pwr)[d+j]);  
	  } );
      }
      else{
	return make_shared<POWERSERIES<n,T>>([pwr, i, d] (unsigned int j){
	    return derivative((*pwr)[j], i-1, d);
	  } );
	
      }
    }
    

    std::shared_ptr<REAL> inverse(const std::shared_ptr<REAL>& x){
      return std::make_shared<REAL>(1/(*x));
    }

    // get coefficient for inverse power series
    template <unsigned int n, class T> 
      std::shared_ptr<typename COEFF_TYPE<n,T>::type> POWERSERIES<n,T>::INVERSION::get_coeff(const unsigned long i) const{
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
	  auto sp = std::make_shared<typename COEFF_TYPE<n,T>::type>(T());
	  cache[j] = sp;
	  for(int k=0; k<=j; k++){
	    cache[j] = add(cache[j], multiply((*pwr)[k],cache[j-k]));
	  }
	  cache[j] = multiply(unary_minus(cache[0]),cache[j]);
	}
      }
      return cache[i];
    }
    
  template <unsigned int n, class T>
  std::shared_ptr<POWERSERIES<n,T>> add(std::shared_ptr<POWERSERIES<n,T>> lhs, std::shared_ptr<POWERSERIES<n,T>> rhs){
    return std::make_shared<POWERSERIES<n,T>>([lhs, rhs] (const unsigned int i) {
	return add((*lhs)[i],(*rhs)[i]);
    });
  }
  
  template <unsigned int n, class T>
  std::shared_ptr<POWERSERIES<n,T>> subtract(std::shared_ptr<POWERSERIES<n,T>> lhs, std::shared_ptr<POWERSERIES<n,T>> rhs){
    return std::make_shared<POWERSERIES<n,T>>([lhs, rhs] (const unsigned int i) {
	return subtract((*lhs)[i],(*rhs)[i]);
    });
  }

  template <unsigned int n, class T>
  std::shared_ptr<POWERSERIES<n,T>> unary_minus(std::shared_ptr<POWERSERIES<n,T>> pwr){
    return std::make_shared<POWERSERIES<n,T>>([pwr] (const unsigned int i) {
	return unary_minus((*pwr)[i]);
    });
  }

  template <unsigned int n, class T>
  std::shared_ptr<POWERSERIES<n,T>> multiply(std::shared_ptr<POWERSERIES<n,T>> lhs, std::shared_ptr<POWERSERIES<n,T>> rhs){
    return std::make_shared<POWERSERIES<n,T>>([lhs, rhs] (const unsigned int i) {
	auto c = std::make_shared<typename COEFF_TYPE<n,T>::type>(T());
	for(int l=0; l<=i;l++){
	  auto coeff = multiply((*lhs)[l],(*rhs)[i-l]);
	  c = add(c,coeff);
	}
	return c;
    });
  }

  std::shared_ptr<REAL> add(std::shared_ptr<REAL> lhs, std::shared_ptr<REAL> rhs){
    return std::make_shared<REAL>(*lhs+*rhs);
  }

  std::shared_ptr<REAL> multiply(std::shared_ptr<REAL> lhs, std::shared_ptr<REAL> rhs){
    return std::make_shared<REAL>(*lhs**rhs);
  }

  std::shared_ptr<REAL> unary_minus(std::shared_ptr<REAL> x){
    return std::make_shared<REAL>(-*x);
  }
    // scalar multiplication
    template <unsigned int n, class T>
      std::shared_ptr<POWERSERIES<n,T>> scalar_multiply(const T& lhs, const std::shared_ptr<POWERSERIES<n,T>>& rhs){
      return make_shared<POWERSERIES<n,T>>([lhs, rhs] (const unsigned int i) {return scalar_multiply(lhs, (*rhs)[i]);});
    }

    std::shared_ptr<REAL> scalar_multiply(const REAL& lhs, const std::shared_ptr<REAL>& rhs){
      return std::make_shared<REAL>(lhs*(*rhs));
    }


    std::shared_ptr<REAL> subtract(std::shared_ptr<REAL> lhs, std::shared_ptr<REAL> rhs){
      return make_shared<REAL>((*lhs)-(*rhs));
    };


  } // PWRSERIES_IMPL namespace
  template <unsigned int n, class T>
  using POWERSERIES = typename PWRSERIES_IMPL::COEFF_TYPE<n+1,T>::type;


} // iRRAM namespace

#endif
