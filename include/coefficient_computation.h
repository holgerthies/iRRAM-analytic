/*-------------------------------------------------
 * Helper class to compute taylor coefficients from function
 ------------------------------------------------*/
#ifndef COEFFICIENT_COMPUTATION_H
#define COEFFICIENT_COMPUTATION_H
#include <vector>
#include <list>
#include "combinatorics.h"

namespace iRRAM{
  template<unsigned long d, class T>
  class coefficient_computation{
  private:
    std::function<coefficient_computation<d-1,T>(const T&)> f;
    REAL M;
    REAL r;
  public:
    template<typename... ARGS>
      coefficient_computation(const std::function<T(const T&, ARGS...)>& f, const REAL& M, const REAL& r) :  M(M), r(r) {
      this->f = [f,M,r] (const T& x) ->coefficient_computation<d-1,T> {
	std::function<T(ARGS...)> g=[f,x] (ARGS... rest) {return f(x,rest...);};
	return coefficient_computation<d-1,T>(g, M,r);
      };
    };
    T get(std::list<unsigned long> ks) const{
      const unsigned long k = ks.front();
      ks.pop_front();
      if(k==0){
	return f(0).get(ks);
      }
      // interpolate function by polynomial 
      int n=2;
      T best;
      sizetype sum_error, interpolation_error, total_error, best_error;
      sizetype_set(best_error, 1 , 1000);
      REAL eps=power(minimum(r/2, 0.5), n);
      do{
	n *= 2;
	eps *= eps;
	T h=eps/int(k);
	T D=1/h;
	T curr;
	for(int i=0; i<=k; i++){
	  T fact = (k-i) % 2 == 0 ? 1 : -1;
	  T p=i*h;
	  curr =curr+ fact*T(choose(k,i))*f(p).get(ks);
	}
	curr *= power(D,k)*inv_factorial(k);
	REAL error = M*power(2,k+1)*eps/power(r,k+1);
	curr.geterror(sum_error);
	sizetype_add(interpolation_error,error.vsize,error.error);
	sizetype_add(total_error, sum_error, interpolation_error);
	curr.seterror(total_error);
	if(sizetype_less(total_error, best_error)){
	  best_error = total_error;
	  best=curr;
	}
      } while(sizetype_less(sum_error, interpolation_error) && interpolation_error.exponent > ACTUAL_STACK.actual_prec);
      return best;
    }

    
  };

  template<class T>
   class coefficient_computation<0,T>{
  private:
    const std::function<T()> f;
    const REAL& M;
    const REAL& r;
  public:
    coefficient_computation(const std::function<T()>& f, const REAL& M, const REAL& r) : f(f), M(M), r(r) {};
    T get(const std::list<unsigned long>& ks ){
      return f();
    }
  };

  template<unsigned long i, unsigned long n, class T>
  struct series_computation_rec{
    std::shared_ptr<coefficient_computation<n,T>> cc;
    std::list<unsigned long> params;
    POWERSERIES<i-1,T> get_series(const unsigned int j) const{
      series_computation_rec<i-1, n, T> next;
      next.cc = cc; 
      next.params = params;
      next.params.push_back(j);
      return POWERSERIES<i-1, T>([next] (const unsigned long k) { return next.get_series(k); });
    }
  };

  template<unsigned long n, class T>
  struct series_computation_rec<1,n,T>{
    std::shared_ptr<coefficient_computation<n,T>> cc;
    std::list<unsigned long> params;
    T get_series(const unsigned long j) const{
      std::list<unsigned long> all_params(params);
      all_params.push_back(j);
      return cc->get(all_params);
    }
  }; 

  template<unsigned long n, class T>
    using series_computation = series_computation_rec<n,n, T>;

  template<unsigned long n, class T, typename... ARGS>
  POWERSERIES<n,T> make_series(const std::function<T(const T&, ARGS...)>& f, const REAL& M, const REAL& r){
    auto sc = std::make_shared<series_computation<n,T>>();
    sc->cc = std::make_shared<coefficient_computation<n,T>>(coefficient_computation<n,T>(f, M, r));
    return POWERSERIES<n,T>([sc] (const unsigned long i) { return sc->get_series(i); });
  }


}
#endif
