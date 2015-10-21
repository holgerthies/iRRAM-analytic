/*---------------------------------------------------------
 * Class template for powerseries (infinite multivariate polynomials)
 * A powerseries in d dimensions is given by a function f:|N^d -> T
 * The result of the function is cached when it is read for the first time
 ----------------------------------------------------------*/
#ifndef POWERSERIES_H
#define POWERSERIES_H
#include "POLYNOMIAL.h"
#include <functional>
#include <cstdarg>
template <unsigned int n, class T>
class POWERSERIES{
  private:
    std::shared_ptr<void> f;
  public:
    template<class... Args,typename std::enable_if<sizeof...(Args) == n, int>::type = 0>
    POWERSERIES(std::shared_ptr<std::function<T(Args... args)>> f): f(f) {};
    // partially evaluate the function as a polynomial
    // of at least the first m coefficients
    T evaluate_partial(std::initializer_list<T> x, const unsigned long start, const unsigned long end){
      for(int m=start; m=end; m++){
      }
    }

    // get the k-th coefficient
    template<class... Args,typename std::enable_if<sizeof...(Args) == n, int>::type = 0>
    T get(Args... args) {
      auto ft =  std::static_pointer_cast<std::function<T(Args...)>>(f);
      return (*ft)(args...);
    }

    int known_coeffs() const{
      return 0;
    }
};
template <class T>
POWERSERIES<1,T> derive(POWERSERIES<1,T>& pwr, const int d){
  return POWERSERIES<1,T>(std::shared_ptr<std::function<T(unsigned long)>>(new std::function<T(unsigned long)>([&pwr,d] (unsigned long i) {
			T prod(1);
				for (unsigned int j = 0; j < d; j++)
					prod *= i+d-j;
				return prod * pwr.get(i+d);
        })));
}
#endif
