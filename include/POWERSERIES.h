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
#include "combinatorics.h"
template <unsigned int n, class T>
class POWERSERIES{
  private:
    std::shared_ptr<std::function<T(const std::vector<unsigned long>&)>> f;
  public:
    POWERSERIES(std::shared_ptr<std::function<T(const std::vector<unsigned long>&)>>  f): f(f) {};
    // partially evaluate the function as a polynomial
    // of at least the first m coefficients
    T evaluate_partial(const std::vector<T>& x, const unsigned long start, const unsigned long end){
      T ans;
      for(unsigned long m=start; m<=end; m++){
        auto P = iRRAM::partitions(m, n);
        for(auto p : P){
          T prod=1;
          for(int i=0; i<n; i++){
            prod *= power(x[i], p[i]);
          }
          prod *= get(p);
          ans += prod;
        }
      }
      return ans;
    }

    // get the coefficient with given index
    T get(const std::vector<unsigned long>& v) {
      return (*f)(v);
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
				return prod * pwr.get({ i+d });
        })));
}
#endif
