/*---------------------------------------------------------
 * Class template for powerseries (infinite polynomials)
 ----------------------------------------------------------*/
#ifndef POWERSERIES_H
#define POWERSERIES_H
#include "POLYNOMIAL.h"
#include <functional>
// TODO this only works for one dimension so far
template <unsigned int n, class T>
class POWERSERIES{
  private:
    std::function<T(unsigned long)> f;
    // the coefficients up to a finite number are cached as polynomial
    POLYNOMIAL<n,T> poly;
  public:
    POWERSERIES(const std::function<T(unsigned long)>& f): f(f) {};
    // partially evaluate the function as a polynomial
    // of at least the first m coefficients
    T evaluate_partial(const T& x, const unsigned long m){
      if(poly.get_degree() <= m){
        int start = poly.get_degree(); 
        for(int i=start; i<=m; i++){
          poly.push(f(i));
        }
      }
      return poly(x);
    }
};
#endif
