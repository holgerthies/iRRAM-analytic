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
// a powerseries of dimension 0 should just be a scalar
// so POWERSERIES<0,T> should be an alias for T
template <unsigned int n, class T>
struct PWRSERIES_TYPE{
typedef POWERSERIES<n,T> type;   
};

template <class T>
struct PWRSERIES_TYPE<0,T>{
typedef T type;
};

template <unsigned int n, class T>
class POWERSERIES{
using seq = std::function<T(const std::vector<unsigned long>&)>;
using sq_ptr = std::shared_ptr<seq>;
private:
sq_ptr f;
public:
    POWERSERIES(sq_ptr f): f(f) {};
    // get the coefficient with given index
    POWERSERIES(seq  f): POWERSERIES(sq_ptr(new seq(f))) {};
    T get(const std::vector<unsigned long>& v) const {
    return (*f)(v);
    }
    typename PWRSERIES_TYPE<n-1,T>::type  operator[](const unsigned long) const;

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

// partially evaluate the function at x as a polynomial
// using coefficients that sum up from start to end
template <unsigned int n, class T>
T evaluate_partial(const std::shared_ptr<POWERSERIES<n, T>> pwr, const std::vector<T>& x, const unsigned long start, const unsigned long end){
  T ans;
  for(unsigned long m=start; m<=end; m++){
    auto P = iRRAM::partitions(m, n);
    for(auto p : P){
      T prod=1;
      for(int i=0; i<n; i++){
        prod *= power(x[i], p[i]);
      }
      prod *= pwr->get(p);
      ans += prod;
    }
  }
  return ans;
}

// add coefficients
template <unsigned int n, class T>
POWERSERIES<n,T> operator+(const POWERSERIES<n,T>& lhs, const POWERSERIES<n,T>& rhs){
  return POWERSERIES<n,T>([lhs, rhs] (const std::vector<unsigned long>& v) {return lhs.get(v)+rhs.get(v);});
}

// subtract coefficients
template <unsigned int n, class T>
POWERSERIES<n,T> operator-(const POWERSERIES<n,T>& lhs, const POWERSERIES<n,T>& rhs){
  return POWERSERIES<n,T>([lhs, rhs] (const std::vector<unsigned long>& v) {return lhs.get(v)-rhs.get(v);});
}

// fix the first parameter and return a POWERSERIES with decreased dimension
template <unsigned int n, class T>
typename PWRSERIES_TYPE<n-1,T>::type POWERSERIES<n,T>::operator[](const unsigned long i) const{
  auto f=this->f;
  return POWERSERIES<n-1,T>([f,i] (const std::vector<unsigned long>& v) {
      std::vector<unsigned long> ind{i};
      ind.insert(ind.end(), v.begin(), v.end());
      return (*f)(ind);
    });
}

// for n=1 just return the scalar given by the coefficient
// todo: since partial template specialization is not allowed
// this only works for REAL at the moment.
template <>
typename PWRSERIES_TYPE<0,REAL>::type POWERSERIES<1,REAL>::operator[](const unsigned long i) const
  {
    return get({i});
  }
// multiply coefficients
template <unsigned int n, class T>
POWERSERIES<n,T> operator*(const POWERSERIES<n,T>& lhs, const POWERSERIES<n,T>& rhs){
  return POWERSERIES<n,T>([lhs, rhs] (const std::vector<unsigned long>& v) {
      T c;
      for(int l=0; l<=v[0];l++){
	auto coeff = lhs[l]*rhs[v[0]-l];
	std::vector<unsigned long> w;
	w.insert(w.begin(), v.begin()+1, v.end());
	c += coeff.get(w);
      }
      return c;
    });
}
template <class T>
POWERSERIES<1,T> operator*(const POWERSERIES<1,T>& lhs, const POWERSERIES<1,T>& rhs){
  return POWERSERIES<1,T>([lhs, rhs] (const std::vector<unsigned long>& v) {
      T c;
      for(int l=0; l<=v[0];l++){
	c += lhs[l]*rhs[v[0]-l];
      }
      return c;
    });
}
} // PWRSERIES_IMPL namespace
template <unsigned int n, class T>
  using POWERSERIES = typename PWRSERIES_IMPL::PWRSERIES_TYPE<n,T>::type;
} // iRRAM namespace

#endif
