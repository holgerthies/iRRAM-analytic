/*-------------------------------------
 * Class for multivariate polynomials
 -------------------------------------*/
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#include <vector>
#include <memory>
#include <string>
namespace iRRAM
{
  // forward declarations
  template<unsigned int n, class T>
  class POLYNOMIAL;

  template<unsigned int n, class T, typename... Ts>
  T evaluate_partial(const POLYNOMIAL<n,T>&, const long, const long, const T&, Ts...);

  template<class T>
  T evaluate_partial(const T&, const long, const long);

  template<unsigned int n, class T>
  int get_total_degree(const POLYNOMIAL<n,T>& P);

  // define coeff type for multivariate polynomial
  template <unsigned int n, class T>
  struct COEFF_TYPE{
    typedef POLYNOMIAL<n-1,T> type;   
  };

  template <class T>
  struct COEFF_TYPE<1,T>{
    typedef T type;   
  };


  template<unsigned int n, class T>
  class POLYNOMIAL{
    using coeff_type = typename COEFF_TYPE<n,T>::type;
  private:
    std::vector<coeff_type> coefficients;
    int total_degree;
  public:
    POLYNOMIAL(std::vector<coeff_type> coeffs) : coefficients(coeffs) {
      total_degree = 0;
      for(int i=0; i<get_degree(); i++){
	total_degree = max(i+get_total_degree(get_coefficient(i)), total_degree);
      }
    };
    friend int get_total_degree<n,T>(const POLYNOMIAL<n,T>&);
    // real degree might be smaller if leading coefficients are 0
    int get_degree() const{
      return coefficients.size();
    }

    // add coefficients to the back
    void push(const coeff_type& x){
      coefficients.push_back(x);
      total_degree = max(total_degree, get_degree()-1+get_total_degree(x));
    }

    // get kth coefficient
    coeff_type get_coefficient(int k) const{
      return coefficients[k];
    }

    template<typename... Ts>
    T operator()(Ts... x) const {
      return evaluate_partial(*this, 0, get_total_degree(*this), x...);
    }
  };
  template<class T>
  int get_total_degree(const T& x){
    return 0;
  }
  template<unsigned int n, class T>
  int get_total_degree(const POLYNOMIAL<n,T>& P){
    return P.total_degree;
  }
  // evaluate only the polynomial given by considering the coefficients
  // between start and end
  template<unsigned int n, class T, typename... Ts>
  T evaluate_partial(const POLYNOMIAL<n,T>& P, const long start, const long end, const T& x, Ts... rest){
    T ans = 0;
    int min_i = (n==1) ? start : 0;
    for(int i=min(end, P.get_degree()-1); i>=min_i; i--){
      ans = evaluate_partial(P.get_coefficient(i), max(0, start-i), end-i, rest...)+ans*x;
    }
    if(start > 0)
      ans *= power(x,start);
    return ans;
  }

  template<class T>
  T evaluate_partial(const T& P, const long start, const long end){
    return P;
  }
}

#endif
