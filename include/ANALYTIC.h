/*-------------------------------------------------
 * Class for real analytic functions 
 ------------------------------------------------*/
#ifndef ANALYTIC_H
#define ANALYTIC_H
#include "POWERSERIES.h"
#include <vector>
namespace iRRAM
{
  // forward declarations
  template <unsigned int n, class T>
  class ANALYTIC;
  template<unsigned int n, class T>
  ANALYTIC<n,T> operator+(const ANALYTIC<n,T>& lhs, const ANALYTIC<n,T>& rhs);
  template<unsigned int n, class T>
  ANALYTIC<n,T> operator-(const ANALYTIC<n,T>& lhs, const ANALYTIC<n,T>& rhs);
  template<unsigned int n, class T>
  ANALYTIC<n,T> operator*(const ANALYTIC<n,T>& lhs, const ANALYTIC<n,T>& rhs);

  template <unsigned int n, class T>
  class ANALYTIC{
    using vec = std::vector<unsigned long>;
    using seq = std::function<T(const vec&)>;
    using sq_ptr = std::shared_ptr<seq>;
    using pwr_ser = POWERSERIES<n,T>;
  private:
    REAL M,r; // maximum of the function and radius of convergence
    std::shared_ptr<pwr_ser> pwr;
  public:

    // constructor from sequence pointer
    ANALYTIC(const sq_ptr f, const REAL& r, const REAL& M ): 
      M(M),
      r(r),
      pwr(new POWERSERIES<n,T>(f))
      {};

    // constructor from sequence
    ANALYTIC(seq f, const REAL& r, const REAL& M): ANALYTIC(sq_ptr(new seq(f)), r, M) {};

    // constructor from T (constant function)
    ANALYTIC(const T& x){
      r = 10000;
      M = abs(x);
      auto series= [x] (const vec& v) -> T{
	for(int i=0; i<n; i++){
	  if(v[i] != 0) return 0;
	}
	return x;
      };      pwr = std::shared_ptr<pwr_ser>(new POWERSERIES<n,T>(sq_ptr(new seq(series))));
    };
    // constructor from powerseries
    ANALYTIC(pwr_ser& P, const REAL& r, const REAL& M):
      M(M),
      r(r),
      pwr(new pwr_ser(P) )
      {};

    // arithmetic operations
    friend ANALYTIC operator+ <>(const ANALYTIC& lhs, const ANALYTIC& rhs);
    friend ANALYTIC operator- <>(const ANALYTIC& lhs, const ANALYTIC& rhs);
    friend ANALYTIC operator* <>(const ANALYTIC& lhs, const ANALYTIC& rhs);

    T operator ()(const std::vector<T>& x) const;

    int get_known_coeffs() const {};

    REAL get_r() const{return r;}
    REAL get_M() const{return M;}

    T get_coeff(const vec& v) const
      {
	return pwr->get(v);
      }

  };

  template<unsigned int n, class T>
  ANALYTIC<n,T> operator+(const ANALYTIC<n,T>& lhs, const ANALYTIC<n,T>& rhs){
    auto spwr=*lhs.pwr+*rhs.pwr;
    return ANALYTIC<n,T>(spwr, minimum(lhs.r, rhs.r), lhs.M+rhs.M);
  }

  template<unsigned int n, class T>
  ANALYTIC<n,T> operator-(const ANALYTIC<n,T>& lhs, const ANALYTIC<n,T>& rhs){
    auto spwr=*lhs.pwr-*rhs.pwr;
    return ANALYTIC<n,T>(spwr, minimum(lhs.r, rhs.r), lhs.M+rhs.M);
  }
  template<unsigned int n, class T>
  ANALYTIC<n,T> operator*(const ANALYTIC<n,T>& lhs, const ANALYTIC<n,T>& rhs){
    auto spwr=*lhs.pwr**rhs.pwr;
    return ANALYTIC<n,T>(spwr, minimum(lhs.r, rhs.r), lhs.M*rhs.M);
  }

  template<unsigned int n, class T>
  ANALYTIC<n,T> power(const ANALYTIC<n,T>& f, unsigned int m){
    // power by repeated squaring
    if(m == 0) return T(1);
    if(m == 1) return ANALYTIC<n,T>(f);
    if(m & 1) return f*power(f, m-1); 
    auto ans = power(f, m/2);
    return ans*ans;
  }

  template<unsigned int n, class T>
  T ANALYTIC<n,T>::operator ()(const std::vector<T>& x) const
  {
    int J = 1;
    T x_max = abs(*std::max_element(x.begin(), x.end(), [] (T x1, T x2) {return abs(x1) < abs(x2);}));
    REAL error = M*power(x_max/r, J+1)*(J+1);
    REAL sum(evaluate_partial(pwr, x,0,J));
    REAL best=sum;
    sizetype best_error, trunc_error, local_error,sum_error;
    sum.geterror(sum_error);
    sizetype_add(trunc_error,error.vsize,error.error);
    sizetype_add(local_error, sum_error, trunc_error);
    best.seterror(local_error);
    best_error = local_error;
    int step_size = 1;
    REAL error_factor = power(x_max/r, step_size);
    while (sizetype_less(sum_error, trunc_error) &&
	   (trunc_error.exponent >= ACTUAL_STACK.actual_prec) ){
      int new_J = J+step_size;
      sum += evaluate_partial(pwr, x, J+1, new_J); 
      error *= error_factor*REAL(new_J)/REAL(J);
      sizetype_add(trunc_error,error.vsize,error.error);
      sum.geterror(sum_error); // get error of partial sum
      J = new_J;
      sizetype_add(local_error, sum_error, trunc_error);
      if (sizetype_less(local_error, best_error)) { 
	best = sum;
	best.seterror(local_error);
	best_error = local_error;
      }
    } 
    return best;
  }
  // some predefined functions
  namespace AnalyticFunction{
    // projection to the i-th coordinate
    template <unsigned int n,unsigned int i, class T>
    ANALYTIC<n,T> projection(){
      static_assert(i <= n && i>0, "projection to invalid parameter");
      using vec = std::vector<unsigned long>;
      // the power series for the projection function
      // returns 1 iff v[i] == 1 and for all i != j v[j] == 0
      // oterwise 0
      auto prj = [] (const vec& v) -> T {
	if(v[i-1] != 1) return 0;
	for(int j = 0; j<n; j++){
	  if(j != i-1 && v[j] != 0) return 0;
	} 
	return 1;
      };  
      return ANALYTIC<n,T>(prj, 1, 1);
      
    }
  }
}


#endif
