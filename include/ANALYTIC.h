/*-------------------------------------------------
 * Class for real analytic functions 
 ------------------------------------------------*/
#ifndef ANALYTIC_H
#define ANALYTIC_H
#include "POWERSERIES.h"
#include <vector>
#include <memory>
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
  template<unsigned int n, class T>
  ANALYTIC<n,T> operator/(const ANALYTIC<n,T>& lhs, const ANALYTIC<n,T>& rhs);
  template<unsigned int n, class T>
  ANALYTIC<n,T> inverse(const ANALYTIC<n,T>&);
  template<unsigned int n, class T>
  ANALYTIC<n,T> derive(const ANALYTIC<n,T>&, unsigned int i, unsigned int d);

  template <unsigned int n, class T>
  class ANALYTIC{
    template<unsigned int d, class F>
    friend class IVP;
    using vec = std::vector<unsigned long>;
    using pwr_ser = POWERSERIES<n,T>;
  private:
    REAL M,r; // maximum of the function and radius of convergence
    std::shared_ptr<pwr_ser> pwr;
  public:
    // constructor from everything that can be used to construct a power series
    template<typename F>
    ANALYTIC(F&& f, const REAL& r, const REAL& M ): 
      M(M),
      r(r),
      pwr(std::make_shared<POWERSERIES<n,T>>(POWERSERIES<n,T>(std::forward<F>(f))))
      {};


    // constructor from T (constant function)
    ANALYTIC(const T& x):
      M(abs(x)),
      r(1000),
      pwr(new POWERSERIES<n,T>(x))
      {};
    // arithmetic operations
    friend ANALYTIC operator+ <>(const ANALYTIC& lhs, const ANALYTIC& rhs);
    friend ANALYTIC operator- <>(const ANALYTIC& lhs, const ANALYTIC& rhs);
    friend ANALYTIC operator* <>(const ANALYTIC& lhs, const ANALYTIC& rhs);
    friend ANALYTIC operator/ <>(const ANALYTIC& lhs, const ANALYTIC& rhs);
    friend ANALYTIC inverse <>(const ANALYTIC&);
    friend ANALYTIC derive <>(const ANALYTIC<n,T>&, unsigned int i, unsigned int d);

    template<unsigned int d, class U>
    friend ANALYTIC<d,U> compose(const ANALYTIC<1,U>& f, const ANALYTIC<d,U>& rhs);

    template<typename... Ts>
    T operator ()(Ts... x) const;

    template<typename... ARGS>
    T get_coeff(ARGS... args);

    int get_known_coeffs() const {};

    REAL get_r() const{return r;}
    REAL get_M() const{return M;}

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

  // d-th derivative w.r.t. to the i-th variable
  template<unsigned int n, class T>
    ANALYTIC<n,T> derive(const ANALYTIC<n,T>& f, const unsigned int i, const unsigned int d){
    auto npwr=derivative(*f.pwr, i, d);
    auto r=f.get_r();
    auto M=f.get_M();
    auto newr = r/2;
    auto p = power(2, d+1);
    T fact=1;
    for(int j=2; j<=d; j++) fact *= d;
    auto newM = M*r*p*fact; 
    return ANALYTIC<n,T>(npwr,newr, newM );
  }

  template<unsigned int n, class T>
  ANALYTIC<n,T> inverse(const ANALYTIC<n,T>& f){
    auto npwr=inverse(*f.pwr);
    auto r=f.get_r();
    auto M=f.get_M();
    auto newr = r;
    // find r' such that the function does not have a zero for all z \in B_r'
    REAL lowerbound=-1;// lower bound for the function
    T f0 = constant_coefficient(*f.pwr);
    while(!positive(lowerbound, -10)){
      newr /= 2;
      lowerbound = abs(f0)-M*r*newr/(r-newr);
    }
    return ANALYTIC<n,T>(npwr, r, 1/lowerbound);
  }

  template<unsigned int n, class T>
  ANALYTIC<n,T> operator/(const ANALYTIC<n,T>& lhs, const ANALYTIC<n,T>& rhs){
    auto inv = inverse(rhs);
    return lhs*inv;
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

  template<unsigned int n, class T, typename... Ts>
  T evaluate(const POWERSERIES<n,T>& pwr, const REAL& M, const REAL& r, const T& x, Ts... rest){
    //iRRAM:: cout<<"evaluating dim" << n << std::endl;
    int stepsize=10;
    int J=0;
    REAL error_factor = abs(x)/r;
    REAL error = M*error_factor/(1-error_factor);
    //iRRAM::cout << "get first coefficient" << std::endl;
    //PWRSERIES_IMPL::print(pwr[0]);
    //iRRAM::cout << "got first coefficient" << std::endl;
    REAL sum(evaluate<n-1,T>(pwr[0], M, r, rest...));
    REAL best=sum;
    sizetype best_error, trunc_error, local_error,sum_error;
    sum.geterror(sum_error);
    sizetype_add(trunc_error,error.vsize,error.error);
    sizetype_add(local_error, sum_error, trunc_error);
    best.seterror(local_error);
    best_error = local_error;
    REAL x0=1; // x[0]^J
    error_factor = power(error_factor, stepsize);
    while (sizetype_less(sum_error, trunc_error) &&
	   (trunc_error.exponent >= ACTUAL_STACK.actual_prec) ){
      J+=stepsize;
      //cout << "J = " << J<< std::endl;
      // horner's method to evaluate polynomial 
      x0 *= x;
      T b = evaluate<n-1,T>(pwr[J], M, r, rest...);
      REAL x1=x0; // x^J
      for(int j=J-1; j>J-stepsize; j--){
	b = evaluate<n-1,T>(pwr[j], M, r, rest...) + b*x;
	x1 *= x;
      }
      sum += b*x0;
      error *= error_factor;
      sizetype_add(trunc_error,error.vsize,error.error);
      sum.geterror(sum_error); // get error of partial sum
      sizetype_add(local_error, sum_error, trunc_error);
      x0 = x1;
      if (sizetype_less(local_error, best_error)) { 
	best = sum;
	best.seterror(local_error);
	best_error = local_error;
      }
    } 
    return best;
    
  }

  template<unsigned int n, class T>
  T evaluate(const T& x, const REAL& M, const REAL& r){
    return x;
  }


  template<unsigned int n, class T>
  template<typename... ARGS>
  T ANALYTIC<n,T>::get_coeff(ARGS... args){
    return pwr->get(args...);
  }
  template<unsigned int n, class T>
  template<typename... Ts>
  T ANALYTIC<n,T>::operator ()(Ts... x) const
  {
    return evaluate<n,T>(*pwr, M,r, x...);
  }

  
    // function composition 
  template<unsigned int m, class T>
  ANALYTIC<m,T> compose(const ANALYTIC<1,T>& f, const ANALYTIC<m,T>& g){
    return ANALYTIC<m,T>(compose(*(f.pwr), *(g.pwr)), f.r, f.M);
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
