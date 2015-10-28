/*-------------------------------------------------
 * Class for real analytic functions on [0,1]
 ------------------------------------------------*/
#ifndef ANALYTIC_H
#define ANALYTIC_H
#include "POWERSERIES.h"
#include <vector>
namespace iRRAM
{
  template <unsigned int n, class T>
  class ANALYTIC{
    using seq = std::function<T(const std::vector<unsigned long>&)>;
    using sq_ptr = std::shared_ptr<seq>;
    private:
      REAL M,r; // maximum of the function and radius of convergence
      std::shared_ptr<POWERSERIES<n,T>> pwr;
    public:
      ANALYTIC(sq_ptr f, const REAL& r, const REAL& M ): 
        M(M),
        r(r),
        pwr(new POWERSERIES<n,T>(f))
        {};
      ANALYTIC(seq f, const REAL& r, const REAL& M): ANALYTIC(sq_ptr(new seq(f)), r, M) {};

      T operator ()(const std::vector<T>& x) const
      {
        int J = -ACTUAL_STACK.actual_prec;
        T x_max = *std::max_element(x.begin(), x.end());
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
             (best_error.exponent >= ACTUAL_STACK.actual_prec) ){
           int new_J = J+step_size;
           sum += evaluate_partial(pwr, x, J+1, new_J); // partial sum
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

      };

      int get_known_coeffs() const {};

      REAL get_r() const{return r;}
      REAL get_M() const{return M;}

      T get_coeff(const std::vector<unsigned long>& v) const
      {
        return pwr->get(v);
      }
  };

  template <class T>
  REAL a(const ANALYTIC<2,T>& F, const std::vector<unsigned long>& v);

  template <class T>
  REAL ai(const ANALYTIC<2,T>& F, const unsigned long n, const unsigned long i){
    if(i == 0 && n==0) return 1;
    if(i==0) return 0;
    REAL ans=0;
    for(unsigned long j=0; j<=n;j++){
        ans+=ai(F,n-j, i-1)*a(F,{j});
    } 
    return ans;
  }

  template <class T>
  REAL a(const ANALYTIC<2,T>& F, const std::vector<unsigned long>& v){
    static std::vector<REAL> cache;
    REAL ans=0;
    auto l = v[0];
    //std::cout << l <<std::endl;
    if(l==0) return REAL(0);
    if(cache.size() > l) return cache[l];
    a(F, {l-1});
    cache.resize(l+1);
    for(auto w : partitions(l-1, 2)){
      auto n = w[0];
      auto k = w[1];
      for(unsigned long i=0; i<=n; i++)
      {
        ans += F.get_coeff({k,i})*ai(F,n,i);
      }
    }
    cache[l] = ans/REAL((int)l);
    return cache[l];
  };

  template <class T>
  ANALYTIC<1,T> solve(const ANALYTIC<2,T>& F){
    using std::vector;
    using seq = std::function<T(const std::vector<unsigned long>&)>;
    seq series = [F] (const vector<unsigned long>& v) {
      return a(F, v);
    };
    return ANALYTIC<1,T>(series, F.get_r(), F.get_M());
  }
}
#endif
