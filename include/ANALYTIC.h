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
      ANALYTIC(const sq_ptr f, const REAL& r, const REAL& M ): 
        M(M),
        r(r),
        pwr(new POWERSERIES<n,T>(f))
        {};
      ANALYTIC(seq f, const REAL& r, const REAL& M): ANALYTIC(sq_ptr(new seq(f)), r, M) {};

      T operator ()(const std::vector<T>& x) const
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

}
#endif
