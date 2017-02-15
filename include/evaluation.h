#ifndef EVALUATION_H
#define EVALUATION_H
#include "Node.h"
#include "SUBSTITUTION.h"
#include "AAREAL.h"
namespace iRRAM{
  
  template<class R, typename... Args>
  struct evaluation_helper
  {
    static R evaluate(const std::shared_ptr<const Node<R,R, Args...>>& f, const R& x, const Args&... xs)
    {
      return evaluation_helper<Args...>::evaluate(substitute(f, x),xs...);
    }
    static R approximate(const std::shared_ptr<const Node<R,R, Args...>>& f, const unsigned int max_coeffs, const R& x, const Args&... xs)
    {
      return evaluation_helper<Args...>::approximate(substitute(f, x),max_coeffs, xs...);
    }
  };
  
    
// evaluate one-dimensional analytic function
  template<class R>
  struct evaluation_helper<R>
  {
    static R evaluate(const std::shared_ptr<const Node<R,R>>& f, const R& x){
      // choose value q
      
      const REAL r = f->get_r_cached();
      
      REAL ax = abs(x);
      REAL q = maximum((ax+r)/2, 0.99*r); //minimum((ax+r)/2, maximum(ax+1, 2*ax));
      REAL B = f->get_M_root(q);
      REAL t = 0.99*r-ax;
      while(positive(log(B)-70, -1)){
        t/=2;
        q = ax+t;
        B = f->get_M_root(q);
      }
      //std::cout << x.as_double() << " " << q.as_double() << " " << B.as_double()<< "\n";
      single_valued code;
      int J=0;
      REAL error_factor = ax/q;
      REAL error = B*error_factor/(1-error_factor);
    
      AAREAL sum(f->get_coefficient_cached(0));
      REAL best=sum.to_real();
      sizetype best_error, trunc_error, local_error,sum_error;
      best.geterror(sum_error);
      
      AAREAL trunc_error_term;
      trunc_error_term.add_error(error);
      trunc_error = real_to_error(error);
      best = (sum+trunc_error_term).to_real();
      best.geterror(best_error);
      AAREAL x0=REAL(1); // x[0]^J
      //REAL error_factor_step = power(error_factor, stepsize);

      while (
             (trunc_error.exponent >= ACTUAL_STACK.actual_prec) ){
        J++;
        //std::cout << J << " " << trunc_error.exponent << "\n";
        REAL b = f->get_coefficient_cached(J);
        sizetype b_error;
        b.geterror(b_error);
        x0 = x0*x;
        sum += b*x0;
        sum.to_real().geterror(sum_error);
        if(J % 4 == 0)
        {
          x0.clean();
          sum.clean();
        }
        error *= error_factor;
        trunc_error_term = AAREAL();
        trunc_error_term.add_error(error);
        trunc_error = real_to_error(error);
        
        REAL local = (sum+trunc_error_term).to_real();
        local.geterror(local_error);
        
        if (sizetype_less(local_error, best_error)) { 
          best = local;
          best_error = local_error;
          // std::cout << sum_error.exponent << " " << trunc_error.exponent << std::endl;
          // std::cout << local_error.exponent  << std::endl;
          
        }
      }
      
      return best;
    }

    static R approximate(const std::shared_ptr<const Node<R,R>>& f, const unsigned int num_coeffs, const R& x){
      REAL sum(0);
      REAL x0=REAL(1); // x[0]^J
      for(int J=0; J<num_coeffs; J++){
        REAL b = f->get_coefficient_cached(J);
        sum += b*x0;
        x0 = x0*x;
      }
      return sum;
    }
  };

  template<class R, class... Args>
  R evaluate_pwr(const std::shared_ptr<const Node<R,Args...>>& f, const Args&... xs)
  {
    return evaluation_helper<Args...>::evaluate(f,xs...);
  }

  template<class R, class... Args>
  R approximate_pwr(const std::shared_ptr<const Node<R,Args...>>& f, const unsigned int num_coeffs, const Args&... xs)
  {
    return evaluation_helper<Args...>::approximate(f,num_coeffs, xs...);
  }
}

#endif
