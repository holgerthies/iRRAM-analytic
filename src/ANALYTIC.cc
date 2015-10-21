#include "iRRAM/core.h"
#include "ANALYTIC.h"
#include "combinatorics.h"
using std::vector;
using std::shared_ptr;
namespace iRRAM
{
  // returns f^(d)(z)/d!
  REAL eval_derivative(shared_ptr<POWERSERIES<1,REAL>> series, const int l, const REAL& z, const int d) {
    // start with a certain number of coefficients and increase in each step
    int J=d+l-ACTUAL_STACK.actual_prec;
    // define powerseries for the derivative
    POWERSERIES<1,REAL> deriv=*series;
    if(d == 0) 
      deriv = *series;
    else
      deriv = shared_ptr<std::function<REAL(unsigned long)>>(new std::function<REAL(unsigned long)>([series, d] (unsigned int j) {
        return series->get(j+d)*choose(j+d, d);
        }));
    REAL sum(deriv.evaluate_partial({z},0, J)); // partial sum
    sizetype sum_error; // and its error
    sum.geterror(sum_error); // get error of first coefficient
    REAL best(sum);
    sizetype best_error, trunc_error, local_error;
    // error bound approximating the error made by considering
    // only J many coefficients (truncation error)
    // it is given by B*l^d*2^(l+1-J)*d*choose(J+d,d)
    REAL error_factor = power(l,d)*power(2, l+1)*(d+1);
    // initial estimate for the number of terms to be considered
    REAL error = error_factor*choose(J+d, d)*power(2, -J);
    sizetype_add(trunc_error,error.vsize,error.error);
    sizetype_add(local_error, sum_error, trunc_error);
    best_error=local_error;
    best.seterror(local_error); // and its error. 
    while (sizetype_less(sum_error, trunc_error) &&
        (best_error.exponent >= ACTUAL_STACK.actual_prec/2) ){
      sum += deriv.evaluate_partial({z}, J+1, 2*J); // partial sum
      J *= 2;
      error = error_factor*choose(J+d, d)*power(2, -J);
      sizetype_add(trunc_error,error.vsize,error.error);
      if(trunc_error.exponent > 1000){
        trunc_error.mantissa = 1;
        trunc_error.exponent = -10000;
      }
      sum.geterror(sum_error); // get error of partial sum
      sizetype_add(local_error, sum_error, trunc_error);
      if (sizetype_less(local_error, best_error)) { 
        best = sum;
        best.seterror(local_error);
        best_error = local_error;
      }
    } 
    return best;
  }
  REAL eval_derivative(const ANALYTIC& f,const int d, const REAL& z) {
     auto df = derive(f, d);
     return df.eval_k(z,0);
  }

  // evaluate a given powerseries with parameter l at point z
  // only works if z <= 1/2l
  REAL eval(shared_ptr<POWERSERIES<1,REAL>> series, const int& l, const REAL& z){
   int J = -ACTUAL_STACK.actual_prec+l+1;
   REAL error_factor = power(z*REAL(l),J);
   REAL error = power(2,l)*error_factor/(1-z*REAL(l));
   REAL sum(series->evaluate_partial({z},0,J));
   REAL best=sum;
   sizetype best_error, trunc_error, local_error,sum_error;
   sum.geterror(sum_error);
   sizetype_add(trunc_error,error.vsize,error.error);
   sizetype_add(local_error, sum_error, trunc_error);
   best.seterror(local_error);
   best_error = local_error;
   while (sizetype_less(sum_error, trunc_error) &&
        (best_error.exponent >= ACTUAL_STACK.actual_prec) ){
      sum += series->evaluate_partial({z}, J+1, 2*J); // partial sum
      J *= 2;
      error *= error_factor;
      error_factor *= error_factor;
      sizetype_add(trunc_error,error.vsize,error.error);
      sum.geterror(sum_error); // get error of partial sum
      sizetype_add(local_error, sum_error, trunc_error);
      if (sizetype_less(local_error, best_error)) { 
        best = sum;
        best.seterror(local_error);
        best_error = local_error;
      }
    } 
    return best;
  }

  ANALYTIC::ANALYTIC(shared_ptr<POWERSERIES<1, REAL>> pwr, const int l){
    // cover the interval with balls of radius 1/(2*l)
    REAL dist = REAL(1)/(2*l);
    REAL center = 0;
    REAL right=1;
    this->l = l;
    shared_ptr<POWERSERIES<1,REAL>> old_series = pwr;
    pwrs.push_back(old_series);
    REAL germ_position = center + dist;
    while (choose(germ_position < right,germ_position > right-dist/REAL(2)) == 1) {
      shared_ptr<POWERSERIES<1,REAL>> new_series(
          new POWERSERIES<1,REAL>(
            std::shared_ptr<std::function<REAL(unsigned long)>>(
              new std::function<REAL(unsigned long)>(
                [dist, old_series, l]
                (unsigned int i) -> REAL { 
                return eval_derivative(ANALYTIC(old_series, l), i,dist);
                })
              )
            )
          );
      
      pwrs.push_back(new_series);
      germ_position += dist;
      old_series = new_series;
    }
  }

  REAL ANALYTIC::operator() (const REAL& x) const{
		REAL center = 0;
    // this is the index of a germ that is close enough to x. 
		int position = ((x-center)*REAL(2*l)).as_INTEGER()-1;
		REAL diff = (x-center)-REAL(position)/REAL(2*l);// this is the distance of x to the center of the germ
    return eval(pwrs[position], l, diff);
  }

  REAL ANALYTIC::eval_k(const REAL& x, const int k){
    return eval(pwrs[k], l, x);
  }

  int ANALYTIC::get_known_coeffs(const int k) const{
    return pwrs[k]->known_coeffs();
  }

  int ANALYTIC::get_l() const{
    return l;
  }
  
  std::shared_ptr<POWERSERIES<1,REAL>> ANALYTIC::get_series(const int k) const{
    return pwrs[k];
  }

  ANALYTIC derive(const ANALYTIC& f, const int d){
   int new_l = f.get_l();
   return ANALYTIC(shared_ptr<POWERSERIES<1,REAL>>(new POWERSERIES<1,REAL>(derive(*f.get_series(0), d))), new_l);
  }

}
