#include "iRRAM/core.h"
#include "ANALYTIC.h"
#include "combinatorics.h"
using std::vector;
using std::shared_ptr;
namespace iRRAM
{

  // evaluate a given powerseries with parameter l at point z
  // only works if z <= 1/2l
  REAL eval(shared_ptr<POWERSERIES<1,REAL>> series, const int& l, const REAL& z){
   using std::endl;
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
   while (false && sizetype_less(sum_error, trunc_error) &&
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
  }

}
