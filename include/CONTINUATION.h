/*-------------------------------------------------
* Class for analytic continuation
 ------------------------------------------------*/
#ifndef CONTINUATION_H
#define CONTINUATION_H
#include "ANALYTIC.h"
namespace iRRAM
{
  // evaluate the series for f^(k)(x)/k!
  template<class T, class... Args>
  class DERIVATIVE_EVALUATOR
  {
  private:
    T x;
    DERIVATIVE_EVALUATOR<Args...> rec_evaluator;
  public:
    REAL get_q(const REAL& r)
    {
      return minimum((r+x)/2, rec_evaluator.get_q(r));
    }

    template<class... Orders>
    static REAL get_B(const REAL& r, const REAL& M,const REAL& q, const int k, const Orders... ks)
    {
      return power(r-q, -k)*get_B(r, M, q, ks...);
    }

    static REAL get_B(const REAL& r, const REAL& M,const REAL& q, const int k)
    {
      return M*power(r-q, -k);
    }

    DERIVATIVE_EVALUATOR(const T x0, const Args... xs):
      x(x0),
      rec_evaluator(DERIVATIVE_EVALUATOR<Args...>(xs...))
    {
    }

    DERIVATIVE_EVALUATOR(const std::vector<T>& xs):
      x(xs[xs.size()-sizeof...(Args)-1]),
      rec_evaluator(DERIVATIVE_EVALUATOR<Args...>(xs))
    {
    }

    REAL get_error_constant(const REAL& q){
      return (1-abs(x)/q)*rec_evaluator.get_error_constant(q);
    }
    
    template<class... Orders>
    T evaluate(const std::shared_ptr<POWERSERIES<sizeof...(Args)+1, T>>& pwr, const REAL& B, const REAL& q, const int k, const Orders... ks)
    {
      int J=0;
      REAL error_factor = abs(x)/q;
      REAL error = B*get_error_constant(q);
      REAL next_B = B;
      REAL sum(rec_evaluator.evaluate((*pwr)[k], next_B, q, ks...));
      REAL best=sum;
      sizetype best_error, trunc_error, local_error,sum_error;
      sum.geterror(sum_error);
      sizetype error_error, error_vsize;
      error.geterror(error_error);
      error.getsize(error_vsize);
      sizetype_add(trunc_error,error_vsize,error_error);
      sizetype_add(local_error, sum_error, trunc_error);
      best.seterror(local_error);
      best_error = local_error;
      REAL x0=1; // (j+k choose k)*x^j
      REAL next_B_factor = 1/q;
      while (sizetype_less(sum_error, trunc_error) &&
             (trunc_error.exponent >= ACTUAL_STACK.actual_prec) ){
        J++;
        next_B *= next_B_factor;
        x0 *= x*REAL(J+k)/J;
        sum += rec_evaluator.evaluate((*pwr)[J+k], next_B, q, ks...)*x0;
        error *= error_factor;
        error.geterror(error_error);
        error.getsize(error_vsize);
        sizetype_add(trunc_error,error_vsize,error_error);
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
  };

  template<class T>
  class DERIVATIVE_EVALUATOR<T>
  {
  private:
    T x;
  public:
    DERIVATIVE_EVALUATOR(const T x):
      x(x)
    {
    }

    DERIVATIVE_EVALUATOR(const std::vector<T>& xs):
      x(xs[xs.size()-1])
    {
    }
  public:
    REAL get_q(const REAL& r)
    {
      return (r+x)/2;
    }

    REAL get_error_constant(const REAL& q)
    {
      return (1-abs(x)/q);
    }

    template<class... Orders>
    static REAL get_B(const REAL& r, const REAL& M,const REAL& q, const int k, const Orders... ks)
    {
      return power(r-q, -k)*get_B(r, M, q, ks...);
    }

    static REAL get_B(const REAL& r, const REAL& M,const REAL& q)
    {
      return M;
    }

    T evaluate(const std::shared_ptr<POWERSERIES<1, T>>& pwr, const REAL& B, const REAL& q, const int k)
    {
      int J=0;
      REAL error_factor = abs(x)/q;
      REAL error = B*get_error_constant(q);
      REAL next_B = B;
      REAL sum(pwr->get(k));
      REAL best=sum;
      sizetype best_error, trunc_error, local_error,sum_error;
      sum.geterror(sum_error);
      sizetype error_error, error_vsize;
      error.geterror(error_error);
      error.getsize(error_vsize);
      sizetype_add(trunc_error,error_vsize,error_error);
      sizetype_add(local_error, sum_error, trunc_error);
      best.seterror(local_error);
      best_error = local_error;
      REAL x0=1; // (j+k choose k)*x^j
      REAL next_B_factor = 1/q;
      while (sizetype_less(sum_error, trunc_error) &&
             (trunc_error.exponent >= ACTUAL_STACK.actual_prec) ){
        J++;
        next_B *= next_B_factor;
        x0 *= x*REAL(J+k)/J;
        sum += pwr->get(J+k)*x0;
        error *= error_factor;
        error.geterror(error_error);
        error.getsize(error_vsize);
        sizetype_add(trunc_error,error_vsize,error_error);
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
  };
  

  template<size_t d, class T, class... Args>
  class CONTINUATION_SERIES
  {
  private:
    std::shared_ptr<POWERSERIES<d,T>> pwr;
  public:
    template<class... Orders>
    CONTINUATION_SERIES(const std::shared_ptr<ANALYTIC<T, Args...>>& node, const Args&... args, const Orders... orders)
    {
      pwr = std::make_shared<POWERSERIES<d,T>>([node, orders..., args...] (unsigned long n) {
          auto continuation = std::make_shared<CONTINUATION_SERIES<d-1, T, Args...>>(node, args..., orders..., n);
          return continuation->get_series();
        });
    };

    template<class... Orders>
    CONTINUATION_SERIES(const std::shared_ptr<ANALYTIC<T, Args...>>& node, const std::vector<T>& xs, const Orders... orders)
    {
      pwr = std::make_shared<POWERSERIES<d,T>>([node, orders..., &xs] (unsigned long n) {
          auto continuation = std::make_shared<CONTINUATION_SERIES<d-1, T, Args...>>(node, xs, orders..., n);
          return continuation->get_series();
        });
    };
    auto get_series() -> decltype(pwr)
    {
      return pwr;
    }

  };
  template<class T, class... Args>
  class CONTINUATION_SERIES<1,T, Args...>
  {
  private:
    std::shared_ptr<POWERSERIES<1,T>> pwr;
  public:
    template<class... Orders>
    CONTINUATION_SERIES(const std::shared_ptr<ANALYTIC<T, Args...>>& node, const Args&... args, const Orders... orders)
    {
      auto series = node->get_series();
      auto evaluator = std::make_shared<DERIVATIVE_EVALUATOR<Args...>>(args...);
      auto r = node->get_r();
      auto M = node->get_M();
      auto q = evaluator->get_q(r);
      pwr = std::make_shared<POWERSERIES<1,T>>([q,r,M, series, orders..., evaluator] (unsigned long n) {
          auto B = evaluator->get_B(r, M, q, orders...);
          return std::make_shared<T>(evaluator->evaluate(series, B, q, orders..., n));
        });
    };

    template<class... Orders>
    CONTINUATION_SERIES(const std::shared_ptr<ANALYTIC<T, Args...>>& node, const std::vector<T>& xs, const Orders... orders)
    {
      auto series = node->get_series();
      auto evaluator = std::make_shared<DERIVATIVE_EVALUATOR<Args...>>(xs);
      auto r = node->get_r();
      auto M = node->get_M();
      auto q = evaluator->get_q(r);
      pwr = std::make_shared<POWERSERIES<1,T>>([q,r,M, series, orders..., evaluator] (unsigned long n) {
          auto B = evaluator->get_B(r, M, q, orders...);
          return std::make_shared<T>(evaluator->evaluate(series, B, q, orders..., n));
        });
    };
    auto get_series() -> decltype(pwr)
    {
      return pwr;
    }

  };

  template <class R, class... Args>
  class CONTINUATION : public Node<R, Args...>{
  using node_ptr = std::shared_ptr<Node<R, Args...>>;
  private:
    node_ptr node;
    std::shared_ptr<ANALYTIC<R,Args...>> analytic;
  public:
    CONTINUATION(const node_ptr& node, const REAL& new_M, const REAL& new_r, Args... args):
      node(node)
    {
      auto f = node->to_analytic();
      auto dpwr= std::make_shared<CONTINUATION_SERIES<sizeof...(args),R,Args...>>(f,args...);
      this->analytic = std::make_shared<ANALYTIC<R,Args...>>(dpwr->get_series(), new_M, new_r);
    }
    CONTINUATION(const node_ptr& node, const REAL& new_M, const REAL& new_r, const std::vector<R>& xs):
      node(node)
    {
      auto f = node->to_analytic();
      auto dpwr= std::make_shared<CONTINUATION_SERIES<sizeof...(Args),R,Args...>>(f,xs);
      this->analytic = std::make_shared<ANALYTIC<R,Args...>>(dpwr->get_series(), new_M, new_r);
    }
    R evaluate(const Args&... args) const override{
      return analytic->evaluate(args...);
    };
    std::shared_ptr<ANALYTIC<R,Args...>> to_analytic() const override{
      return analytic;
    };

    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::CONTINUATION;
    }
  };

  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> continuation(const std::shared_ptr<Node<R,Args...>>& node, const REAL& new_M, const REAL& new_r, Args... args)
  {
    return std::make_shared<CONTINUATION<R, Args...>>(node, new_M, new_r, args...);
  }

  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> continuation(const std::shared_ptr<Node<R,Args...>>& node, const REAL& new_M, const REAL& new_r, const std::vector<R>& xs)
  {
    return std::make_shared<CONTINUATION<R, Args...>>(node, new_M, new_r, xs);
  }

} // namespace iRRAM


#endif
