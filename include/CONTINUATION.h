/*-------------------------------------------------
* Class for analytic continuation
 ------------------------------------------------*/
#ifndef CONTINUATION_H
#define CONTINUATION_H
#include "ANALYTIC.h"
#include "POLYNOMIAL.h"
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

    REAL get_error_constant(const REAL& q) const{
      return (1-abs(x)/q)*rec_evaluator.get_error_constant(q);
    }
    
    T evaluate(const std::shared_ptr<POWERSERIES<sizeof...(Args)+1, T>>& pwr, const REAL& B, const REAL& q, const tutil::n_tuple<sizeof...(Args)+1, int>& ks) const
    {
      auto k = std::get<0>(ks);
      int J=0;
      REAL error_factor = abs(x)/q;
      REAL error = B*get_error_constant(q);
      REAL next_B = B;
      REAL sum(rec_evaluator.evaluate((*pwr)[k], next_B, q, tutil::tail(ks)));
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
        sum += rec_evaluator.evaluate((*pwr)[J+k], next_B, q, tutil::tail(ks))*x0;
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

    REAL get_error_constant(const REAL& q) const
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

    T evaluate(const std::shared_ptr<POWERSERIES<1, T>>& pwr, const REAL& B, const REAL& q, const std::tuple<int> ks) const
    {
      auto k = std::get<0>(ks);
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
  


  template <class R, class... Args>
  class CONTINUATION : public Node<R, Args...>{
  using node_ptr = std::shared_ptr<Node<R, Args...>>;
  public:
    node_ptr node;
    std::vector<R> center;
  private:
    REAL r,M, max_x;
    DERIVATIVE_EVALUATOR<Args...> evaluator;
  public:

    CONTINUATION(const node_ptr& node, const REAL& new_M, const REAL& new_r, const std::vector<R>& xs):
      node(node), center(xs), r(new_r), M(new_M), evaluator(xs)
    {
      max_x=0;
      for(auto x : xs) max_x = maximum(max_x, abs(x));
    }

    CONTINUATION(const node_ptr& node, const REAL& new_M, const REAL& new_r, Args... args):
      CONTINUATION(node, new_M, new_r, std::vector<R>{args...})
    {
    }


    std::string to_string() const override
    {
      return "CONTINUATION("+node->to_string()+")";
      
    }

    virtual size_t get_size() const override{
      return 1+node->get_size();
    }

    REAL get_r() const override {
      return r;
    };

    REAL get_M(const REAL& r) const override {
      return M;
    };

    void reset_visited() const override
    {
      if(this->visited){
        this->visited = false;
        node->reset_visited();
      }
    }

    int count_nodes() const override
    {
      if(!this->visited){
        this->visited = true;
        int n=1+node->count_nodes();
        return n;
      }
      return 0;
    }

    R get_coefficient(const tutil::n_tuple<sizeof...(Args),size_t>& idx) const override
    {
      return evaluator.evaluate(this->node->pwr, 2*max_x, this->node->get_M_root(2*max_x), idx );
    }

    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::CONTINUATION;
    }
  };

  template <class R, class... Args>
  class TRANSPOSITION : public Node<R, Args...>{
  using node_ptr = std::shared_ptr<Node<R, Args...>>;
  public:
    node_ptr node;
    std::vector<R> center;
  private:
    REAL max_x;
    DERIVATIVE_EVALUATOR<Args...> evaluator;
  public:

    TRANSPOSITION(const node_ptr& node, const std::vector<R>& xs):
      node(node), center(xs),  evaluator(xs)
    {
      max_x=0;
      for(auto x : xs) max_x = maximum(max_x, abs(x));
    }

    TRANSPOSITION(const node_ptr& node, Args... args):
      TRANSPOSITION(node, std::vector<R>{args...})
    {
    }

    void reset_visited() const override
    {
      if(this->visited){
        this->visited = false;
        node->reset_visited();
      }
    }

    int count_nodes() const override
    {
      if(!this->visited){
        this->visited = true;
        int n=1+node->count_nodes();
        return n;
      }
      return 0;
    }
    std::string to_string() const override
    {
      std::string ans = "TRANSPOSE("+node->to_string()+";";
      for(int i=0; i<center.size(); i++){
        ans += std::to_string(center[i].as_double());
        if(i != center.size()-1) ans += ", ";
      }
      ans += ")";
      return ans;
    }


    REAL get_r() const override {
      return node->get_r_cached()-max_x;
    };

    REAL get_M(const REAL& r) const override {
      return node->get_M_root(r+max_x);
    };

    R get_coefficient(const tutil::n_tuple<sizeof...(Args),size_t>& idx) const override
    {
      
      REAL r = minimum(0.9*this->node->get_r_cached(), maximum(max_x+1, 2*max_x));
      REAL M =node->get_M_root(r);
      
      return evaluator.evaluate(this->node->get_pwr(), r, M, idx );
    }

    virtual size_t get_size() const override{
      return 1+node->get_size();
    }

    R evaluate(const Args&... args) const override{
      std::vector<R> x(sizeof...(args));
      for(int i=0; i<x.size(); i++){
        x[i] = center[i]+tutil::get(i, args...);
      }
      return this->node->evaluate(x);
    };

    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::TRANSPOSITION;
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

  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> transpose(const std::shared_ptr<Node<R,Args...>>& node,  const std::vector<R>& xs)
  {
    return std::make_shared<TRANSPOSITION<R, Args...>>(node,  xs);
  }

} // namespace iRRAM


#endif
