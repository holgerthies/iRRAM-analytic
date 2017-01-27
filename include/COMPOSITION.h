/*-------------------------------------------------
* Class for composition
 ------------------------------------------------*/
#ifndef COMPOSITION_H
#define COMPOSITION_H
#include "ANALYTIC.h"
namespace iRRAM
{

  template<int v, class R, class... Args>
  class deriv_cache
  {
  private:
    using node_ptr = std::shared_ptr<Node<R, Args...>>;
    mutable std::vector<deriv_cache<v+1, R, Args...>> cache;
    node_ptr f;
  public:
    deriv_cache(const node_ptr f) : f(f) {}

    node_ptr get_f() const{
      return f;
    }

    node_ptr get_derivative(const tutil::n_tuple<sizeof...(Args)-v,size_t>& idx) const{
      int n = std::get<0>(idx);
      if(cache.size() == 0){
        cache.push_back(deriv_cache<v+1, R, Args...>(f));
      }
      int sz=cache.size();
      for(int i=sz; i<=n; i++)
      {
        auto fs=REAL(1)/REAL(i)*pderive(cache[i-1].get_f(), v, 1);
        simplify(fs);
        auto deriv = deriv_cache<v+1, R, Args...>(fs);
        cache.push_back(deriv);
      }
      return cache[n].get_derivative(tutil::tail(idx));
    }
  };

  template<class R, class... Args>
  class deriv_cache<sizeof...(Args),R,Args...> 
  {
  private:
    using node_ptr = std::shared_ptr<Node<R, Args...>>;
    node_ptr f;
  public:
    deriv_cache(const node_ptr& f) : f(f) {}
    
    node_ptr get_f() const{
      return f;
    }

    node_ptr get_derivative(const std::tuple<>& idx) const{
      return f;
    }
  };

  template <class R, class... Args>
  class COMPOSITION : public Node<R, Args...>{
  private:
    using node_ptr = std::shared_ptr<Node<R, Args...>>;
    using node_ptr1d = std::shared_ptr<Node<R,R>>;
    mutable std::unique_ptr<deriv_cache<0, R, Args...>> cache;
  public:
    node_ptr1d lhs;
    node_ptr rhs;
    COMPOSITION(const node_ptr1d& lhs, const node_ptr& rhs):
      lhs(lhs), rhs(rhs)
    {
    }

    R evaluate(const Args&... args) const override{
      return lhs->evaluate_cached(rhs->evaluate_cached(args...));
    }

    REAL get_r() const override {
      REAL r = rhs->get_r_cached();
      while(!positive(lhs->get_r_cached()-rhs->get_M_root(r), -10)){
        r /= 2;
      }
      
      return r;
    }
    
    REAL get_M(const REAL& r) const override {
      return lhs->get_M_root(rhs->get_M_cached(r));
    }

    std::string to_string() const override
    {
      return this->lhs->to_string()+"("+this->rhs->to_string()+")";
      
    }

    virtual size_t get_size() const override{
      return 1+lhs->get_size()+rhs->get_size();
    }


    void reset_visited() const override
    {
      if(this->visited){
        this->visited = false;
        lhs->reset_visited();
        rhs->reset_visited();
      }
    }

    int count_nodes() const override
    {
      if(!this->visited){
        this->visited = true;
        int n=1+lhs->count_nodes()+rhs->count_nodes();
        return n;
      }
      return 0;
    }

    R get_coefficient(const tutil::n_tuple<sizeof...(Args),size_t>& idx) const override
    {
      if(!cache)
        cache.reset(new deriv_cache<0,R,Args...>(compose(lhs, rhs)));
      auto d= cache->get_derivative(idx);
      std::vector<R> Z(sizeof...(Args));
      return d->evaluate(Z);
    }

    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::COMPOSITION;
    }
  };

  // composition operators
  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> compose(const std::shared_ptr<Node<R,R>>& lhs,const std::shared_ptr<Node<R,Args...>>& rhs)
  {
    return std::make_shared<COMPOSITION<R, Args...>>(lhs, rhs);
  }


} // namespace iRRAM


#endif
