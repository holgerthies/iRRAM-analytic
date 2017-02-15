/*-------------------------------------------------
* Class for variable substitution
 ------------------------------------------------*/
#ifndef SUBSTITUTION_H
#define SUBSTITUTION_H
#include "ANALYTIC.h"
#include "tutil.h"
namespace iRRAM
{

  template<class R, class... Args>
  std::shared_ptr<Node<R, R>> fix_rest(const std::shared_ptr<Node<R, Args...>>& f, const tutil::n_tuple<sizeof...(Args)-1, size_t>& params);

  
  template <class R, class... Args>
  class FIXFIRST : public Node<R, Args...>{
    using node_ptr = std::shared_ptr<Node<R, R, Args...>>;
  public:
    node_ptr node;
    size_t i;

    FIXFIRST(const node_ptr& node, const size_t& i):
      node(node), i(i)
    {
      
    }

    REAL get_r() const override
    {
      return node->get_r_cached();
    }

    REAL get_M(const REAL& r) const override
    {
      REAL r2 = minimum(2*r, (r+node->get_r_cached())/2);
      return node->get_M_root(r2)/(power(r2, i)*power(1-r/r2, sizeof...(Args)));
    }

    R get_coefficient(const tutil::n_tuple<sizeof...(Args),size_t>& idx) const override
    {
      return node->get_coefficient_cached(std::tuple_cat(std::make_tuple(i), idx));
    }

    virtual size_t get_size() const override{
      return 1+node->get_size();
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
      std::string ans="FIX_FIRST("+this->node->to_string()+";"+std::to_string(this->i)+")";
      return ans;
    }

    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::FIX;
    }
  };

  template <class R, class... Args>
  class FIXREST : public Node<R, R>{
    using node_ptr = std::shared_ptr<const Node<R, Args...>>;
    using ntuple =tutil::n_tuple<sizeof...(Args)-1,size_t>;
  public:
    node_ptr node;
    ntuple params;

    FIXREST(const node_ptr& node, const ntuple& params):
      node(node), params(params)
    {
      
    }

    REAL get_r() const override
    {
      return node->get_r_cached();
    }

    REAL get_M(const REAL& r) const override
    {
      
      auto sum = tutil::accumulate_tuple(params);
      REAL r2 = minimum(2*r, (r+node->get_r_cached())/2);
      return node->get_M_root(r2)/(power(r/r2, sum)*(1-r/r2));
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

    R get_coefficient(const std::tuple<size_t>& idx) const override
    {
      return node->get_coefficient_cached(std::tuple_cat(idx, params));
    }

    virtual size_t get_size() const override{
      return 1+node->get_size();
    }

    std::string to_string() const override
    {
      std::string ans="FIX_REST("+this->node->to_string()+";"+tutil::to_string(this->params)+")";
      return ans;
    }

    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::FIXREST;
    }
  };
  template <class R, class... Args>
  class SUBSTITUTION : public Node<R, Args...>{
    using node_ptr = std::shared_ptr<const Node<R, R, Args...>>;
  public:
    node_ptr node;
    R x0;

    template<class... orders>
    SUBSTITUTION(const node_ptr& node, const R& x0):
      node(node), x0(x0)
    {
      
    }

    REAL get_r() const override
    {
      return node->get_r_cached();
    }

    REAL get_M(const REAL& r) const override
    {
      return node->get_M_root(maximum(r, x0));
    }

    R get_coefficient(const tutil::n_tuple<sizeof...(Args),size_t>& idx) const override
    {
      auto f_idx = fix_rest(node, idx);
      return f_idx->evaluate_root(x0);
    }

    R evaluate(const Args&... args) const override
    {
      return node->evaluate_cached(x0, args...);
    }

    virtual size_t get_size() const override{
      return 1+node->get_size();
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
      std::string ans="SUBSTITUTE("+this->node->to_string()+";"+std::to_string(this->x0.as_double())+")";
      return ans;
    }

    

    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::SUBSTITUTION;
    }
  };

  template<class R, class... Args>
  std::shared_ptr<Node<R, Args...>> fix_first(const std::shared_ptr<Node<R, R, Args...>>& f, const size_t i)
  {
    return std::make_shared<FIXFIRST<R,Args...>>(f,i);
  }

  template<class R, class... Args>
  std::shared_ptr<Node<R, R>> fix_rest(const std::shared_ptr<const Node<R, Args...>>& f, const tutil::n_tuple<sizeof...(Args)-1, size_t>& params)
  {
    return std::make_shared<FIXREST<R,Args...>>(f,params);
  }
  
  template<class R, class... Args>
  std::shared_ptr<Node<R, Args...>> substitute(const std::shared_ptr<const Node<R, R, Args...>>& f, const R& x)
  {
    return std::make_shared<SUBSTITUTION<R,Args...>>(f,x);
  }
  
} // namespace iRRAM


#endif
