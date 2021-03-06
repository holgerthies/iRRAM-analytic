/*-------------------------------------------------
* Class for derivative
 ------------------------------------------------*/
#ifndef DERIVATIVE_H
#define DERIVATIVE_H
#include "ANALYTIC.h"
#include "COMPOSITION.h"
#include "tutil.h"
namespace iRRAM
{

  template <class R, class... Args>
  class DERIVATIVE;

  template<class... Args>
  REAL factorial(const std::tuple<Args...>& orders)
  {
    return REAL(factorial(std::get<0>(orders)))*factorial(tutil::tail(orders));
  }

  template<>
  REAL factorial(const std::tuple<>& orders)
  {
    return 1;
  }
    
  
  template<class... Params>
  REAL get_coefficient_factor(const std::tuple<Params...>& orders, const std::tuple<Params...>& idx)
  {
    REAL ans=1;
    for(int j=std::get<0>(idx)+1; j<=std::get<0>(idx)+std::get<0>(orders); j++)
    {
      ans *= j;
    }
    return ans*get_coefficient_factor(tutil::tail(orders), tutil::tail(idx));
  }

  template<>
  REAL get_coefficient_factor(const std::tuple<>& orders, const std::tuple<>& idx)
  {
    return 1;
  }

  template<class... Params>
  REAL get_M_factor(const REAL& r, const REAL& new_r, const std::tuple<Params...>& orders)
  {
    auto n = std::get<0>(orders);
    return REAL(factorial(n))*get_M_factor(r, new_r, tutil::tail(orders))*power(r-new_r,-(n+1));
  }

  template<>
  REAL get_M_factor(const REAL& r, const REAL& new_r,const std::tuple<>& orders)
  {
    return 1;
  }

  template <class R, class... Args>
  class DERIVATIVE : public Node<R, Args...>{
    using node_ptr = std::shared_ptr<Node<R, Args...>>;
  public:
    node_ptr node;
    tutil::repeat<sizeof...(Args), std::tuple, size_t> orders_t;

    template<class... orders>
    DERIVATIVE(const node_ptr& node, orders... ods):
      node(node), orders_t(std::tuple<orders...>(ods...))
    {
      
    }

    REAL get_r() const override
    {
      
      return node->get_r_cached();
    }

    REAL get_M(const REAL& r) const override
    {
      
      auto sum = tutil::accumulate_tuple(orders_t);
      REAL fact = factorial(orders_t);
      
      std::vector<REAL> test_rs = {(r+node->get_r_cached())/2, r+2, 2*r, 1.1*r};
      REAL M =node->get_M_root(test_rs[0])/power(test_rs[0]-r, sum);
      if(positive(M-10000, 0))
      {
        for(int i=1; i<test_rs.size(); i++){
          REAL r1 = minimum(test_rs[i], test_rs[0]);
          M = minimum(M, node->get_M_root(r1)/power(r1-r, sum));
        }
      }
      return fact*M;
      
    }

    R get_coefficient(const tutil::n_tuple<sizeof...(Args),size_t>& idx) const override
    {
      auto t = tutil::sum_tuples(idx, orders_t);
      return get_coefficient_factor(orders_t, idx)*node->get_coefficient_cached(t);
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
      std::string ans="DERIVE("+this->node->to_string()+";"+tutil::to_string(this->orders_t)+")";
      return ans;
    }

    

    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::DERIVATIVE;
    }
  };




  // derivative operators
  template <class R, class... Args, class... orders>
  std::shared_ptr<Node<R, Args...>> derive(const std::shared_ptr<Node<R,Args...>>& node,const orders... ods)
  {
    
    return std::make_shared<DERIVATIVE<R, Args...>>(node, ods...);
  }

  template<int dim>
  struct vderive
  {
    template<class R, class... Args, class... orders>
    static std::shared_ptr<Node<R, Args...>> get(const std::shared_ptr<Node<R,Args...>>& node,const int variable, const int order, const orders... ods)
    {
      if(sizeof...(Args) - dim == variable)
        return vderive<dim-1>::get(node, variable, order, ods..., order);
      return vderive<dim-1>::get(node, variable, order, ods..., 0);
    }
  };

  template<>
  struct vderive<0>
  {
    template<class R, class... Args, class... orders>
    static std::shared_ptr<Node<R, Args...>> get(const std::shared_ptr<Node<R,Args...>>& node,const int variable, const int order, const orders... ods)
    {
      return derive(node, ods...);
    }
  };
  
    
  template<class R, class... Args>
  std::shared_ptr<Node<R, Args...>> pderive(const std::shared_ptr<Node<R,Args...>>& node,const int variable, const int order)
  {
    if(order == 0) return node;
    return vderive<sizeof...(Args)>::get(node, variable, order);
    
  }

  template<size_t n, class R, class... Args>
  struct tderive_helper
  {
    template<class... orders>
    static std::shared_ptr<Node<R,Args...>> get(const std::shared_ptr<Node<R, Args...>>& node, const tutil::n_tuple<sizeof...(Args), int> tuple, const orders... ods)
    {
      return tderive_helper<n-1, R, Args...>::get(node,tuple, ods...,std::get<sizeof...(Args)-n>(tuple));
    }
  };

  template<class R, class... Args>
  struct tderive_helper<0, R, Args...>
  {
    template<class... orders>
    static std::shared_ptr<Node<R,Args...>> get(const std::shared_ptr<Node<R, Args...>>& node, const tutil::n_tuple<sizeof...(Args), int> tuple, const orders... ods)
    {
      return std::make_shared<DERIVATIVE<R, Args...>>(node, ods...);
    }
  };

  // derivative from std::tuple
  template<class R, class... Args>
  std::shared_ptr<Node<R,Args...>> tderive(const std::shared_ptr<Node<R, Args...>>& node, const tutil::n_tuple<sizeof...(Args),int>& tuple)
  {
    return tderive_helper<sizeof...(Args), R, Args...>::get(node, tuple);
    
  }
  
} // namespace iRRAM


#endif
