/*-------------------------------------------------
* Class for derivative
 ------------------------------------------------*/
#ifndef DERIVATIVE_H
#define DERIVATIVE_H
#include "ANALYTIC.h"
#include "tutil.h"
namespace iRRAM
{

  template <class R, class... Args>
  class DERIVATIVE;

    
  template<size_t v, class R, class... Args>
    struct simplify_multiplication_helper
    {
      static std::shared_ptr<Node<R,Args...>> get(const std::shared_ptr<DERIVATIVE<R,Args...>>& f)
      {
        auto orders_t = f->orders_t;
        const int order = std::get<v>(orders_t);
        if(order == 0)
          return simplify_multiplication_helper<v+1,R,Args...>::get(f);
        auto child = std::dynamic_pointer_cast<MULTIPLICATION<R,Args...>>(f->node);
        auto ans = pderive(child->lhs, v, order)*child->rhs;
        for(int i=1; i<=order; i++){
          ans = ans+REAL(choose(order,i))*pderive(child->lhs, v, order-i)*pderive(child->rhs, v, i);
        }
        //ans = ans->simplify();
        return tderive(ans, tutil::tuple_replace<v>(orders_t, 0))->simplify();
      }
    };

  template<class R, class... Args>
  struct simplify_multiplication_helper<sizeof...(Args),R,Args...>
    {
      static std::shared_ptr<Node<R,Args...>> get(const std::shared_ptr<DERIVATIVE<R,Args...>>& f)
      {
        return f->node->simplify();
      }
      
    };
    
  
  
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
  private:
    node_ptr node;
    
    template<size_t I, typename D, typename... Ts>
    friend struct simplify_multiplication_helper;
    
    tutil::repeat<sizeof...(Args), std::tuple, size_t> orders_t;

    
  public:
    template<class... orders>
    DERIVATIVE(const node_ptr& node, orders... ods):
      node(node), orders_t(std::tuple<orders...>(ods...))
    {
      
    }

    REAL get_r() const override
    {
      return node->get_r();
    }

    REAL get_M(const REAL& r) const override
    {
      REAL ans=get_M_factor(node->get_r(), r, orders_t);
      int d =sizeof...(Args);
      ans *= power(node->get_r()*node->get_M(node->get_r()), d);
      return ans;
    }

    R get_coefficient(const tutil::n_tuple<sizeof...(Args),size_t>& idx) const override
    {
      auto t = tutil::sum_tuples(idx, orders_t);
      return get_coefficient_factor(orders_t, idx)*node->get_coefficient(t);
    }


    std::shared_ptr<Node<R,Args...>> simplify() const override
    {
      switch(this->node->get_type()){
      case ANALYTIC_OPERATION::ADDITION:
        {

          auto child = std::dynamic_pointer_cast<ADDITION<R,Args...>>(this->node);
          auto new_node = std::make_shared<ADDITION<R,Args...>>(tderive(child->lhs,orders_t), tderive(child->rhs, orders_t));
         return new_node->simplify();
        }
      case ANALYTIC_OPERATION::SUBTRACTION:
        {

          auto child = std::dynamic_pointer_cast<SUBTRACTION<R,Args...>>(this->node);
          auto new_node = std::make_shared<SUBTRACTION<R,Args...>>(tderive(child->lhs,orders_t), tderive(child->rhs, orders_t));
         return new_node->simplify();
        }
      case ANALYTIC_OPERATION::SCALAR_ADDITION:
        {

          auto child = std::dynamic_pointer_cast<SCALAR_ADDITION<R,Args...>>(this->node);
          auto new_node = tderive(child->node->simplify(),orders_t);
          return new_node->simplify();
        }
      case ANALYTIC_OPERATION::MULTIPLICATION:
        {
          
          return simplify_multiplication_helper<0, R, Args...>::get(std::make_shared<DERIVATIVE>(*this));
        }

      case ANALYTIC_OPERATION::SCALAR_MULTIPLICATION:
        {

          auto child = std::dynamic_pointer_cast<SCALAR_MULTIPLICATION<R,Args...>>(this->node);
          auto new_node = std::make_shared<SCALAR_MULTIPLICATION<R,Args...>>(tderive(child->node,orders_t), child->scalar);
         return new_node->simplify();
        } 

      // case ANALYTIC_OPERATION::COMPOSITION:
      //   {

      //     auto child = std::dynamic_pointer_cast<COMPOSITION>(this->node);
      //     auto new_node = std::make_shared<SCALAR_MULTIPLICATION>(std::make_shared<DERIVATIVE>(child->node->simplify(),ods..), this->scalar);
      //    return new_node->simplify();
      //   }
      case ANALYTIC_OPERATION::DERIVATIVE:{
          auto child = std::dynamic_pointer_cast<DERIVATIVE>(this->node);
          auto new_node = tderive(child->node, tutil::sum_tuples(child->orders_t, this->orders_t));
          return new_node->simplify();
          
      }
      case ANALYTIC_OPERATION::POLYNOMIAL:{
        auto child = std::dynamic_pointer_cast<poly_impl::POLY<R, Args...>>(this->node);
        auto new_node = std::make_shared<poly_impl::POLY<R, Args...>>(derive(*child, this->orders_t));
        return new_node;
          
          
      }
      case ANALYTIC_OPERATION::ANALYTIC:
        return std::make_shared<DERIVATIVE>(*this);
      default:
        return std::make_shared<DERIVATIVE>(*this);
      }

    };


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
