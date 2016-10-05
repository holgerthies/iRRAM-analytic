/*-------------------------------------------------
* Class for derivative
 ------------------------------------------------*/
#ifndef DERIVATIVE_H
#define DERIVATIVE_H
#include "ANALYTIC.h"
#include "tutil.h"
namespace iRRAM
{
  template <size_t d, class T, class... orders>
  std::shared_ptr<POWERSERIES<d,T>> get_derivative_series(const std::shared_ptr<POWERSERIES<d,T>>&, const int, const int, const orders...);

  template <size_t d, class T>
  std::shared_ptr<POWERSERIES<d,T>> get_derivative_series(const std::shared_ptr<POWERSERIES<d,T>>&, const int);

  template <size_t d, class T>
  class SERIES_DERIVATIVE : public SERIES_OPERATOR<d,T>
  {
  private:
    std::shared_ptr<POWERSERIES<d,T>> pwr;
    unsigned int variable;
    int order;
  public:
    SERIES_DERIVATIVE(const std::shared_ptr<POWERSERIES<d,T>>& pwr, const unsigned int variable, const int order ):
      pwr(pwr), variable(variable), order(order) 
    {
    }

    std::shared_ptr<POWERSERIES<d-1,T>> get_coeff(const unsigned long n) const override 
    {
      if(variable >= d || order==0)
        return (*pwr)[n];
      if(variable == 0){
        T fact=1;
        for(int k=n+1; k<=n+order; k++) fact *= T(k);
        std::shared_ptr<SERIES_OPERATOR<d-1,T>> multiplication = std::make_shared<SERIES_SCALAR_MULTIPLICATION<d-1,T>>((*pwr)[order+n], fact);
        return get_series(multiplication);
      }
      std::shared_ptr<SERIES_OPERATOR<d-1,T>> ds = std::make_shared<SERIES_DERIVATIVE<d-1,T>>((*pwr)[n], variable-1, order);
      return get_series(ds);
    }

  };


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
    
  
  
    
  template <class T>
  class SERIES_DERIVATIVE<1,T> : public SERIES_OPERATOR<1,T>
  {
  private:
    std::shared_ptr<POWERSERIES<1,T>> pwr;
    unsigned int variable;
    int order;
  public:
    SERIES_DERIVATIVE(const std::shared_ptr<POWERSERIES<1,T>>& pwr, const unsigned int variable, const int order ):
      pwr(pwr), variable(variable), order(order) 
    {
      
    }

    std::shared_ptr<T> get_coeff(const unsigned long n) const override 
    {
      if(variable != 0 || order==0)
        return (*pwr)[n];
      T fact=1;
      for(int k=n+1; k<=n+order; k++) fact *= T(k);
      T multiplication = pwr->get(order+n)*fact;
      return std::make_shared<T>(multiplication);
    }

  };

  REAL get_deriv_M_factor(const REAL& r, const REAL& r0)
  {
    return 1;
  }

  template<class... orders>
  REAL get_deriv_M_factor(const REAL& r,const REAL& r0, int d, orders... rest)
  {
    auto p = power(2/(r-r0), d);
    REAL fact=1;
    for(int j=2; j<=d; j++) fact *= j;
    return fact*p*get_deriv_M_factor(r, r0, rest...);
  }


  template <class R, class... Args>
  class DERIVATIVE : public Node<R, Args...>{
    using node_ptr = std::shared_ptr<Node<R, Args...>>;
  private:
    node_ptr node;
    template<size_t I, typename D, typename... Ts>
    friend struct simplify_multiplication_helper;
    
    std::function<REAL(const REAL&, const Args&...)> M_function;
    std::shared_ptr<ANALYTIC<R,Args...>> analytic;
    tutil::repeat<sizeof...(Args), std::tuple, int> orders_t;
    
    
  public:
    template<class... orders>
    DERIVATIVE(const node_ptr& node, orders... ods):
      node(node)
    {
      orders_t = std::tuple<orders...>(ods...);
      auto f = node->to_analytic();
      auto dpwr= get_derivative_series<sizeof...(Args), R>(f->get_series(),0, ods...);
      auto new_r = f->get_r()/2;
      auto new_M =  f->get_M()*get_deriv_M_factor(f->get_r(),0, ods...);
      this->M_function = [node, ods...] (const REAL& r, const Args&... args) {
        return node->get_M(r+node->get_r()/2,args...)*get_deriv_M_factor(node->get_r(), r, ods...);
      };
      this->analytic = std::make_shared<ANALYTIC<R,Args...>>(dpwr, new_M, new_r);
      
    }

    REAL get_r() const override {
      return node->get_r();
    };
    REAL get_M(const REAL& r, const Args&... args) const override {
      return M_function(r,args...);
    };

    R evaluate(const Args&... args) const override;
    std::shared_ptr<ANALYTIC<R,Args...>> to_analytic() const override;

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

  // member definitions 
  template <size_t d, class T, class... orders>
  std::shared_ptr<POWERSERIES<d,T>> get_derivative_series(const std::shared_ptr<POWERSERIES<d,T>>& pwr,const int variable, const int order, const orders... rest) 
  {
    std::shared_ptr<SERIES_OPERATOR<d, T>> derivative = std::make_shared<SERIES_DERIVATIVE<d,T>>(pwr, variable, order);
    auto dpwr = get_series(derivative);
    return get_derivative_series<d,T>(dpwr,variable+1, rest...);
  }

  template <size_t d, class T>
  std::shared_ptr<POWERSERIES<d,T>> get_derivative_series(const std::shared_ptr<POWERSERIES<d,T>>& pwr, const int variable)
  {
    return pwr;
  }

  template <class R, class... Args>
  R DERIVATIVE<R, Args...>::evaluate(const Args&... args) const
  {
    return this->analytic->evaluate(args...);
  }


  template <class R, class... Args>
  std::shared_ptr<ANALYTIC<R,Args...>> DERIVATIVE<R,Args...>::to_analytic() const
    {
      return analytic;
    }

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
