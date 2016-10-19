/*-------------------------------------------------
* Class for multiplication
 ------------------------------------------------*/
#ifndef MULTIPLICATION_H
#define MULTIPLICATION_H
#include "ANALYTIC.h"
namespace iRRAM
{

  
  template<size_t n, class R,class... Args>
  struct cauchy_product
  {
    using nodeptr = std::shared_ptr<Node<R, Args...>>;
    
    nodeptr lhs,rhs;
    cauchy_product(const nodeptr& lhs, const nodeptr& rhs) :
      lhs(lhs), rhs(rhs) 
    {
    };
    
    R get(const tutil::n_tuple<n,size_t>& idx, const tutil::n_tuple<sizeof...(Args)-n,size_t>& idxl, const tutil::n_tuple<sizeof...(Args)-n,size_t>& idxr)
    {
      R ans=0;
      auto p = cauchy_product<n-1, R, Args...>(lhs,rhs);
      
      for(int i=0; i<=std::get<0>(idx); i++){
        ans += p.get(tutil::tail(idx), std::tuple_cat(idxl,std::make_tuple(i)), std::tuple_cat(idxr, std::make_tuple(std::get<0>(idx)-i)));
      }
      return ans;
    }
  };
  
  template<class R,class... Args>
  struct cauchy_product<0,R, Args...>
  {
    using nodeptr = std::shared_ptr<Node<R, Args...>>;
    
    nodeptr lhs,rhs;
    cauchy_product(const nodeptr& lhs, const nodeptr& rhs) :
      lhs(lhs), rhs(rhs) 
    {
    };
    
    R get(const std::tuple<>& idx, const tutil::n_tuple<sizeof...(Args),size_t>& idxl, const tutil::n_tuple<sizeof...(Args),size_t>& idxr)
    {
      return lhs->get_coefficient(idxl)*rhs->get_coefficient(idxr);
    }
  };

  template <class R, class... Args>
  class MULTIPLICATION : public BinaryNode<R, Args...>{
    using BinaryNode<R,Args...>::BinaryNode;
  public:
    R evaluate(const Args&... args) const override{
      return this->lhs->evaluate(args...)*this->rhs->evaluate(args...);
    };

    REAL get_r() const override {
      return minimum(this->lhs->get_r(), this->rhs->get_r());
    };

    REAL get_M(const REAL& r) const override {
      return this->lhs->get_M(r)*this->rhs->get_M(r);
    };

    R get_coefficient(const tutil::n_tuple<sizeof...(Args),size_t>& idx) const override
    {
      return cauchy_product<sizeof...(Args), R, Args...>(this->lhs, this->rhs).get(idx, std::tuple<>(), std::tuple<>());
    }


    std::string to_string() const override
    {
      return "("+this->lhs->to_string()+"x"+this->rhs->to_string()+")";
    }

    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::MULTIPLICATION;
    }
  };

  template <class R, class... Args>
  class SCALAR_MULTIPLICATION : public ScalarNode<R, Args...>{
    using ScalarNode<R,Args...>::ScalarNode;
  public:
    R evaluate(const Args&... args) const override
    {
      return this->scalar*this->node->evaluate(args...);
    }
    
    virtual std::string to_string() const override
    {
      return "("+this->node->to_string()+"_*_"+std::to_string(this->scalar.as_double())+")";
      
    }
    REAL get_r() const override {
      return this->node->get_r();
    };
    REAL get_M(const REAL& r) const override {
      return this->node->get_M(r)*abs(this->scalar);
    };

    R get_coefficient(const tutil::n_tuple<sizeof...(Args),size_t>& idx) const override
    {
      return this->scalar * this->node->get_coefficient(idx);
    }

    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::SCALAR_MULTIPLICATION;
    }
  };


  // multiplication operators
  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> operator*(const std::shared_ptr<Node<R,Args...>>& lhs,const std::shared_ptr<Node<R,Args...>>& rhs)
  {
    
    return std::make_shared<MULTIPLICATION<R, Args...>>(lhs, rhs);
  }

  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> operator*(const std::shared_ptr<Node<R,Args...>>& lhs,const R& rhs)
  {
    return std::make_shared<SCALAR_MULTIPLICATION<R, Args...>>(lhs, rhs);
  }

  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> operator*(const R& lhs, const std::shared_ptr<Node<R,Args...>>& rhs)
  {
    return rhs*lhs;
    
  }
  

} // namespace iRRAM


#endif
