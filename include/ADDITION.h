/*-------------------------------------------------
* Class for addition
 ------------------------------------------------*/
#ifndef ADDITION_H
#define ADDITION_H
#include "ANALYTIC.h"
#include "POLYNOMIAL.h"
namespace iRRAM
{

  template <class R, class... Args>
  class ADDITION : public BinaryNode<R, Args...>{
    using BinaryNode<R,Args...>::BinaryNode;
  public:

    R evaluate(const Args&... args) const override;

    REAL get_r() const override {
      return minimum(this->lhs->get_r(), this->rhs->get_r());
    };

    REAL get_M(const REAL& r) const override {
      return this->lhs->get_M(r)+this->rhs->get_M(r);
    };
    
    R get_coefficient(const tutil::n_tuple<sizeof...(Args),size_t>& idx) const override
    {
      return this->lhs->get_coefficient(idx)+this->rhs->get_coefficient(idx);
    }

    std::string to_string() const override
    {
      return "("+this->lhs->to_string()+"+"+this->rhs->to_string()+")";
    }

    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::ADDITION;
    }

  };

  template <class R, class... Args>
  class SCALAR_ADDITION : public ScalarNode<R, Args...>{
    using ScalarNode<R,Args...>::ScalarNode;
  public:
    R evaluate(const Args&... args) const override;

    REAL get_r() const override {
      return this->node->get_r();
    };
    REAL get_M(const REAL& r) const override {
      return this->node->get_M(r)+abs(this->scalar);
    };

    R get_coefficient(const tutil::n_tuple<sizeof...(Args), size_t>& idx ) const override
    {
      auto ans = this->node->get_coefficient(idx);
      if(tutil::all_zero(idx))
        ans += this->scalar;
      return ans;
    }

    std::string to_string() const override
    {
      return "("+this->node->to_string()+"_+_"+std::to_string(this->scalar.as_double())+")";
      
    }

    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::SCALAR_ADDITION;
    }
  };

  // member definitions 
  template <class R, class... Args>
  R ADDITION<R, Args...>::evaluate(const Args&... args) const
  {
    return this->rhs->evaluate(args...)+this->lhs->evaluate(args...);
  }

  template <class R, class... Args>
  R SCALAR_ADDITION<R, Args...>::evaluate(const Args&... args) const
  {
    return this->node->evaluate(args...)+this->scalar;
  }

  // addition operators
  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> operator+(const std::shared_ptr<Node<R,Args...>>& lhs,const std::shared_ptr<Node<R,Args...>>& rhs)
  {
    return std::make_shared<ADDITION<R, Args...>>(lhs, rhs);
  }


  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> operator+(const std::shared_ptr<Node<R,Args...>>& lhs,const R& rhs)
  {
    return std::make_shared<SCALAR_ADDITION<R, Args...>>(lhs, rhs);
  }

  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> operator+(const R& lhs, const std::shared_ptr<Node<R,Args...>>& rhs)
  {
    return rhs+lhs;
  }

} // namespace iRRAM


#endif
