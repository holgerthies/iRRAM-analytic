/*-------------------------------------------------
* Class for subtraction
 ------------------------------------------------*/
#ifndef SUBTRACTION_H
#define SUBTRACTION_H
#include "ANALYTIC.h"
namespace iRRAM
{

  template <class R, class... Args>
  class SUBTRACTION : public BinaryNode<R, Args...>{
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
      return this->lhs->get_coefficient(idx)-this->rhs->get_coefficient(idx);
    }
 
    
    std::shared_ptr<Node<R,Args...>> simplify() const override
    {
      return std::make_shared<SUBTRACTION>(this->lhs->simplify(), this->rhs->simplify());
    };

    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::SUBTRACTION;
    }
  };

  // member definitions 
  template <class R, class... Args>
  R SUBTRACTION<R, Args...>::evaluate(const Args&... args) const
  {
    return this->lhs->evaluate(args...)-this->rhs->evaluate(args...);
  }

  // subtraction operators
  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> operator-(const std::shared_ptr<Node<R,Args...>>& lhs,const std::shared_ptr<Node<R,Args...>>& rhs)
  {
    
    return std::make_shared<SUBTRACTION<R, Args...>>(lhs, rhs);
  }


  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> operator-(const std::shared_ptr<Node<R,Args...>>& lhs,const R& rhs)
  {
    
    return std::make_shared<SCALAR_ADDITION<R, Args...>>(lhs, -1*rhs);
  }

  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> operator-(const R& lhs, const std::shared_ptr<Node<R,Args...>>& rhs)
  {
    return -1*rhs+lhs;
  }

  // unary minus
  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> operator-(const std::shared_ptr<Node<R,Args...>>& node)
  {
    return -1*node;
  }

  

} // namespace iRRAM


#endif
