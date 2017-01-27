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
      return minimum(this->lhs->get_r_cached(), this->rhs->get_r_cached());
    };

    REAL get_M(const REAL& r) const override {
      return this->lhs->get_M_cached(r)+this->rhs->get_M_cached(r);
    };

    R get_coefficient(const tutil::n_tuple<sizeof...(Args),size_t>& idx) const override
    {
      return this->lhs->get_coefficient_cached(idx)-this->rhs->get_coefficient_cached(idx);
    }
 
    
    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::SUBTRACTION;
    }
  };

  // member definitions 
  template <class R, class... Args>
  R SUBTRACTION<R, Args...>::evaluate(const Args&... args) const
  {
    return this->lhs->evaluate_cached(args...)-this->rhs->evaluate_cached(args...);
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
    return REAL(-1)*rhs+lhs;
  }

  // unary minus
  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> operator-(const std::shared_ptr<Node<R,Args...>>& node)
  {
    return REAL(-1)*node;
  }

  

} // namespace iRRAM


#endif
