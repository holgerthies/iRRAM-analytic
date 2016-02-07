/*-------------------------------------------------
* Class for subtraction
 ------------------------------------------------*/
#ifndef SUBTRACTION_H
#define SUBTRACTION_H
#include "ANALYTIC.h"
namespace iRRAM
{
  template <size_t d, class T>
  class SERIES_SUBTRACTION : public SERIES_BINARY_OPERATOR<d,T>
  {
    using SERIES_BINARY_OPERATOR<d,T>::SERIES_BINARY_OPERATOR;
  public:
    std::shared_ptr<POWERSERIES<d-1,T>> get_coeff(const unsigned long n) const override 
    {
      std::shared_ptr<SERIES_OPERATOR<d-1,T>> subtraction = std::make_shared<SERIES_SUBTRACTION<d-1, T>>((*this->lhs)[n], (*this->rhs)[n]);
      return get_series(subtraction);
    }

  };

  template<class T>
  class SERIES_SUBTRACTION<1,T> : public SERIES_BINARY_OPERATOR<1,T>
  {
    using SERIES_BINARY_OPERATOR<1,T>::SERIES_BINARY_OPERATOR;
  public:
    std::shared_ptr<T> get_coeff(const unsigned long n) const override
    {
      return std::make_shared<T>(this->lhs->get(n)-this->rhs->get(n));
    }
  };
  

  template <class R, class... Args>
  class SUBTRACTION : public BinaryNode<R, Args...>{
    using BinaryNode<R,Args...>::BinaryNode;
  public:
    R evaluate(const Args&... args) const override;
    std::shared_ptr<ANALYTIC<R,Args...>> to_analytic() const override;
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

  template <class R, class... Args>
  std::shared_ptr<ANALYTIC<R,Args...>> SUBTRACTION<R,Args...>::to_analytic() const
    {
      auto l = this->lhs->to_analytic();
      auto r = this->rhs->to_analytic();
      std::shared_ptr<SERIES_OPERATOR<sizeof...(Args), R>> subtraction= std::make_shared<SERIES_SUBTRACTION<sizeof...(Args), R>>(l->get_series(), r->get_series());
      auto add_pwr = get_series(subtraction);
      auto add_M = l->get_M()+r->get_M();
      auto add_r = minimum(l->get_r(), r->get_r());
      return std::make_shared<ANALYTIC<R, Args...>> (add_pwr, add_M, add_r);
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
