/*-------------------------------------------------
* Class for addition
 ------------------------------------------------*/
#ifndef ADDITION_H
#define ADDITION_H
#include "ANALYTIC.h"
namespace iRRAM
{
  template <size_t d, class T>
  class SERIES_ADDITION : public SERIES_BINARY_OPERATOR<d,T>
  {
    using SERIES_BINARY_OPERATOR<d,T>::SERIES_BINARY_OPERATOR;
  public:
    std::shared_ptr<POWERSERIES<d-1,T>> get_coeff(const unsigned long n) const override 
    {
      std::shared_ptr<SERIES_OPERATOR<d-1,T>> addition = std::make_shared<SERIES_ADDITION<d-1, T>>((*this->lhs)[n], (*this->rhs)[n]);
      return get_series(addition);
    }

  };

  template<class T>
  class SERIES_ADDITION<1,T> : public SERIES_BINARY_OPERATOR<1,T>
  {
    using SERIES_BINARY_OPERATOR<1,T>::SERIES_BINARY_OPERATOR;
  public:
    std::shared_ptr<T> get_coeff(const unsigned long n) const override
    {
      return std::make_shared<T>(this->lhs->get(n)+this->rhs->get(n));
    }
  };
  

  template <class R, class... Args>
  class ADDITION : public BinaryNode<R, Args...>{
    using BinaryNode<R,Args...>::BinaryNode;
  public:
    R evaluate(const Args&... args) const;
    std::shared_ptr<ANALYTIC<R,Args...>> to_analytic() const;
  };

  template <class R, class... Args>
  class SCALAR_ADDITION : public Node<R, Args...>{
    using node_ptr = std::shared_ptr<Node<R, Args...>>;
  private:
    R scalar;
    node_ptr node;
  public:
    SCALAR_ADDITION(const node_ptr& node, const R& scalar):
      scalar(scalar), node(node) {};
    R evaluate(const Args&... args) const;
    std::shared_ptr<ANALYTIC<R,Args...>> to_analytic() const;
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
    return node->evaluate(args...)+scalar;
  }

  template <class R, class... Args>
  std::shared_ptr<ANALYTIC<R,Args...>> ADDITION<R,Args...>::to_analytic() const
    {
      auto l = this->lhs->to_analytic();
      auto r = this->rhs->to_analytic();
      std::shared_ptr<SERIES_OPERATOR<sizeof...(Args), R>> addition= std::make_shared<SERIES_ADDITION<sizeof...(Args), R>>(l->get_series(), r->get_series());
      auto add_pwr = get_series(addition);
      auto add_M = l->get_M()+r->get_M();
      auto add_r = minimum(l->get_r(), r->get_r());
      return std::make_shared<ANALYTIC<R, Args...>> (add_pwr, add_M, add_r);
    }

  // addition operator
  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> operator+(const std::shared_ptr<Node<R,Args...>>& lhs,const std::shared_ptr<Node<R,Args...>>& rhs)
  {
    return std::make_shared<ADDITION<R, Args...>>(lhs, rhs);
  }


} // namespace iRRAM


#endif
