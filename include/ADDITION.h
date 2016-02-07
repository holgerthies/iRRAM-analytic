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
  
  template <size_t d, class T>
  class SERIES_SCALAR_ADDITION : public SERIES_SCALAR_OPERATOR<d,T>
  {
    using SERIES_SCALAR_OPERATOR<d,T>::SERIES_SCALAR_OPERATOR;
  public:
    std::shared_ptr<POWERSERIES<d-1,T>> get_coeff(const unsigned long n) const override 
    {
      if(n == 0)
      {
        std::shared_ptr<SERIES_OPERATOR<d-1,T>> addition = std::make_shared<SERIES_SCALAR_ADDITION<d-1, T>>((*this->series)[n], this->scalar);
      return get_series(addition);
      }
      return (*this->series)[n];
    }

  };

  template<class T>
  class SERIES_SCALAR_ADDITION<1,T> : public SERIES_SCALAR_OPERATOR<1,T>
  {
    using SERIES_SCALAR_OPERATOR<1,T>::SERIES_SCALAR_OPERATOR;
  public:
    std::shared_ptr<T> get_coeff(const unsigned long n) const override
    {
      if(n==0) return std::make_shared<T>(this->series->get(0)+this->scalar);
      return (*this->series)[n];
    }
  };
  

  template <class R, class... Args>
  class ADDITION : public BinaryNode<R, Args...>{
    using BinaryNode<R,Args...>::BinaryNode;
  public:
    R evaluate(const Args&... args) const override;
    std::shared_ptr<ANALYTIC<R,Args...>> to_analytic() const override;
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
    std::shared_ptr<ANALYTIC<R,Args...>> to_analytic() const override;
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

  template <class R, class... Args>
  std::shared_ptr<ANALYTIC<R,Args...>> SCALAR_ADDITION<R,Args...>::to_analytic() const
    {
      auto f = this->node->to_analytic();
      std::shared_ptr<SERIES_OPERATOR<sizeof...(Args), R>> addition= std::make_shared<SERIES_SCALAR_ADDITION<sizeof...(Args), R>>(f->get_series(), this->scalar);
      auto add_pwr = get_series(addition);
      auto add_M = f->get_M()+this->scalar;
      auto add_r = f->get_r();
      return std::make_shared<ANALYTIC<R, Args...>> (add_pwr, add_M, add_r);
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
