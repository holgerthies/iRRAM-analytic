/*-------------------------------------------------
* Class for multiplication
 ------------------------------------------------*/
#ifndef MULTIPLICATION_H
#define MULTIPLICATION_H
#include "ANALYTIC.h"
namespace iRRAM
{
  template <size_t d, class T>
  class SERIES_MULTIPLICATION : public SERIES_BINARY_OPERATOR<d,T>
  {
    using SERIES_BINARY_OPERATOR<d,T>::SERIES_BINARY_OPERATOR;
  public:
    std::shared_ptr<POWERSERIES<d-1,T>> get_coeff(const unsigned long n) const override 
    {
      auto c = std::make_shared<POWERSERIES<d-1,T>>(T());
      for(int l=0; l<=n;l++){
        std::shared_ptr<SERIES_OPERATOR<d-1,T>> multiplication = std::make_shared<SERIES_MULTIPLICATION<d-1, T>>((*this->lhs)[l], (*this->rhs)[n-l]);
        std::shared_ptr<SERIES_OPERATOR<d-1,T>> add = std::make_shared<SERIES_ADDITION<d-1,T>>(c, get_series(multiplication));
        c = get_series(add);
      }
      return c;
    }

  };

  template<class T>
  class SERIES_MULTIPLICATION<1,T> : public SERIES_BINARY_OPERATOR<1,T>
  {
    using SERIES_BINARY_OPERATOR<1,T>::SERIES_BINARY_OPERATOR;
  public:
    std::shared_ptr<T> get_coeff(const unsigned long n) const override
    {
      auto c = T();
      for(int l=0; l<=n;l++){
        auto coeff = this->lhs->get(l)*this->rhs->get(n-l);
        c = c+coeff;
      }
      return std::make_shared<T>(c);
    }
  };
  

  template <size_t d, class T>
  class SERIES_SCALAR_MULTIPLICATION : public SERIES_SCALAR_OPERATOR<d,T>
  {
    using SERIES_SCALAR_OPERATOR<d,T>::SERIES_SCALAR_OPERATOR;
  public:
    std::shared_ptr<POWERSERIES<d-1,T>> get_coeff(const unsigned long n) const override 
    {
      std::shared_ptr<SERIES_OPERATOR<d-1,T>> multiplication = std::make_shared<SERIES_SCALAR_MULTIPLICATION<d-1, T>>((*this->series)[n], this->scalar);
      return get_series(multiplication);
    }

  };

  template<class T>
  class SERIES_SCALAR_MULTIPLICATION<1,T> : public SERIES_SCALAR_OPERATOR<1,T>
  {
    using SERIES_SCALAR_OPERATOR<1,T>::SERIES_SCALAR_OPERATOR;
  public:
    std::shared_ptr<T> get_coeff(const unsigned long n) const override
    {
      return std::make_shared<T>(this->series->get(n)*this->scalar);
    }
  };

  template <class R, class... Args>
  class MULTIPLICATION : public BinaryNode<R, Args...>{
    using BinaryNode<R,Args...>::BinaryNode;
  public:
    R evaluate(const Args&... args) const override;
    std::shared_ptr<ANALYTIC<R,Args...>> to_analytic() const override;

    REAL get_r() const override {
      return minimum(this->lhs->get_r(), this->rhs->get_r());
    };
    REAL get_M(const REAL& r, const Args&... args) const override {
      return this->lhs->get_M(r, args...)*this->rhs->get_M(r, args...);
    };

    std::shared_ptr<Node<R,Args...>> simplify() const override
    {
      
      return std::make_shared<MULTIPLICATION>(this->lhs->simplify(), this->rhs->simplify());
    };

    virtual std::string to_string() const override
    {
      return "x("+this->lhs->to_string()+","+this->rhs->to_string()+")";
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
    R evaluate(const Args&... args) const override; 
    std::shared_ptr<ANALYTIC<R,Args...>> to_analytic() const override;
    
    virtual std::string to_string() const override
    {
      return "x("+this->node->to_string()+","+std::to_string(this->scalar.as_double())+")";
      
    }
    REAL get_r() const override {
      return this->node->get_r();
    };
    REAL get_M(const REAL& r, const Args&... args) const override {
      return this->node->get_M(r, args...)*this->scalar;
    };

    std::shared_ptr<Node<R,Args...>> simplify() const override
    {
      
      switch(this->node->get_type()){
      case ANALYTIC_OPERATION::SCALAR_MULTIPLICATION:
        {

          auto child = std::dynamic_pointer_cast<SCALAR_MULTIPLICATION>(this->node);
          auto new_node = std::make_shared<SCALAR_MULTIPLICATION>(child->node->simplify(), this->scalar*child->scalar);
         return new_node->simplify();
        }
      default:
        return std::make_shared<SCALAR_MULTIPLICATION>(this->node->simplify(), this->scalar);
      }
      
    }
    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::SCALAR_MULTIPLICATION;
    }
  };

  // member definitions 
  template <class R, class... Args>
  R MULTIPLICATION<R, Args...>::evaluate(const Args&... args) const
  {
    return this->lhs->evaluate(args...)*this->rhs->evaluate(args...);
  }

  template <class R, class... Args>
  R SCALAR_MULTIPLICATION<R, Args...>::evaluate(const Args&... args) const
  {
    return this->node->evaluate(args...)*this->scalar;
  }

  template <class R, class... Args>
  std::shared_ptr<ANALYTIC<R,Args...>> MULTIPLICATION<R,Args...>::to_analytic() const
    {
      auto l = this->lhs->to_analytic();
      auto r = this->rhs->to_analytic();
      std::shared_ptr<SERIES_OPERATOR<sizeof...(Args), R>> multiplication= std::make_shared<SERIES_MULTIPLICATION<sizeof...(Args), R>>(l->get_series(), r->get_series());
      auto new_pwr = get_series(multiplication);
      auto new_M = l->get_M()*r->get_M();
      auto new_r = minimum(l->get_r(), r->get_r());
      return std::make_shared<ANALYTIC<R, Args...>> (new_pwr, new_M, new_r);
    }

  template <class R, class... Args>
  std::shared_ptr<ANALYTIC<R,Args...>> SCALAR_MULTIPLICATION<R,Args...>::to_analytic() const
    {
      auto f = this->node->to_analytic();
      std::shared_ptr<SERIES_OPERATOR<sizeof...(Args), R>> multiplication= std::make_shared<SERIES_SCALAR_MULTIPLICATION<sizeof...(Args), R>>(f->get_series(), this->scalar);
      auto new_pwr = get_series(multiplication);
      auto new_M = f->get_M()*this->scalar;
      auto new_r = f->get_r();
      return std::make_shared<ANALYTIC<R, Args...>> (new_pwr, new_M, new_r);
    }

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
