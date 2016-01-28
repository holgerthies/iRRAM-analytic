/*-------------------------------------------------
* Class for derivative
 ------------------------------------------------*/
#ifndef DERIVATIVE_H
#define DERIVATIVE_H
#include "ANALYTIC.h"
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

  REAL get_deriv_M_factor()
  {
    return 1;
  }

  template<class... orders>
  REAL get_deriv_M_factor(int d, orders... rest)
  {
    auto p = power(2, d+1);
    REAL fact=1;
    for(int j=2; j<=d; j++) fact *= j;
    return fact*p*get_deriv_M_factor(rest...);
  }


  template <class R, class... Args>
  class DERIVATIVE : public Node<R, Args...>{
    using node_ptr = std::shared_ptr<Node<R, Args...>>;
  private:
    node_ptr node;
    std::shared_ptr<ANALYTIC<R,Args...>> analytic;
    
  public:
    template<class... orders>
    DERIVATIVE(const node_ptr& node, orders... ods):
      node(node)
    {
      auto f = node->to_analytic();
      auto dpwr= get_derivative_series<sizeof...(Args), R>(f->get_series(),0, ods...);
      auto new_r = f->get_r()/2;
      auto new_M =  f->get_M()/f->get_r()*get_deriv_M_factor(ods...);
      this->analytic = std::make_shared<ANALYTIC<R,Args...>>(dpwr, new_M, new_r);
    }

    R evaluate(const Args&... args) const;
    std::shared_ptr<ANALYTIC<R,Args...>> to_analytic() const;
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


} // namespace iRRAM


#endif
