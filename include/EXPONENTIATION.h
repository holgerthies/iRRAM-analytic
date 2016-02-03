/*-------------------------------------------------
* Class for exponentiation
 ------------------------------------------------*/
#ifndef EXPONENTIATION_H
#define EXPONENTIATION_H
#include "ANALYTIC.h"
namespace iRRAM
{


  template <class R, class... Args>
  class EXPONENTIATION : public Node<R, Args...>{
    using node_ptr = std::shared_ptr<Node<R, Args...>>;
  private:
    node_ptr node;
  public:
    EXPONENTIATION(node_ptr node):
      node(node)
    {
    }

    R evaluate(const Args&... args) const;
    std::shared_ptr<ANALYTIC<R,Args...>> to_analytic() const;
  };

  // member definitions 

  template <class R, class... Args>
  R EXPONENTIATION<R, Args...>::evaluate(const Args&... args) const
  {
    return exp(node->evaluate(args...));
  }


  template <class R, class... Args>
  std::shared_ptr<ANALYTIC<R,Args...>> EXPONENTIATION<R,Args...>::to_analytic() const
    {
      std::function<R(unsigned long)> exp_series = [] (unsigned long n) -> R{
        return inv_factorial(n);
      };
      auto exppwr = std::make_shared<POWERSERIES<1,R>>(exp_series);
      auto f = node->to_analytic();
      std::shared_ptr<SERIES_OPERATOR<sizeof...(Args), R>> exponentiation= std::make_shared<SERIES_COMPOSITION<sizeof...(Args), R>>(exppwr, f->get_series());
      REAL r = f->get_r();
      REAL M = power(euler(), f->get_M());
      return std::make_shared<ANALYTIC<R, Args...>>(get_series(exponentiation), M, r);
    }

  // exponentiation operators
  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> exp(const std::shared_ptr<Node<R,Args...>>& node)
  {
    
    return std::make_shared<EXPONENTIATION<R, Args...>>(node);
  }

} // namespace iRRAM


#endif
