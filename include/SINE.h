/*-------------------------------------------------
* Class for sine
 ------------------------------------------------*/
#ifndef SINE_H
#define SINE_H
#include "ANALYTIC.h"
namespace iRRAM
{


  template <class R, class... Args>
  class SINE : public Node<R, Args...>{
    using node_ptr = std::shared_ptr<Node<R, Args...>>;
  private:
    node_ptr node;
  public:
    SINE(node_ptr node):
      node(node)
    {
    }

    R evaluate(const Args&... args) const;
    std::shared_ptr<ANALYTIC<R,Args...>> to_analytic() const;
  };

  template <class R, class... Args>
  class COSINE : public Node<R, Args...>{
    using node_ptr = std::shared_ptr<Node<R, Args...>>;
  private:
    node_ptr node;
  public:
    COSINE(node_ptr node):
      node(node)
    {
    }

    R evaluate(const Args&... args) const;
    std::shared_ptr<ANALYTIC<R,Args...>> to_analytic() const;
  };
  // member definitions 

  template <class R, class... Args>
  R SINE<R, Args...>::evaluate(const Args&... args) const
  {
    return sin(node->evaluate(args...));
  }

  template <class R, class... Args>
  R COSINE<R, Args...>::evaluate(const Args&... args) const
  {
    return cos(node->evaluate(args...));
  }

  template <class R, class... Args>
  std::shared_ptr<ANALYTIC<R,Args...>> SINE<R,Args...>::to_analytic() const
    {
      std::function<R(unsigned long)> sin_series = [] (unsigned long n) -> R{
        if(n % 2 == 0) return 0;
        if((n-1) % 4 == 0) return inv_factorial(n);
        return -inv_factorial(n);
      };
      auto sinpwr = std::make_shared<POWERSERIES<1,R>>(sin_series);
      auto f = node->to_analytic();
      std::shared_ptr<SERIES_OPERATOR<sizeof...(Args), R>> sine= std::make_shared<SERIES_COMPOSITION<sizeof...(Args), R>>(sinpwr, f->get_series());
      REAL r = f->get_r();
      REAL M = power(euler(), f->get_M());
      return std::make_shared<ANALYTIC<R, Args...>>(get_series(sine), M, r);
    }

  template <class R, class... Args>
  std::shared_ptr<ANALYTIC<R,Args...>> COSINE<R,Args...>::to_analytic() const
    {
      std::function<R(unsigned long)> cos_series = [] (unsigned long n) -> R{
        if(n % 2 == 1) return 0;
        if(n % 4 == 0) return inv_factorial(n);
        return -inv_factorial(n);
      };
      auto cospwr = std::make_shared<POWERSERIES<1,R>>(cos_series);
      auto f = node->to_analytic();
      std::shared_ptr<SERIES_OPERATOR<sizeof...(Args), R>> cosine= std::make_shared<SERIES_COMPOSITION<sizeof...(Args), R>>(cospwr, f->get_series());
      REAL r = f->get_r();
      REAL M = power(euler(), f->get_M());
      return std::make_shared<ANALYTIC<R, Args...>>(get_series(cosine), M, r);
    }
  // sine operators
  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> sin(const std::shared_ptr<Node<R,Args...>>& node)
  {
    
    return std::make_shared<SINE<R, Args...>>(node);
  }

  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> cos(const std::shared_ptr<Node<R,Args...>>& node)
  {
    
    return std::make_shared<COSINE<R, Args...>>(node);
  }

} // namespace iRRAM


#endif
