/*-------------------------------------------------
* Class for analytic continuation
 ------------------------------------------------*/
#ifndef CONTINUATION_H
#define CONTINUATION_H
#include "ANALYTIC.h"
namespace iRRAM
{
  template<size_t d, class T, class... Args>
  class CONTINUATION_SERIES
  {
  private:
    std::shared_ptr<POWERSERIES<d,T>> pwr;
  public:
    template<class... Orders>
    CONTINUATION_SERIES(const std::shared_ptr<Node<T, Args...>>& node, const Args&... args, const Orders... orders)
    {
      pwr = std::make_shared<POWERSERIES<d,T>>([node, orders..., args...] (unsigned long n) {
          auto continuation = std::make_shared<CONTINUATION_SERIES<d-1, T, Args...>>(node, args..., orders..., n);
          return continuation->get_series();
        });
    };

    auto get_series() -> decltype(pwr)
    {
      return pwr;
    }

  };
  
  template<class T, class... Args>
  class CONTINUATION_SERIES<1,T, Args...>
  {
  private:
    std::shared_ptr<POWERSERIES<1,T>> pwr;
  public:
    template<class... Orders>
    CONTINUATION_SERIES(const std::shared_ptr<Node<T, Args...>>& node, const Args&... args, const Orders... orders)
    {
      pwr = std::make_shared<POWERSERIES<1,T>>([node, orders..., args...] (unsigned long n) {
          
          return std::make_shared<T>(inv_factorial(n, orders...)*derive(node,n, orders...)->evaluate(args...));
        });
    };

    auto get_series() -> decltype(pwr)
    {
      return pwr;
    }

  };

  template <class R, class... Args>
  class CONTINUATION : public Node<R, Args...>{
  using node_ptr = std::shared_ptr<Node<R, Args...>>;
  private:
    node_ptr node;
    std::shared_ptr<ANALYTIC<R,Args...>> analytic;
  public:
    CONTINUATION(const node_ptr& node, const REAL& new_M, const REAL& new_r, Args... args):
      node(node)
    {
      auto dpwr= std::make_shared<CONTINUATION_SERIES<sizeof...(args),R,Args...>>(node,args...);
      this->analytic = std::make_shared<ANALYTIC<R,Args...>>(dpwr->get_series(), new_M, new_r);
    }
    R evaluate(const Args&... args) const override{
      return analytic->evaluate(args...);
    };
    std::shared_ptr<ANALYTIC<R,Args...>> to_analytic() const override{
      return analytic;
    };
  };

  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> continue_around(const std::shared_ptr<Node<R,Args...>>& node, const REAL& new_M, const REAL& new_r, Args... args)
  {
    return std::make_shared<CONTINUATION<R, Args...>>(node, new_M, new_r, args...);
  }

} // namespace iRRAM


#endif
