/*-------------------------------------------------
 * Class for real analytic functions 
 ------------------------------------------------*/
#ifndef ANALYTIC_H
#define ANALYTIC_H
#include "POWERSERIES.h"
#include "Node.h"
#include "tutil.h"
#include <vector>
#include <memory>
namespace iRRAM
{

  template <class R, class... Args>
  class ANALYTIC : public Node<R,Args...> {
    using pwr_series = POWERSERIES<sizeof...(Args),R>;
    using pwr_series_ptr = std::shared_ptr<pwr_series>;
  public:
    virtual ~ANALYTIC() = default;

    // get cached coefficient
    template<class... Indices>
    R get_coeff(const Indices&... ind)
    {
      return pwr->get_coeff(ind...);
    }
    
    virtual std::shared_ptr<Node<R,Args...>> simplify() const override
    {
    };

    virtual std::string to_string() const override
    {
      auto addr = reinterpret_cast<std::uintptr_t>(pwr.get()) % 100;
      return "f"+std::to_string(addr);
    };

    ANALYTIC_OPERATION get_type() const override 
    {
      return ANALYTIC_OPERATION::ANALYTIC;
    }

    pwr_series_ptr get_series()
    {
      return pwr;
    }

  };
  template<size_t n>
  using REAL_ANALYTIC = tutil::repeat<n+1, ANALYTIC, REAL>;

  template<size_t n, class T>
  std::shared_ptr<tutil::repeat<n+1, Node, REAL>> make_analytic()
  {
    return std::make_shared<T>();
  }

} // namespace iRRAM

#endif
