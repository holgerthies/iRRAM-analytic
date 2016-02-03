/*-------------------------------------------------
 * Class for real analytic functions 
 ------------------------------------------------*/
#ifndef ANALYTIC_H
#define ANALYTIC_H
#include "POWERSERIES.h"
#include "Node.h"
#include <vector>
#include <memory>
namespace iRRAM
{

  
  template <class R, class... Args>
  class ANALYTIC : public Node<R,Args...> {
    using pwr_series = POWERSERIES<sizeof...(Args),R>;
    using pwr_series_ptr = std::shared_ptr<pwr_series>;
  private:
    pwr_series_ptr pwr;
    REAL M,r; // maximum of the function and radius of convergence
  public:
    ANALYTIC(const pwr_series_ptr& pwr, const REAL M, const REAL r):
      pwr(pwr), M(M), r(r) 
    {
    };
    
    // template<typename... ARGS>
    // T get_coeff(ARGS... args);

    R evaluate(const Args&... args) const;

    template<class... Indices>
    R get_coeff(const Indices&... ind)
    {
      return pwr->get_coeff(ind...);
    }

    std::shared_ptr<ANALYTIC> to_analytic() const
    {
      return std::make_shared<ANALYTIC>(*this);
    }

    REAL get_r() const{return r;}

    REAL get_M() const{return M;}

    pwr_series_ptr get_series()
    {
      return pwr;
    }

  };

  template <class R, class... Args>
  R ANALYTIC<R, Args...>::evaluate(const Args&... args) const
  {
    return PWRSERIES_IMPL::evaluate(pwr, M, r, args...);
  }

  
  

  // construct from anything that can be used to construct a powerseries
  template <class R, class... Args, class D>
  std::shared_ptr<Node<R,Args...>> make_analytic(const D&& pwr, const REAL& M, const REAL& r)
  {
    auto series = std::make_shared<POWERSERIES<sizeof...(Args),REAL>>(pwr);
     return std::make_shared<ANALYTIC<R,Args...>>(series, M, r);
  }

  
} // namespace iRRAM

#include "ADDITION.h"
#include "SUBTRACTION.h"
#include "MULTIPLICATION.h"
#include "DIVISION.h"
#include "DERIVATIVE.h"
#include "COMPOSITION.h"
#include "IVPSOLVER.h"
#include "CONTINUATION.h"
#include "SINE.h"
#include "EXPONENTIATION.h"


#endif
