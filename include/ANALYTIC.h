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
    std::function<R(const Args&...)> algorithm;
    bool has_algorithm=false;
  public:
    ANALYTIC(const pwr_series_ptr& pwr, const REAL M, const REAL r):
      pwr(pwr), M(M), r(r) 
    {
    };
    
    // template<typename... ARGS>
    // T get_coeff(ARGS... args);

    R evaluate(const Args&... args) const override;

    template<class... Indices>
    R get_coeff(const Indices&... ind)
    {
      return pwr->get_coeff(ind...);
    }

    std::shared_ptr<ANALYTIC> to_analytic() const override
    {
      return std::make_shared<ANALYTIC>(*this);
    }

    REAL get_r() const override {return r;}

    REAL get_M() const {
      return M;
    };
    REAL get_M(const REAL& r, const Args&... args) const override {
      return M;
    };

    std::shared_ptr<Node<R,Args...>> simplify() const override
    {
      return std::make_shared<ANALYTIC>(*this);
    };

    std::string to_string() const override
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

    void add_algorithm(std::function<R(const Args&...)> f)
    {
      algorithm = f;
      has_algorithm = true;
    }
  };
    

  template <class R, class... Args>
  R ANALYTIC<R, Args...>::evaluate(const Args&... args) const 
  {
    if(has_algorithm)
      return algorithm(args...);
    return PWRSERIES_IMPL::evaluate(pwr, M, r, args...);
  }

  
  // construct from anything that can be used to construct a powerseries
  template <class R, class... Args, class D>
  std::shared_ptr<Node<R,Args...>> make_analytic(const D&& pwr, const REAL& M, const REAL& r)
  {
    auto series = std::make_shared<POWERSERIES<sizeof...(Args),REAL>>(pwr);
     return std::make_shared<ANALYTIC<R,Args...>>(series, M, r);
  }
  
  template<int dim, class... Args>
  struct real_analytic
  {
    using type = typename real_analytic<dim-1, REAL, Args...>::type;
    template<class D>
    static type make(const D&& pwr, const REAL& M, const REAL& r) 
    {
      return real_analytic<dim-1, REAL, Args...>::make(std::move(pwr), M, r);
    }
  };
  

  template<class... Args>
  struct real_analytic<0, Args...>
  {
    using type = std::shared_ptr<Node<REAL, Args...>>;
    
    template<class D>
    static type make(const D&& pwr, const REAL&M, const REAL& r)
    {
      return make_analytic<REAL, Args...>(std::move(pwr), M, r);


    }
  };

  // template<class dim, class D>
  // auto make_real_analytic(const D&& pwr, const REAL& M, const REAL& r) -> decltype(real_analytic<dim>::make(pwr, M, r))
  // {
  //   return real_analytic<dim>::make(pwr, M, r);
    
  // }
  

  
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
