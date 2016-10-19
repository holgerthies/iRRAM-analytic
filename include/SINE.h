/*-------------------------------------------------
* Class for sine
 ------------------------------------------------*/
#ifndef SINE_H
#define SINE_H
#include "ANALYTIC.h"
#include "COMPOSITION.h"
namespace iRRAM
{


  struct SINE : REAL_ANALYTIC<1>{
  public:
    REAL evaluate(const REAL& x) const override
    {
      return sin(x);
    }

    REAL get_r() const override
    {
      return 100;
    }

    REAL get_M(const REAL& r) const override
    {
      return exp(r);
    }

    std::string to_string() const override
    {
      return "sin";
    }

    
    REAL get_coeff(const size_t n) const override
    {
      if(n % 2 == 0) return 0;
      if((n-1) % 4 == 0) return inv_factorial(n);
      return -inv_factorial(n);
    }

    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::SINE;
    }

  };

  struct COSINE : REAL_ANALYTIC<1>{
  public:
    REAL evaluate(const REAL& x) const override
    {
      return cos(x);
    }

    REAL get_r() const override
    {
      return 100;
    }

    REAL get_M(const REAL& r) const override
    {
      return exp(r);
    }

    std::string to_string() const override
    {
      return "cos";
    }

    
    REAL get_coeff(const size_t n) const override
    {
      if(n % 2 == 1) return 0;
      if(n % 4 == 0) return inv_factorial(n);
      return -inv_factorial(n);
    }

    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::COSINE;
    }

  };

  auto sine_function =make_analytic<1,SINE>(); 
  auto cosine_function =make_analytic<1,COSINE>(); 
  // sine operators
  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> sin(const std::shared_ptr<Node<R,Args...>>& node)
  {
    return compose(sine_function, node);
  }

  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> cos(const std::shared_ptr<Node<R,Args...>>& node)
  {
    auto cosine = make_analytic<1,COSINE>();
    return compose(cosine_function, node);
  }

} // namespace iRRAM


#endif
