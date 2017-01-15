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

  template <class... Params>
  struct coefficient_function
  {
    virtual REAL get_coeff(const Params... params) const = 0;
  };

  template<size_t n, class R, class... Args>
  struct coefficient_getter
  {
    template<class... Params>
    static R get(const ANALYTIC<R,Args...>& f, const tutil::n_tuple<n, size_t>& idx, const Params&... params)
    {
      return coefficient_getter<n-1, R, Args...>::get(f, tutil::tail(idx), params..., std::get<0>(idx));
    }
  };
  

  template<class R, class... Args>
  struct coefficient_getter<1, R, Args...>
  {
    template<class... Params>
    static R get(const ANALYTIC<R,Args...>& f, const std::tuple<size_t>& idx, const Params&... params)
    {
      return f.get_coeff(params..., std::get<0>(idx));
    }
  };
  

  template <class R, class... Args>
  class ANALYTIC : public Node<R,Args...>, public tutil::repeat<sizeof...(Args), coefficient_function, size_t> {
    using pwr_series = POWERSERIES<sizeof...(Args),R>;
    using pwr_series_ptr = std::shared_ptr<pwr_series>;
  public:
    virtual ~ANALYTIC() = default;


    virtual R get_coefficient(const tutil::n_tuple<sizeof...(Args), size_t>& t) const override
    {
      return coefficient_getter<sizeof...(Args), R, Args...>::get(*this, t);
    }

    virtual std::string to_string() const override
    {
      auto addr = reinterpret_cast<std::uintptr_t>(this->get_pwr().get()) % 100;
      return "f"+std::to_string(addr);
    };


    virtual ANALYTIC_OPERATION get_type() const override 
    {
      return ANALYTIC_OPERATION::ANALYTIC;
    }
   
    void reset_visited() const override
    {
      this->visited = false;
    }

    int count_nodes() const override
    {
      if(!this->visited){
        this->visited = true;
        return 1;
      }
      return 0;
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
