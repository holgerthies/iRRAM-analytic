/*-------------------------------------------------
* Class for composition
 ------------------------------------------------*/
#ifndef COMPOSITION_H
#define COMPOSITION_H
#include "ANALYTIC.h"
namespace iRRAM
{

  template <size_t d, class R, class... Args>
  struct coefficient_getter_cache
  {
    using node_ptr = std::shared_ptr<Node<R, Args...>>;
    node_ptr rhs;
    std::shared_ptr<coefficient_getter_cache<d-1, R, Args...>> next;
    coefficient_getter_cache(const node_ptr& rhs): rhs(rhs), next(std::make_shared<coefficient_getter_cache<d-1,R,Args...>>(rhs)) {}

    R get(const size_t k, const tutil::n_tuple<d,size_t>& idx, const tutil::n_tuple<sizeof...(Args)-d,size_t>& r_ind, const tutil::n_tuple<sizeof...(Args)-d,size_t>& B_ind)
    {
      R ans=0;
      if(k == 0 && std::get<0>(idx) == 0) return 1;
      if(k == 0) return 0;
      for(int i=0; i<=std::get<0>(idx); i++){
        ans += next->get(k, tutil::tail(idx), std::tuple_cat(r_ind,std::make_tuple(i)), std::tuple_cat(B_ind, std::make_tuple(std::get<0>(idx)-i)));
      }
      return ans;
    }
  };

  template <class R, class... Args>
  struct coefficient_getter_cache<0, R, Args...>
  {
    using node_ptr = std::shared_ptr<Node<R, Args...>>;
    node_ptr rhs;
    coefficient_getter_cache(const node_ptr& rhs): rhs(rhs) {}

    R get(const size_t k, const std::tuple<>& idx, const tutil::n_tuple<sizeof...(Args),size_t>& r_ind, const tutil::n_tuple<sizeof...(Args),size_t>& B_ind)
    {
      auto root = std::make_shared<coefficient_getter_cache<sizeof...(Args), R, Args...>>(rhs);
      return rhs->get_coefficient(r_ind)*root->get(k-1, B_ind, std::tuple<>(), std::tuple<>());
    }
  };

  template <class R, class... Args>
  class COMPOSITION : public Node<R, Args...>{
    using node_ptr = std::shared_ptr<Node<R, Args...>>;
    using node_ptr1d = std::shared_ptr<Node<R,R>>;
  private:
    std::shared_ptr<coefficient_getter_cache<sizeof...(Args), R, Args...>> cache;
    mutable std::vector<std::vector<R>> B; // cache
  public:
    node_ptr1d lhs;
    node_ptr rhs;
    COMPOSITION(const node_ptr1d& lhs, const node_ptr& rhs):
      cache(std::make_shared<coefficient_getter_cache<sizeof...(Args), R, Args...>>(rhs)),lhs(lhs), rhs(rhs) 
    {
    }

    R evaluate(const Args&... args) const override{
      return lhs->evaluate(rhs->evaluate(args...));
    }

    REAL get_r() const override {
      return rhs->get_r();
    }
    
    REAL get_M(const REAL& r) const override {
      return lhs->get_M(r);
    }

    std::string to_string() const override
    {
      return this->lhs->to_string()+"("+this->rhs->to_string()+")";
      
    }

    R get_coefficient(const tutil::n_tuple<sizeof...(Args),size_t>& idx) const override
    {
      R ans=0;
      for(int i=0;i<=std::get<0>(idx); i++)
      {
        ans += lhs->get_coefficient(i)*cache->get(i, idx, std::tuple<>(), std::tuple<>());
      }
      return ans;
    }

    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::COMPOSITION;
    }
  };

  // composition operators
  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> compose(const std::shared_ptr<Node<R,R>>& lhs,const std::shared_ptr<Node<R,Args...>>& rhs)
  {
    
    return std::make_shared<COMPOSITION<R, Args...>>(lhs, rhs);
  }


} // namespace iRRAM


#endif
