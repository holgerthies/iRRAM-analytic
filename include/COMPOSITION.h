/*-------------------------------------------------
* Class for composition
 ------------------------------------------------*/
#ifndef COMPOSITION_H
#define COMPOSITION_H
#include "ANALYTIC.h"
namespace iRRAM
{

  template <size_t d, class T>
  class SERIES_COMPOSITION : public SERIES_OPERATOR<d,T>
  {
  private:
    std::shared_ptr<POWERSERIES<1,T>> lhs;
    std::shared_ptr<POWERSERIES<d,T>> rhs;
    mutable std::vector<std::vector<std::shared_ptr<POWERSERIES<d-1,T>>>> B; // cache
  public:
    SERIES_COMPOSITION(const std::shared_ptr<POWERSERIES<1,T>>& lhs, const std::shared_ptr<POWERSERIES<d,T>>& rhs):
      lhs(lhs), rhs(rhs)
    {
    }

    std::shared_ptr<POWERSERIES<d-1,T>> get_coeff(const unsigned long n) const override 
    {
      
      if(B.size() == 0){
        auto one = std::make_shared<POWERSERIES<d-1,T>>(1);
        B.push_back(std::vector<std::shared_ptr<POWERSERIES<d-1,T>>>{one});
      }
      auto zero = std::make_shared<POWERSERIES<d-1,T>>(0);
      std::fill_n(std::back_inserter(B), n-B.size()+1, std::vector<std::shared_ptr<POWERSERIES<d-1,T>>>());
      std::fill_n(std::back_inserter(B[0]), (n+1)-B[0].size(), zero);
      for(int k=1; k<=n; k++){
        auto Bksize = B[k].size();
        std::fill_n(std::back_inserter(B[k-1]), n-B[k-1].size()+1, zero);
        for(int j=Bksize; j<=n; j++){
          std::shared_ptr<SERIES_OPERATOR<d-1,T>> bkj = std::make_shared<SERIES_MULTIPLICATION<d-1,T>>(B[k-1][j], (*rhs)[0]);
          for(int t=1; t<=j; t++){
            std::shared_ptr<SERIES_OPERATOR<d-1,T>> multiplication = std::make_shared<SERIES_MULTIPLICATION<d-1,T>>(B[k-1][j-t], (*rhs)[t]);
            bkj = std::make_shared<SERIES_ADDITION<d-1,T>>(get_series(bkj), get_series(multiplication));
          }
          B[k].push_back(get_series(bkj));
        }
      }
      auto ans = std::make_shared<POWERSERIES<d-1,T>>(0);
      for(int k=0; k<=n; k++)
      {
        std::shared_ptr<SERIES_OPERATOR<d-1,T>> multiplication = std::make_shared<SERIES_SCALAR_MULTIPLICATION<d-1,T>>(B[k][n],lhs->get(k) );
        std::shared_ptr<SERIES_OPERATOR<d-1,T>> addition = std::make_shared<SERIES_ADDITION<d-1,T>>(ans, get_series(multiplication));
        ans = get_series(addition);
      }
      return ans;
    }

  };

  template <class T>
    class SERIES_COMPOSITION<1,T> : public SERIES_OPERATOR<1,T>
  {
  private:
    std::shared_ptr<POWERSERIES<1,T>> lhs;
    std::shared_ptr<POWERSERIES<1,T>> rhs;
    mutable std::vector<std::vector<T>> B; // cache
  public:
    SERIES_COMPOSITION(const std::shared_ptr<POWERSERIES<1,T>>& lhs, const std::shared_ptr<POWERSERIES<1,T>>& rhs):
      lhs(lhs), rhs(rhs)
    {
    }

    std::shared_ptr<T> get_coeff(const unsigned long n) const override 
    {
      if(B.size() == 0){
        B.push_back(std::vector<T>{1});
      }
      std::fill_n(std::back_inserter(B), n-B.size()+1, std::vector<T>());
      std::fill_n(std::back_inserter(B[0]), n-B[0].size()+1, 0);
      for(int k=1; k<=n; k++){
        std::fill_n(std::back_inserter(B[k]), n-B[k].size()+1, 0);
        for(int j=0; j<=n; j++){
          for(int t=0; t<=j; t++){
            B[k][j] += rhs->get(t)*B[k-1][j-t];
          }
        }
      }
      T ans = 0;
      for(int k=0; k<=n; k++)
      {
        ans += lhs->get(k)*B[k][n];
      }
      return std::make_shared<T>(ans);
    }

  };


  template <class R, class... Args>
  class COMPOSITION : public Node<R, Args...>{
    using node_ptr = std::shared_ptr<Node<R, Args...>>;
    using node_ptr1d = std::shared_ptr<Node<R,R>>;
  private:
    node_ptr1d lhs;
    node_ptr rhs;
  public:
    COMPOSITION(const node_ptr1d& lhs, node_ptr rhs):
      lhs(lhs), rhs(rhs)
    {
    }

    R evaluate(const Args&... args) const;
    std::shared_ptr<ANALYTIC<R,Args...>> to_analytic() const;
  };

  // member definitions 

  template <class R, class... Args>
  R COMPOSITION<R, Args...>::evaluate(const Args&... args) const
  {
    return lhs->evaluate(rhs->evaluate(args...));
  }


  template <class R, class... Args>
  std::shared_ptr<ANALYTIC<R,Args...>> COMPOSITION<R,Args...>::to_analytic() const
    {
      auto l = lhs->to_analytic();
      auto r = rhs->to_analytic();
      std::shared_ptr<SERIES_OPERATOR<sizeof...(Args), R>> composition= std::make_shared<SERIES_COMPOSITION<sizeof...(Args), R>>(l->get_series(), r->get_series());
      return std::make_shared<ANALYTIC<R, Args...>>(get_series(composition), l->get_M(), r->get_r());
    }

  // composition operators
  template <class R, class... Args, class... orders>
  std::shared_ptr<Node<R, Args...>> compose(const std::shared_ptr<Node<R,R>>& lhs,const std::shared_ptr<Node<R,Args...>>& rhs)
  {
    
    return std::make_shared<COMPOSITION<R, Args...>>(lhs, rhs);
  }


} // namespace iRRAM


#endif
