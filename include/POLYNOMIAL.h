/*-------------------------------------
 * Class for multivariate polynomials
 -------------------------------------*/
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#include <vector>
#include <memory>
#include <string>
#include "node.h"
#include "tutil.h"
namespace iRRAM
{
  // forward declarations
  namespace poly_impl{
    template<class T, class... Ts>
    class POLY;
  }

  template<size_t n, class T>
  using POLYNOMIAL = tutil::repeat<n+1,poly_impl::POLY, T>;


  template<size_t n, class T, typename... Ts>
  T evaluate_polynomial(const POLYNOMIAL<n,T>&, const T&, Ts...);

  template<class T>
  T evaluate_polynomial(const T&);

  template<size_t n, class T>
  int get_total_degree(const POLYNOMIAL<n,T>& P);

  // define coeff type for multivariate polynomial
  template <size_t n, class T>
  struct COEFF_TYPE{
    typedef POLYNOMIAL<n-1,T> type;   
  };

  template <class T>
  struct COEFF_TYPE<1,T>{
    typedef T type;   
  };

  namespace poly_impl{
    std::string append_xs(int v0)
    {
      return "";
    }

    template<typename... Ds>
    std::string append_xs(int v0, int d, Ds... rest)
    {
      if(d == 0) return append_xs(v0+1, rest...);
      std::string ans = "x"+std::to_string(v0);
      if(d > 1) ans += "^"+std::to_string(d);
      return ans+append_xs(v0+1, rest...);
    }

    template<typename... Ds>
    std::string poly_to_string(const REAL& P, Ds... degs) 
    {
      auto ans = std::to_string(P.as_double());
      ans += append_xs(0,degs...);
      return ans;
    }

    template<class T, class... Ts, typename... Ds>
    std::string poly_to_string(const POLY<T,Ts...>& P, Ds... degs)
    {
      std::string ans;
      ans += poly_to_string(P.get_coefficient(0),degs...,0);
      for(int i=1; i<P.get_degree(); i++){
        ans += " + "+poly_to_string(P.get_coefficient(i),degs..., i);
      }
      return ans;
    }

    template<class T, class... Ts>
    class POLY : public Node<T, Ts...>
    {
    public:
      static const int n = sizeof...(Ts);
      using coeff_type = typename COEFF_TYPE<n,T>::type;
    private:
      std::vector<coeff_type> coefficients;
      int total_degree;
    public:
      POLY() {};
      POLY(std::vector<coeff_type> coeffs) : coefficients(coeffs) {
        total_degree = 0;
        for(int i=0; i<get_degree(); i++){
          total_degree = max(i+get_total_degree(get_coefficient(i)), total_degree);
        }
      };
      friend int get_total_degree<n,T>(const POLY<T,Ts...>&);
      // real degree might be smaller if leading coefficients are 0
      int get_degree() const{
        return coefficients.size();
      }

      T evaluate(const Ts&... args) const override
      {
        return evaluate_polynomial(*this, args...);
      }
      
      std::string to_string() const override
      {
        return poly_to_string(*this);
      }
      REAL get_r() const override
      {
        return 100;
      }

      REAL get_M(const REAL& r, const Ts&... args) const override
      {
        return 1000;
      }

      ANALYTIC_OPERATION get_type() const override
      {
        return ANALYTIC_OPERATION::POLYNOMIAL;
      }


      std::shared_ptr<Node<T, Ts...>> simplify() const override
      {
        return std::make_shared<POLY>(*this);
      }

      std::shared_ptr<ANALYTIC<T, Ts...>> to_analytic() const override
      {
        return std::make_shared<ANALYTIC<T,Ts...>>(std::make_shared<POWERSERIES<sizeof...(Ts),T>>(), 100, 1000);
      }

      // add coefficients to the back
      void push(const coeff_type& x){
        coefficients.push_back(x);
        total_degree = max(total_degree, get_degree()-1+get_total_degree(x));
      }

      // get kth coefficient
      coeff_type get_coefficient(int k) const{
        if(k >= get_degree()) return coeff_type();
        return coefficients[k];
      }

      T operator()(const Ts&... args)
      {
        return evaluate(args...);
      }

    };
  }
  
  template<class T>
  int get_total_degree(const T& x){
    return 0;
  }
  template<size_t n, class T>
  int get_total_degree(const POLYNOMIAL<n,T>& P){
    return P.total_degree;
  }
  // evaluate only the polynomial given by considering the coefficients
  // between start and end
  template<class T, typename... Ts>
  T evaluate_polynomial(const poly_impl::POLY<T,T,Ts...>& P, const T& x, Ts... rest){
    T ans = 0;
    for(int i=P.get_degree()-1; i>=0; i--){
      ans = evaluate_polynomial(P.get_coefficient(i), rest...)+ans*x;
    }
    return ans;
  }

  template<class T>
  T evaluate_polynomial(const T& P){
    return P;
  }

  template<class T,class... Ts>
  poly_impl::POLY<T,Ts...> operator+(const poly_impl::POLY<T,Ts...>& lhs, const poly_impl::POLY<T,Ts...>& rhs )
  {
    if(lhs.get_degree() == 0) return rhs;
    if(rhs.get_degree() == 0) return lhs;
    
    using coeff_type = typename poly_impl::POLY<T,Ts...>::coeff_type;
    std::vector<coeff_type> coeffs(max(lhs.get_degree(),rhs.get_degree()));
    for(int i=0; i<coeffs.size(); i++){
      coeffs[i] = lhs.get_coefficient(i) + rhs.get_coefficient(i);
    }
    return poly_impl::POLY<T,Ts...>(coeffs);
  }

  template<class T,class... Ts>
  poly_impl::POLY<T,Ts...> operator*(const REAL& lhs, const poly_impl::POLY<T,Ts...>& rhs )
  {
    if(rhs.get_degree() == 0) return poly_impl::POLY<T,Ts...>();
    using coeff_type = typename poly_impl::POLY<T,Ts...>::coeff_type;
    std::vector<coeff_type> coeffs(rhs.get_degree());
    for(int i=0; i<coeffs.size(); i++){
      coeffs[i] = lhs*rhs.get_coefficient(i);
    }
    return poly_impl::POLY<T,Ts...>(coeffs);
  }

  template<size_t n, class T>
  POLYNOMIAL<n,T> operator*(const POLYNOMIAL<n,T>& lhs, const REAL& rhs )
  {
    return rhs*lhs;
  }

  template<class T,class... Ts, class... Ds>
  poly_impl::POLY<T,Ts...> operator*(const poly_impl::POLY<T,Ds...>& lhs, const poly_impl::POLY<T,Ts...>& rhs )
  {
    if(rhs.get_degree() == 0 || lhs.get_degree() == 0) return poly_impl::POLY<T,Ts...>();
    using coeff_type = typename poly_impl::POLY<T,Ts...>::coeff_type;
    std::vector<coeff_type> coeffs(lhs.get_degree()+rhs.get_degree());
    for(int i=0; i<lhs.get_degree(); i++){
      auto P = lhs.get_coefficient(i)*rhs;
      if(i+P.get_degree() > coeffs.size())
        coeffs.resize(i+P.get_degree());
      for(int j=0; j<P.get_degree(); j++){
        coeffs[i+j] = coeffs[i+j]+P.get_coefficient(j);
      }
    }
    return poly_impl::POLY<T,Ts...>(coeffs);
  }

  template <size_t n>
  struct node_type 
  {
    using type = tutil::repeat<n+1,Node, REAL>;
    
  };
  
  template<>
  struct node_type<0>
  {
    typedef REAL type;
  };
  
  template<size_t n>
  std::shared_ptr<typename node_type<n>::type> make_polynomial(const tutil::concat<n, std::initializer_list, REAL>& coeff_list )
  {
    using coeff_type = typename COEFF_TYPE<n,REAL>::type;
    std::vector<coeff_type> coeffs(coeff_list.size());
    int i=0;
    for(auto& L : coeff_list){
      coeffs[i] = *std::dynamic_pointer_cast<coeff_type>(make_polynomial<n-1>(L));
      ++i;
    }
    return std::make_shared<POLYNOMIAL<n,REAL>>(coeffs);
  }

  template<>
  std::shared_ptr<REAL> make_polynomial<0>(const REAL& coeff)
  {
    return std::make_shared<REAL>(coeff);
  }
}

#endif
