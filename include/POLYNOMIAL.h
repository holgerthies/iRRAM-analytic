/*-------------------------------------
 * Class for multivariate polynomials
 -------------------------------------*/
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#include <vector>
#include <numeric>
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

  template<size_t d>
  struct continuation_helper;

  template<size_t n, class T, typename... Ts>
  T evaluate_polynomial(const POLYNOMIAL<n,T>&, const T&, Ts...);

  template<class T>
  T evaluate_polynomial(const T&);

  template<size_t n, class T, typename... Ts>
  T upper_bound(const POLYNOMIAL<n,T>&, const T&, Ts...);

  template<class T>
  T upper_bound(const T&);


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
      ans += poly_to_string(P.get_coeff(0),degs...,0);
      for(int i=1; i<P.get_degree(); i++){
        ans += " + "+poly_to_string(P.get_coeff(i),degs..., i);
      }
      return ans;
    }

    REAL get_coefficient_tuple(const REAL& p, const std::tuple<>& t)
    {
      return p;
    }
      
    template<class T, class... Ts>
    T get_coefficient_tuple(const POLY<T, Ts...>& p, const tutil::n_tuple<sizeof...(Ts),size_t>& t)
    {
      return get_coefficient_tuple(p.get_coeff(std::get<0>(t)), tutil::tail(t));
    }

    template<class T, class... Ts>
    class POLY : public Node<T, Ts...>
    {
    public:
      using Node<T,Ts...>::evaluate;
      static const int n = sizeof...(Ts);
      using coeff_type = typename COEFF_TYPE<n,T>::type;
      int max_degree;
    private:
      std::vector<coeff_type> coefficients;
      REAL M, r;
    public:
      POLY() = default;
      POLY(const std::vector<coeff_type>& coeffs) : coefficients(coeffs) {
        r=1;
        M=upper_bound(*this, r);
       max_degree = get_degree();
    for(int i=0; i<get_degree(); i++){
      max_degree = std::max(max_degree, get_max_degree(get_coeff(i)));
    }
 
      };
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

      void set_r(const REAL& r)
      {
        this->r = r;
      }

      void set_M(const REAL& M)
      {
        this->M = M;
      }

      void set_rM(const REAL& r)
      {
        this->set_r(r);
        
        this->set_M(upper_bound(*this, this->get_r()));
      }


      REAL get_r() const override
      {
        return r;
      }

      REAL get_M(const REAL& r) const override
      {
        return M;
      }


      ANALYTIC_OPERATION get_type() const override
      {
        return ANALYTIC_OPERATION::POLYNOMIAL;
      }

      std::shared_ptr<POLY> continuation(const std::vector<T>& center)
      {
        std::vector<coeff_type> coeffs(get_degree());
        auto ans = std::make_shared<POLY>(continuation_helper<n>::get(*this,center));
        return ans;
      }
      // add coefficients to the back
      void push(const coeff_type& x){
        coefficients.push_back(x);
      }

      // get kth coefficient
      coeff_type get_coeff(int k) const{
        if(k >= get_degree()) return coeff_type();
        return coefficients[k];
      }

      T get_coefficient(const tutil::n_tuple<sizeof...(Ts),size_t>& idx) const override
      {
        return get_coefficient_tuple(*this, idx);
      }
 


      T operator()(const Ts&... args)
      {
        return evaluate(args...);
      }

    };
  }
  
  int get_max_degree(const REAL& x){
    return 0;
  }

  template<class T, class... Ts>
  int get_max_degree(const poly_impl::POLY<T,Ts...>& P){
    return P.max_degree;
  }

  template<class T, typename... Ts>
  T evaluate_polynomial(const poly_impl::POLY<T,T,Ts...>& P, const T& x, Ts... rest){
    T ans = 0;
    for(int i=P.get_degree()-1; i>=0; i--){
      ans = evaluate_polynomial(P.get_coeff(i), rest...)+ans*x;
    }
    return ans;
  }

  template<class T>
  T evaluate_polynomial(const T& P){
    return P;
  }


  template<class T, typename... Ts>
  T upper_bound(const poly_impl::POLY<T,T,Ts...>& P, const T& r){
    T ans = 0;
    for(int i=P.get_degree()-1; i>=0; i--){
      ans = upper_bound(P.get_coeff(i), r)+abs(ans)*r;
    }
    return ans;
  }

  template<class T>
  T upper_bound(const T& P, const T& r){
    return abs(P);
  }
  
  template<class T,class... Ts>
  poly_impl::POLY<T,Ts...> operator+(const poly_impl::POLY<T,Ts...>& lhs, const poly_impl::POLY<T,Ts...>& rhs )
  {
    if(lhs.get_degree() == 0) return rhs;
    if(rhs.get_degree() == 0) return lhs;
    
    using coeff_type = typename poly_impl::POLY<T,Ts...>::coeff_type;
    std::vector<coeff_type> coeffs(max(lhs.get_degree(),rhs.get_degree()));
    for(int i=0; i<coeffs.size(); i++){
      coeffs[i] = lhs.get_coeff(i) + rhs.get_coeff(i);
    }
    return poly_impl::POLY<T,Ts...>(coeffs);
  }

  template<class T,class... Ts>
  poly_impl::POLY<T,Ts...> operator+(const REAL& lhs, const poly_impl::POLY<T,Ts...>& rhs )
  {
    using coeff_type = typename poly_impl::POLY<T,Ts...>::coeff_type;
    std::vector<coeff_type> coeffs(max(1,rhs.get_degree()));
    coeffs[0] = lhs+rhs.get_coeff(0);
    for(int i=1; i<coeffs.size(); i++){
      coeffs[i] = rhs.get_coeff(i);
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
      coeffs[i] = lhs*rhs.get_coeff(i);
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
    std::vector<coeff_type> coeffs(lhs.get_degree()+rhs.get_degree()-1);
    for(int i=0; i<coeffs.size(); i++){
      for(int j=0; j<=i; j++){
        coeffs[i] = coeffs[i]+lhs.get_coeff(j)*rhs.get_coeff(i-j);
      }
    }
    return poly_impl::POLY<T,Ts...>(coeffs);
  }

  REAL derive(const REAL& P)
  {
    return P;
  }

  template<class T,class... Ts, class... IDX>
  poly_impl::POLY<T,Ts...> derive(const poly_impl::POLY<T,Ts...>& P, const int order, const IDX... orders)
  {
    if(order >= P.get_degree()) return poly_impl::POLY<T,Ts...>();
    using coeff_type = typename poly_impl::POLY<T,Ts...>::coeff_type;
    std::vector<coeff_type> coeffs(P.get_degree()-order);
    T fact=factorial(order);
    for(int n=0; n<coeffs.size(); n++){
      if(n > 0){
        fact *= T(order+n);
        fact /= T(n);
      }
      coeffs[n] = fact*derive(P.get_coeff(n+order), orders...);
    }
    return poly_impl::POLY<T,Ts...>(coeffs);
  }

  REAL continuation_derivative(const REAL& P)
  {
    return P;
  }

  template<class T,class... Ts, class... IDX>
  poly_impl::POLY<T,Ts...> continuation_derivative(const poly_impl::POLY<T,Ts...>& P, const int order, const IDX... orders)
  {
    if(P.get_degree() <= order) return poly_impl::POLY<T,Ts...>();
    using coeff_type = typename poly_impl::POLY<T,Ts...>::coeff_type;
    std::vector<coeff_type> coeffs(P.get_degree()-order);
    T fact=1;
    for(int n=0; n<coeffs.size(); n++){
      if(n > 0){
        fact *= T(order+n);
        fact /= T(n);
      }
      coeffs[n] = fact*continuation_derivative(P.get_coeff(n+order), orders...);
    }
    return poly_impl::POLY<T,Ts...>(coeffs);
  }

  template<size_t d>
  struct continuation_helper
  {

    template<class T,class... Ts, class... IDX>
    static typename COEFF_TYPE<d+1,T>::type get(const typename COEFF_TYPE<d+1,T>::type& P,  const std::vector<T>& center, const IDX... orders)
    {
      using coeff_type = typename COEFF_TYPE<d,T>::type;
      
      std::vector<coeff_type> coeffs(P.get_degree());
      for(int n=0; n<coeffs.size(); n++){
      REAL fact=1;
        REAL x0 = 1;
        for(int k = 0; k < P.get_degree()-n; k++)
        {
          if(k > 0)
            fact *= T(n+k)/T(k);
          coeffs[n] = coeffs[n]+fact*x0*continuation_helper<d-1>::get(P.get_coeff(n+k), center, orders..., n);
          x0 *= center[center.size()-d];
        }
      }
      return POLYNOMIAL<d,T>(coeffs);
    }
  };
  

  template<>
  struct continuation_helper<0>
  {

    template<class T, class... IDX>
    static T get(const REAL& P, const std::vector<T>& center, const IDX... orders)
    {
      return P;
    }
  };

  template<size_t n>
  struct vderive_helper
  {
    template<class T,class... Ts, class... IDX>
    static poly_impl::POLY<T,Ts...> get(const poly_impl::POLY<T,Ts...>& P, const std::vector<int>& orders, const IDX... idx)
    {
      return vderive_helper<n-1>::get(P, orders, orders[n-1], idx...);
    }
  };
  

  template<>
  struct vderive_helper<0>
  {
    template<class T,class... Ts, class... IDX>
    static poly_impl::POLY<T,Ts...> get(const poly_impl::POLY<T,Ts...>& P, const std::vector<int>& orders, const IDX... idx)
    {
      return derive(P, idx...);
    }
  };
  
  template<size_t n, class T, class... Ts>
  struct poly_tderive_helper
  {
    template<class... orders>
    static poly_impl::POLY<T,Ts...> get(const poly_impl::POLY<T,Ts...>& node, const tutil::n_tuple<sizeof...(Ts), int> tuple, const orders... ods)
    {
      return poly_tderive_helper<n-1, T, Ts...>::get(node,tuple, ods...,std::get<sizeof...(Ts)-n>(tuple));
    }
  };

  template<class T, class... Ts>
  struct poly_tderive_helper<0, T, Ts...>
  {
    template<class... orders>
    static poly_impl::POLY<T,Ts...> get(const poly_impl::POLY<T,Ts...>& node, const tutil::n_tuple<sizeof...(Ts), int> tuple, const orders... ods)
    {
      return derive(node, ods...);
    }
  };

  // derivative from std::tuple
  template<class T, class... Ts>
  poly_impl::POLY<T,Ts...> derive(const poly_impl::POLY<T,Ts...>& node, const tutil::n_tuple<sizeof...(Ts),int>& tuple)
  {
    return poly_tderive_helper<sizeof...(Ts), T, Ts...>::get(node, tuple);
    
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
  
  template<class T,class... Ts>
  poly_impl::POLY<T,Ts...> derive(const poly_impl::POLY<T,Ts...>& P, const std::vector<int>& orders)
  {
    return vderive_helper<sizeof...(Ts)>::get(P, orders);
  }
  
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

  template<size_t n, size_t m>
  struct variable_symbol_helper
  {
    static POLYNOMIAL<n, REAL> get()
    {
      using coeff_type = typename COEFF_TYPE<n,REAL>::type;
      std::vector<coeff_type> coeffs(1);
      coeffs[0] = variable_symbol_helper<n-1, m-1>::get();
      return POLYNOMIAL<n,REAL>(coeffs);
    }
  };
  
  template<size_t n>
  struct variable_symbol_helper<n,0>
  {
    static POLYNOMIAL<n,REAL> get()
    {
      using coeff_type = typename COEFF_TYPE<n,REAL>::type;
      std::vector<coeff_type> coeffs(2, coeff_type());
      coeffs[1] = 1+coeffs[1];
      return POLYNOMIAL<n,REAL>(coeffs);
    }
  };

  template<size_t n, size_t m>
  std::shared_ptr<typename node_type<n>::type> variable_symbol()
  {
    static_assert(m < n, "variable symbol exceeds number of variables.");
    return std::make_shared<POLYNOMIAL<n,REAL>>(variable_symbol_helper<n,m>::get());
  }
}

#endif
