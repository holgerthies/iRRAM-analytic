/*-------------------------------------------------
* Affine Arithmetic Implementation
 ------------------------------------------------*/
#ifndef AAREAL_H
#define AAREAL_H
#include <unordered_map>
#include <map>
#include <iRRAM/core.h>
#include <chrono>
using namespace std::chrono;

namespace iRRAM{
  struct error_map
  {
    std::vector<std::pair<int, REAL>> rep;
    auto begin() -> decltype(rep.begin())
    {
      return rep.begin();
    }
    auto end() -> decltype(rep.end())
    {
      return rep.end();
    }
    auto begin() const -> decltype(rep.begin())
    {
      return rep.begin();
    }
    auto end() const -> decltype(rep.end())
    {
      return rep.end();
    }

    auto find(int id) -> decltype(rep.begin())
    {
      for(auto it = rep.begin(); it != rep.end(); it++){
        if(it->first == id)
          return it;
      }
      return rep.end();
    }

    void clear()
    {
      rep.clear();
    }

    void add(const int id, const REAL& x)
    {
      rep.push_back({id, x});
    }

    int size() const
    {
      return rep.size();
    }

    REAL& operator[](const int& id) 
    {
      auto it = find(id);
      if(it != rep.end())
        return it->second;
      else{
        rep.push_back({id, 0});
        return rep[size()-1].second;
      }
        
      
    }
  };
  
  REAL error_to_real(const REAL& x)
  {
    sizetype r;
    x.geterror(r);
    REAL err = scale(REAL(signed(r.mantissa)),r.exponent);
    sizetype_exact(err.error);
    return err;
  }
  sizetype real_to_error(const REAL& error)
  {
    sizetype ans;
    REAL zero_one=0;
    sizetype l;
    sizetype_set(l,1,0);
    zero_one.seterror(l);
    REAL e = zero_one*error;
    e.geterror(ans);
    return ans;
    
  }

  class AAREAL
  {
  private:
    REAL center;
    error_map errors;
    static unsigned int errorsymbol;
    static unsigned int next_error_symbol(){ return errorsymbol++; }
    static double duration;
    
    
  public:
    AAREAL();
    AAREAL(const REAL&);
    AAREAL& operator+=(const AAREAL&);
    AAREAL& operator-=(const AAREAL&);
    AAREAL& operator-=(const REAL&);
    AAREAL& operator*=(const REAL&);
    void clean();
    void add_error(const REAL& e);
    REAL error_radius() const;
    friend AAREAL operator+(const AAREAL&, const AAREAL&);
    friend AAREAL operator-(const AAREAL&, const AAREAL&);
    friend AAREAL operator*(const AAREAL&, const AAREAL&);
    REAL to_real() const;
    void print_error();
  };
  unsigned int AAREAL::errorsymbol = 0;
  double AAREAL::duration = 0;
// move to cc file later
  AAREAL::AAREAL()
  {
    center=0;
  }

  AAREAL::AAREAL(const REAL& x)
  {
    center = x;
    sizetype center_error;
    sizetype_exact(center_error);
    center.seterror(center_error);
    REAL error = error_to_real(x);
    auto s = AAREAL::next_error_symbol();
    errors.add(s,error);
  }
  REAL AAREAL::error_radius() const
  {
    REAL ans=0;
    for (const auto &i : errors) {
      ans += i.second;
    }
    return ans;
  }

  void AAREAL::add_error(const REAL& e)
  {
    auto s = AAREAL::next_error_symbol();
    errors.add(s,e);
  }

  REAL AAREAL::to_real() const
  {
    
    REAL total_error=error_radius();
    //std::cout << "error symbols: " << errors.size() << std::endl;
    
    REAL result=center;
    REAL zero_one=0;
    sizetype l;
    sizetype_set(l,1,0);
    zero_one.seterror(l);
    result += zero_one*total_error;
    return result;
  }

  AAREAL& AAREAL::operator+=(const AAREAL& r)
  {
    center += r.center;
    for(const auto& e : r.errors){
      errors[e.first] += e.second;
    }
    return *this;
  }

  AAREAL& AAREAL::operator-=(const REAL& r)
  {
    center -= r;
    return *this;
  }

  AAREAL& AAREAL::operator-=(const AAREAL& r)
  {
    center -= r.center;
    for(const auto& e : r.errors){
      errors[e.first] -= e.second;
    }
    return *this;
  }

  AAREAL operator+(const AAREAL& lhs, const AAREAL& rhs)
  {
    AAREAL ans=lhs;
    ans += rhs;
    return ans;
  }

  AAREAL operator-(const AAREAL& lhs, const AAREAL& rhs)
  {
    AAREAL ans=lhs;
    ans -= rhs;
    return ans;
  }


  AAREAL& AAREAL::operator*=(const REAL& x)
  {
    center *= x;
    for(auto& e : errors){
      e.second *= x;
    }
    return *this;
  }

  AAREAL operator*(const REAL& x, const AAREAL& y)
  {
    AAREAL ans=y;
    ans *= x;
    return ans;
  }

  AAREAL operator*(const AAREAL& lhs, const AAREAL& rhs)
  {
    REAL m = lhs.center*rhs.center;
    AAREAL ans(m);
    for(const auto& e : lhs.errors){
      ans.errors.add(e.first,rhs.center*e.second);
      
    }
    for(const auto& e : rhs.errors){
      ans.errors[e.first] += lhs.center*e.second;
    }
    REAL na_error=lhs.error_radius()*rhs.error_radius();
    auto index = AAREAL::next_error_symbol();
    ans.errors.add(index,na_error);
    return ans;
  }
  // collect small errors in new error
  void AAREAL::clean()
  {
    if(errors.size() >= 2){
      REAL zero_one=0;
      sizetype l;
      
      sizetype_set(l,1,0);
      zero_one.seterror(l);
    
      REAL s=center,e;
      error_map next_errors;
      REAL th = 1;//power(2, ACTUAL_STACK.actual_prec/2);
      for (const auto &i : errors) {
        if(i.second > th)
          next_errors.add(i.first, i.second);
        else
          s += i.second * zero_one;
      }
      e = error_to_real(s);
      sizetype_exact(s.error);
      center = s;
      errors = next_errors;
      errors.add(next_error_symbol(),e);
    }
  
  }
  void AAREAL::print_error()
  {
  //   sizetype total_error;
  //   sizetype_exact(total_error);
  //   for(auto& e : errors){
  //     std::cout << e.first << ": " << e.second.mantissa << ", " << e.second.exponent << std::endl;
  //     sizetype es;
  //     sizetype_add(es, total_error, e.second);
  //     total_error = es;
  //   }
  //     std::cout << "SUM: " << total_error.mantissa << ", " << total_error.exponent << std::endl;
   }
}
#endif /* AAREAL_H */
