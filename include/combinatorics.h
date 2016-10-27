/*------------------------------------------------------------
 * collection of combinatorical functions
-------------------------------------------------------------*/

#ifndef COMBINATORICS_h
#define COMBINATORICS_h
#include <vector>
namespace iRRAM{
  std::vector<std::vector<unsigned long>> partitions(const unsigned long n, const unsigned long k);
  std::vector<std::vector<unsigned long>> bounded_count(const std::vector<unsigned long>& bound);
  INTEGER factorial(int n);
  INTEGER choose(int n, int k);
  REAL inv_factorial(const int n);
  REAL inv_factorial();
  template<class... Args>
  REAL inv_factorial(const int n, Args... rest)
  {
    return inv_factorial(n)*inv_factorial(rest...);
  }
  // get all partitions of size k of the number n
  std::vector<std::vector<unsigned long>> partitions(const unsigned long n, const unsigned long k){
    if(k == 1) return std::vector<std::vector<unsigned long>>{{n}};
    std::vector<std::vector<unsigned long>> ans;
    for(int i=0; i<=n; i++){
      for(std::vector<unsigned long> p : partitions(n-i, k-1)){
        p.push_back(i);
        ans.push_back(p);
      }
    }
    return ans;
  }

  INTEGER factorial(const unsigned long n){
    static std::vector<INTEGER> ans={1};
    for(int j=ans.size(); j<=n; j++){
      ans.push_back(ans.back()*j);
    }
    return ans[n];
  }

  INTEGER choose(int n, int k){
    static std::vector<std::vector<INTEGER>> mem={{1}};
    if(k > n) return 0;
    if(mem.size() > n && mem[n].size() > k)
      return mem[n][k];
    if(mem.size() <= n){
      // guarantee all needed coefficients known
      choose(n-1, k);
      mem.resize(n+1);
      mem[n] = std::vector<INTEGER>{1};
    }
    if(mem[n].size() <= k){
      choose(n, k-1);
      mem[n].resize(k+1);
    }
    if(k == 0 || n==k) mem[n][k] = 1;
    else mem[n][k] = choose(n-1,k-1)+choose(n-1,k);
    return mem[n][k];
  }

  std::vector<std::vector<unsigned long>> bounded_count(const std::vector<unsigned long>& bound, const int size){
    if(size == 0) return std::vector<std::vector<unsigned long>>{std::vector<unsigned long>()};
    std::vector<std::vector<unsigned long>> ans; 
    auto rest=bounded_count(bound, size-1);
    for(auto& v : rest){
      v.resize(size);
      for(int i=0; i<=bound[size-1]; i++){
          v[size-1] = i;
          ans.push_back(v);
      }
    }
    return ans;
  }
  std::vector<std::vector<unsigned long>> bounded_count(const std::vector<unsigned long>& bound){
    return bounded_count(bound, bound.size());
  } 
  REAL inv_factorial()
  {
    return 1;
    
  }
REAL inv_factorial(const int n){
  using std::log;
  if ((n!=0)&&(n*log(n)-n > 2*-ACTUAL_STACK.actual_prec)){
    REAL return_value(0);
    sizetype error;
    sizetype_set(error,1,ACTUAL_STACK.actual_prec);
    return_value.seterror(error);
    return return_value;
  }
  if (n==0)
    return REAL(1);
  REAL inv_fact=inv_factorial(n-1)/REAL(n);
  return inv_fact;
}
}
#endif
