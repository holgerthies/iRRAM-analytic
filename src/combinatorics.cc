#include "iRRAM/core.h"
#include "combinatorics.h"
#include <utility>
#include <unordered_map>
#include <map>
using std::vector;
using std::map;
using std::pair;
using std::make_pair;
namespace iRRAM
{
  // get all partitions of size k of the number n
  vector<vector<unsigned long>> partitions(const unsigned long n, const unsigned long k){
    if(k == 1) return vector<vector<unsigned long>>{{n}};
    vector<vector<unsigned long>> ans;
    for(int i=0; i<=n; i++){
      for(vector<unsigned long> p : partitions(n-i, k-1)){
        p.push_back(i);
        ans.push_back(p);
      }
    }
    return ans;
  }

  INTEGER factorial(int n){
    static vector<INTEGER> ans={1};
    for(int j=ans.size(); j<=n; j++){
      ans.push_back(ans.back()*j);
    }
    return ans[n];
  }

  REAL choose(int n, int k){
    static map<pair<int,int>, REAL> mem;
    if(mem.find(make_pair(n,k)) != mem.end())
      return mem[make_pair(n,k)];
    if(k == 0 || n==k) return 1;
    if(k > n) return 0;
    mem[make_pair(n,k)] = choose(n-1,k)*REAL(n)/REAL(n-k); 
    return mem[make_pair(n,k)];
  }
}
