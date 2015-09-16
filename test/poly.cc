#include "POLYNOMIAL.h"
#include <iostream>
using namespace std;
int main(){
  POLYNOMIAL<1,int> test({10,1,2});
  cout << test(3) << endl;
  return 0;
}
