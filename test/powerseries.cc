#include "POWERSERIES.h"
#include <iostream>
using namespace std;
int main(){
  std::function<int(unsigned long)> series = [] (int n) {return n;};
  POWERSERIES<1,int> test(series);
  cout << test.evaluate_partial(2, 5) << endl;
  return 0;
}

