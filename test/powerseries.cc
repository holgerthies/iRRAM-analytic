#include "iRRAM.h"
#include "POWERSERIES.h"
#include <iostream>
using namespace std;
using namespace iRRAM;
void compute(){
  using seq = std::function<REAL(const vector<unsigned long>&)>;
  using seq_ptr = shared_ptr<seq>;
  seq series = [] (const vector<unsigned long>& v) {int n=v[0]; int m=v[1]; return REAL(n+m);};
  shared_ptr<POWERSERIES<2,REAL>> test(new POWERSERIES<2,REAL>(seq_ptr(new seq(series))));
  iRRAM::cout << evaluate_partial(test, {3,5}, 2, 5) << endl;
  iRRAM::cout << (*test)[3][1] << endl;
  iRRAM::cout << ((*test)*(*test))[2][0] << endl;
}

