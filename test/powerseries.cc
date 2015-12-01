#include "iRRAM.h"
#include "POWERSERIES.h"
#include <iostream>
using namespace std;
using namespace iRRAM;
void compute(){
  using seq = std::function<POWERSERIES<1,REAL>(const unsigned long)>;
  using seq_ptr = shared_ptr<seq>;
  function<REAL(const unsigned long)> series0 = [] (const unsigned long n) {return REAL(int(n+1));};
  seq series = [] (const unsigned long n) {
    auto subseq = [n] (unsigned long m) {return REAL(int(n+m+1));};
    return POWERSERIES<1,REAL>(subseq);
	  
  };

  POWERSERIES<1, REAL> pwr1(series0);

  function<REAL(unsigned long, unsigned long, unsigned long)> series2 = [] (unsigned long n, unsigned long m, const unsigned long k) {return REAL(int(n+m-k));};

  shared_ptr<POWERSERIES<2,REAL>> test(new POWERSERIES<2,REAL>(seq_ptr(new seq(series))));
  POWERSERIES<3,REAL> test2(series2);
  auto comp = compose(pwr1, *test);
  auto inv1 = inverse(pwr1);
  iRRAM::cout << comp[2][1] << endl;
  iRRAM::cout << inv1[2] << endl;
  auto one=pwr1*inv1;
  iRRAM::cout << one[0] << "," << one[1] << ", " << one[2] << endl;
  auto inv2 = inverse(*test);
  auto one2=(*test)*inv2;
  iRRAM::cout << one2[0][0] << "," << one2[0][1] << ", " << one2[1][0] << endl;
  //iRRAM::cout << evaluate_partial(test, {3,5}, 2, 5) << endl;
  iRRAM::cout << (*test)[300][245] << endl;
  iRRAM::cout << test2[300][245][1000] << endl;
  iRRAM::cout << ((*test)*(*test))[1][1] << endl;
}

