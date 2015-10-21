/*------------------------------------------------------------
 * collection of combinatorical functions
-------------------------------------------------------------*/

#ifndef COMBINATORICS_h
#define COMBINATORICS_h
#include <vector>
namespace iRRAM{
  std::vector<std::vector<int>> partitions(const int n, const int k);
  INTEGER factorial(int n);
  REAL choose(int n, int k);
}
#endif
