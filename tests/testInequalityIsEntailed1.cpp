/* Test program for the planar entailment.
 */

#include <iostream>
#include "planar.hh"

using namespace std;
using namespace Tvpi;

  
int main() {
  Inequality i = Inequality(-1, 1, 0);
  Inequality j = Inequality(-1, 1, 1);
  //  cerr << j << " includes " << i << ": " << j.includes(i) << endl;
  bool correct= j.includes(i) && !i.includes(j);

  // Return 1 on error.
  return !correct;
};
