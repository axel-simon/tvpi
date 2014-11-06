/* Test program for the planar entailment.
 */

#include <iostream>
#include "planar.hh"

using namespace std;
using namespace Tvpi;

  
int main() {
  // These are inequalities form tests/testPolyhedronHull1
  Inequality inp1 = Inequality(1, 3, 6);
  Inequality inp2 = Inequality(-1, 1, 10);
  Inequality resM = Inequality(-4, 11, 67);
  Inequality resN = Inequality(-4, 11, 68);
  Inequality resP = Inequality(-4, 11, 69);
  
  bool correct= 
    !resM.includes(inp1,inp2) &&
    resN.includes(inp1,inp2) &&
    resP.includes(inp1,inp2);

  // Return 1 on error.
  return !correct;
};
