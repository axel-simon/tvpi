/* Test program for the planar entailment.
 */

#include <iostream>
#include "planar.hh"

using namespace std;
using namespace Tvpi;

  
int main() {
  Inequality inpPos = Inequality(1, 3, 6);
  Inequality inpNeg = Inequality(-1, -3, -6);
  Inequality inpLessPos = Inequality(1, 3, 4);
  Inequality test1 = Inequality(1, 3, 8);
  
  bool correct=
    test1.includes(inpPos,inpNeg) &&
    test1.includes(inpPos, inpLessPos);

  // Return 1 on error.
  return !correct;
};
