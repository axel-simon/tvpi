/* Test program for the planar entailment.
 */

#include <iostream>
#include "planar.hh"
#include "polyhedron.hh"

using namespace std;
using namespace Tvpi;

  
int main() {
  Interval<false> xBound(6);
  Interval<false> yBound(3);
  Polyhedron p;

  vector<Inequality*> eqs;
  eqs.push_back(new Inequality(-1,-2,-12));

  p.addInequalitySet<false>(xBound, yBound, eqs);
  //  cout << "x in " << xBound << ", y in " << yBound << endl;

  bool correct= !xBound.isEmpty() && !yBound.isEmpty();

  // Return 1 on error.
  return !correct;
};
