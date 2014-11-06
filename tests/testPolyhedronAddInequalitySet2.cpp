/* Test program to check redundancy removal.
 */

#include <iostream>
#include "interval.hh"
#include "planar.hh"
#include "polyhedron.hh"

using namespace std;
using namespace Tvpi;

  
int main() {
  Polyhedron<true> p;

  vector<Inequality*> eqs;
  eqs.push_back(new Inequality(Point(12,2), Point(6,8)));
  eqs.push_back(new Inequality(-1, -8, -56));

  Interval<true> xBound;
  Interval<true> yBound;
  xBound.updateUpper(mpq_class(10));
  p.addInequalitySet(xBound, yBound, eqs);
  //cout << p << "x: " << xBound << " y: " << yBound << endl;

  // Check that the two inequalities have changed the upper bound on x
  // to 8 by testing that 9 is not in xBound anymore.
  bool correct=!xBound.includes(Interval<true>(9));
  // Return 1 on error.
  return !correct;
};
