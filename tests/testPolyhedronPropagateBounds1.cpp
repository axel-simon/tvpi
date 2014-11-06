/* Test program for propagation of updated bounds.
 */

#include <iostream>
#include "planar.hh"
#include "polyhedron.hh"

using namespace std;
using namespace Tvpi;

  
int main() {
  Polyhedron<true> p;

  vector<Inequality*> eqs;
  eqs.push_back(new Inequality(Point(13,-6), Point(12,-1)));
  eqs.push_back(new Inequality(Point(12,-1), Point(11,2)));
  eqs.push_back(new Inequality(Point(11,2), Point(9,5)));
  eqs.push_back(new Inequality(Point(9,5), Point(2,8)));
  eqs.push_back(new Inequality(Point(0,8), Point(-7,7)));
  eqs.push_back(new Inequality(Point(-7,7), Point(-10,6)));
  eqs.push_back(new Inequality(Point(-10,6), Point(-12,3)));
  eqs.push_back(new Inequality(Point(-12,3), Point(-13,-4)));
  eqs.push_back(new Inequality(Point(-13,-4), Point(-8,-7)));

  Interval<true> xBound;
  Interval<true> yBound;
  p.addInequalitySet(xBound, yBound, eqs);
  //  cout << "before:\n" << p
  //       << "x in " << xBound << ", y in " << yBound << endl;

  bool correct=true;

  Interval<true> xBoundA = Interval<true>(xBound);
  Interval<true> yBoundA = Interval<true>(yBound);
  // if x>=5 then y<7
  xBoundA.updateLower(5);
  p.propagateXBounds(xBoundA, yBoundA);
  correct = !yBoundA.updateUpper(7) && correct;

  Interval<true> xBoundB = Interval<true>(xBound);
  Interval<true> yBoundB = Interval<true>(yBound);
  // if x<=10 then -6 < y <= 6
  xBoundB.updateUpper(-10);
  p.propagateXBounds(xBoundB, yBoundB);
  correct = !yBoundB.updateUpper(7) && !yBoundB.updateLower(-6) && correct;

  Interval<true> xBoundC = Interval<true>(xBound);
  Interval<true> yBoundC = Interval<true>(yBound);
  // if y>=0 then -13 < x < 12
  yBoundC.updateLower(0);
  p.propagateYBounds(xBoundC, yBoundC);
  correct = !xBoundC.updateUpper(12) && !xBound.updateLower(-13) && correct;

  Interval<true> xBoundD = Interval<true>(xBound);
  Interval<true> yBoundD = Interval<true>(yBound);
  // if y<=-5 then -12  < y
  yBoundD.updateUpper(-5);
  p.propagateYBounds(xBoundD, yBoundD);
  correct = !xBoundD.updateLower(-12) && correct;

  // Return 1 on error.
  return !correct;
};
