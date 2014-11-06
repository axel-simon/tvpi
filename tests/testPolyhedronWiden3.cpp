/* Test program to check widening.
 */

#include <iostream>
#include "interval.hh"
#include "planar.hh"
#include "polyhedron.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;
  Polyhedron<true> p1;
  Interval<true> xBound1;
  Interval<true> yBound1;
  xBound1.updateLower(mpq_class(-77));
  xBound1.updateUpper(mpq_class(39));
  yBound1.updateLower(mpq_class(1));
  yBound1.updateUpper(mpq_class(39));
  vector<Inequality*> e1;
  e1.push_back(new Inequality(-1,1,78));
  e1.push_back(new Inequality(1,-1,0));
  p1.addInequalitySet(xBound1, yBound1, e1);

  if (chat)
    cerr << "p1: x = " << xBound1 << ", y = " << yBound1 <<endl << p1;

  Polyhedron<true> p2;
  Interval<true> xBound2;
  Interval<true> yBound2;
  vector<Inequality*> e2;
  xBound2.updateLower(mpq_class(-116));
  xBound2.updateUpper(mpq_class(39));
  yBound2.updateLower(mpq_class(1));
  yBound2.updateUpper(mpq_class(39));
  e2.push_back(new Inequality(-1,1,117));
  e2.push_back(new Inequality(1,-1,0));
  p2.addInequalitySet(xBound2, yBound2, e2);

  if (chat)
    cerr << "p2: x = " << xBound2 << ", y = " << yBound2 <<endl << p2;

  Polyhedron<true> p3=p2;
  Interval<true> xBound3 = xBound2;
  Interval<true> yBound3 = yBound2;
  
  mpz_class extrapolate = 3;
  xBound3.widen(xBound1, extrapolate);
  yBound3.widen(yBound1, extrapolate);
  if (chat)
    cerr << "extrapolated bounds: x = " << xBound3 << ", y = "
	 << yBound3 << endl;
  p3.widen(xBound3, yBound3, xBound1, yBound1, p1, extrapolate);

  extrapolate = 0;
  xBound2.widen(xBound1, extrapolate);
  yBound2.widen(yBound1, extrapolate);
  p2.widen(xBound2, yBound2, xBound1, yBound1, p1, extrapolate);

  if (chat)
    cerr << "extrapolated: x = " << xBound3 << ", y = " << yBound3 << endl
	 << p3;
  if (chat)
    cerr << "widened: x = " << xBound2 << ", y = " << yBound2 << endl << p2;
  
  // Full widening needs to be coarser than extrapolation.
  bool correct =
    p2.includes(xBound3, yBound3, p3) &&
    !p3.includes(xBound2, yBound2, p2);

  // Return 1 on error.
  return !correct;
};
