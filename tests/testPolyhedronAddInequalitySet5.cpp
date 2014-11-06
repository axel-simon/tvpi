/* Test program for removal of inequalities when tightening bounds.
 */

#include <iostream>
#include "planar.hh"
#include "polyhedron.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;
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
  xBound.updateUpper(14);
  yBound.updateUpper(8);
  p.addInequalitySet(xBound, yBound, eqs);
  if (chat) cout << "before:\n" << p
		 << "x in " << xBound << ", y in " << yBound << endl;

  // -20 is the lower bound if x is at most 14. Make sure setting this
  // lower bound explicitly doesn't change anything.
  bool correct=!yBound.updateLower(mpq_class(-20));

  // Return 1 on error.
  return !correct;
};
