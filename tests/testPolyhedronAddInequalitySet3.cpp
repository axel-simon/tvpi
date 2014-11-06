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
  eqs.push_back(new Inequality(Point(13,-1), Point(11,4)));
  eqs.push_back(new Inequality(Point(12,2), Point(6,8)));
  eqs.push_back(new Inequality(Point(0,-7), Point(5,-6)));
  eqs.push_back(new Inequality(Point(2,-7), Point(12,-3)));
  eqs.push_back(new Inequality(Point(6,-6), Point(13,-2)));
  eqs.push_back(new Inequality(Point(9,-5), Point(12,-1)));
  eqs.push_back(new Inequality(Point(11,-4), Point(13,1)));

  Interval<true> xBound;
  Interval<true> yBound;
  yBound.updateLower(mpq_class(2));
  yBound.updateUpper(mpq_class(2));
  xBound.updateLower(mpq_class(5));
  xBound.updateUpper(mpq_class(5));
  p.addInequalitySet(xBound, yBound, eqs);
  //cout << p << "x: " << xBound << " y: " << yBound << endl;

  // Check that all inequalities were removed.
  bool correct=p.getNoOfInequalities()==0;
  // Return 1 on error.
  return !correct;
};
