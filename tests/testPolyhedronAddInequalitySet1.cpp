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
  //  yBound.updateUpper(mpq_class(1));
  xBound.updateUpper(mpq_class(10));
  p.addInequalitySet(xBound, yBound, eqs);
  //cout << eqs.size() << " new inequalities:" << endl
  //     << p << "x: " << xBound << " y: " << yBound << endl;

  bool correct=true;

  // Try to insert 2x - 3y <= 30. This is redundant in the integral
  // polyhedron.
  eqs.clear();
  eqs.push_back(new Inequality(2, -3, 30));
  p.addInequalitySet(xBound, yBound, eqs);
  if (eqs.size()!=0) correct=false;

  // Try to insert 2x - 3y <= 29 . This is already a proper integral cut.
  eqs.clear();
  eqs.push_back(new Inequality(2, -3, 29));
  p.addInequalitySet(xBound, yBound, eqs);
  if (eqs.size()!=1) correct=false;

  // Return 1 on error.
  return !correct;
};
