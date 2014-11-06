/* Test program for checking sorting and ordering.
 */

#include <iostream>
#include "polyhedron.hh"

using namespace std;
using namespace Tvpi;

  
int main() {
  Polyhedron p;
  Point p0 = Point(6,0);
  Point p1 = Point(4,2);
  Point p2 = Point(1,3);
  Point p3 = Point(-1,3);
  Point p4 = Point(-4,2);
  Point p5 = Point(-6,0);
  Point p6 = Point(-4,-2);
  Point p7 = Point(-1,-3);
  Point p8 = Point(1,-3);
  Point p9 = Point(4,-2);

  vector<Inequality*> eqs;
  eqs.push_back(new Inequality(p0,p1));
  eqs.push_back(new Inequality(p1,p2));
  eqs.push_back(new Inequality(p3,p4));
  eqs.push_back(new Inequality(p4,p5));
  eqs.push_back(new Inequality(p5,p6));
  eqs.push_back(new Inequality(p6,p7));
  eqs.push_back(new Inequality(p8,p9));
  eqs.push_back(new Inequality(p9,p0));

  Interval<false> xBound;
  Interval<false> yBound;
  yBound.updateUpper(mpq_class(3));
  yBound.updateLower(mpq_class(-3));
  p.addInequalitySet<false>(xBound, yBound, eqs, 0);

  Polyhedron q=p;
  //cout << "before:\n" << p;
  p.swapVars();
  //cout << "swapped:\n" << p;
  p.swapVars();
  //cout << "swapped again:\n" << p;

  bool correct= p.includes(xBound, yBound, q) && q.includes(xBound, yBound, p);

  // Return 1 on error.
  return !correct;
};
