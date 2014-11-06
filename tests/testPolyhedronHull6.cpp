/* Test program to check entailment between polyhedra.
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
  xBound1.updateLower(mpq_class(1));
  xBound1.updateUpper(mpq_class(60));
  yBound1.updateLower(mpq_class(33));
  yBound1.updateUpper(mpq_class(65));
  vector<Inequality*> e1;
  e1.push_back(new Inequality(-1,1,63));
  e1.push_back(new Inequality(-1,-1,-65));
  e1.push_back(new Inequality(1,-1,-1));
  p1.addInequalitySet(xBound1, yBound1, e1);

  if (chat)
    cerr << "p1: x = " << xBound1 << ", y = " << yBound1 <<endl << p1;

  Polyhedron<true> p2;
  Interval<true> xBound2;
  Interval<true> yBound2;
  xBound2.updateLower(mpq_class(1));
  xBound2.updateUpper(mpq_class(56));
  yBound2.updateLower(mpq_class(29));
  yBound2.updateUpper(mpq_class(61));

  vector<Inequality*> e2;
  e2.push_back(new Inequality(-1,1,55));
  e2.push_back(new Inequality(-1,-1,-61));
  e2.push_back(new Inequality(1,-1,3));
  p2.addInequalitySet(xBound2, yBound2, e2);

  if (chat)
    cerr << "p2: x = " << xBound2 << ", y = " << yBound2 <<endl << p2;

  Polyhedron<true> p3 = Polyhedron<true>(p1, xBound1, yBound1,
			     p2, xBound2, yBound2);
  Interval<true> xBound3(xBound1, xBound2);
  Interval<true> yBound3(yBound1, yBound2);

  if (chat)
    cerr << "hull: x = " << xBound3 << ", y = " << yBound3 <<endl << p3;

  bool p3p1 = p3.includes(xBound1, yBound1, p1);
  if (chat) cerr << "hull includes p1 : " << (p3p1 ? "yes" : "no") << endl;
  bool p3p2 = p3.includes(xBound2, yBound2, p2);
  if (chat) cerr << "hull includes p2 : " << (p3p2 ? "yes" : "no") << endl;

  xBound3.widen(xBound1);
  yBound3.widen(yBound1);
  p3.widen(xBound3, yBound3, xBound1, yBound1, p1);

  if (chat)
    cerr << "widened: x = " << xBound3 << ", y = " << yBound3 <<endl << p3;
  
  bool p3p1_ = p3.includes(xBound1, yBound1, p1);
  if (chat) cerr << "widened includes p1 : " << (p3p1 ? "yes" : "no") << endl;
  bool p3p2_ = p3.includes(xBound2, yBound2, p2);
  if (chat) cerr << "widened includes p2 : " << (p3p2 ? "yes" : "no") <<endl;

  bool correct= p3p1 && p3p2 && p3p1_ && p3p2_;

  // Return 1 on error.
  return !correct;
};
