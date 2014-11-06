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
  yBound1.updateLower(mpq_class(-1));

  vector<Inequality*> eqs1;
  eqs1.push_back(new Inequality(Point(-1,1), Point(-3,3)));
  eqs1.push_back(new Inequality(Point(-3,3), Point(-6,4)));
  eqs1.push_back(new Inequality(Point(-6,4), Point(-8,2)));
  eqs1.push_back(new Inequality(Point(-8,2), Point(-7,-1)));
  eqs1.push_back(new Inequality(Point(-4,-1), Point(-1,1)));

  p1.addInequalitySet(xBound1, yBound1, eqs1);
  if (argc>1) cerr << "first polyhedron: x in " << xBound1 << ", y in " << yBound1
                   << endl << p1;

  Polyhedron<true> p2;
  Interval<true> xBound2;
  xBound2.updateUpper(mpq_class(14));
  Interval<true> yBound2;
  vector<Inequality*> eqs2;
  eqs2.push_back(new Inequality(Point(14,5), Point(12,8)));
  eqs2.push_back(new Inequality(Point(12,8), Point(9,9)));
  eqs2.push_back(new Inequality(Point(9,9), Point(5,8)));
  eqs2.push_back(new Inequality(Point(5,8), Point(3,6)));
  eqs2.push_back(new Inequality(Point(3,6), Point(2,2)));
  eqs2.push_back(new Inequality(Point(2,2), Point(4,-2)));
  eqs2.push_back(new Inequality(Point(4,-2), Point(8,-3)));
  eqs2.push_back(new Inequality(Point(8,-3), Point(12,-2)));
  eqs2.push_back(new Inequality(Point(12,-2), Point(14,1)));

  p2.addInequalitySet(xBound2, yBound2, eqs2);
  if (argc>1) cerr << "second polyhedron: x in " << xBound2 << ", y in " << yBound2
                   << endl << p2;

  Interval<true> xBound3 = Interval<true>(xBound1, xBound2);
  Interval<true> yBound3 = Interval<true>(yBound1, yBound2);
  Polyhedron<true> p3 = Polyhedron<true>(p1, xBound1, yBound1,
			     p2, xBound2, yBound2);

  if (argc>1) cerr << "hull polyhedron: x in " << xBound3 << ", y in " << yBound3
                   << endl << p3;

  bool p3p1 = p3.includes(xBound1, yBound1, p1);
  if (argc>1) cerr << "p3 includes p1: " << p3p1 << endl;
  bool p3p2 = p3.includes(xBound2, yBound2, p2);
  if (argc>1) cerr << "p3 includes p2: " << p3p2 << endl;
  bool correct=p3p1 && p3p2;
    
  // Return 1 on error.
  return !correct;
};
