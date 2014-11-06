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
  vector<Inequality*> first;
  first.push_back(new Inequality(Point(26,3), Point(22,11)));
  first.push_back(new Inequality(Point(22,11), Point(20,13)));
  first.push_back(new Inequality(Point(20,13), Point(17,15)));
  first.push_back(new Inequality(Point(17,15), Point(11,14)));
  first.push_back(new Inequality(Point(11,14), Point(1,6)));
  first.push_back(new Inequality(Point(1,6), Point(2,4)));
  first.push_back(new Inequality(Point(2,4), Point(4,3)));
  first.push_back(new Inequality(Point(4,3), Point(13,2)));
  first.push_back(new Inequality(Point(13,2), Point(26,3)));
  
  Interval<true> xBound1;
  Interval<true> yBound1;
  p1.addInequalitySet(xBound1, yBound1, first);

  if (chat)
    cerr << "first polyhedron p1: x=" << xBound1 << ", y=" << yBound1
	 << " ===========================" << endl << p1;

  Polyhedron<true> p2;
  vector<Inequality*> second;
  second.push_back(new Inequality(Point(26,3), Point(21,13)));
  second.push_back(new Inequality(Point(21,13), Point(17,15)));
  second.push_back(new Inequality(Point(17,15), Point(11,14)));
  second.push_back(new Inequality(Point(11,14), Point(6,12)));
  second.push_back(new Inequality(Point(6,12), Point(3,10)));
  second.push_back(new Inequality(Point(3,10), Point(1,6)));
  second.push_back(new Inequality(Point(1,6), Point(0,1)));
  second.push_back(new Inequality(Point(0,1), Point(26,3)));
  
  Interval<true> xBound2;
  Interval<true> yBound2;
  p2.addInequalitySet(xBound2, yBound2, second);

  if (chat) {
    cerr << "second polyhedron p2: x=" << xBound2 << ", y=" << yBound2
	 << " ===========================" << endl << p2;
    cerr << "p1 includes p2:" << 
      (p1.includes(xBound2, yBound2, p2) ? "yes" : "no") << endl;
    cerr << "p2 includes p1:" << 
      (p2.includes(xBound1, yBound1, p1) ? "yes" : "no") << endl;
  };


  Interval<true> xBound3(yBound1);
  Interval<true> yBound3(xBound1);
  Polyhedron<true> p3=p1;
  p3.swapVars();
  Interval<true> xBound4(yBound2);
  Interval<true> yBound4(xBound2);
  Polyhedron<true> p4=p2;
  p4.swapVars();

  mpz_class extrapolate = 0;
  xBound2.widen(xBound1, extrapolate);
  yBound2.widen(yBound1, extrapolate);
  p2.widen(xBound2, yBound2, xBound1, yBound1, p1, extrapolate);
  xBound4.widen(xBound3, extrapolate);
  yBound4.widen(yBound3, extrapolate);
  p4.widen(xBound4, yBound4, xBound3, yBound3, p3, extrapolate);

  p4.swapVars();

  if (chat) {
    cerr << "p2 widen p1: x = " << xBound2 << ", y = " << yBound2
	 << endl << p2;
    cerr << "p4 widen p3:  x = " << yBound4 << ", y = " << xBound4
	 << endl << p4;
  };

  bool correct=
    p4.includes(yBound2, xBound2, p2) &&
    p2.includes(yBound4, xBound4, p4);

  // Return 1 on error.
  return !correct;
};
