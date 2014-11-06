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
  Interval<true> xBound1= Interval<true>(7);
  Interval<true> yBound1= Interval<true>(8);
  Polyhedron<true> p2;
  Interval<true> xBound2= Interval<true>(3);
  Interval<true> yBound2= Interval<true>(4);
  Polyhedron<true> p3 = Polyhedron<true>(p1, xBound1, yBound1,
			     p2, xBound2, yBound2);
  Interval<true> xBound3 = Interval<true>(xBound1, xBound2);
  Interval<true> yBound3 = Interval<true>(yBound1, yBound2);

  if (chat)
    cerr << "hull polyhedron: x in " << xBound3 << ", y in " << yBound3
         << endl << p3;

  bool p3p1 = p3.includes(xBound1, yBound1, p1);
  if (chat) cerr << "p3 includes p1: " << p3p1 << endl;
  bool p3p2 = p3.includes(xBound2, yBound2, p2);
  if (chat) cerr << "p3 includes p2: " << p3p2 << endl;
  bool correct=p3p1 && p3p2;
    
  // Return 1 on error.
  return !correct;
};
