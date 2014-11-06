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

  Polyhedron p1;
  vector<Inequality*> eqs1;
  eqs1.push_back(new Inequality(Point(13,-1), Point(11,4)));
  eqs1.push_back(new Inequality(Point(12,2), Point(6,8)));
  eqs1.push_back(new Inequality(Point(0,-7), Point(5,-6)));
  eqs1.push_back(new Inequality(Point(2,-7), Point(12,-3)));
  eqs1.push_back(new Inequality(Point(6,-6), Point(13,-2)));
  eqs1.push_back(new Inequality(Point(9,-5), Point(12,-1)));
  eqs1.push_back(new Inequality(Point(11,-4), Point(13,1)));

  Interval<true> xBound1;
  Interval<true> yBound1;
  p1.addInequalitySet(xBound1, yBound1, eqs1);

  if (chat)
    cerr << "first polyhedron p1: ===========================" << endl << p1;
  Polyhedron p2=p1;
  vector<Inequality*> eqs2;
  eqs2.push_back(new Inequality(Point(0,-8), Point(13,-2)));

  if (chat) {
    cerr << "x in " << xBound1 << ", y in " << yBound1 << endl;
    cerr << "adding " << *eqs2[0] << endl;
  };
  p2.addInequalitySet(xBound1, yBound1, eqs2);

  if (chat) {
    cerr << "second polyhedron p2: ===========================" << endl << p2;

    cerr << "p1 includes p2:" << 
      (p1.includes(xBound1, yBound1, p2) ? "yes" : "no") << endl;
    cerr << "p2 includes p1:" << 
      (p2.includes(xBound1, yBound1, p1) ? "yes" : "no") << endl;
  };
  bool correct=
    p1.includes(xBound1, yBound1, p2) && 
    !p2.includes(xBound1, yBound1, p1);
  // Return 1 on error.
  return !correct;
};
