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

  Interval<true> xBound1;
  Interval<true> yBound1;
  Polyhedron p1;
  vector<Inequality*> eqs1;
  eqs1.push_back(new Inequality(-10, 3, -40623));
  eqs1.push_back(new Inequality(-7, 2, -28437));
  p1.addInequalitySet(xBound1, yBound1, eqs1);

  // Create an intersection point that's not on the Z grid.
  yBound1.updateLower(3);
  yBound1.updateUpper(59);
  xBound1.updateLower(4064);
  xBound1.updateUpper(4080);


  if (chat)
    cerr << "polyhedron p1: ===========================" << endl
	 << "x : " << xBound1 << ", y : " << yBound1 << endl << p1;

  Interval<true> xBound2;
  xBound2.updateLower(4063);
  xBound2.updateUpper(4080);
  Interval<true> yBound2;
  Polyhedron p2;
  vector<Inequality*> eqs2;
  eqs2.push_back(new Inequality(-10, 3, -40623));
  eqs2.push_back(new Inequality(-7, 2, -28437));
  p2.addInequalitySet(xBound2, yBound2, eqs2);

  // Force the polyhedron to have meaningless inequalities.
  yBound2.updateLower(1);
  yBound2.updateUpper(1);

  if (chat)
    cerr << "polyhedron p2: ===========================" << endl
	 << "x : " << xBound2 << ", y : " << yBound2 << endl << p2;

  Interval<true> xBoundHull(xBound1, xBound2);
  Interval<true> yBoundHull(yBound1, yBound2);
  Polyhedron hull(p1, xBound1, yBound1,
		  p2, xBound2, yBound2);

  
  if (chat) {
    cerr << "convex hull of p1 and p2: ===========================" << endl
	 << "x : " << xBoundHull << ", y : " << yBoundHull << endl << hull;

    cerr << "hull includes p1:" << 
      (hull.includes(xBound1, yBound1, p1) ? "yes" : "no") << endl;
    cerr << "hull includes p2:" << 
      (hull.includes(xBound2, yBound2, p2) ? "yes" : "no") << endl;
  };

  bool correct =
    hull.includes(xBound1, yBound1, p1) && 
    hull.includes(xBound2, yBound2, p2);

  // Return 1 on error.
  return !correct;
};
