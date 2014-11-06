/* Test program for the closure and resultants calculation.

  In particular, this test represents the integral version of the
  closure/integral hull counter example.
 */

#include <iostream>
#include "planar.hh"
#include "polyhedron.hh"
#include "tvpi.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;
  TvpiVar ident[] = { 0, 1, 2, 3, 4, 5 };
  DenseTvpi<false> dom;
  DomVar x = dom.createVariable();
  DomVar y = dom.createVariable();
  DomVar z = dom.createVariable();
  Point xzB = Point(8,1);
  Point xzCC = Point(7,4);
  Point xzDD = Point(5,6);
  Point xzD = Point(6,7);
  Point xzE = Point(5,7);
  Point xzEE = Point(3,7);
  Point xzG = Point(0,8);
  Point yzA = Point(0,-7);
  Point yzC = Point(4,-7);
  Point yzBB = Point(4,-6);
  Point yzD = Point(6,-6);
  Point yzF = Point(8,-2);
  Inequality* eqsXZ[] = {
    new Inequality(xzB, xzD),
    new Inequality(xzCC, xzDD),
    new Inequality(xzDD, xzEE),
    new Inequality(xzEE, xzG),
    new Inequality(xzE, xzG)
  };
  Inequality* eqsYZ[] = {
    new Inequality(yzA, yzD),
    new Inequality(yzA, yzBB),
    new Inequality(yzBB, yzF),
    new Inequality(yzC, yzF)
  };
  dom.addInequalities(x,z, 5, eqsXZ);
  dom.addInequalities(y,z, 4, eqsYZ);
  if (chat) cout << dom;
  if (chat) cout << "Intersection points:" << endl;
  Polyhedron xyP = dom.getProjection(x,y);
  if (chat) for(size_t i=0; i<xyP.getNoOfInequalities()-1; i++)
    cout << Point(*xyP[i], *xyP[i+1]) << endl; 

  // The point <4,15,5.5> is just about not entailed in the closure of
  // the integral hull.
  DenseTvpi<false> dom2;
  DomVar x2=dom2.createVariable(4);
  dom2.createVariable(15);
  // A clumsy way to create z=5.5:
  DomVar z2=dom2.createVariable();
  Inequality* eqsZ2[] = {
    new Inequality(1, 2, 15),
    new Inequality(-1, -2, -15)
  };
  dom2.addInequalities(x2, z2, 2, eqsZ2);
  bool correct=!dom.includes(dom2, ident);
  return !correct;

};
