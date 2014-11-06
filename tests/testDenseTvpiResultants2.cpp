/* Test program for the closure and resultants calculation.

  In particular, this test represents the rational version of the
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
  TvpiVar ident[] = {0,1,2 };
  DenseTvpi<true> dom;
  DomVar x = dom.createVariable();
  DomVar y = dom.createVariable();
  DomVar z = dom.createVariable();
  Point xzA = Point(8,-2);
  Point xzB = Point(8,1);
  Point xzC = Point(7,5);
  Point xzD = Point(6,7);
  Point xzE = Point(5,7);
  Point xzF = Point(2,8);
  Point xzG = Point(0,8);
  Point yzA = Point(0,-7);
  Point yzB = Point(2,-7);
  Point yzC = Point(4,-7);
  Point yzD = Point(6,-6);
  Point yzE = Point(7,-5);
  Point yzF = Point(8,-2);
  Point yzG = Point(8,2);
  Inequality* eqsXZ[] = {
    //new Inequality(xzA, xzB),
    new Inequality(xzB, xzD),
    new Inequality(xzC, xzF),
    new Inequality(xzE, xzG),
    0
  };
  Inequality* eqsYZ[] = {
    new Inequality(yzA, yzD),
    new Inequality(yzB, yzE),
    new Inequality(yzC, yzF),
    0,//new Inequality(yzF, yzG),
  };
  dom.addInequalities(x,z, 3, eqsXZ);
  if (chat) cout << "Before:" << endl << dom;
  dom.addInequalities(y,z, 3, eqsYZ);
  if (chat) cout << "After:" << endl << dom;
  if (chat) cout << "Intersection points:" << endl;
  Polyhedron<true> xyP = dom.getProjection(x,y);
  if (chat) for(size_t i=0; i+1<xyP.getNoOfInequalities(); i++)
    cout << Point(*xyP[i], *xyP[i+1]) << endl;

  // The point <8,10,-> is not entailed in the closure.
  DenseTvpi<true> dom2;
  dom2.createVariable(8);
  dom2.createVariable(10);
  dom2.createVariable();
  if (chat) cout << "Testing:" << endl << dom2;
  bool correct=!dom.includes(dom2, ident) && xyP.getNoOfInequalities()==8;
  return !correct;
};
