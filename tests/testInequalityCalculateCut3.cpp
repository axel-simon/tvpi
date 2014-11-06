/* Test program for shrinking around integral grid.
 */

#include <iostream>
#include "planar.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;
  Inequality e1 = Inequality(Point(14,2),Point(-4,1));
  Inequality e2 = Inequality(Point(12,2),Point(0,1));
  Inequality e3 = Inequality(Point(6,2),Point(0,0));
  assert(e1.lessThanPi(e2));
  if (chat)
    cerr << "calculating cuts between e1=" << e1 << " and e2=" << e2 << endl;
  Inequality* cut1 = e1.calculateCut(e2);
  Inequality* lastCut1=0;
  while (cut1) {
    if (chat)
      cerr << "created cut " << *cut1 << endl;
    delete lastCut1;
    lastCut1=cut1;
    cut1 = cut1->calculateCut(e2);
  };
  Inequality* cut2 = e2.calculateCut(e3);
  if (chat)
    cerr << "calculating cuts between e2=" << e2 << " and e3=" << e3 << endl;
  Inequality* lastCut2=0;
  while (cut2) {
    if (chat)
	cerr << "created cut " << *cut2 << endl;
    delete lastCut2;
    lastCut2=cut2;
    cut2 = cut2->calculateCut(e3);
  };
  Point p =Point(*lastCut1, *lastCut2);
  if (chat)
    cerr << "intersection point of cuts is " << p << endl;
  bool correct=p.getX().get_den()==5 && p.getY().get_den()==5;

  // Return 1 on error.
  return !correct;
};
