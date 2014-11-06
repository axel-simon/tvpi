/* Test program for shrinking around integral grid.
 * Test for correct rounding of division results on negative numbers.
 */

#include <iostream>
#include "planar.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;
  Inequality e1 = Inequality(-4,-1,-11);
  Inequality e2 = Inequality(0,-1,-1);
  if (chat)
    cerr << "calculating cuts between e1=" << e1 << " and e2=" << e2 << endl;
  Inequality* cut1 = e1.calculateCut(e2);
  if (chat)
    if (cut1) cerr << "created cut " << *cut1 << endl;
  bool correct=
    cut1->getA()==mpz_class(-2) &&
    cut1->getB()==mpz_class(-1) &&
    cut1->getC()==mpz_class(-7);

  // Return 1 on error.
  return !correct;
};
