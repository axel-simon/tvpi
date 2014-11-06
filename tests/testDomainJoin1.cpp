/* Test program that joins two simple domains and checks multiplicity.
 */

#define IN_TVPI_CHECK

#include <iostream>
#include "affine.hh"
#include "common.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;
  Domain domA, domB;
  DomVar x,y;
  x=Domain::createVariable();
  y=Domain::createVariable();

  Interval<true> a = Interval<true>(4096);
  Interval<true> b = Interval<true>(4080);
  Interval<true> zero = Interval<true>(0l);
  Interval<true> one = Interval<true>(1);

  bool res;
  res = domA.setVariable(x, 12, a);
  assert(res);
  res = domA.setVariable(y,  255, zero);
  assert(res);
  res = domB.setVariable(x, 4, b);
  assert(res);
  res = domB.setVariable(y, 0, one);
  assert(res);

  domA.sane();
  domB.sane();
  if (chat) cerr << "Domains after setting two ranges:"
		 << endl << domA << endl << domB;

  Domain domAA = domA;
  Domain domBB = domB;

  domA.joinWiden(domB, mpz_class(-1));
  domA.sane();

  if (chat)
    cerr << "Domain after join first way:" << endl << domA;

  domBB.joinWiden(domAA, mpz_class(-1));
  domBB.sane();

  if (chat)
    cerr << "Domain after join second way:" << endl << domBB;

  if (chat)
    cerr << "Calculating entailment between both polyhedra." << endl;

  bool correct = domB.entails(domA) && domAA.entails(domBB) &&
		 domA.entails(domBB) && domBB.entails(domA);

  // Return 1 on error.
  return !correct;
};
