/* Test program to perform an if-test with splits and rejoins the state space.
 */

#include <iostream>
#include "interval.hh"
#include "affine.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;

  Domain dOrig;

  dOrig.createVariable();
  dOrig.createVariable();
  dOrig.createVariable();
  DomVar y = dOrig.createVariable();
  DomVar x = dOrig.createVariable();
  // A bit of fiddleing to make this match the funnycount.c testcase.
  dOrig.update(x,x-2);
  dOrig.setVariableToTop(x);
  x=x-2;
  Result sat;
  vector<LinComponent> ll;
  ll.push_back(LinComponent(-1, x));
  ll.push_back(LinComponent(6, y));
  sat = dOrig.inequality(ll, -175);
  ll.clear();
  assert(sat==resChanged);
  ll.push_back(LinComponent(-1, x));
  ll.push_back(LinComponent(1, y));
  sat = dOrig.inequality(ll, 0);
  ll.clear();
  assert(sat==resChanged);
  ll.push_back(LinComponent(36, x));
  ll.push_back(LinComponent(-41, y));
  sat = dOrig.inequality(ll, 0);
  ll.clear();
  assert(sat==resChanged);

  if (chat) cerr << "domain before split: " << endl << dOrig;

  // Simulate an "if (x!=35)" statement
  Domain dLess(dOrig);
  Domain dGreater(dOrig);
  ll.push_back(LinComponent(1,x));
  sat = dLess.inequality(ll, -34);
  assert(sat==resChanged);
  if (chat) cerr << "x<35 domain: " << endl << dLess;

  ll.push_back(LinComponent(-1,x));
  sat = dGreater.inequality(ll, 36);
  assert(sat==resChanged);
  if (chat) cerr << "x>35 domain: " << endl << dGreater;
  
  Domain dRes(dLess);
  dRes.joinWiden(dGreater, -1);

  if (chat) cerr << "hull of both: " << endl << dRes;
  bool resLess = dLess.entails(dRes);
  if (chat) cerr << "hull includes x<35 domain: " << resLess << endl;
  bool resGreater = dGreater.entails(dRes);
  if (chat) cerr << "hull includes x>35 domain: " << resGreater << endl;
  bool origRes = dRes.entails(dOrig);
  if (chat) cerr << "original includes hull: " << origRes << endl;

  bool correct= resLess && resGreater && origRes;

  // Return 1 on error.
  return !correct;
};
