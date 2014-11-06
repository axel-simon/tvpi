/* Test program for inserting constraints when variables have
   non-trivial multiplicity.
 */

#include <iostream>
#include "affine.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;

  Domain dom;
  DomVar x2 = dom.createVariable();
  DomVar x3 = dom.createVariable();
  DomVar x4 = dom.createVariable();
  bool res;
  res = dom.setVariable(x2, 0, Interval<true>(0,65535));
  assert(res);

  if (chat) cerr << "initial domain:" << endl << dom;

  // add x3=4x2+4
  vector<LinComponent> ll;
  ll.push_back(LinComponent(4,x2));
  ll.push_back(LinComponent(-1,x3));
  res = dom.inequality(ll, 4, NULL, false);
  if (chat) cerr << "domain after adding 4x_2-x_3+4=0: res="
		   << res << endl << dom;
  assert(res==resChanged);

  // add x4=4x2
  ll.clear();
  ll.push_back(LinComponent(4,x2));
  ll.push_back(LinComponent(-1,x4));
  res = dom.inequality(ll, 0, NULL, false);
  if (chat) cerr << "domain after adding 4x_2-x_4=0: res="
		   << res << endl << dom;
  assert(res==resChanged);

  // Evaluate x2!=2
  Domain dom1(dom);
  ll.clear();
  ll.push_back(LinComponent(1,x2));
  res = dom1.inequality(ll, -1, NULL, false);
  if (chat) cerr << "domain after x_2<2: res=" << res << endl << dom1;
  assert(res==resChanged);

  Domain dom2(dom);
  ll.clear();
  ll.push_back(LinComponent(-1,x2));
  res = dom2.inequality(ll, 3, NULL, false);
  if (chat) cerr << "domain after x_2>2: res=" << res << endl << dom2;
  assert(res==resChanged);

  dom2.joinWiden(dom1, mpz_class(-1));
  if (chat) cerr << "join of both:" << endl << dom2;

  bool correct= dom2.entails(dom) && dom.entails(dom2);
  return !correct;
};
