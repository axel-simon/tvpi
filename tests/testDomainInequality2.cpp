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
  DomVar x5 = dom.createVariable();
  DomVar x6 = dom.createVariable();
  bool res;
  res = dom.setVariable(x2, 2, Interval<true>(3,19));
  assert(res);
  res = dom.setVariable(x3, 4, Interval<true>(-7,42));
  assert(res);
  res = dom.setVariable(x4, 3, Interval<true>(-8,43));
  assert(res);
  res = dom.setVariable(x5, 1, Interval<true>(1,5));
  assert(res);
  res = dom.setVariable(x6, 3, Interval<true>(2,69));
  assert(res);

  if (chat) cerr << "initial domain:" << endl << dom;

  {
    Domain dom2 = dom;
    vector<LinComponent> ll;
    ll.push_back(LinComponent(-2,x2)); // -32
    ll.push_back(LinComponent( 2,x3)); //  64
    ll.push_back(LinComponent(-1,x4)); // -32
    ll.push_back(LinComponent(-4,x5)); // -16
    ll.push_back(LinComponent( 1,x6)); //  56
    res = dom2.inequality(ll, -48, NULL, true);
    if (chat) cerr << "domain after adding -2x_2+2x_3-x_4-4x_5+x_6-48=0: res="
		   << res << endl << dom2;
    assert(res==resChanged);
  };

  {
    Domain dom2 = dom;
    vector<LinComponent> ll;
    ll.push_back(LinComponent(-2,x2)); // -32
    ll.push_back(LinComponent( 2,x3)); //  64
    ll.push_back(LinComponent(-1,x4)); // -32
    ll.push_back(LinComponent(-4,x5)); // -16
    ll.push_back(LinComponent( 1,x6)); //  56
    res = dom2.inequality(ll, -44, NULL, true);
    assert(res==resUnsatisfiable);
  };

  {
    Domain dom2 = dom;
    vector<LinComponent> ll;
    ll.push_back(LinComponent(-2,x2)); // -16
    ll.push_back(LinComponent( 1,x3)); //  32
    ll.push_back(LinComponent(-1,x5)); // -16
    res = dom2.inequality(ll, 0, NULL, true);
    if (chat) cerr << "domain after adding -2x_2-x_4-4x_5=0: res="
		   << res << endl << dom2;
    assert(res==resChanged);
  };

  bool correct=true;
  return !correct;
};
