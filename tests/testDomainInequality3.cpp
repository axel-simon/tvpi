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
  res = dom.inequality(ll, 4, NULL, true);
  if (chat) cerr << "domain after adding 4x_2-x_3+4=0: res="
		   << res << endl << dom;
  assert(res==resChanged);

  // add x4=4x2
  ll.clear();
  ll.push_back(LinComponent(4,x2));
  ll.push_back(LinComponent(-1,x4));
  res = dom.inequality(ll, 0, NULL, true);
  if (chat) cerr << "domain after adding 4x_2-x_4=0: res="
		   << res << endl << dom;
  assert(res==resChanged);

  signed long lower, upper;
  Mult m;
  int bits;
  bits = dom.queryValue(&lower, &upper, &m, x2);
  if (chat) cerr << "bounds of x2:  lower = " << lower << ", upper = " << upper
		 << ", mult = " << (int) m << endl;
  bool correct= bits==0x7 && lower==0 && upper==65535 && m==0;
  bits = dom.queryValue(&lower, &upper, &m, x3);
  if (chat) cerr << "bounds of x3:  lower = " << lower << ", upper = " << upper
		 << ", mult = " << (int) m << endl;
  correct= correct && bits==0x7 && lower==4 && upper==262144 && m==2;
  bits = dom.queryValue(&lower, &upper, &m, x4);
  if (chat) cerr << "bounds of x4:  lower = " << lower << ", upper = " << upper
		 << ", mult = " << (int) m << endl;
  correct= correct && bits==0x7 && lower==0 && upper==262140 && m==2;
  return !correct;
};
