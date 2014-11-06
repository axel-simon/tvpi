/* Test program for replacing one variable with another.
 */

#define IN_TVPI_CHECK

#include <iostream>
#include "affine.hh"
#include "common.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;
  Domain dom;
  DomVar vars[7];
  for (int i=0; i<7; i++) vars[i]=dom.createVariable();

  Interval<true> mone  = Interval<true>(-1);
  Interval<true> zero  = Interval<true>(0l);
  Interval<true> one   = Interval<true>(1);
  Interval<true> two   = Interval<true>(2);
  Interval<true> three = Interval<true>(3);
  Interval<true> four  = Interval<true>(4);

  dom.setVariable(vars[1], 0, Interval<true>(zero,one));
  dom.setVariable(vars[2], 0, Interval<true>(mone,zero));
  dom.setVariable(vars[3], 0, Interval<true>(two, three));
  dom.sane();
  if (chat) cerr << "Domain after setting three ranges:" << endl << dom;

  vector<LinComponent> ll;
  ll.push_back(LinComponent(3,vars[1]));
  ll.push_back(LinComponent(2,vars[3]));
  ll.push_back(LinComponent(5,vars[2]));
  ll.push_back(LinComponent(-1,vars[4]));
  ll.push_back(LinComponent(7,vars[5]));
  ll.push_back(LinComponent(-7,vars[5]));
  ll.push_back(LinComponent(2,vars[3]));
  ll.push_back(LinComponent(2,vars[3]));
  bool res1 = dom.inequality(ll, 0, NULL, true);
  assert(res1);
  dom.sane();
  if (chat)
    cerr << "Domain after inserting first linear equality:" << endl << dom;

  dom.setVariable(vars[5], 0, four);
  dom.sane();

  ll.clear();
  ll.push_back(LinComponent(3,vars[1]));
  ll.push_back(LinComponent(9,vars[2]));
  ll.push_back(LinComponent(6,vars[3]));
  ll.push_back(LinComponent(4,vars[4]));
  ll.push_back(LinComponent(2,vars[5]));
  ll.push_back(LinComponent(-1,vars[6]));
  bool res2 = dom.inequality(ll, 0, NULL, true);
  dom.sane();

  Domain dom2= dom;
  if (chat) 
    cerr << "Domain after inserting second linear equality:" << endl << dom;

  dom.update(vars[5],vars[4]);
  dom.sane();

  if (chat) cerr << "setting x_" << vars[4] << " := x_" << vars[5]
		 << " yields " << endl << dom;

  ll.clear();
  bool correct = res1 && res2;
  if (chat) cerr << "again " << dom;

  // Return 1 on error.
  return !correct;
};
