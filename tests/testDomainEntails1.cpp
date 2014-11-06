/* Test program for checking entailment in wonderful ways.
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
  DomVar vars[4];
  for (int i=0; i<4; i++) vars[i]=dom.createVariable();

  dom.setVariable(vars[0], 0, Interval<true>(0,10));
  dom.setVariable(vars[1], 0, Interval<true>(0,20));
  dom.setVariable(vars[2], 0, Interval<true>(0,30));
  dom.setVariable(vars[3], 0, Interval<true>(0,40));
  dom.sane();
  if (chat) cerr << "Domain after setting four ranges:" << endl << dom;

  Domain dom2 = dom;

  {
    vector<LinComponent> ll;
    ll.push_back(LinComponent(1,vars[0]));
    bool res1 = dom2.inequality(ll, -5, NULL, false);
    assert(res1);
    dom2.demote();
    dom2.sane();
  };
  if (chat)
    cerr << "Restriciting range of v_0 and v_1 of copy:" << endl << dom2;

  // Force the variables into the relational domain.
  {
    vector<LinComponent> ll;
    ll.push_back(LinComponent(1,vars[0]));
    ll.push_back(LinComponent(1,vars[1]));
    ll.push_back(LinComponent(1,vars[2]));
    ll.push_back(LinComponent(1,vars[3]));
    Result res1 = dom.inequality(ll, -100, NULL, false);
    assert(res1==resRedundant);
    dom.sane();
  };

  {
    vector<LinComponent> ll;
    ll.push_back(LinComponent(1,vars[0]));
    ll.push_back(LinComponent(-1,vars[1]));
    bool res1 = dom.inequality(ll, -5, NULL, false);
    assert(res1);
    dom.sane();
  };
  if (chat)
    cerr << "Adding v_0 - v_1 <= 5 in original:" << endl << dom;

  bool correct = dom2.entails(dom);

  // Return 1 on error.
  return !correct;
};
