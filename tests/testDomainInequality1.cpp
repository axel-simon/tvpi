/* Test program for checking adding an inequality where variables have
   a non-zero multiplicity.
 */

#define IN_TVPI_CHECK

#include <iostream>
#include "affine.hh"
#include "common.hh"
#include <assert.h>

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;
  Domain dom;
  DomVar x = dom.createVariable();
  DomVar y = dom.createVariable();
  DomVar z = dom.createVariable();

  dom.setVariable(x, 2, Interval<true>(8,64));
  dom.setVariable(y, 0, Interval<true>(1,200));
  dom.setVariable(z, 0, Interval<true>(1,200));
  dom.sane();

  Result res;
  if (chat) cerr << "Domain initially:" << endl << dom;

  {
    vector<LinComponent> ll;
    ll.push_back(LinComponent(1,z));
    ll.push_back(LinComponent(-1,x));
    res = dom.inequality(ll, 5, NULL, true);
  };

  if (chat) cerr << "Domain after adding z=x-5:" << endl << dom;
  assert(res==resChanged);

  {
    vector<LinComponent> ll;
    ll.push_back(LinComponent(1,y));
    ll.push_back(LinComponent(-1,z));
    res = dom.inequality(ll, -5, NULL, true);
  };

  if (chat) cerr << "Domain after adding y=z+5:" << endl << dom;
  assert(res==resChanged);

  {
    vector<LinComponent> ll;
    ll.push_back(LinComponent(-1,z));
    res = dom.inequality(ll, 0, NULL, false);
  };

  if (chat) cerr << "Domain after adding z>=0" << endl << dom;
  assert(res==resRedundant);

  {
    vector<LinComponent> ll;
    ll.push_back(LinComponent(-1,z));
    res = dom.inequality(ll, 4, NULL, false);
  };

  if (chat) cerr << "Domain after adding z>=4" << endl << dom;
  assert(res==resChanged);

  signed long lower, upper;
  Mult m;
  dom.queryValue(&lower, &upper, &m, y);
  if (chat) cerr << "lower = " << lower << endl;
  bool correct = lower==12 && m==0;
  // used to be 2, but propagation of multiplicity doesn't happen anymore

  // Return 1 on error.
  return !correct;
};
