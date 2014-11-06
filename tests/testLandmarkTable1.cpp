/* Test program for adding redundant inequalities and widen until they are
   non-redundant.
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

  LandmarkTable lm;

  DomVar x_i, x_o;
  dom.createVariable();
  x_i=dom.createVariable();
  x_o=dom.createVariable();

  dom.setVariable(x_i, 0, Interval<true>(32));
  dom.setVariable(x_o, 0, Interval<true>(128));

  DomVar x_iTemp = dom.createVariable();
  DomVar x_oTemp = dom.createVariable();

  if (chat) cerr << "Domain before loop:" << endl << dom;

  Domain* prev = NULL;

  while (true) {

    if (prev) {
      // Ask for the widening mode.
      WideningInfo wi;
      lm.calcNoOfIterations(wi);
      mpz_class steps = wi.getIterations();
      if (chat) cerr << lm;
      switch (sgn(steps)) {
        case -1: {
	  if (chat) cerr << "No widening, just join." << endl;
	  lm.promote();
        }; break;
        case 0: {
	  if (chat) cerr << "Imprecise widening." << endl;
        }; break;
        case 1: {
	  steps--;
	  if (chat) cerr << "Widening by " << steps << " steps" << endl;
	  lm.clear();
        }; break;
      };
      if (dom.entails(*prev)) {
	if (chat) cerr << "stable:" << endl << *prev << endl;
	break;
      };

      dom.joinWiden(*prev, steps);
    
      if (chat) cerr << "Domain after join:" << endl << dom;
    } else if (chat) cerr << "No join after first iteration" << endl;

    delete prev;
    prev = new Domain(dom);

    // intersect with loop invariant
    {
      Domain domExit = dom;
      vector<LinComponent> ll;
      ll.push_back(LinComponent(1,x_i));
      assert(ll.size()==1);
      domExit.inequality(ll, -128, &lm, true);
    };
    vector<LinComponent> ll;
    ll.push_back(LinComponent(1,x_i));
    assert(ll.size()==1);
    Result res = dom.inequality(ll, -127, &lm, false);
    assert(res!=resUnsatisfiable);

    // increment i by one
    ll.clear();
    ll.push_back(LinComponent(-1,x_iTemp));
    ll.push_back(LinComponent(1,x_i));
    res = dom.inequality(ll, 1, &lm, true);
    assert(res==resChanged);
    dom.update(x_iTemp, x_i);
    dom.setVariableToTop(x_iTemp);
    dom.sane();

    // increment p by four
    ll.clear();
    ll.push_back(LinComponent(-1,x_oTemp));
    ll.push_back(LinComponent(1,x_o));
    res = dom.inequality(ll, 4, &lm, true);
    assert(res==resChanged);
    dom.update(x_oTemp, x_o);
    dom.setVariableToTop(x_oTemp);
    dom.sane();

    if (chat) cerr << "Domain after loop body:" << endl << dom;

  };

  signed long up;
  signed long low;
  Mult m;
  assert(prev);
  assert(prev->queryValue(&low, &up, &m, x_o)==7);

  bool correct = low==128 && up==512 && m==2 &&
    prev->queryPolyhedron(x_i, x_o)->getNoOfInequalities()==2;

  // Return 1 on error.
  return !correct;
};
