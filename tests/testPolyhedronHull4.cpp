/* This program tests, amongst others, that an inequality is dropped
   even if both end points lie outside the square.
 */

#include <iostream>
#include "interval.hh"
#include "planar.hh"
#include "polyhedron.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;
  Polyhedron<true> p1;
  Interval<true> xBound1;
  Interval<true> yBound1;
  xBound1.updateUpper(mpq_class(-1));

  vector<Inequality*> eqs1;
  eqs1.push_back(new Inequality(1,2,0));
  eqs1.push_back(new Inequality(1,-2,0));

  p1.addInequalitySet(xBound1, yBound1, eqs1);
  p1.enforceBounds(xBound1, yBound1);

  if (chat)
    cerr << "first polyhedron: x in " << xBound1 << ", y in " << yBound1
	 << endl << p1;

  Polyhedron<true> p2;
  Interval<true> xBound2;
  Interval<true> yBound2;


  vector<Inequality*> eqs2;
  eqs2.push_back(new Inequality(1,1,0));

  p2.addInequalitySet(xBound2, yBound2, eqs2);
  p2.enforceBounds(xBound2, yBound2);

  if (chat)
    cerr << "second polyhedron: x in " << xBound2 << ", y in " << yBound2
	 << endl << p2;

  Polyhedron<true> p3 = Polyhedron<true>(p1, xBound1, yBound1,
			     p2, xBound2, yBound2);

  Interval<true> xBound3 = Interval<true>(xBound1, xBound2);
  Interval<true> yBound3 = Interval<true>(yBound1, yBound2);

  if (chat)
    cerr << "hull polyhedron: x in " << xBound3 << ", y in " << yBound3
	 << endl << p3;

  bool p3p1 = p3.includes(xBound1, yBound1, p1);
  if (chat) cerr << "p3 includes p1: " << p3p1 << endl;
  bool p3p2 = p3.includes(xBound2, yBound2, p2);
  if (chat) cerr << "p3 includes p2: " << p3p2 << endl;
  bool p2p3 = p2.includes(xBound3, yBound3, p3);
  if (chat) cerr << "p2 includes p3: " << p2p3 << endl;
  if (chat) cerr << "p1 after enforceBounds: x in " << xBound1 << ", y in " 
		 << yBound1 << endl << p1 << endl;
  bool correct=p3p1 && p3p2 && p2p3 && p3.getNoOfInequalities()==1;
    
  // Return 1 on error.
  return !correct;
};
