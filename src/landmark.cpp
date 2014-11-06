// landmark.cpp: storing inequalities that are currently unsatisfiable but may
// enable a new branch in the program later iterations

#include "landmark.hh"
#include <iostream>
#include <assert.h>

#undef DEBUG_LANDMARKS

#ifdef __GNUC__
namespace stdext =__gnu_cxx;
#endif

using namespace std;

#ifdef NDEBUG
#undef DEBUG_LANDMARKS
#endif

namespace Tvpi {

  // Add a new landmark.
  void Landmarks::addLandmark(const vector<LinComponent>& ts,
			      mpz_class& val) {
    assert(sgn(val)>0);
    LmList::iterator iter = lmlist.begin();
    while (true) {
      if (iter==lmlist.end()) {
	lmlist.push_front(Landmark(ts,val));
#ifdef DEBUG_LANDMARKS
	cerr << "landmark: added new " << lmlist.front() << endl;
#endif // DEBUG_LANDMARKS
	return;
      };
      if (iter->terms==ts) {
	// Replace the current value if the new value is smaller (or
	// the current value is unset (zero)).
	if (sgn(iter->curVal)==0 || val<iter->curVal) {
	  iter->curVal=val;
#ifdef DEBUG_LANDMARKS
	  cerr << "landmark: augmented existing " << *iter << endl;
#endif // DEBUG_LANDMARKS
	};
	return;
      };
      iter++;
    }
  }

  // Calculate the number of iterations until any of the redundant
  // inequalities would become non-redundant.
  bool Landmarks::calcNoOfIterations(WideningInfo& wi) {
    bool changed=false;
    for (LmList::iterator iter = lmlist.begin();
	 iter!=lmlist.end(); iter++) {
      assert(sgn(iter->curVal)>=0);
      if (sgn(iter->prevVal)==0) {
	// The basic block that created this landmark was probably
	// executed only once. Hence, we need another iteration to
	// infer a delta.

	// Whenever the analyser calls this function without promoting
	// the landmarks at least once, then the current value is zero.
	if (sgn(iter->curVal)!=0) wi.hasPending=true;

	// However, the non-zero value in curVal might have been set
	// in the second iteration around the SCC if a specific branch
	// became active after the first iteration finished. 
	// Furthermore, the resulting state might be stable
	// immediately, such that no new landmarks will be added for
	// this statement. We need to change prevVal and curVal such
	// that we don't force another iteration (via some other path
	// in the SCC) indefinitely.
	iter->prevVal = iter->curVal;
      } else {
	mpz_class diff = iter->prevVal - iter->curVal;
#ifdef DEBUG_LANDMARKS
	cerr << "landmark: diff of terms " << iter->terms
	     << " is " << diff << "= " << iter->prevVal
	     << " - " << iter->curVal << endl;
#endif // DEBUG_LANDMARKS
	switch (sgn(diff)) {
	  // The difference should not be negative. This case is
	  // analogous to the one mentioned above, except that two
	  // inequalities where added in the last iteration of which the
	  // smaller one now is non-redundant. We better do another
	  // iteration to get a difference of the remaining inequality.
	case -1: {
	  wi.hasPending=true;
	  // Again, ensure that the difference is zero next time
	  // round. If it isn't and this basic block happens to be
	  // stable after this iteration, this landmark would never be
	  // updated again and we're in for an infinite loop.
	  iter->prevVal = iter->curVal;
	}; break;
	  // The difference is positive. Estimate how many more
	  // iterations are needed until the feasible space reaches
	  // the as-of-yet unsatisfiable inequality.
	case 1: {
	  mpz_fdiv_q(diff.get_mpz_t(), iter->curVal.get_mpz_t(),
		     diff.get_mpz_t());
	  // Always set iterations to the new count if iterations is
	  // has not been set so far. Also, replace the current number
	  // of iterations if the new estimate in diff is smaller.
	  if ((!wi.hasIterations) || (wi.iterations>diff)) {
	    mpz_swap(wi.iterations.get_mpz_t(), diff.get_mpz_t());
	    wi.hasIterations=true;
	    changed=true; 
	  }
	}; break;
	  // A difference of zero means that inequality will stay
	  // infeasible forever (or at least until some dead branch of
	  // the current SCC becomes live).
	default: break;
	}
      }
    }
    return changed;
  }

  ostream& operator<<(ostream& s, const Landmarks::Landmark& lm) {
    s << lm.terms << " cur=" << lm.curVal << ", prev=" << lm.prevVal;
    return s;
  }

  ostream& operator<<(ostream& s, const Landmarks& lm) {
    for (Landmarks::LmList::const_iterator iter = lm.lmlist.begin();
	 iter!=lm.lmlist.end(); iter++) {
      if (iter!=lm.lmlist.begin()) s << ", ";
      s << *iter;
    };
    return s;
  }
  
  // Calculate the number of iterations until any of the redundant
  // inequalities would become non-redundant.
  bool LandmarkTable::calcNoOfIterations(WideningInfo& wi) {
    bool changed = false;
    for (vector<Landmarks>::iterator iter = landmarks.begin();
	 iter != landmarks.end(); iter++)
      if (iter->calcNoOfIterations(wi)) {
	wi.statementNo = (iter-landmarks.begin());
	changed = true;
      };
    return changed;
  }

  ostream& operator<<(ostream& s, const LandmarkTable& lm) {
    s << "Landmark table:" << endl;
    for (vector<Landmarks>::const_iterator iter = lm.landmarks.begin();
	 iter != lm.landmarks.end(); iter++)
      s << (iter-lm.landmarks.begin()) << ": " << *iter << endl;
    return s;
  }

}

