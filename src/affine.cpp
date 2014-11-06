// affine.cpp - Reduced domain including affine and TVPI linear relations.
//

#include "affine.hh"
#include "interval.hh"
#include "tvpi.hh"
#include "polyhedron.hh"
#include "memory.hh"
#include <algorithm>

using namespace std;

#undef DEBUG_INTERSECT
#undef DEBUG_LINOPT
#undef DEBUG_SUBSTITUTE
#undef DEBUG_DIVIDEMULT
#undef DEBUG_LANDMARKS
#undef DEBUG_DEMOTE
#undef DEBUG_ENTAIL
#undef DEBUG_ENTAIL_FAIL
#undef DEBUG_JOIN
#undef DEBUG_PROMOTE

#undef JOIN_CHECK_ENTAILMENT

#define HAVE_EQUALITYSUBST


namespace Tvpi {

  // FIXME: this must be a template
  inline void
  swap(LinComponent& l1, LinComponent& l2) {
    l1.swap(l2);
  }

  DomVar Domain::lastTempVar = -1;
  Domain::VarTable Domain::typedVars;
  Domain::TypeTable Domain::registeredTypes;
  const DomVar Domain::invalidVar = 0;

  // Dump the Domain dude.
  std::ostream& operator<<(std::ostream& s, const Domain& d) {
    d.sane();
    vector<Domain::RangeInfo>::const_iterator iter = d.ranges.begin();
    while (iter!=d.ranges.end()) {
      s << "r" << (iter-d.ranges.begin());
      if (iter->xRef<0)
	s << " t_" << (-iter->xRef); else
	s << " x_" << iter->xRef;
      s << " = " << iter->interval;
      Domain::ConsTable::const_iterator viIter = d.cons.find(iter->xRef);
      assert(viIter!=d.cons.end());
      assert(viIter->second.isRange);
      assert(viIter->second.var==(size_t) (iter-d.ranges.begin()));
      if (viIter->second.mult==invalidMult) s << ", mult n/a";
	else s << ", mult " << (int) viIter->second.mult;
      if (iter->xRef>0)
	s << ", nom range " << Domain::getVariableTop(iter->xRef)
	  << ", type " << Domain::getVariableType(iter->xRef);
      s << endl;
      
      iter++;
    };
    d.relDomain.showDist(s);
    {
      Mult mults[d.relDomain.size()];
      for(size_t i=0; i<d.relDomain.size(); i++) {
	Domain::ConsTable::const_iterator viIter =
	  d.cons.find(d.relDomain.xRef(i));
	assert(viIter!=d.cons.end());
	assert(!viIter->second.isRange);
	mults[i]=viIter->second.mult;
      };
      d.relDomain.output(s, mults);
    };
    return s;
  };

  // Retrieve the multiplicity and/or the value of a variable.
  int Domain::queryValue(signed long* low, signed long* upp,
			 Mult* m, DomVar v) const {
    if (v>=(DomVar) typedVars.size()) v=invalidVar;
    int res=0;
    ConsTable::const_iterator iter = cons.find(v);
    Interval<isIntegral> i;
    if (iter==cons.end()) {
      // The variable is not tracked and hence takes on the range of
      // its type.
      if (m) *m=0;
      i = getVariableTop(v);
    } else {
      res|= 4;
      if (m) {
	*m=iter->second.mult;
      };
      i = (iter->second.isRange ? ranges[iter->second.var].interval :
	   relDomain.getInterval(iter->second.var));
      i = i.stretch(iter->second.mult);
    };
    if (!i.lowerIsInfinite()) {
      mpz_class l = i.getLower();
      if (low) *low = mpz_get_si(l.get_mpz_t());
      if (mpz_fits_slong_p(l.get_mpz_t())) res|= 1;
    };
    if (!i.upperIsInfinite()) {
      mpz_class u = i.getUpper();
      if (upp) {
	if (mpz_fits_slong_p(u.get_mpz_t())) {
	  res|= 2;
	  *upp = mpz_get_si(u.get_mpz_t());
	} else if (mpz_fits_ulong_p(u.get_mpz_t())) {
	  res|= 8;
	  *upp = (signed long) mpz_get_ui(u.get_mpz_t());
	}
      }
    }
    return res;
  }

  // Set the multiplicity of a variable. Returns false if the domain
  // became unsatisfiable.
  bool Domain::enforceMult(DomVar var, Mult m) {
    checkDomVar(var,true);
    if (m==0) return true;
    ConsTable::iterator iter = cons.find(var);
    if (iter==cons.end()) iter = createRange(var, getVariableTop(var));
    // For the sake of efficiency, don't multiply by 2^255.
    if (m==invalidMult) {
      if (!iter->second.isRange) {
        // The variable is in the relational domain. Remove it and insert it
        // as range.
        redirectToRange(iter->second.var);
        iter = cons.find(var);
        if (iter==cons.end()) {
          iter = createRange(var, getVariableTop(var));
        }
      };
      assert(iter->second.isRange);
      // Dividing the interval by 2^255 means it's only going to be 
      // satisfiable if it contains zero. If it does, it will be zero after
      // the division. Implement this logic by restricting the range to zero.
      ranges[iter->second.var].interval.updateUpper(0);
      ranges[iter->second.var].interval.updateLower(0);
      iter->second.mult=invalidMult;
      bool empty = ranges[iter->second.var].interval.isEmpty();
      if (empty) return false;
      // Remove the range if the type of the variable happens to have the 
      // default range of zero.
      if (ranges[iter->second.var].interval==getVariableTop(var)) {
        removeRange(iter->second.var);
        cons.erase(iter);
      };
      assert(sane());
      return true;
    };
    if (iter->second.mult<m) {
      Mult diff = m-iter->second.mult;
      iter->second.mult=m;
      if (iter->second.isRange) {
	return ranges[iter->second.var].interval.enforceMult(diff);
      } else {
	vector<Inequality*> set;
	mpz_class f = 0;
	mpz_setbit(f.get_mpz_t(), (unsigned long int) diff);
	set.reserve(2);
	set.push_back(new Inequality(-1, f, 0));
	set.push_back(new Inequality(1, -f, 0));
	TvpiVar temp = relDomain.createVariable();
	Result res = relDomain.addInequalitySet(iter->second.var,
						temp, set, true);
	relDomain.update(iter->second.var);
	relDomain.setXRef(iter->second.var, var);
	iter->second.mult=m;	
        assert(sane());
      	return (res!=resUnsatisfiable);
      }
    }
    return true;
  }


  // Project out a variable. The variable will be unbounded after
  // this call.
  void Domain::projectOut(DomVar var) {
    checkDomVar(var,);
    ConsTable::iterator iter = cons.find(var);
    if (iter==cons.end()) {
      // If a typed variable is not in this domain, it implicitly has
      // the bounds of its type. In order to ensure that the variable
      // is truly unbounded, insert it as unbounded interval.
      if (var>invalidVar) createRange(var, Interval<isIntegral>());
      assert(sane());
      return;
    };
    iter->second.mult=0;
    if (iter->second.isRange) {
      if (var>invalidVar)
	ranges[iter->second.var].interval=Interval<isIntegral>();
      else {
	removeRange(iter->second.var);
	cons.erase(iter);
      }
    } else {
      removeTvpi(iter->second.var);
      cons.erase(iter);
      if (var>invalidVar) createRange(var, Interval<isIntegral>());
    }
  }

  // Remove a range.
  void Domain::removeRange(TvpiVar var) {
    assert(var<ranges.size());
    if (var<ranges.size()-1) {
      // Replace this dying variable with the last one.
      ConsTable::iterator iter = cons.find(ranges.back().xRef);
      assert(iter!=cons.end());
      assert(iter->second.isRange);
      iter->second.var = var;
      ranges[var].interval.swap(ranges.back().interval);
      ranges[var].xRef=ranges.back().xRef;
    };
    ranges.pop_back();
  }

  // Remove a tvpi variable.
  void Domain::removeTvpi(TvpiVar var) {
    assert(var<relDomain.size());
    if (var==relDomain.size()-1) relDomain.removeLast(); else {
      ConsTable::iterator iter = cons.find(relDomain.xRef(relDomain.size()-1));
      assert(iter!=cons.end());
      assert(!iter->second.isRange);
      iter->second.var = var;
      relDomain.update(var);
    }
  }

  // Promote a variable into the relational domain.
  TvpiVar Domain::promoteRange(ConsTable::iterator& iter) {
    assert(iter->second.isRange);
    TvpiVar res = relDomain.createVariable(ranges[iter->second.var].interval,
					   iter->first);
    removeRange(iter->second.var);
    iter->second.isRange=false;
    iter->second.var=res;
    return res;
  }

  // Change an entry in the cons table from refering to a relational
  // variable to refering to a range. If the range is unbounded, v is
  // deleted completely.
  void Domain::redirectToRange(TvpiVar v) {
    DomVar dVar = relDomain.xRef(v);
    ConsTable::iterator iter = cons.find(dVar);
    assert(iter!=cons.end());
    // Remove the variable from the hashtable if the variable takes on the
    // value of it's type.
    if (relDomain.getInterval(v) == getVariableTop(dVar)) {
      cons.erase(iter);
      return;
    };
    iter->second.isRange=true;
    iter->second.var=ranges.size();
    ranges.push_back(RangeInfo(relDomain.getInterval(v), dVar));
  };

  // Demote variables from the relational domain.
  void Domain::demote() {
    // Remove all temporary variables.
    {
      DomVar kill[cons.size()];
      size_t lastKill = 0;
      for (ConsTable::iterator iter = cons.begin(); iter!=cons.end(); iter++)
	if (iter->first<invalidVar) kill[lastKill++]=iter->first;
	else if (iter->second.isRange) {
	  // Remove variables that can take on exactly the values in
	  // of their type.
	  if (ranges[iter->second.var].interval==getVariableTop(iter->first) &&
	      (iter->second.mult==0 || iter->second.mult==invalidMult))
	    kill[lastKill++]=iter->first;
	};
      for (size_t i=0; i<lastKill; i++) {
#ifdef DEBUG_DEMOTE
	if (i==0) cerr << "Demote: removing temporary vars: "; else
	  cerr << ", ";
	showVar(cerr, kill[i]);
#endif // DEBUG_DEMOTE
	setVariableToTop(kill[i]);
      }
    };

    size_t size = relDomain.size();
    bool demote[size];

    for (TvpiVar b=0; b<size; b++) demote[b]=true;

    for (TvpiVar y=0; y<size; y++) {
      if (!demote[y]) break;
      for (TvpiVar x=0; x<y; x++)
	// Reset the demote flag of x and y if there are inequalities
	// in the (x,y) projection. However, if the range of either x
	// and y are singleton ranges, then the inequalities we see
	// are all redundant, in which case we need to demote.
	if (relDomain.getProjection(x,y).getNoOfInequalities()>0 &&
	    !relDomain.getInterval(x).isSingleton() &&
	    !relDomain.getInterval(y).isSingleton()) {
	  demote[x]=false;
	  demote[y]=false;
	}
    };

#ifdef DEBUG_DEMOTE
    {
      bool isFirst = true;
      for (TvpiVar v=0; v<size; v++) if (demote[v]) {
	if (isFirst) {
	  cerr << "Demoting this=" << hex << this << dec << " ";
	  isFirst=false;
	} else cerr << ", ";
	cerr << "p_" << v;
      };
      if (!isFirst) cerr << endl;
    };
#endif // DEBUG_DEMOTE

    // Remove all variables that have their flag set.
    TvpiVar v=0;
    while (v<size) {
      // Ensure that the last variable in the TVPI system is not
      // supposed to be demoted.
      do {
	TvpiVar lastVar=size-1;
	if (!demote[lastVar]) break;
#ifdef DEBUG_DEMOTE
	cerr << "demoting last TVPI var p_" << lastVar
	     << ", that is x_" << relDomain.xRef(lastVar) << endl;
#endif // DEBUG_DEMOTE
	redirectToRange(lastVar);
	relDomain.removeLast();
	size--;
	assert(size==relDomain.size());
      } while (v<size);

      if (v>=size) break;

      // Demote a variable that is not the last variable of the
      // relational domain.
      if (demote[v]) {
	// At this point, the last variable is one that remains in the
	// relational domain.
	assert(v+1<size);
	assert(!demote[size-1]);
	// Demote v.
#ifdef DEBUG_DEMOTE
	cerr << "demoting inner TVPI var p_" << v << endl;
#endif // DEBUG_DEMOTE
	redirectToRange(v);
	removeTvpi(v);
	size--;
	assert(size==relDomain.size());
      };
      v++;
    };
    assert(sane());
  }

  // Let the given target variable contain all the values it had before and
  // all values of the source variable.
  void Domain::augment(DomVar source, DomVar target) {
    checkDomVar(source,);
    checkDomVar(target,);
    ConsTable::iterator sIter = cons.find(source);
    if (sIter==cons.end()) {
      Interval<isIntegral> st = getVariableTop(source);
      ConsTable::iterator tIter = cons.find(target);
      if (tIter!=cons.end())
	if (tIter->second.isRange)
	  st.join(ranges[tIter->second.var].interval);
	else
	  st.join(relDomain.getInterval(tIter->second.var));
      else
	st.join(getVariableTop(target));
      setVariable(target,0,st);
      return;
    };
    ConsTable::iterator tIter = cons.find(target);
    if (tIter==cons.end()) return;
    // This function cannot deal with multiplicity correctly.
    assert(tIter->second.mult==sIter->second.mult);
    // Force both variables into the relational domain and update one
    // with the other.
    if (sIter->second.isRange) promoteRange(sIter);
    if (tIter->second.isRange) promoteRange(tIter);
    relDomain.augment(sIter->second.var, tIter->second.var);
    assert(sane());
  }

  // Replace the given target variable with the valuation of another
  // variable.
  void Domain::update(DomVar source, DomVar target) {
    checkDomVar(source,);
    checkDomVar(target,);
    // Force both variables into the relational domain and update one
    // with the other. The result is an equality relationship between
    // the two variables.
    ConsTable::iterator sIter = cons.find(source);
    if (sIter==cons.end())
      sIter = createRange(source, Interval<isIntegral>());

    ConsTable::iterator tIter = cons.find(target);
    if (tIter==cons.end()) {
      tIter = createRange(target, Interval<isIntegral>());
      // Insertion might have invalidated the sIter iterator.
      sIter = cons.find(source);
      assert(sIter!=cons.end());
    };
    if (sIter->second.isRange) promoteRange(sIter);
    if (tIter->second.isRange) promoteRange(tIter);
    relDomain.update(sIter->second.var, tIter->second.var);
    tIter->second.mult=sIter->second.mult;
    assert(sane());
  }

  // Check if this domain describes a sub-space of the prev
  // domain. The result is true if this domain entails the prev
  // domain (i.e. if the prev domain includes this domain)..
  bool Domain::entails(Domain& prev) {
    // The relation between indices of the current and previous domain.
    TvpiVar prevToCur[prev.relDomain.size()];
    for(size_t p=0; p<prev.relDomain.size(); p++) prevToCur[p]=invalidTvpiVar;

    for(ConsTable::const_iterator viPrevIter = prev.cons.begin();
	viPrevIter!=prev.cons.end(); viPrevIter++) {
      ConsTable::iterator viThisIter = cons.find(viPrevIter->first);

      if (viThisIter==cons.end()) {
	// The variable is not in this domain. Perform special tests
	// that compare the variable with the range implied by the
	// type of the variable.
	if (viPrevIter->second.isRange) {
	  if (!prev.ranges[viPrevIter->second.var].interval.
	      stretch(viPrevIter->second.mult).
	      includes(getVariableTop(viPrevIter->first))) return false;
	} else {
	  // !viPrevIter->second.isRange

	  // Insert the top range of the variable into the relational
	  // domain, and then use the relational entailment test.
	  viThisIter = createTvpiRange(viPrevIter->first,
				       getVariableTop(viPrevIter->first));
	  prevToCur[viPrevIter->second.var] = viThisIter->second.var;
	}
	continue;
      };

      // Check that multiplicity is matching. Invalid multiplicities
      // (i.e. the multiplicity of zero) are a bit tricky, as it is
      // matching any other multiplicity. So deal with invalidMult
      // specially.
      Mult thisMult = viThisIter->second.mult;
      Mult prevMult = viPrevIter->second.mult;
      if (thisMult==invalidMult)
	if (prevMult==invalidMult) thisMult=prevMult=0;
	else thisMult = prevMult;
      else  if (prevMult==invalidMult) prevMult = thisMult;
      // The previous iterate had a greater multiplicity, thus contains fewer points
      // than this domain.
      if (thisMult<prevMult) {
#ifdef DEBUG_ENTAIL_FAIL
        cerr << "entail(x_" << viPrevIter->first << ") mult: cur: 2^"
	     << (int) thisMult << " |=  prev: 2^"
	     << (int) prevMult << endl;
#endif // DEBUG_ENTAIL_FAIL
	return false;
      };
      if (thisMult>prevMult) {
	// Weaken the multiplicity of this domain until it
	// reaches that of the prev.
	Mult diff = thisMult - prevMult;
	if (viThisIter->second.isRange)
	  ranges[viThisIter->second.var].interval=
	    ranges[viThisIter->second.var].interval.stretch(diff);
	else
	  relDomain.stretch(viThisIter->second.var, diff);
	viThisIter->second.mult=prevMult;
      };

      if (viThisIter->second.isRange) {
	if (viPrevIter->second.isRange) {
	  if (!prev.ranges[viPrevIter->second.var].interval.
	      includes(ranges[viThisIter->second.var].interval)) {
#ifdef DEBUG_ENTAIL_FAIL
	    cerr << "entail(x_" << viPrevIter->first << ") ranges: cur: "
		 << ranges[viThisIter->second.var].interval << " |=  prev: "
		 << prev.ranges[viPrevIter->second.var].interval << endl;
#endif // DEBUG_ENTAIL_FAIL
	    return false;
	  }
	} else {
	  // !viPrevIter->second.isRange

	  // The variable in this domain is a range while it is in a
	  // polyhedron in the previous domain. Insert this range into
	  // the relational part of this domain, thereby pushing the
	  // entailment check into the relational domain.
	  prevToCur[viPrevIter->second.var] = promoteRange(viThisIter);
	}
      } else {
	// !viThisIter->second.isRange
	if (viPrevIter->second.isRange) {
	  // The previous variable is only a bound and the current is in the
	  // relational domain. A simple test determines entailment.
	  if (!prev.ranges[viPrevIter->second.var].interval.
	      includes(relDomain.getInterval(viThisIter->second.var))) {
#ifdef DEBUG_ENTAIL_FAIL
	    cerr << "entail(x_" << viPrevIter->first << ") poly/range: cur: "
		 << relDomain.getInterval(viThisIter->second.var)
		 << " |=  prev: "
		 << prev.ranges[viPrevIter->second.var].interval << endl;
#endif // DEBUG_ENTAIL_FAIL
	    return false;
	  }
	} else {
	  // !viPrevIter->second.isRange
	  prevToCur[viPrevIter->second.var] = viThisIter->second.var;
	}
      }
    };
    assert(sane());
    return prev.relDomain.includes(relDomain, prevToCur);
  }

  // Calculate the join of both and store the result in this domain. 
  // Note the weird argument order: The previous domain that is
  // already stored at a program location remains untouched. It it the
  // new, changed state that is joined with the previous state.
  void Domain::joinWiden(Domain& prev, mpz_class extrapolate) {

#ifdef DEBUG_PROMOTE
    cerr << "join: *this:\n" << *this << "prev:\n" << prev;
#endif // DEBUG_PROMOTE

    // Remove temporary variable in this domain if they are absent
    // from the prev domain. For typed variables, insert the range of
    // their type into the prev domain (which leaves prev semantically
    // unchanged).
    {
      TvpiVar kill[cons.size()];
      size_t lastKill = 0;
      for(ConsTable::iterator thisIter = cons.begin();
	  thisIter!=cons.end(); thisIter++) {
	ConsTable::iterator prevIter = prev.cons.find(thisIter->first);
	if (prevIter==prev.cons.end()) {
	  if (thisIter->second.isRange &&
	      ranges[thisIter->second.var].interval==
	      getVariableTop(thisIter->first)) {
	    // Remove variables in this domain that shouldn't be here
	    // in the first place.
	    kill[lastKill++]=thisIter->first;
	    removeRange(thisIter->second.var);
	  } else if (thisIter->first<invalidVar) {
	    // Kill any temporary variable that is not in the
	    // previous domain.
	    kill[lastKill++]=thisIter->first;
	    if (thisIter->second.isRange) removeRange(thisIter->second.var);
	    else removeTvpi(thisIter->second.var);
	  } else {
	    // For normal variables, insert their maximum range as
	    // determined by their type into prev. Create a range in
	    // the relational domain as the variable is bound to get
	    // involved in a linear relationship. The only case in
	    // which no linear relationship would ensue is when the
	    // variable in this domain had exactly the range of the
	    // variable's type. These, however, were removed in the
	    // first test of this loop.
#ifdef DEBUG_PROMOTE
	    cerr << "var "; showVar(cerr, thisIter->first);
	    cerr << " not in prev, add with type "
		 << typedVars[thisIter->first] << ", interval "
		 << registeredTypes[typedVars[thisIter->first]] <<endl;
#endif // DEBUG_PROMOTE
	    prev.createTvpiRange(thisIter->first,
				 getVariableTop(thisIter->first));
	  }
	}
      };
      for (size_t i=0; i<lastKill; i++) {
#ifdef DEBUG_PROMOTE
	cerr << "temp var "; showVar(cerr, kill[i]);
	cerr << " not in prev, killed." << endl;
#endif // DEBUG_PROMOTE
	cons.erase(kill[i]);
      };
    };

    assert(cons.size()<=prev.cons.size());
    assert(sane());

    // Add top values into this domain for all typed variables that are in
    // the prev domain, but not in this.
    for (ConsTable::iterator prevIter = prev.cons.begin();
	 prevIter != prev.cons.end(); prevIter++) 
      if (prevIter->first>invalidVar &&
	  cons.find(prevIter->first)==cons.end()) {
#ifdef DEBUG_PROMOTE
	cerr << "var "; showVar(cerr, prevIter->first);
	cerr << " not in cur, add with type "
	     << typedVars[prevIter->first] << ", interval "
	     << registeredTypes[typedVars[prevIter->first]] <<endl;
#endif // DEBUG_PROMOTE
	createTvpiRange(prevIter->first,
			getVariableTop(prevIter->first));
    };

    // If it weren't for temporary variables in prev, this would be an
    // equality.
    assert(cons.size()<=prev.cons.size());

    // From this point onwards, no new variables are added.
    TvpiVar curToPrev[relDomain.size()+ranges.size()];
    Mult multPrev[prev.relDomain.size()+prev.ranges.size()];

#ifndef NDEBUG
    // Make DenseTvpi::join() fail if prev.cons does not map to all
    // variables in prev.relDomain.
    for(size_t t=0; t<relDomain.size()+ranges.size(); t++)
      curToPrev[t]=invalidTvpiVar;

    // Make DenseTvpi::join() fail if not all multiplicities are set.
    for(size_t p=0; p<prev.relDomain.size()+prev.ranges.size(); p++)
      multPrev[p]=invalidMult;
#endif

    // Promote all variables that appear as range in one and as
    // TVPI variables in the other domain.
    for(ConsTable::iterator viThisIter = cons.begin();
	viThisIter!=cons.end(); viThisIter++) {
      ConsTable::iterator viPrevIter = prev.cons.find(viThisIter->first);
      assert(viPrevIter!=prev.cons.end());
#ifdef DEBUG_PROMOTE
      cerr << (viThisIter->second.isRange ? "r_" : "p_")
	   << viThisIter->second.var << ", prev "
	   << (viPrevIter->second.isRange ? "r_" : "p_")
	   << viPrevIter->second.var << " ";
#endif // DEBUG_PROMOTE

      // Calculate sane values for multiplicity. See entails().
      Mult thisMult = viThisIter->second.mult;
      Mult prevMult = viPrevIter->second.mult;
      if (thisMult==invalidMult)
	if (prevMult==invalidMult) thisMult=prevMult=0;
	else thisMult = prevMult;
      else if (prevMult==invalidMult) prevMult = thisMult;

      // No adjustment of multiplicity necessary when...
      Mult diff=0;

      if (thisMult>prevMult) {
	// ...we can weaken the multiplicity of this domain until it
	// reaches that of the prev.
	Mult diff = thisMult - prevMult;
	if (viThisIter->second.isRange)
	  ranges[viThisIter->second.var].interval=
	    ranges[viThisIter->second.var].interval.stretch(diff);
	else
	  relDomain.stretch(viThisIter->second.var, diff);
      } else if (thisMult<prevMult) {
	// We need to scale the value in the previous domain if the variable
	// has a larger multiplicity. A task we pass on to DenseTvpi::join().
	diff = prevMult - thisMult;
      };

      // The joined variable will have the minium of both
      // multiplicities, except if the variable is zero in both
      // domains. This case is dealt with below.
      viThisIter->second.mult = min(prevMult, thisMult);

      if (viThisIter->second.isRange) {
	if (viPrevIter->second.isRange) {
	  if (ranges[viThisIter->second.var].interval==
	      prev.ranges[viPrevIter->second.var].interval.stretch(diff)) {
	    // If both ranges are the same, no relational information
	    // will be created. Reset the multiplicity in case the
	    // interval is zero.
	    if (ranges[viThisIter->second.var].interval.isZero()) {
	      viThisIter->second.mult = invalidMult;
	    }; 
#ifdef DEBUG_PROMOTE
	    cerr << "keeping both ranges" << endl;
#endif // DEBUG_PROMOTE
	  } else {
	    // Promote the two ranges.
	    TvpiVar domVarThis = promoteRange(viThisIter);
	    TvpiVar domVarPrev = prev.promoteRange(viPrevIter);
	    multPrev[domVarPrev]=diff;
#ifdef DEBUG_PROMOTE
	    cerr << "promoting both ranges to (cur) p_" << domVarThis
		 << " and (prev) p_" << domVarPrev << endl;
#endif // DEBUG_PROMOTE
	    curToPrev[domVarThis]=domVarPrev;
	  }
	} else { // !viPrevIter->second.isRange
	  TvpiVar domVarThis = promoteRange(viThisIter);
#ifdef DEBUG_PROMOTE
	  cerr << "promoting this range to p_" << domVarThis << endl;
#endif // DEBUG_PROMOTE
	  multPrev[viPrevIter->second.var]=diff;
	  curToPrev[domVarThis]=viPrevIter->second.var;
	}
      } else { // !viThisIter->second.isRange
	if (viPrevIter->second.isRange) {
	  TvpiVar domVarPrev = prev.promoteRange(viPrevIter);
#ifdef DEBUG_PROMOTE
	  cerr << "promoting prev range to p_" << domVarPrev << endl;
#endif // DEBUG_PROMOTE
	  multPrev[domVarPrev]=diff;
	  curToPrev[viThisIter->second.var]=domVarPrev;
	} else { // !viPrevIter->second.isRange
#ifdef DEBUG_PROMOTE
	  cerr << "both already in tvpi" << endl;
#endif // DEBUG_PROMOTE
	  multPrev[viPrevIter->second.var]=diff;
	  curToPrev[viThisIter->second.var]=viPrevIter->second.var;
	}
      }
    };

#ifdef DEBUG_JOIN
    cerr << "var mapping cur->prev: ";
    for (size_t i=0; i<relDomain.size(); i++) {
      if (i) cerr << ", ";
      cerr << "p_" << i << " -> p_" <<curToPrev[i];
    };
    cerr << endl << "before join: ";
    cerr << "this->relDomain:\n" << this->relDomain
	 << "prev.relDomain:\n" << prev.relDomain;
#endif // DEBUG_JOIN

#ifdef JOIN_CHECK_ENTAILMENT
    DenseTvpi<isIntegral> current(relDomain);
#endif

    // Only join the relational domains since variables that are still
    // stored as ranges have the same range (otherwise they would have
    // been promoted).
    relDomain.join(prev.relDomain, curToPrev, multPrev);

#ifdef JOIN_CHECK_ENTAILMENT
    bool correct = prev.entails(*this);
    if (!correct) {
      cerr << "hull: " << relDomain << endl;
      cerr << "current: " << current << endl;
      cerr << "prev: " << prev.relDomain << endl;
      cerr << "var mapping cur->prev: ";
      for (size_t i=0; i<relDomain.size(); i++) {
	if (i) cerr << ", ";
	cerr << "p_" << i << " -> p_" << curToPrev[i];
      };
      cerr << endl;
    };
    assert(correct);
#endif

#ifdef DEBUG_JOIN
    cerr << "after join: ";
    cerr << "this->relDomain:\n" << this->relDomain;
#endif // DEBUG_JOIN

    // Widen if asked to.
    if (sgn(extrapolate)>=0) 
      relDomain.widen(prev.relDomain, curToPrev, extrapolate);
    
    assert(sane());
  }

  // Calculate the maximum of the given expression and return it
  // in value. Returns false if the expression is unbounded.
  bool Domain::linOpt(vector<LinComponent>& comps, mpz_class& value) {
    // Calculate the maximum of all non-relational terms.
    mpz_class constant = 0;

    vector<TVPIComponent> newComps;
    newComps.reserve(comps.size());
    for (vector<LinComponent>::iterator iter = comps.begin();
	 iter!=comps.end(); iter++) {
      assert(sgn(iter->coefficient)!=0);
      ConsTable::iterator vIter = cons.find(iter->variable);
      if (vIter==cons.end()) {
	const Interval<isIntegral> &i=getVariableTop(iter->variable);
        if (sgn(iter->coefficient)>0) {
          if (i.upperIsInfinite()) return false;
	  constant+=i.getUpper()*iter->coefficient;
	} else {
          assert(sgn(iter->coefficient)<0);
	  if (i.lowerIsInfinite()) return false;
	  constant+=i.getLower()*iter->coefficient;
	}
      } else if (vIter->second.isRange) {
	const Interval<isIntegral> &i=ranges[vIter->second.var].interval;
	mpz_class bound;
        if (sgn(iter->coefficient)>0) {
          if (i.upperIsInfinite()) return false;
	  bound=i.getUpper()*iter->coefficient;
	} else {
          assert(sgn(iter->coefficient)<0);
	  if (i.lowerIsInfinite()) return false;
	  bound=i.getLower()*iter->coefficient;
	}
	if (vIter->second.mult!=invalidMult)
	  mpz_mul_2exp(bound.get_mpz_t(), bound.get_mpz_t(),
		       vIter->second.mult);
	constant+=bound;
      } else {
	newComps.push_back(TVPIComponent(iter->coefficient,
					 vIter->second.var,
					 vIter->second.mult));
      }
    };
    // Query the maximum of the relational part.
    sort(newComps.begin(), newComps.end());
    bool isFinite = relDomain.linOpt(newComps, value);
    value+=constant;
    return isFinite;
  }


  // Set a new variable to a specific value.
  bool Domain::setVariable(DomVar var, Mult m, Interval<isIntegral> c) {
    checkDomVar(var,true);

    Mult m_ = calcMult(c);
    if (m_>m) m=m_;

    bool isSat = c.enforceMult(m);
    if (!isSat) return false;

    if (c==getVariableTop(var)) {
	setVariableToTop(var);
	return true;
    };
    ConsTable::iterator iter = cons.find(var);
    if (iter==cons.end()) {
      TvpiVar idx = ranges.size();
      ranges.push_back(RangeInfo(c, var));
      cons.insert(ConsTable::value_type(var,VarInfo(idx, true, m)));
      assert(sane());
      return true;
    };
    VarInfo& vi=iter->second;
    vi.mult=m;
    if (vi.isRange) {
      ranges[vi.var].interval=c;
      return true;
    };
    removeTvpi(vi.var);
    vi.isRange=true;
    vi.var=ranges.size();
    ranges.push_back(RangeInfo(c, var));
    assert(sane());
    return true;
  }

  // Intersect with inequality a_1 x_1 + ... + a_n x_n + c <= 0. If
  // isEquality is true, replace "<=" with "=".
  Result Domain::inequality(vector<LinComponent>& comps, mpz_class c,
			    LandmarkTable* lm,
			    bool isEquality) {
    // Join duplicate variables and remove those with zero
    // coefficient.
    canonicalizeTerms(comps);

    // No variables means no landmarks and no intersection.
    if (comps.size()==0) 
      return (isEquality ?
	      (sgn(c)==0 ? resRedundant : resUnsatisfiable) :
	      (sgn(c)<=0 ? resRedundant : resUnsatisfiable));

#ifdef DEBUG_INTERSECT
    cerr << "Intersecting with ";
    showTerms(cerr, comps, c);
    if (isEquality) cerr << "=0"; else cerr << "<=0";
    cerr << endl;
#endif // DEBUG_INTERSECT

    // Check if there are temporary variables in the expression. Since the
    // varaibles are sorted, any temporary variable will be in the first slot.
    bool hasTemps = comps[0].variable<invalidVar;

    if (lm && (!hasTemps)) {
      mpz_class extreme;

      if (isEquality) {
	// Check whether the expression with >= 0 is unsatisfiable.
	bool isFinite = linOpt(comps, extreme);
#ifdef DEBUG_LANDMARKS
	cerr << "maximum of terms ";
        showTerms(cerr, comps, c);
	if (isFinite) cerr << " is " << extreme;
	else cerr << "does not exist";
#endif // DEBUG_LANDMARKS

	extreme=-extreme;
	extreme=extreme-c;

#ifdef DEBUG_LANDMARKS
	if (isFinite) cerr << " distance is " << extreme << endl;
	else cerr << endl;
#endif // DEBUG_LANDMARKS

	if (isFinite && sgn(extreme)>0) {
	  // It is, insert this inequality as landmark. 
	  lm->addLandmark(comps, extreme);
          return resUnsatisfiable;
	}
      };

      for ( vector<LinComponent>::iterator iter = comps.begin();
	    iter!=comps.end(); iter++) iter->coefficient=-iter->coefficient;

      // Check whether the inequality is unsatisfiable.
      bool isFinite = linOpt(comps, extreme);

#ifdef DEBUG_LANDMARKS
      cerr << "minimum of terms ";
      showTerms(cerr, comps, c);
      if (isFinite) cerr << " is " << (-extreme);
      else cerr << "does not exist";
#endif // DEBUG_LANDMARKS

      extreme=c-extreme;

#ifdef DEBUG_LANDMARKS
      if (isFinite) cerr << " distance is " << extreme << endl;
      else cerr << endl;
#endif // DEBUG_LANDMARKS

      if (isFinite && sgn(extreme)>0) {
	lm->addLandmark(comps, extreme);
	return resUnsatisfiable;
      }

      for ( vector<LinComponent>::iterator iter = comps.begin();
	    iter!=comps.end(); iter++) iter->coefficient=-iter->coefficient;

    };

    // Multiply each coefficient of this inequality by the
    // multiplicity of the variable in the domain. After this call,
    // all coefficients are cast in terms of the underlying domain
    // where variables are tracked at 1/2^m of their normal range.
    multiplyByMultiplicity(comps, c, isEquality);

    // Refine multiplicity information.
    if (isEquality) {
      Mult newMult[comps.size()];
      Result res = calcNewMult(comps, c, newMult);

      // If the terms of the equality imply a higher multiplicity
      // than the constant, the whole domain becomes unsatisfiable. 
      // No landmarks will be gathered for this case. This is ok,
      // since unless the multiplicity of the expression changes, no
      // matter how much we extrapolate, the multiplicity
      // information will always mean that this constraint is
      // unsatisfiable.
      if (res==resUnsatisfiable) {
#ifdef DEBUG_INTERSECT
	cerr << "domain unsatisfiable due to multiplicity of constant" << endl;
#endif // DEBUG_INTERSECT
	return resUnsatisfiable;
      };
      if (res==resChanged) {
#ifdef DEBUG_INTERSECT
	cerr << "Force multiplicity of variables to be: ";
	for (size_t i=0; i<comps.size(); i++) {
	  if (i>0) cerr << ", ";
	  cerr << (int) newMult[i];
	};
	cerr << endl;
#endif // DEBUG_INTERSECT
	// The multiplicity information of some of the variables can
	// be changed.
	/* makeUnique(); -- if we ever reference-count the domain */
	for (size_t i=0; i<comps.size(); i++) {
	  bool isSat = enforceMult(comps[i].variable, newMult[i]);
	  if (!isSat) return resUnsatisfiable;
	}
      }
    };

    // Turn domain variables into TVPI variables. This function may
    // promote ranges to the relational domain, but does not change
    // the feasible state space in any way. Hence, it is ok to call
    // this function even if this domain is shared (which isn't
    // implemented yet).
    bool canApprox = renameVariables(comps, c);

    // Bail out here if renameVariables has inlined all variables.
    if (comps.size()==0) 
      return (isEquality ?
	      (sgn(c)==0 ? resRedundant : resUnsatisfiable) :
	      (sgn(c)<=0 ? resRedundant : resUnsatisfiable));

    // If more than two variables are unbounded, the given constraint
    // is non-redundant, but cannot be added to the relational domain
    // since it cannot be expressed as a set of TVPI inequalities.
    if (!canApprox) {
#ifdef DEBUG_INTERSECT
      cerr << "No TVPI approximation possible!" << endl;
#endif // DEBUG_INTERSECT
      return resChanged;
    };

    sort(comps.begin(), comps.end());

    Result res = relDomain.approximateInequality(comps, c, isEquality);

    // For each variable in the relational domain, check if its range is
    // zero and force multiplicity to invalidMult if that is the case.
    for (TvpiVar i = 0; i < relDomain.size(); i++)
      if (relDomain.getInterval(i).isZero()) {
	ConsTable::iterator iter = cons.find(relDomain.xRef(i));
	assert(iter!=cons.end());
	iter->second.mult=invalidMult;
      };
    
    // For every variable in comps, check if there is another variable
    // that is in an equality relationship and has a higher multiplicity.
    // If so, try to increase the multiplicity of this variable.
    
   /*
    // The function seekMultFromEquality is not yet implemented.
    if (res!=resUnsatisfiable) {
      Mult mults[relDomain.size()];
      for (TvpiVar x = 0; x<relDomain.size(); x++) {
        ConsTable::const_iterator vIter = cons.find(relDomain.xRef(x));
        assert(vIter!=cons.end());
        assert(!vIter->second.isRange);
        mults[x]=vIter->second.mult;
      };
      for ( vector<LinComponent>::iterator cIter = comps.begin();
	    cIter!=comps.end(); cIter++) {
	ConsTable::iterator vIter = cons.find(cIter->variable);
	assert(vIter!=cons.end());
	assert(!cIter->second.isRange);
	Mult m = relDomain.seekMultFromEquality(vIter->second.var, mults);
	if (m>vIter->second.mult) {
	  unsigned long diff = m-vIter->second.mult;
	  vector<Inequality*> set(2);
	  mpz_class a(0);
	  mpz_setbit(a.get_mpz_t(), diff);
	  set.push_back(new Inequality(-1, a, 0));
	  set.push_back(new Inequality( 1,-a, 0));
	  TvpiVar temp=createVariable();
	  Result r = addInequalitySet(vIter->second.var, temp, set, true);
	  if (r==resUnsatisfiable) return resUnsatisfiable;
	  if (r==resChanged) res=resChanged;
	  update(vIter->second.var)
	  vIter->second.mult=m;
	}
      };
      */
    return res;
  }

  // Scale the coefficients and the constant.
  void Domain::multiplyByMultiplicity(std::vector<LinComponent>& comps,
				      mpz_class& c, bool isEquality) const {

#ifdef DEBUG_DIVIDEMULT
    cerr << "before multiplying with var's multiplicity: ";
#endif // DEBUG_DIVIDEMULT

    mpz_class gcd = 0;
    if (isEquality) gcd = c;

    for (vector<LinComponent>::iterator cIter = comps.begin();
	 cIter != comps.end(); cIter++) {
      ConsTable::const_iterator vIter = cons.find(cIter->variable);
#ifdef DEBUG_DIVIDEMULT
      if (cIter!=comps.begin()) cerr << " + ";
      cerr << cIter->coefficient << " x_" << cIter->variable;
      cerr << "/2^";
      if (vIter==cons.end()) cerr << "0";
      else cerr << (int) vIter->second.mult;
#endif // DEBUG_DIVIDEMULT
      // Variables that are not in the domain implicitly take on the range
      // of their type which usually has non-zero multiplicity.
      if (vIter==cons.end() &&
          calcMult(getVariableTop(cIter->variable))==invalidMult) {
        cIter->coefficient=0;
        continue;
      };
      // If the variable is not in the domain, the default range has a
      // multiplicity of zero (unless the default range is 0).
      // Multiplying with 2^0=1 is a no-op, so we don't do that.
      if (vIter==cons.end() || vIter->second.mult==0) {
        // Adjust the gcd.
        if (sgn(gcd)==0) gcd=cIter->coefficient; else
	  mpz_gcd(gcd.get_mpz_t(),
		  gcd.get_mpz_t(),
		  cIter->coefficient.get_mpz_t());
        continue;
      };
      // The variable is constant zero, remove it from the terms.
      if (vIter->second.mult==invalidMult) {
        assert((vIter->second.isRange ?
	     	ranges[vIter->second.var].interval.isZero() :
	     	relDomain.getInterval(vIter->second.var).isZero()));
	cIter->coefficient=0;
	continue;
      };
      mpz_mul_2exp(cIter->coefficient.get_mpz_t(),
		   cIter->coefficient.get_mpz_t(),
		   vIter->second.mult);
      if (sgn(gcd)==0) gcd=cIter->coefficient; else
	mpz_gcd(gcd.get_mpz_t(),
		gcd.get_mpz_t(),
		cIter->coefficient.get_mpz_t());
    };

    gcd=abs(gcd);
        
#ifdef DEBUG_DIVIDEMULT
    cerr << " + " << c;
    if (isEquality) cerr << " = 0" << endl; else cerr << " <= 0" << endl;
    cerr << "before dividing by gcd = " << gcd << ": ";
    for (vector<LinComponent>::const_iterator iter = comps.begin();
	 iter!=comps.end(); iter++) {
      if (iter!=comps.begin()) cerr << " + ";
      cerr << iter->coefficient << " x_" << iter->variable;
    };
    cerr << " + " << c;
    if (isEquality) cerr << " = 0" << endl; else cerr << " <= 0" << endl;
#endif // DEBUG_DIVIDEMULT

    if (gcd>1) {
      if (isEquality)
	mpz_divexact(c.get_mpz_t(),
		     c.get_mpz_t(),
		     gcd.get_mpz_t());
      else
	mpz_cdiv_q(c.get_mpz_t(),
		   c.get_mpz_t(),
		   gcd.get_mpz_t());
      for (vector<LinComponent>::iterator cIter = comps.begin();
	   cIter != comps.end(); cIter++)
	mpz_divexact(cIter->coefficient.get_mpz_t(),
		     cIter->coefficient.get_mpz_t(),
		     gcd.get_mpz_t());
    };

#ifdef DEBUG_DIVIDEMULT
    cerr << "after multiplying with var's multiplicity: ";
    for (vector<LinComponent>::const_iterator iter = comps.begin();
	 iter!=comps.end(); iter++) {
      if (iter!=comps.begin()) cerr << " + ";
      cerr << iter->coefficient << " x_" << iter->variable;
    };
    cerr << " + " << c;
    if (isEquality) cerr << " = 0" << endl; else cerr << " <= 0" << endl;
#endif // DEBUG_DIVIDEMULT
  }

  // Calculate a minimum multiplicity of every variable.
  Result Domain::calcNewMult(vector<LinComponent>& comps,
			     const mpz_class& c,
			     Mult newMult[]) const {
    const size_t compsSize = comps.size();
    Mult termMult[compsSize];
    Mult multRhs = calcMult(c);

    for (size_t i=0; i<compsSize; i++)
      termMult[i] = calcMult(comps[i].coefficient);

    // Refine the constant. This is only meaningful if the constant
    // isn't zero.
    if (multRhs!=invalidMult) {
      // Calculate a minimal multiplicity for the rhs. If it turns out
      // to be larger than multRhs, then the constraint turns the
      // domain unsatisfiable.
      unsigned int newMultRhs = invalidMult;
      for (size_t i=0; i<compsSize && newMultRhs>0; i++)
	if (termMult[i]<newMultRhs) newMultRhs=termMult[i];

#ifdef DEBUG_DIVIDEMULT
      cerr << "min mult of constant is at least " << (int) newMultRhs
	   << " and was " << (int) multRhs << endl;
#endif // DEBUG_DIVIDEMULT
      // If the terms imply that the rhs should have a higher
      // multiplicity then this constraint is not satisfiable at all.
      if (newMultRhs>multRhs) return resUnsatisfiable;
      // We can only increase the multiplicity of a term if the
      // right-hand side is larger than the minimum multiplicity of the
      // terms.
      if (multRhs==newMultRhs) return resRedundant;
    };

    bool changed=false;

    // Refine each term in sequence.
    for (size_t term=0; term<compsSize; term++) {
      unsigned int newMultTerm = multRhs;
      // Find a higher multiplicity of the variable.
      for (size_t cur=0; cur<compsSize && newMultTerm>termMult[term]; cur++) {
	if (cur==term) continue;
	if (termMult[cur]<newMultTerm) newMultTerm=termMult[cur];
      };
      // Check if the multiplicity of a variable has changed. Ignore if
      // the multiplicity has changed to invalidMult, since that's the
      // special case when the constraint x=0 is added.
      if (newMultTerm!=invalidMult && newMultTerm>termMult[term]) {
#ifdef DEBUG_DIVIDEMULT
	cerr << " mult of x_" << comps[term].variable
	     << " must change from " << (int) termMult[term]
	     << " to " << (int) newMultTerm << endl;
#endif // DEBUG_DIVIDEMULT
	Mult m = newMultTerm-termMult[term];
	// Multiply the coefficient of the variable that has just
	// changed in the domain. This might incur a gcd of the
	// coefficients that is greater than one, but that doesn't
	// really matter.
	mpz_mul_2exp(comps[term].coefficient.get_mpz_t(),
		     comps[term].coefficient.get_mpz_t(), m);
	changed=true;
      };
      // The new multiplicity of the variable is the calculated multiplicity minus
      // the multiplicity of the coefficient.
      newMult[term]=(newMultTerm>termMult[term] ? newMultTerm-termMult[term] : 0);
    };

    return (changed ? resChanged : resRedundant);
  }


  // Ensure that each variable occurs at most once and has
  // a non-zero coefficient.
  void Domain::canonicalizeTerms(std::vector<LinComponent>& comps) {
    if (comps.size()==0) return;

    // Join variables that are mentioned more than once.
    vector<LinComponent>::iterator cur=comps.begin();
    vector<LinComponent>::iterator next=cur+1;
    // Join variables that are mentioned more than once.
    while (next!=comps.end()) {
      if (cur->variable==next->variable) {
#ifdef DEBUG_INTERSECT
	cerr << "joining vars at pos " << cur-comps.begin() << " and "
	     << next-comps.begin() << endl;
#endif // DEBUG_INTERSECT
	next->coefficient+=cur->coefficient;
	cur->coefficient=0;
      };	
      cur=next++;
    };

    cur=comps.begin();
    while (cur!=comps.end()) {
      if (sgn(cur->coefficient)==0) {
#ifdef DEBUG_INTERSECT
	cerr << "nuking zero var in pos " << cur-comps.begin() << endl;
#endif // DEBUG_INTERSECT
	vector<LinComponent>::iterator repl=comps.end()-1;
	if (repl!=cur) cur->swap(*repl);
	comps.pop_back();
      } else cur++;
    }
  }

  // Change variable indices from DomVar to TvpiVar in the relational
  // domain. This function may promote ranges to the relational domain
  // and is therefore not const.
  bool Domain::renameVariables(vector<LinComponent>& comps, mpz_class& c) {
    if (comps.size()==0) return false;

#ifdef DEBUG_SUBSTITUTE
    cerr << "before: ";
    for (vector<LinComponent>::const_iterator iter = comps.begin();
	 iter!=comps.end(); iter++) {
      if (iter!=comps.begin()) cerr << " + ";
      cerr << iter->coefficient << " x_" << iter->variable;
    };
    cerr << " + " << c << " = 0" << endl;
#endif // DEBUG_SUBSTITUTE

#ifdef HAVE_EQUALITYSUBST
    // Replace variables with other variables that have a smaller
    // index and are in an equality relationship.
    for (vector<LinComponent>::iterator cIter = comps.begin();
	 cIter != comps.end(); cIter++) {
      if (sgn(cIter->coefficient)==0) continue;
      // Don't try to find an equality relationship if the variable is
      // not in the relational domain.
      ConsTable::const_iterator vIter = cons.find(cIter->variable);
      if (vIter==cons.end()) continue;
      if (vIter->second.isRange) continue;
#ifdef DEBUG_SUBSTITUTE
	cerr << "looking for substitute for " << cIter->coefficient
	     << "x_" << cIter->variable << endl;
#endif // DEBUG_SUBSTITUTE
      mpz_class m;
      TvpiVar oldVar = vIter->second.var;
      TvpiVar subVar = relDomain.substEquality(vIter->second.var,
					       cIter->coefficient, c, m);
      if (subVar!=oldVar) {
	// A substitution occurred.
	if (m!=1) {
	  // Multiply the coefficients by m.
	  for (vector<LinComponent>::iterator iter = comps.begin();
	       iter != comps.end(); iter++)
	    if (iter!=cIter) iter->coefficient*=m;
	};
	cIter->variable=relDomain.xRef(subVar);
#ifdef DEBUG_SUBSTITUTE
	cerr << "replacing with " << cIter->coefficient
	     << "x_" << relDomain.xRef(subVar) << endl;
	if (m!=1) {
    	  for (vector<LinComponent>::const_iterator iter = comps.begin();
	       iter!=comps.end(); iter++) {
	    if (iter!=comps.begin()) cerr << " + ";
	    cerr << iter->coefficient << " x_" << iter->variable;
    	  };
	  cerr << " + " << c << " = 0" << endl;
	};
#endif // DEBUG_SUBSTITUTE
      };
    };

#ifdef DEBUG_SUBSTITUTE
    cerr << "after: ";
    for (vector<LinComponent>::const_iterator iter = comps.begin();
	 iter!=comps.end(); iter++) {
      if (iter!=comps.begin()) cerr << " + ";
      cerr << iter->coefficient << " x_" << iter->variable;
    };
    cerr << " + " << c << " = 0" << endl;
#endif // DEBUG_SUBSTITUTE

    // Remove duplicate variables and zero entries.
    sort(comps.begin(), comps.end());
    canonicalizeTerms(comps);

#endif // HAVE_EQUALITYSUBST

    // Inline those variables that have a constant value. Promote
    // those to the relational domain that don't. Track up to two
    // variables that are not in the domain.
    size_t fstVar = (size_t) -1;
    size_t sndVar = (size_t) -1;

    size_t cur=0;
    while (cur<comps.size()) {
      ConsTable::iterator vIter = cons.find(comps[cur].variable);
      if ((vIter==cons.end() && comps[cur].variable<invalidVar) ||
	  (vIter!=cons.end() && vIter->second.isRange &&
	   ranges[vIter->second.var].interval.upperIsInfinite() &&
	   ranges[vIter->second.var].interval.lowerIsInfinite())) {
	// This variable is unbounded. Track the index of this
	// variable.
#ifdef DEBUG_SUBSTITUTE
	cerr << "variable at pos " << cur << " is unbounded" << endl;
#endif // DEBUG_SUBSTITUTE
	if (sndVar!= (size_t) -1) return false;
	sndVar=fstVar;
	fstVar=cur;
	cur++;
      } else {
	bool isConst = 
	  (vIter==cons.end() ?
	   getVariableTop(comps[cur].variable).isSingleton() :
	   (vIter->second.isRange ?
	    ranges[vIter->second.var].interval.isSingleton() :
	    relDomain.getInterval(vIter->second.var).isSingleton())
	   );
	if (isConst) {
#ifdef DEBUG_SUBSTITUTE
	  cerr << "inlining constant variable at " << cur <<endl;
#endif // DEBUG_SUBSTITUTE

	  // Incorporate constant-valued variable into the constant c.
	  if (vIter==cons.end())
	    comps[cur].coefficient*=
	      getVariableTop(comps[cur].variable).getUpper();
	  else if (vIter->second.isRange)
	    comps[cur].coefficient*=
	      ranges[vIter->second.var].interval.getUpper();
	  else
	    comps[cur].coefficient*=
	      relDomain.getInterval(vIter->second.var).getUpper();
	  c+=comps[cur].coefficient;
	  size_t repl=comps.size()-1;
	  if (repl!=cur) comps[cur].swap(comps[repl]);
	  comps.pop_back();
	} else {
	  // Replace the DomVar with a TvpiVar of the relational
	  // domain. Promote DomVar to relational domain if necessary.
	  if (vIter==cons.end())
	    vIter=createTvpiRange(comps[cur].variable,
				  getVariableTop(comps[cur].variable));
	  if (vIter->second.isRange) promoteRange(vIter);
#ifdef DEBUG_SUBSTITUTE
	  cerr << "translating DomVar ";
	  showVar(cerr, comps[cur].variable);
	  cerr << " to TvpiVar " << vIter->second.var << endl;
#endif // DEBUG_SUBSTITUTE
	  comps[cur].variable=vIter->second.var;
	  cur++;
	}
      }
    };

    assert(sane());

    // Deal with up to two terms whose variables are unbounded.
    while (fstVar!= (size_t) -1) {
      TvpiVar var;
      ConsTable::iterator vIter = cons.find(comps[fstVar].variable);
      if (vIter==cons.end()) {
	vIter = createTvpiRange(comps[fstVar].variable,
				getVariableTop(comps[fstVar].variable));
	var = vIter->second.var;
      } else if (vIter->second.isRange)
	var = promoteRange(vIter);
      else var = vIter->second.var;

#ifdef DEBUG_SUBSTITUTE
      cerr << "patching unrestricted var ";
      showVar(cerr, comps[fstVar].variable);
      cerr << " at pos " << fstVar << " mapped to p_" << var << endl;
#endif // DEBUG_SUBSTITUTE
      assert(!vIter->second.isRange);
      comps[fstVar].variable = var;
      fstVar=sndVar;
      sndVar=(size_t) -1;
    };

    assert(sane());

#ifdef DEBUG_SUBSTITUTE
    cerr << "lin expr over TVPI:";
    for(vector<LinComponent>::const_iterator termIter=comps.begin();
	termIter!=comps.end(); termIter++) {
      cerr << " ";
      if (termIter!=comps.begin() &&
	  termIter->coefficient>=0) cerr << "+";
      cerr << termIter->coefficient << " p_" << termIter->variable;
    };
    switch (sgn(c)) {
    case -1: cerr << c; break;
    case 0: break;
    case 1: cerr << " + " << c;
    };
    cerr << endl;
 #endif // DEBUG_SUBSTITUTE

    return true;
  }

  // Query the relational information between two variables.
  Polyhedron<isIntegral>* Domain::queryPolyhedron(DomVar vx, DomVar vy) const {
    Domain::ConsTable::const_iterator xPtr = cons.find(vx);
    Domain::ConsTable::const_iterator yPtr = cons.find(vy);
    if (xPtr==cons.end()) return 0;
    if (yPtr==cons.end()) return 0;
    if (xPtr->second.isRange) return 0;
    if (yPtr->second.isRange) return 0;
    if (xPtr->second.var==yPtr->second.var) return 0;
    if (xPtr->second.var<yPtr->second.var) return
      new Polyhedron<isIntegral>(relDomain.
				 getProjection(xPtr->second.var,
					       yPtr->second.var).
				 stretch(xPtr->second.mult,
					 yPtr->second.mult));
    Polyhedron<isIntegral>* p = 
      new Polyhedron<isIntegral>(relDomain.getProjection(yPtr->second.var,
							 xPtr->second.var).
				 stretch(yPtr->second.mult,
					 xPtr->second.mult));
    p->swapVars();
    return p;
  }
  
  void Domain::showVar(std::ostream& s, DomVar v) const {
    if (v==invalidVar) s << "invalid"; else 
      if (v>0) s << "x_" << v; else s << "t_" << (-v);
  }
  
  void Domain::showTerms(ostream& s,
			 const std::vector<LinComponent>& terms,
			 const mpz_class& c) const {
    for(std::vector<LinComponent>::const_iterator termIter=terms.begin();
	termIter!=terms.end(); termIter++) {
      s << " ";
      if (termIter!=terms.begin() &&
	  termIter->coefficient>=0) s << "+";
      s << termIter->coefficient << " ";
      showVar(s, termIter->variable);
    };
    switch (sgn(c)) {
    case -1: s << c; break;
    case 0: break;
    case 1: s << " + " << c;
    };
  }

  // Check if the relationship between DomVars and TvpiVars is consistent.
  bool Domain::sane() const {
    relDomain.sane();
    bool isSane=true;
    for (ConsTable::const_iterator iter = cons.begin();
	 iter!=cons.end(); iter++) {
      if (iter->second.isRange) {
	if (ranges[iter->second.var].xRef!=iter->first) {
	  cerr << "DomVar ";
	  showVar(cerr, iter->first);
	  cerr << " maps to range no. " << iter->second.var
	       << " but the range maps back to ";
	  showVar(cerr, ranges[iter->second.var].xRef);
	  cerr << endl;
	  isSane=false;
	}
      } else {
	if (relDomain.xRef(iter->second.var)!=iter->first) {
	  cerr << "DomVar ";
	  showVar(cerr, iter->first);
	  cerr << " maps to relational var no. " << iter->second.var
	       << " but the TvpiVar maps back to ";
	  showVar(cerr, relDomain.xRef(iter->second.var));
	  cerr << endl;
	  isSane=false;
	}
      }
    };
    for (vector<RangeInfo>::const_iterator iter = ranges.begin();
	 iter!=ranges.end(); iter++) {
      ConsTable::const_iterator cIter = cons.find(iter->xRef);
      if (cIter==cons.end()) {
	cerr << "Range no " << iter-ranges.begin()
	     << " maps to non-existent DomVar ";
	showVar(cerr, iter->xRef);
	cerr << endl;
	isSane=false;
      }
    };
    for (size_t i = 0; i<relDomain.size(); i++) {
      DomVar xRef = relDomain.xRef(i);
      ConsTable::const_iterator cIter = cons.find(xRef);
      if (cIter==cons.end()) {
	cerr << "Relational var no " << i
	     << " maps to non-existent DomVar ";
	showVar(cerr, xRef);
	cerr << endl;
	isSane=false;
      }
    };
    return isSane;
  };

  defMemDbg(,Domain,d,D)

};
