// tvpi_c.cpp - C interface to the C++ Tvpi class.
//

#include "tvpi_c.h"
#include "lincomponent.hh"
#include "landmark.hh"
#include "tvpi.hh"
#include "planar.hh"
#include "affine.hh"
#include <iostream>
#include <vector>
#include <assert.h>

using namespace Tvpi;

#define isZ true


template<class ToType, class FromType>
const ToType* to_const(const FromType* x) {
  return reinterpret_cast<const ToType*>(x);
}

template<class ToType, class FromType>
ToType* to_nonconst(FromType* x) {
  return reinterpret_cast<ToType*>(x);
}


LandmarkTable_p tvpiLandmarkTableNew() {
  return to_nonconst<LandmarkTable_t, LandmarkTable>(new LandmarkTable());
}

void tvpiLandmarkTableReserve(LandmarkTable_p lm, size_t statements) {
  to_nonconst<LandmarkTable, LandmarkTable_t>(lm)->reserve(statements);
}

void tvpiLandmarkTableFree(LandmarkTable_p lm) {
  delete to_nonconst<LandmarkTable, LandmarkTable_t>(lm);
}

void tvpiLandmarkTableSetStatementNo(LandmarkTable_p lm,
				     size_t statement) {
  to_nonconst<LandmarkTable, LandmarkTable_t>(lm)->setStatementNo(statement);
}

void tvpiLandmarkTableClear(LandmarkTable_p lm) {
  to_nonconst<LandmarkTable, LandmarkTable_t>(lm)->clear();
}

void tvpiLandmarkTablePromote(LandmarkTable_p lm) {
  to_nonconst<LandmarkTable, LandmarkTable_t>(lm)->promote();
}

size_t tvpiLandmarkCalcNoOfIterations(LandmarkTable_p *lm,
				      size_t lmCount,
				      unsigned long* iters,
				      int* stmt) {
  size_t lmIdx=lmCount;
  WideningInfo wi;
  if (lm==NULL) return lmCount;
  for (size_t i=0; i<lmCount; i++)
    if (to_nonconst<LandmarkTable, LandmarkTable_t>(lm[i])->
	calcNoOfIterations(wi)) lmIdx = i;
  mpz_class iterations = wi.getIterations();
  if (iters) *iters=0;
  if (sgn(iterations)<0) lmIdx=lmCount+1; // Don't apply widening.
  if (iters && mpz_fits_ulong_p(iterations.get_mpz_t()))
    *iters=mpz_get_ui(iterations.get_mpz_t());
  if (stmt) *stmt = wi.getStatementNo();
  return lmIdx;
}

void tvpiLandmarkTableDump(LandmarkTable_p lm) {
  std::cerr << *to_const<LandmarkTable, LandmarkTable_t>(lm);
}

Domain_p tvpiNew() {
  return to_nonconst<Domain_t, Domain>(new Domain());
}

void tvpiFree(Domain_p d) {
  delete to_nonconst<Domain, Domain_t>(d);
}

Domain_p tvpiCopy(Domain_p d) {
  return to_nonconst<Domain_t, Domain>
    (new Domain(*to_const<Domain, Domain_t>(d)));
}

// Fill the given cells with the bounds of the inteval of v in d. If
// upper was set then the 0th result bit is set, if lower was set the
// 1st result bit is set.
int tvpiQueryValue(Domain_p d,
		   signed long* lowerP, signed long* upperP,
		   Mult* multP, DomVar v) {
  return to_const<Domain, Domain_t>(d)->queryValue(lowerP, upperP, multP, v);
}

void tvpiProjectOut(Domain_p d, DomVar v) {
  to_nonconst<Domain, Domain_t>(d)->projectOut(v);
}

int tvpiVarsInDomain(Domain_p d) {
  return to_const<Domain, Domain_t>(d)->varsInDomain();
}

int tvpiRelVarsInDomain(Domain_p d) {
  return to_const<Domain, Domain_t>(d)->relVarsInDomain();
}

int tvpiVarInDomain(Domain_p d, DomVar v) {
  return (to_const<Domain, Domain_t>(d)->
	  queryValue(NULL, NULL, NULL, v)&4)!=0;
}

DomVar tvpiGetNextTypedVarId(TypeId ty) {
  return Domain::createVariable(ty);
}

DomVar tvpiGetNextTempVarId() {
  return Domain::createTempVariable();
}

TypeId tvpiRegisterType(TypeId ty,
			signed long lower,
			unsigned long upper) {
  Interval<isZ> c;
  c.updateLower(mpq_class(lower));
  c.updateUpper(mpq_class(upper));
  return Domain::registerType(c, ty);
}

void tvpiResetVariables() {
  Domain::resetVariables();
}

void tvpiSetVariableToTop(Domain_p d, DomVar v) {
  to_nonconst<Domain, Domain_t>(d)->setVariableToTop(v);
}

Result tvpiIntersect(Domain_p d,
		     struct Term lc[], size_t lc_size,
		     int isEquality,
		     signed long val,
		     LandmarkTable_p lm) {
  std::vector<LinComponent> terms;
  terms.reserve(lc_size);
  for (struct Term* lcp = &lc[0]; lcp<&lc[lc_size]; lcp++)
    terms.push_back(LinComponent(mpz_class(lcp->coeff), lcp->var));
  return to_nonconst<Domain, Domain_t>(d)->
    inequality(terms, mpz_class(val),
	       to_nonconst<LandmarkTable, LandmarkTable_t>(lm), isEquality);
}

Result tvpiIntersectU(Domain_p d,
		      struct Term lc[], size_t lc_size,
		      int isEquality,
		      unsigned long val,
		      LandmarkTable_p lm) {
  std::vector<LinComponent> terms;
  terms.reserve(lc_size);
  for (struct Term* lcp = &lc[0]; lcp<&lc[lc_size]; lcp++)
    terms.push_back(LinComponent(mpz_class(lcp->coeff), lcp->var));
  return to_nonconst<Domain, Domain_t>(d)->
    inequality(terms, mpz_class(val),
	       to_nonconst<LandmarkTable, LandmarkTable_t>(lm), isEquality);
}

int tvpiSetMultiplicity(Domain_p t1,
                        DomVar var,
                        unsigned int mult) {
  return to_nonconst<Domain, Domain_t>(t1)->enforceMult(var, mult);
}

void tvpiAugment(Domain_p d, DomVar source, DomVar target) {
  to_nonconst<Domain, Domain_t>(d)->augment(source, target);
}

void tvpiUpdate(Domain_p d, DomVar source, DomVar target) {
  to_nonconst<Domain, Domain_t>(d)->update(source, target);
}

void tvpiSetVariable(Domain_p d, DomVar v,
		     signed long* lowerP, signed long* upperP, Mult mult) {
  Interval<isZ> i;
  if (lowerP) i.updateLower(*lowerP);
  if (upperP) i.updateUpper(*upperP);
  to_nonconst<Domain, Domain_t>(d)->setVariable(v, mult, i);
}

void tvpiSetVariableU(Domain_p d, DomVar v,
		      signed long* lowerP, unsigned long* upperP, Mult mult) {
  Interval<isZ> i;
  if (lowerP) i.updateLower(*lowerP);
  if (upperP) i.updateUpper(*upperP);
  to_nonconst<Domain, Domain_t>(d)->setVariable(v, mult, i);
}

void tvpiDemote(Domain_p d1) {
  to_nonconst<Domain, Domain_t>(d1)->demote();
}

int tvpiEntails(Domain_p d1, Domain_p d2) {
  return to_nonconst<Domain, Domain_t>(d1)->entails(*to_nonconst<Domain, Domain_t>(d2));
}

// Note: the arguments get swapped here. The first domain is joined
// into the second from the caller's point of view.
void tvpiWiden(Domain_p d1, Domain_p d2, unsigned long steps) {
  mpz_class w(steps);
  to_nonconst<Domain, Domain_t>(d2)->joinWiden(*to_nonconst<Domain, Domain_t>(d1), w);
}

void tvpiJoin(Domain_p d1, Domain_p d2) {
  mpz_class w(-1);
  to_nonconst<Domain, Domain_t>(d2)->joinWiden(*to_nonconst<Domain, Domain_t>(d1), w);
}

void tvpiDump(const_Domain_p d) {
  std::cerr << *to_const<Domain, Domain_t>(d);
}

void tvpiSane(const_Domain_p d) {
  to_const<Domain, Domain_t>(d)->sane();
}

void tvpiQueryXRefs(const_Domain_p d, DomVar* xRefPtr, int* isRelPtr) {
  to_const<Domain, Domain_t>(d)->queryXRefs(xRefPtr, isRelPtr);
}

const_Polyhedron_p tvpiQueryPolyhedron(const_Domain_p d, DomVar x, DomVar y) {
  return to_nonconst<Polyhedron_t, Polyhedron<isIntegral> >
    (to_const<Domain, Domain_t>(d)->queryPolyhedron(x,y));
}

void tvpiFreePolyhedron(Polyhedron_p p) {
  delete to_nonconst<Polyhedron<isIntegral> , Polyhedron_t>(p);
}

int tvpiPolyhedronGetNoOfInequalities(const_Polyhedron_p p) {
  return to_const<Polyhedron<isIntegral> , Polyhedron_t>(p)->
    getNoOfInequalities();
}

const_Inequality_p tvpiPolyhedronGetNth(const_Polyhedron_p p, int idx) {
  return to_const<Inequality_t, Inequality>
    ((*to_const<Polyhedron<isIntegral>, Polyhedron_t>(p))[idx]);
}

int tvpiInequalityGetLength(const_Inequality_p e) {
  int s=6;
  s+=mpz_sizeinbase(to_const<Inequality, Inequality_t>(e)->getA().get_mpz_t(),
		    10);
  s+=mpz_sizeinbase(to_const<Inequality, Inequality_t>(e)->getB().get_mpz_t(),
		    10);
  s+=mpz_sizeinbase(to_const<Inequality, Inequality_t>(e)->getC().get_mpz_t(),
		    10);
  return s;
}

void tvpiInequalityWrite(const_Inequality_p e, char* s) {
  mpz_get_str(s, 10, to_const<Inequality, Inequality_t>(e)->getA().get_mpz_t());
  while (*s) s++;
  *s++='x';
  mpz_get_str(s, 10, to_const<Inequality, Inequality_t>(e)->getB().get_mpz_t());
  while (*s) s++;
  *s++='y';
  mpz_get_str(s, 10, to_const<Inequality, Inequality_t>(e)->getC().get_mpz_t());
}
