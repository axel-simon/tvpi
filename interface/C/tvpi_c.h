/* tvpi_c.h - C header file that reexports the TVPI classes as C
 *   functions and types.
 */

#ifndef TVPI_C_H_
#define TVPI_C_H_

#include "common.hh"

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#undef TVPI_TYPE_DECLARATION
#define TVPI_TYPE_DECLARATION(Type) \
typedef struct Type ## _tag Type ## _t; \
typedef Type ## _t * Type ## _p; \
typedef Type ## _t const* const_ ## Type ## _p

  TVPI_TYPE_DECLARATION(Polyhedron);
  TVPI_TYPE_DECLARATION(Inequality);
  TVPI_TYPE_DECLARATION(Domain);
  TVPI_TYPE_DECLARATION(LandmarkTable);

  LandmarkTable_p tvpiLandmarkTableNew();

  void tvpiLandmarkTableReserve(LandmarkTable_p lm, size_t statements);

  void tvpiLandmarkTableFree(LandmarkTable_p lm);

  void tvpiLandmarkTableSetStatementNo(LandmarkTable_p lm,
				       size_t statement);

  void tvpiLandmarkTableClear(LandmarkTable_p lm);

  void tvpiLandmarkTablePromote(LandmarkTable_p lm);

  // Given a set of lmCount tables stored at lm, calculate the
  // smallest number of iterations such that one of the unsatisfiable
  // inequalities added will be satisfiable.
  /* ret val      | iters | stmt | condition
     -------------|-------|------|--------------------------------------
     lmCount+1    |     0 |  (-1)| perform a normal join
     lmCount      |     0 |  (-1)| perform full widening
     0..lmCount-1 |     i |    s | extrapolate i steps
  */
  size_t tvpiLandmarkCalcNoOfIterations(LandmarkTable_p *lm,
					size_t lmCount,
					unsigned long* iters,
					int* stmt);

  void tvpiLandmarkTableDump(LandmarkTable_p lm);

  Domain_p tvpiNew(void);

  void tvpiFree(Domain_p d);

  Domain_p tvpiCopy(Domain_p d);

  // Fill the given cells with the bounds of the inteval of v in d. If
  // lower was set then the 0th result bit is set, if upper was set
  // the 1st result bit is set. The 2nd bit is set if the variables
  // was in the domain in the first place (in which case multP was
  // set). If upper was set with an unsigned value, the 3rd bit is
  // set.
  int tvpiQueryValue(Domain_p d,
		     signed long* lowerP, signed long* upperP,
		     Mult* multP, DomVar v);

  void tvpiProjectOut(Domain_p d, DomVar v);

  int tvpiVarsInDomain(Domain_p d);

  int tvpiRelVarsInDomain(Domain_p d);

  int tvpiVarInDomain(Domain_p d, DomVar v);

  DomVar tvpiGetNextTypedVarId(TypeId ty);

  DomVar tvpiGetNextTempVarId(void);

  TypeId tvpiRegisterType(TypeId ty,
			  signed long lower,
			  unsigned long upper);

  void tvpiResetVariables();

  void tvpiSetVariableToTop(Domain_p d, DomVar v);

  struct Term {
    signed long coeff;
    DomVar var;
  };

  // Intersect the domain with the given equality.
  enum Result tvpiIntersect(Domain_p d,
			    struct Term lc[], size_t lc_size,
			    int isEquality,
			    signed long val,
			    LandmarkTable_p lm);

  enum Result tvpiIntersectU(Domain_p d,
			     struct Term lc[], size_t lc_size,
			     int isEquality,
			     unsigned long val,
			     LandmarkTable_p lm);

  int tvpiSetMultiplicity(Domain_p t1,
			  DomVar var,
			  unsigned int mult);

  void tvpiAugment(Domain_p d, DomVar source, DomVar target);

  void tvpiUpdate(Domain_p d, DomVar source, DomVar target);

  void tvpiSetVariable(Domain_p d, DomVar v,
		       signed long* lowerP, signed long* upperP, Mult mult);

  void tvpiSetVariableU(Domain_p d, DomVar v,
			signed long* lowerP, unsigned long* upperP, Mult mult);

  void tvpiDemote(Domain_p d1);

  int tvpiEntails(Domain_p d1, Domain_p d2);

  void tvpiWiden(Domain_p d1, Domain_p d2, unsigned long steps);

  void tvpiJoin(Domain_p d1, Domain_p d2);

  void tvpiDump(const_Domain_p d);

  void tvpiSane(const_Domain_p d);

  void tvpiQueryXRefs(const_Domain_p d, DomVar* xRefPtr, int* isRelPtr);

  const_Polyhedron_p tvpiQueryPolyhedron(const_Domain_p d, DomVar x, DomVar y);

  void tvpiFreePolyhedron(Polyhedron_p p);

  int tvpiPolyhedronGetNoOfInequalities(const_Polyhedron_p p);

  const_Inequality_p tvpiPolyhedronGetNth(const_Polyhedron_p p, int idx);

  int tvpiInequalityGetLength(const_Inequality_p e);

  void tvpiInequalityWrite(const_Inequality_p e, char* s);
  
#ifdef __cplusplus
}
#endif // __cplusplus

#endif // TVPI_C_H_

