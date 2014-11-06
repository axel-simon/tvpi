/* Affine relationships.
 *
 * This module implements lightweight affine relationships between
 * variables. The equations represented are generalized to
 * inequalities, i.e. a typical constraint is [3..5] = x + 3 y + 2 z,
 * which describes two inequalities bounding the linear terms from
 * below and above. y and z are expressed with either simple intervals
 * or in a TVPI domain.
 *
 * The description above is bullocks. This module doesn't do any of that.
 * It implements a reduced product between the multiplicity domain and
 * the TVPI/interval domain.
 */

#ifndef __AFFINE_H
#define __AFFINE_H

#include "common.hh"
#include "interval.hh"
#include "tvpi.hh"
#include "landmark.hh"
#include "memory.hh"
#include <ext/hash_map>
#include <vector>

#undef DEBUG_REFCOUNT

#ifdef __GNUC__
namespace stdext =__gnu_cxx;
#endif

#define isIntegral true

// Declaration needed to make main a friend. This is used by the test
// programs to verify private information within the classes.
#ifdef IN_TVPI_CHECK
int main(int);
#endif // IN_TVPI_CHECK

#define checkDomVar(var,retVal) \
  if (var>=(DomVar) typedVars.size()) return retVal

namespace Tvpi {

  class Domain;

  class VarInfo {

    friend class Domain;
    friend std::ostream& operator<<(std::ostream& s, const Domain& d);

  public:

    // Constructor taking the Variable number and an isRange flag.
    inline VarInfo(TvpiVar v, bool r, Mult m) : 
      mult(m), isRange(r), var(v) {};

    // No of lower zero bits of this variable. The domain value must be
    // multiplied by 2^multiplicity to get the actual value of this
    // variable. The value in this variable ranges between 0 and
    // invalidMult, the latter is used to represent the vaulue zero.
    Mult mult;

    // This flag is true if the variable index refers to an index into
    // the ranges vector, if it is false, it refers to the relational
    // domain.
    bool isRange;

    // This variable can be found at this particular index in the relational
    // domain or in the ranges vector of the domain.
    TvpiVar var;

  };

  typedef std::vector<DomVar> VarSet;

  class Domain {

  private:
    typedef stdext::hash_map<DomVar, VarInfo> ConsTable;

    // This table maps domain variables to their descriptions. If a
    // variable is not contained in here, it takes on the range that
    // is given by its type. Temporary variable (which don't have a
    // type) are unbounded.
    ConsTable cons;

    // Counters of temporary domain variables. This counter starts
    // negative and always decreases.
    static DomVar lastTempVar;

    // A table of all non-temporary variables and their type.
    typedef std::vector<TypeId> VarTable;
    static VarTable typedVars;

    // A table of all types.
    typedef std::vector<Interval<isIntegral> > TypeTable;
    static TypeTable registeredTypes;

    // The illegal variable, never created but used internally as a
    // placeholder.
    static const DomVar invalidVar;

    // The underlying relational domain.
    DenseTvpi<isIntegral> relDomain;

    struct RangeInfo {
      RangeInfo(const Interval<isIntegral> i, DomVar v) :
	interval(i), xRef(v) {};
      Interval<isIntegral> interval;
      DomVar xRef;
    };

    typedef std::vector<RangeInfo> Ranges;

    // This vector contains variables that are not part of the
    // relational domain.
    Ranges ranges;

  public:
    // Construct a new domain object.
    Domain() {
      if (registeredTypes.empty()) resetVariables();
    };

    // Retrieve the multiplicity and/or the value of a variable. The
    // 0th bit in the return value is set if the lower bound was set,
    // the 1st bit in the return value is set if the upper bound was
    // set. If the 3rd bit is unset, then the
    // variable can take on exactly those values that are implied by
    // its type. The 4th bit is set if the upper bound was set to an
    // unsigned value.
    int queryValue(signed long* low,
		   signed long* upp,
		   Mult* m, DomVar v) const;

    // Set the minimum multiplicity of a variable. Returns false if
    // the domain became unsatisfiable.
    bool enforceMult(DomVar var, Mult m);

  public:
    // Project out a variable. The variable will be unbounded after
    // this call.
    void projectOut(DomVar var);

  private:

    // Remove a range.
    void removeRange(TvpiVar var);

    // Remove a tvpi variable.
    void removeTvpi(TvpiVar var);

    // Promote a variable into the relational domain.
    TvpiVar promoteRange(ConsTable::iterator& iter);

    // Change an entry in the cons table from refering to a relational
    // variable to refering to a range. If the range is equivalent to
    // the type of the variable, v is deleted completely.
    void redirectToRange(TvpiVar v);

  public:

    // Move all variables that out of the relational domain that have
    // no relational information associated with them.
    void demote();

  private:

    // Create a new range.
    inline ConsTable::iterator createRange(DomVar var,
					   const Interval<isIntegral>& i) {
      assert(cons.find(var)==cons.end());
      TvpiVar v = ranges.size();
      Mult m = calcMult(i);
      ranges.push_back(RangeInfo(i, var));
      return cons.insert(ConsTable::value_type(var,
					       VarInfo(v, true, m))).first;
    };

    // Create a range in the polyhedral domain.
    inline 
    ConsTable::iterator createTvpiRange(DomVar var,
					const Interval<isIntegral>& i) {
      assert(cons.find(var)==cons.end());
      TvpiVar res = relDomain.createVariable(i, var);
      Mult m = calcMult(i);
      return cons.insert(ConsTable::value_type(var,
					       VarInfo(res, false, m))).first;
    };

  public:
    // Create a new variable. The resulting variable is restricted to
    // the given type, in all domains in which it is mentioned.
    inline static DomVar createVariable(TypeId ty=0) {
      if (registeredTypes.size()<=ty) registeredTypes.resize(ty+1);
      DomVar var = typedVars.size();
      typedVars.push_back(ty);
      return var;
    };

    inline static DomVar createTempVariable() { return lastTempVar--; };

    // Reset all type information, the temporary variable counter, etc.
    inline static void resetVariables() { 
      lastTempVar=-1;
      registeredTypes.clear();
      typedVars.clear();
#ifndef NDEBUG
      DomVar first =
#endif
	createVariable();
      assert(first==invalidVar);
      createVariable();
    };

    // Register a new type. The given interval determines upper and
    // lower bounds that are conservative assumptions on possible
    // values of variables of that type. If the given type argument is
    // (-1), a new type id will be created. The a non-negative type id
    // is given, the bounds for that type id are replace. This feature
    // is useful to register types with fixed ids and should probably
    // not be used once the analysis runs.
    TypeId static registerType(const Interval<isIntegral>& c,
			       TypeId ty = (TypeId) -1) {
      if (ty != (TypeId) -1) {
	if (registeredTypes.size()<=ty) registeredTypes.resize(ty+1);
	registeredTypes[ty]=c;
      } else {
	ty = registeredTypes.size();
	registeredTypes.push_back(c);
      };
      return ty;
    };
      
    // Let the given target variable contain all the values it had before and
    // all values of the source variable.
    void augment(DomVar source, DomVar target);

    // Replace the given target variable with the valuation of another
    // variable.
    void update(DomVar source, DomVar target);

    // Set a variable to its maximum bounds as given by its registered
    // type. Set temporary variables to an unbounded range.
    void setVariableToTop(DomVar var) {
      ConsTable::iterator iter = cons.find(var);
      if (iter==cons.end()) return;
      if (iter->second.isRange)
	removeRange(iter->second.var);
      else
	removeTvpi(iter->second.var);
      cons.erase(iter);
    };

    // Return the maximum range for a given variable. Using this
    // function ensures that the ranges in the registeredType array
    // are not changed.
    inline static const Interval<isIntegral> getVariableTop(DomVar var) {
      if (var<=invalidVar)
	return Interval<isIntegral>();
      else
	return registeredTypes[typedVars[var]];
    };

    // Query the variable type.
    inline static const TypeId getVariableType(DomVar var) {
      assert(var>=0);
      assert((size_t) var<typedVars.size());
      return typedVars[var];
    };

    // Calculate the maximum of the given expression and return it
    // in value. Returns false if the expression is unbounded.
    bool linOpt(std::vector<LinComponent>& comps, mpz_class& value);

    // Set a given variable to a specific value. The variable is
    // projected out first if it is still in the domain. If the given
    // interval an the multiplicity information imply an empty domain,
    // false is returned.
    bool setVariable(DomVar var, Mult m, Interval<isIntegral> c);

    // Intersect with inequality a_1 x_1 + ... + a_n x_n + c <= 0. If
    // isEquality is true, replace "<=" with "=". The multTgt denotes
    // the variable that is defined by this equation, it may refer to
    // a non-existent variable. The function takes an optional
    // parameter to a LandmarkTable class which tracks redundant
    // inequalities.
    Result inequality(std::vector<LinComponent>& comps, mpz_class c,
		      LandmarkTable* lm = NULL,
		      bool isEquality = false);

  private:

    // Descale the coefficients and the constant.
    void multiplyByMultiplicity(std::vector<LinComponent>& comps,
				mpz_class& c, bool isEquality) const;

    // Calculate a minimum multiplicity of every variable. Change the
    // coefficients and the constant if necessary. Note that the
    // function can only refine the multiplicity in an equality
    // constraint. If such an equality constraint implies a higher
    // multiplicity of the constant, the return value will inevitably
    // be resUnsatisfiable, hence the constness of c. A return value
    // of resChanged implies that the newMult array reflects the
    // minimum multiplicity of each variable.
    Result calcNewMult(std::vector<LinComponent>& comps,
		       const mpz_class& c,
		       Mult newMult[]) const;

    // Sort terms, ensure that each variable occurs at most once and
    // has a non-zero coefficient.
    void canonicalizeTerms(std::vector<LinComponent>& comps);

    // Inline constant variables and rename variables to the
    // relational domain. Returns false if the linear expression
    // ranges over more than two variables that are not in the domain,
    // to indicate that approximation in the TVPI domain is not
    // possible. If the result is true, all variables in the linear
    // expressions are in the TVPI domain. In this successful case,
    // the variables are changed from DomVars to TvpiVars.
    bool renameVariables(std::vector<LinComponent>& comps,
			 mpz_class& c);

  public:
    // Check if this domain describes a sub-space of the prev domain.
    bool entails(Domain& prev);

    // Calculate the join of both and store the result in this domain. 
    // If extrapolate is greater than zero, extrapolate the change
    // between the previous and the join 'extrapolate' times. If
    // extrapolate is 0, widen the join with respect to the previous
    // domain. If extrapolate is negative, do nothing after
    // calculating the join. Note the weird argument order: The
    // previous domain that is already stored at a program location
    // remains untouched. It it the new, changed state that is joined
    // with the previous state.
    void joinWiden(Domain& prev, mpz_class extrapolate);

    // Ask for the number of non-top variables in this instance.
    int varsInDomain() const { return cons.size(); };

    // Ask for the number of relational variables in this instance.
    int relVarsInDomain() const { return relDomain.size(); }

    // Ask for the xRefs of all relational variables. The first
    // passed-in array must have n elements, the second one (n*n-n)/2
    // elements where n is the number returned by relVarsInDomain.
    void queryXRefs(DomVar* xRefPtr, int* relPtr) const {
      relDomain.queryXRefs(xRefPtr, relPtr);
    }

    // Check if a given variable is in the domain.
    bool isInDomain(DomVar v) const { return cons.find(v)!=cons.end(); };

    // Query the relational information between two variables. The
    // function returns a planar polyhedron that contains
    // inequalities. The polyhedron has to be freed. The result is
    // NULL if there is no relational information between the given
    // variables.
    Polyhedron<isIntegral>* queryPolyhedron(DomVar vx, DomVar vy) const;

  private:

    // Show the name of a variable.
    void showVar(std::ostream& s, DomVar v) const;

    void showTerms(std::ostream& s,
		   const std::vector<LinComponent>& terms,
		   const mpz_class& c) const;

  public:

    declMemDbg;

    // Dump the Domain dude.
    friend std::ostream& operator<<(std::ostream& s, const Domain& d);

    friend inline bool 
    compareVariable(Domain::ConsTable::const_iterator iterA,
		    Domain::ConsTable::const_iterator iterB);
#ifdef IN_TVPI_CHECK
    friend int ::main(int);
#endif // IN_TVPI_CHECK

    // Check if the relationship between DomVars and TvpiVars is consistent.
    bool sane() const;

  };

  std::ostream& operator<<(std::ostream& s, const Domain& d);

};

#endif
