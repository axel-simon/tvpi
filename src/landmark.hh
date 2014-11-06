// landmark.hh: storing inequalities that are currently redundant but may
// affect the state space in later iterations

#ifndef __LANDMARK_H
#define __LANDMARK_H

#include "common.hh"
#include "lincomponent.hh"
#include <gmpxx.h>
#include <vector>
#include <list>

#ifdef __GNUC__
namespace stdext =__gnu_cxx;
#endif

namespace Tvpi {

  class Landmarks;
  class LandmarkTable;

  // Information on the landmark that is reached next. This structure
  // is passed to Landmarks::calcNoOfIterations when then fills this
  // with the appropriate information.
  class WideningInfo {
    // The number of iterations until the next landmark becomes a
    // non-redundant inequality.
    mpz_class iterations;

    // This flag is true if iterations contains a valid value.
    bool hasIterations;

    // This flag is true if any of the landmarks had only one value
    // set which means that no estimate of iteration numbers could be
    // inferred. If this flag is set, it might be worthwhile to do
    // another iteration through the SCC to infer a second bound on
    // the landmarks.
    bool hasPending;

    // The statement number of the landmark that has last fired.
    int statementNo;

  public:
    // Constructor that initializes this structure so that no
    // elements are set.
    WideningInfo() : hasIterations(false),
		     hasPending(false),
		     statementNo(-1) {};

    // Ask for the number of iterations. This call resets the landmark
    // that was found when this structure was passed to
    // LandmarkTable::calcNoOfIterations. The return value is -1 if no
    // widening should be done, it is zero if standard widening should
    // be applied and it is positive for extrapolation of the given
    // number of steps.
    inline mpz_class getIterations() {
      if (hasPending) {
	hasPending=false;
	hasIterations=false;
	return -1;
      };
      if (hasIterations) {
	hasIterations=false;
	return iterations;
      };
      return 0;
    };

    // Ask for the statement number that determined the number of
    // iterations. The return value is indeterminate if the last call
    // to getIterations did not return a positive value.
    inline int getStatementNo() { return statementNo; };

    friend class Landmarks;
    friend class LandmarkTable;
  };

  // A set of landmarks. Each landmark is a TVPI expression with two
  // rhs values. The rhs values denote the distance of the halfspace
  // that the linear expression describes from the feasible space. Two
  // rhs values are stored to express the distance from the feasible
  // space in the last iteration and in the current iteration.
  class Landmarks {

    struct Landmark {
      std::vector<LinComponent> terms;
      mpz_class curVal, prevVal;
      Landmark(const std::vector<LinComponent>& ts,
	       mpz_class& val) : terms(ts), curVal(val), prevVal(0) {};

      friend std::ostream& operator<<(std::ostream& s, const Landmark& lm);

    };

    // The list of landmarks.
    typedef std::list<Landmark> LmList;
    LmList lmlist;

  public:
    // Add a new landmark. Save the unsatisfiable inequality by storing
    // its normal vector. The value denotes the distance
    // between the half-space and the currently feasible space. It may
    // not be negative (since the inequality would otherwise not be
    // feasible).
    void addLandmark(const std::vector<LinComponent>& ts, mpz_class& val);

    // Promote the current values to previous value and mark the
    // current values as infinite.
    inline void promote() {
      for(LmList::iterator iter = lmlist.begin(); iter!=lmlist.end(); iter++) {
	mpz_swap(iter->prevVal.get_mpz_t(), iter->curVal.get_mpz_t());
	if (sgn(iter->curVal)!=0) iter->curVal=0;
      }
    };

    // Empty the set of landmarks. This function must be called after
    // each widening since the widened state does not preserve the
    // distances between the redundant inequalities and the feasible
    // space.
    inline void clear() {
      lmlist.clear();
    };

    // Calculate the number of iterations until any of the redundant
    // inequalities would become non-redundant. Returns true if an
    // extrapolation statement was found.
    bool calcNoOfIterations(WideningInfo& li);

    // Print the set of landmarks.
    friend std::ostream& operator<<(std::ostream& s, const Landmarks& lm);
    friend std::ostream& operator<<(std::ostream& s, const Landmark& lm);

  };

  // A LandmarkTable contains as-of-now redundant inequalities for a
  // whole basic block. A basic block can contain several statements. 
  // The LandmarkTable contains a separate list of redundant
  // inequalities for each of statement.
  class LandmarkTable {

    // An array of lists of landmarks. Each element of this array
    // contains the landmarks of one statement.
    std::vector<Landmarks> landmarks;

    // The current statement number.
    size_t statement;

  public:
    LandmarkTable() : statement(0) {};

    // Reserve space for the given number of statements.
    inline void reserve(size_t statements) {
      landmarks.reserve(statements);
    };

    // Set the statement number.
    inline void setStatementNo(size_t stmtNo) { statement = stmtNo; };

    // Add a new landmark.
    void addLandmark(std::vector<LinComponent>& ts,
		     mpz_class& val) {
      if (statement>=landmarks.size()) landmarks.resize(statement+1);
      landmarks[statement].addLandmark(ts,val);
    };

    // Clear all landmarks. This function must be called after each
    // widening since the widened state does not preserve the
    // distances between the redundant inequalities and the feasible
    // space.
    inline void clear() {
      for (std::vector<Landmarks>::iterator iter = landmarks.begin();
	   iter != landmarks.end(); iter++) iter->clear();
    };

    inline void promote() {
      for (std::vector<Landmarks>::iterator iter = landmarks.begin();
	   iter != landmarks.end(); iter++) iter->promote();
    };

    // Calculate the number of iterations until any of the redundant
    // inequalities would become non-redundant. Returns true if an
    // statement was found were extrapolation is possible.
    bool calcNoOfIterations(WideningInfo& li);

    // Print the table.
    friend std::ostream& operator<<(std::ostream& s, const LandmarkTable& lm);

  };

}

#endif // __LANDMARK_H
