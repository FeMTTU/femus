#ifndef __linear_solverM_h__
#define __linear_solverM_h__

#include <memory>
#include <cstdio>

#include "FEMTTUConfig.h"
#include "SolverPackageEnum.hpp"
#include "PrecondtypeEnum.hpp"
#include "SolvertypeEnum.hpp"
#include "LinSysPde.hpp"

#include "vector"
using std::vector;

// forward declarations
template <typename T> class AutoPtr;  //
class SparseMatrix;
class NumericVector;
template <typename T> class ShellMatrix;  //
class Preconditioner;

// ================================================
// This class provides a uniform interface for linear solvers.  This base
// class is overloaded to provide linear solvers from different packages
// like PETSC or LASPACK.
// ===============================================
class LinearSolverM : public lsysPde {
  // ================================
  // DATA
  // ==============================
protected:
  /// Enum stating which type of iterative solver to use.
  SolverType _solver_type;

  /// Enum statitng with type of preconditioner to use.
  PreconditionerType _preconditioner_type;

  /// Holds the Preconditioner object to be used for the linear solves.
  Preconditioner *_preconditioner;

  /// Flag indicating if the data structures have been initialized.
  bool _is_initialized;

  /// The number of elements in a Vanka Block
  //unsigned _num_elem_vanka_block;

public:

  /// Boolean flag to indicate whether we want to use an identical preconditioner to the previous solve.
  bool same_preconditioner;

  // =================================
  // CONSTR/DESTR
  // ================================
  ///  Constructor. Initializes Solver data structure
  LinearSolverM (const unsigned &igrid, mesh* other_msh);
  /// Destructor.
  virtual ~LinearSolverM();

  /// Release all memory and clear data structures.
  virtual void clear() {}

  /// Initialize data structures if not done so already.
  virtual void init() = 0;

  /// Builds a \p LinearSolverM using the linear solver in \p solver_package
  
  static std::auto_ptr<LinearSolverM > build(const unsigned &igrid, mesh* other_msh,const SolverPackage solver_package =LSOLVER);
  
  // =================================
  // SETTING FUNCTIONS
  // ================================
  ///Set the tolerance for the solver
  virtual void set_tolerances(const double &rtol, const double &atol,
                              const double &divtol, const unsigned &maxits, const unsigned &index) = 0;

  /// Sets the type of solver to use.
  void set_solver_type (const SolverType st)  {
    _solver_type = st;
  }

  /// Sets the type of preconditioner to use.
  void set_preconditioner_type (const PreconditionerType pct);

  /// Attaches a Preconditioner object to be used
  void attach_preconditioner(Preconditioner * preconditioner);

  // =================================
  // RETURN FUNCTIONS
  // ================================

  /// @returns true if the data structures are
  bool initialized () const {
    return _is_initialized;
  }

  /// Returns the type of solver to use.
  SolverType solver_type () const {
    return _solver_type;
  }

  /// Returns the type of preconditioner to use.
  PreconditionerType preconditioner_type () const;

  // =================================
  // OUR virtual FUNCTIONS
  // ================================

  /// Set the number of elements of the Vanka Block
  virtual void set_num_elem_vanka_block(const unsigned num_elem_vanka_block)=0;
  
  /// Call the Vanka(Schur) smoother-solver using the PetscLibrary.
  virtual std::pair< int, double> solve(const vector <unsigned> &VankaIndex, const short unsigned &NSchurVar,const bool &Schur) = 0;
  /// Call the Gmres smoother-solver
  virtual std::pair< int, double> solve() = 0;
  
//   /// This function calls the solver
//   virtual std::pair< int, double> solve (SparseMatrix&,  // System Matrix
//                                          NumericVector&, // Solution vector
//                                          NumericVector&, // RHS vector
//                                          const double,      // Stopping tolerance
//                                          const  int) = 0; // N. Iterations
// 
//   ///   *This function calls the solver with preconditioner
//   virtual std::pair< int, double> solve (SparseMatrix&,  // System Matrix
//                                          SparseMatrix&,  // Preconditioning Matrix
//                                          NumericVector&, // Solution vector
//                                          NumericVector&, // RHS vector
//                                          const double,      // Stopping tolerance
//                                          const  int) = 0; // N. Iteration

  /// Prints a useful message about why the latest linear solve con(di)verged.
  virtual void print_converged_reason() = 0;
};

// -------------------- inline functions ---------

// =============================================
inline LinearSolverM::LinearSolverM(const unsigned &igrid, mesh* other_msh) :
  lsysPde(other_msh),
  _solver_type(GMRES),
  _preconditioner(NULL),
  _is_initialized(false),
  same_preconditioner(false) {

  if(igrid==0){
    _preconditioner_type=LU_PRECOND;
  }
  else{
    _preconditioner_type=ILU_PRECOND;
  }
}

inline LinearSolverM::~LinearSolverM() {
  this->clear();
}

#endif // #ifdef __solver_h__
