#ifndef __linear_solverM_h__
#define __linear_solverM_h__

#include <memory>

#include "FemusExtLib_conf.hpp"

#include "Typedefs_conf.hpp" 

#include "SolverPackageEnum.hpp"      // #include "enum_solver_package.h"
#include "Precondtype_enum.hpp"        // #include "enum_preconditioner_type.h"
#include "Solvertype_enum.hpp"         // #include "enum_solver_type.h"



// forward declarations
template <typename T> class AutoPtr;
class SparseMatrixM;
class NumericVectorM;
template <typename T> class ShellMatrix;
 class PreconditionerM;
 

// ================================================
// This class provides a uniform interface for linear solvers.  This base
// class is overloaded to provide linear solvers from different packages
// like PETSC.
// ===============================================
class LinearSolverM
#ifdef LM_REFCOUNT
: public ReferenceCountedObject<LinearSolverM > 
#endif
{
  // ================================
  // DATA
  // ==============================
protected:
  /// Enum stating which type of iterative solver to use.
  SolverTypeM _solver_type;

  /// Enum statitng with type of preconditioner to use.
  PreconditionerTypeM _preconditioner_type;
  /// Holds the Preconditioner object to be used for the linear solves.
  PreconditionerM * _preconditioner;
  
   /// Flag indicating if the data structures have been initialized.
  bool _is_initialized;
  
public:
  
   /// Boolean flag to indicate whether we want to use an identical preconditioner to the previous solve. 
  bool same_preconditioner;
  
  // =================================
  // CONSTR/DESTR
  // ================================
  ///  Constructor. Initializes Solver data structure
  LinearSolverM ();
    
  /// Destructor.
  virtual ~LinearSolverM();
  /// Release all memory and clear data structures.
  virtual void clear() {}

  /// Initialize data structures if not done so already.
  virtual void init() = 0;
  /// Builds a \p LinearSolverM using the linear solver in \p solver_package
  static std::auto_ptr<LinearSolverM > build(const SolverPackage solver_package = LSOLVER);

  // =================================
  // SETTING FUNCTIONS 
  // ================================
  /// Sets the type of solver to use.
  void set_solver_type (const SolverTypeM st)  { _solver_type = st; }
  
  /// Sets the type of preconditioner to use.
  void set_preconditioner_type (const PreconditionerTypeM pct);
  /// Attaches a Preconditioner object to be used
  void attach_preconditioner(PreconditionerM * preconditioner);
  
  // =================================
  // RETURN FUNCTIONS 
  // ================================
  
  /// @returns true if the data structures are
  bool initialized () const { return _is_initialized; }
  
  /// Returns the type of solver to use.
  SolverTypeM solver_type () const { return _solver_type; }

  /// Returns the type of preconditioner to use.
  PreconditionerTypeM preconditioner_type () const;

  // =================================
  // SOLVE FUNCTIONS 
  // ================================
  
//   //sandro
//   /// This function calls the solver
//   virtual std::pair<unsigned int, Real> solve(
//     SparseMatrixM&/*Sys Mtrx*/,NumericVectorM&/*sol_vec*/, NumericVectorM&/*rhs*/,
//     const double /*tol*/,const unsigned int/*NIter*/) = 0;
// 
//   ///   *This function calls the solver with preconditioner
//   virtual std::pair<unsigned int, Real> solve(
//     SparseMatrixM&/*Sys Mtrx*/,SparseMatrixM&/*Prec Mtrx*/,
//     NumericVectorM&/*sol_vec*/,NumericVectorM&/*rhs*/,
//     const double /*tol*/,const unsigned int/*NIter*/) = 0; 
  
 
  /// This function calls the solver
  virtual std::pair<unsigned int, Real> solve (SparseMatrixM&,  // System Matrix
					       NumericVectorM&, // Solution vector
					       NumericVectorM&, // RHS vector
					       const double,      // Stopping tolerance
					       const unsigned int) = 0; // N. Iterations
  
  ///   *This function calls the solver with preconditioner 
  virtual std::pair<unsigned int, Real> solve (SparseMatrixM&,  // System Matrix
					       SparseMatrixM&,  // Preconditioning Matrix
					       NumericVectorM&, // Solution vector
					       NumericVectorM&, // RHS vector
					       const double,      // Stopping tolerance
					       const unsigned int) = 0; // N. Iteration
 /// This function calls the solver	       
  /// Prints a useful message about why the latest linear solve con(di)verged.
  virtual void print_converged_reason() = 0;
};

// -------------------- inline functions ---------

// =============================================
inline LinearSolverM::LinearSolverM() :
    _solver_type(GMRESM),
    _preconditioner_type(ILU_PRECONDM),
    _preconditioner(NULL),
    _is_initialized(false),
    same_preconditioner(false) {}
// ========================================================
inline LinearSolverM::~LinearSolverM() {  this->clear();}

#endif // #ifdef __solver_h__
