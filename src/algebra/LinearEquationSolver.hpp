/*=========================================================================

 Program: FEMuS
 Module: LinearEquationSolver
 Authors: Eugenio Aulisa, Simone Bnà, Giorgio Bornia

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_algebra_LinearEquationSolver_hpp__
#define __femus_algebra_LinearEquationSolver_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include <vector>
#include <memory>
#include <cstdio>

#include "FemusConfig.hpp"
#include "SolverPackageEnum.hpp"
#include "PrecondtypeEnum.hpp"
#include "SolvertypeEnum.hpp"
#include "LinearEquation.hpp"
#include "LinearEquationSolverEnum.hpp"
#include "MgTypeEnum.hpp"
#include "FieldSplitTree.hpp"

#include <petscksp.h>

namespace femus {


  using std::vector;

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
  class SparseMatrix;
  class NumericVector;
  class Preconditioner;

  /**
  * This class provides a uniform interface for linear solvers.  This base
  * class is overloaded to provide linear solvers from different packages
  * like PETSC or LASPACK.
  */

  class LinearEquationSolver : public LinearEquation {

    public:

      /**  Constructor. Initializes Solver data structure */
      LinearEquationSolver(const unsigned &igrid, Solution *other_solution);

      /** Destructor. */
      virtual ~LinearEquationSolver();

      /** Release all memory and clear data structures. */
      virtual void Clear() {}

      /** Builds a \p LinearEquationSolver using the linear solver in \p solver_package */
      static std::unique_ptr<LinearEquationSolver> build(const unsigned &igrid, Solution* other_solution,
          const LinearEquationSolverType & smoother_type, const SolverPackage solver_package = LSOLVER);

      /** Set the tolerance for the solver */
      virtual void SetTolerances(const double& rtol, const double& atol,
                                 const double& divtol, const unsigned& maxits,
                                 const unsigned& restart) = 0;

      virtual KSP* GetKSP() {
        std::cout << "Warning GetKSP() is not available for this smoother\n";
        abort();
      }

      virtual void MGInit(const MgSmootherType & mg_smoother_type, const unsigned &levelMax, const SolverType & mgSolverType) {
        std::cout << "Warning InitMG(...) is not available for this smoother\n";
        abort();
      }

      virtual void MGClear() {
        std::cout << "Warning ClearMG() is not available for this smoother\n";
        abort();
      };

      virtual void MGSetLevel(LinearEquationSolver *LinSolver, const unsigned &levelMax,
                              const vector <unsigned> &variable_to_be_solved,
                              SparseMatrix* PP, SparseMatrix* RR,
                              const unsigned &npre, const unsigned &npost
                              ) = 0;

      virtual void MGSolve(const bool ksp_clean) = 0;
      
      virtual void SetRichardsonScaleFactor(const double & richardsonScaleFactor) = 0; 

      /** Sets the type of solver to use. */
      void set_solver_type (const SolverType st)  {
        _levelSolverType = st;
      }

      /** Sets the type of preconditioner to use. */
      void set_preconditioner_type(const PreconditionerType pct);

      /** Attaches a Preconditioner object to be used */
      void attach_preconditioner(Preconditioner * preconditioner);

      /** @returns true if the data structures are */
      bool initialized() const {
        return _is_initialized;
      }

      /** Returns the type of solver to use. */
      SolverType solver_type() const {
        return _levelSolverType;
      }

      /** Returns the type of preconditioner to use. */
      PreconditionerType preconditioner_type() const;

      void PrintSolverInfo(const bool & printInfo) {
        _printSolverInfo = printInfo;
      }

      /** Set the number of elements of the Vanka Block */
      virtual void SetElementBlockNumber(const unsigned & block_elemet_number) {
        std::cout << "Warning SetElementBlockNumber(const unsigned &) is not available for this smoother\n";
      };
      virtual void SetElementBlockNumberFluid(const unsigned & block_elemet_number, const unsigned &overlap) {
        std::cout << "Warning SetElementBlockNumberFluid(const unsigned &) is not available for this smoother\n";
      };
      virtual void SetElementBlockNumberSolid(const unsigned & block_elemet_number, const unsigned &overlap) {
        std::cout << "Warning SetElementBlockNumberSolid(const unsigned &) is not available for this smoother\n";
      };
      virtual void SetElementBlockNumberPorous(const unsigned & block_elemet_number, const unsigned &overlap) {
        std::cout << "Warning SetElementBlockNumberPorous(const unsigned &) is not available for this smoother\n";
      };

      virtual void SetFieldSplitTree(FieldSplitTree * fieldSplitTree) {
        std::cout << "SetFieldSplitTree(const FieldSpliTreeStructure & fieldSplitTree) is not available for this smoother\n";
      };

      /** To be Added */
      virtual void SetElementBlockNumber(const char all[], const unsigned & overlap = 1) {
        std::cout << "Warning SetElementBlockNumber(const char [], const unsigned & ) is not available for this smoother\n";
      };

      /** To be Added */
      virtual void SetNumberOfSchurVariables(const unsigned short & NSchurVar) {
        std::cout << "Warning SetNumberOfSchurVariables(const unsigned short &) is not available for this smoother\n";
      };

      /** Call the smoother-solver using the PetscLibrary. */
      virtual void Solve(const vector <unsigned> &VariableTobeSolved, const bool &ksp_clean) = 0;


      /** @deprecated Old solver with algebra objects passed as arguments TODO think of removing */
      virtual std::pair<unsigned int, double> solve(SparseMatrix&,   // System Matrix
          SparseMatrix&,  // prec
          NumericVector&, // Solution vector
          NumericVector&, // RHS vector
          const double,      // Stopping tolerance
          const unsigned int) {
        std::cout << "If I call this it's wrong" << std::endl;
        abort();
      }

    protected:

      /** Enum stating which type of iterative solver to use. */
      SolverType _levelSolverType;
      SolverType _mgSolverType;

      /** Enum statitng with type of preconditioner to use. */
      PreconditionerType _preconditioner_type;

      /** Holds the Preconditioner object to be used for the linear solves. */
      Preconditioner *_preconditioner;

      /** Flag indicating if the data structures have been initialized. */
      bool _is_initialized;

      /// Boolean flag to indicate whether we want to use an identical preconditioner to the previous solve.
      bool same_preconditioner;

      bool _printSolverInfo;

  };

  /**
   * -------------------- inline functions ---------------------
   */

  inline LinearEquationSolver::LinearEquationSolver(const unsigned &igrid, Solution *other_solution) :
    LinearEquation(other_solution),
    _levelSolverType(GMRES),
    _mgSolverType(GMRES),
    _preconditioner(NULL),
    _is_initialized(false),
    same_preconditioner(false) {

    if(igrid == 0) {
      _preconditioner_type = LU_PRECOND;
    }
    else {
      _preconditioner_type = ILU_PRECOND;
    }
  }

  inline LinearEquationSolver::~LinearEquationSolver() {
    this->Clear();
  }


} //end namespace femus



#endif
