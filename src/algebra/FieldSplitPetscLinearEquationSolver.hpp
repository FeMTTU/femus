/*=========================================================================

  Program: FEMUS
  Module: FieldSplitPetscLinearEquationSolver
  Authors: Eugenio Aulisa, Guoyi Ke

  Copyright (c) FEMTTU
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/

#ifndef __femus_algebra_FieldSplitPetscLinearEquationSolver_hpp__
#define __femus_algebra_FieldSplitPetscLinearEquationSolver_hpp__

#include "FemusConfig.hpp"

#ifdef HAVE_PETSC

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "GmresPetscLinearEquationSolver.hpp"
#include "FieldSplitTree.hpp"

namespace femus {

  /**
   * This class inherits the class GmresPetscLinearEquationSolver. In this class the solver is implemented using the PETSc package
   **/

  class FieldSplitPetscLinearEquationSolver : public GmresPetscLinearEquationSolver {

    public:

      /**  Constructor. Initializes Petsc data structures */
      FieldSplitPetscLinearEquationSolver(const unsigned& igrid, Solution *other_solution);

      /** Destructor */
      ~FieldSplitPetscLinearEquationSolver() {};

    private:

      void SetFieldSplitTree(FieldSplitTree* fieldSplitTree);

      /** To be Added */
      void BuildBdcIndex(const vector <unsigned>& variable_to_be_solved);

      void SetPreconditioner(KSP& subksp, PC& subpc);

      // member data
      FieldSplitTree* _fieldSplitTree;

  };

// =================================================

  inline FieldSplitPetscLinearEquationSolver::FieldSplitPetscLinearEquationSolver(const unsigned& igrid, Solution *other_solution)
    : GmresPetscLinearEquationSolver(igrid, other_solution) {
  }

} //end namespace femus


#endif
#endif
