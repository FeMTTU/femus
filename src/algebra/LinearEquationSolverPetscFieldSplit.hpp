/*=========================================================================

  Program: FEMUS
  Module: LinearEquationSolverPetscFieldSplit
  Authors: Eugenio Aulisa, Guoyi Ke

  Copyright (c) FEMTTU
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/

#ifndef __femus_algebra_LinearEquationSolverPetscFieldSplit_hpp__
#define __femus_algebra_LinearEquationSolverPetscFieldSplit_hpp__

#include "FemusConfig.hpp"

#ifdef HAVE_PETSC

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "LinearEquationSolverPetsc.hpp"
#include "FieldSplitTree.hpp"

namespace femus {

  /**
   * This class inherits the class LinearEquationSolverPetsc. In this class the solver is implemented using the PETSc package
   **/

  class LinearEquationSolverPetscFieldSplit : public LinearEquationSolverPetsc {

    public:

      /**  Constructor. Initializes Petsc data structures */
      LinearEquationSolverPetscFieldSplit(const unsigned& igrid, Solution *other_solution);

      /** Destructor */
      ~LinearEquationSolverPetscFieldSplit() {};

    private:

      void SetFieldSplitTree(FieldSplitTree* fieldSplitTree);

      /** To be Added */
      void BuildBdcIndex(const vector <unsigned>& variable_to_be_solved);

      void SetPreconditioner(KSP& subksp, PC& subpc);

      void FielSlipTreeIsNotDefined();
      // member data
      FieldSplitTree* _fieldSplitTree;

  };

// =================================================

  inline LinearEquationSolverPetscFieldSplit::LinearEquationSolverPetscFieldSplit(const unsigned& igrid, Solution *other_solution)
    : LinearEquationSolverPetsc(igrid, other_solution),
    _fieldSplitTree(NULL){
  }

} //end namespace femus


#endif
#endif
