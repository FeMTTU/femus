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

#ifndef __femus_enums_FieldSplitTree_hpp__
#define __femus_enums_FieldSplitTree_hpp__


#include <cstdlib>
#include <iostream>

#include <vector>
#include <map>
#include <string>

#include "SolvertypeEnum.hpp"
#include "PrecondtypeEnum.hpp"
#include "PetscVector.hpp"
#include "PetscPreconditioner.hpp"


namespace femus {
  class FieldSplitTree {
    public:
      //single split constructor
      FieldSplitTree( const SolverType& solver, const PreconditionerType& preconditioner, const std::vector < unsigned >& fields, std::string name );

      //multiple split constructor
      FieldSplitTree( const SolverType& solver, const PreconditionerType& preconditioner, std::vector < FieldSplitTree*> childBranch, std::string name );

      ~FieldSplitTree();

      void PrintFieldSplitTree( const unsigned& counter = 0 );

      void BuildIndexSet( const std::vector< std::vector < unsigned > >& KKoffset, const unsigned& iproc, const unsigned& nprocs, const unsigned& level );

      void SetPC( KSP& ksp, const unsigned& level) ; 

      const unsigned& GetNumberOfSplits() {
        return _numberOfSplits;
      }

      std::string GetName() const {
        return _name;
      }

      const std::vector < unsigned >& GetFieldsInSplit( const unsigned& i ) {
        return _fieldsSplit[i];
      }

      const std::vector < unsigned >& GetAllFields() {
        return _fieldsAll;
      }

      const SolverType& GetSolver() {
        return _solver;
      }

      const PreconditionerType& GetPreconditioner() {
        return _preconditioner;
      }

      FieldSplitTree* GetFather() const;

      FieldSplitTree* GetChild( const unsigned& i );
      
    private:

      void SetPetscSolverType(KSP& ksp);
      
      SolverType _solver;
      PreconditionerType _preconditioner;
      unsigned _numberOfSplits;
      FieldSplitTree* _father;
      std::vector < FieldSplitTree* > _child;
      std::vector < std::vector < unsigned > > _fieldsSplit;
      std::vector < unsigned > _fieldsAll;
      std::string _name;
      std::vector < PetscInt* > _isSplitIndexPt;
      std::vector < std::vector < IS > > _isSplit;
      double _rtol;
      double _abstol;
      double _dtol;
      unsigned _maxits;

  };


}


#endif
