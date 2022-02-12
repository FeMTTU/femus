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

#include "SchurFactTypeEnum.hpp"
#include "SchurPreType.hpp"

#include "Mesh.hpp"

namespace femus {

  class LinearEquationSolverPetscFieldSplit;

  class FieldSplitTree {
    public:
      //single split constructor
      FieldSplitTree (const SolverType& solver, const PreconditionerType& preconditioner, const std::vector < unsigned >& fields, std::string name);

      FieldSplitTree (const SolverType& solver, const PreconditionerType& preconditioner, const std::vector < unsigned >& fields, const std::vector < unsigned >& solutionType, std::string name);

      void FieldSplitTreeBuild (const SolverType& solver, const PreconditionerType& preconditioner, const std::vector < unsigned >& fields, std::string name);
      //multiple split constructor
      FieldSplitTree (const SolverType& solver, const PreconditionerType& preconditioner, std::vector < FieldSplitTree*> childBranch, std::string name);

      ~FieldSplitTree();

      void PrintFieldSplitTree (const unsigned& counter = 0);

      void BuildIndexSet (const std::vector< std::vector < unsigned > >& KKoffset, const unsigned& iproc,
                          const unsigned& nprocs, const unsigned& level, const LinearEquationSolverPetscFieldSplit *solver);

      void BuildASMIndexSet (const unsigned& level, const LinearEquationSolverPetscFieldSplit *solver);

      void SetPC (KSP& ksp, const unsigned& level) ;

      void SetSolver (KSP &ksp, const SolverType &solver);      

      void SetTolerances (const double& rtol, const double& abstol, const double& dtol, 
                               const unsigned& maxits, const unsigned &restart = 30);
      void SetSchurFactorizationType (const SchurFactType& schurFactType);
      void SetSchurPreType (const SchurPreType& schurPreType);

      void SetRichardsonScaleFactor (const double & richardsonScaleFactor) {
        _richardsonScaleFactor = richardsonScaleFactor;
      }

      const unsigned& GetNumberOfSplits() {
        return _numberOfSplits;
      }

      std::string GetName() const {
        return _name;
      }

      const std::vector < unsigned >& GetFieldsInSplit (const unsigned& i) {
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

      FieldSplitTree* GetChild (const unsigned& i);

    private:

      //void SetPetscSolverType (KSP& ksp);
      void SetSchurFactorizationType (PC &pc);
      void SetSchurPreType (PC &pc);

      SolverType _solver;
      PreconditionerType _preconditioner;
      unsigned _numberOfSplits;
      FieldSplitTree* _father;
      std::vector < FieldSplitTree* > _child;
      std::vector < std::vector < unsigned > > _fieldsSplit;
      std::vector < unsigned > _fieldsAll;
      std::vector < unsigned > _solutionType;
      std::string _name;
      std::vector < PetscInt* > _isSplitIndexPt;
      std::vector < std::vector < IS > > _isSplit;
      double _rtol;
      double _abstol;
      double _dtol;
      unsigned _maxits;
      unsigned _restart;
      double _richardsonScaleFactor;
      
      SchurFactType _schurFactType;
      SchurPreType _schurPreType;

      std::vector < std::vector< std::vector < unsigned > > >_MatrixOffset;


      //for ASM pourposes
    public:

      void SetAsmStandard (const bool &standard) {
        _asmStandard = standard;
        if (standard) _asmOverlapping = 1;
        else _asmOverlapping = 0;
      };

      void SetAsmBlockSize (const unsigned &BlockSize) {
        _asmBlockSize[0] = BlockSize;
        _asmBlockSize[1] = BlockSize;
        _asmStandard = false;
        _asmOverlapping = 0;
      };

      void SetAsmBlockSizeSolid (const unsigned &BlockSize) {
        _asmBlockSize[0] = BlockSize;
        _asmStandard = false;
        _asmOverlapping = 0;
      };

      void SetAsmBlockSizeFluid (const unsigned &BlockSize) {
        _asmBlockSize[1] = BlockSize;
        _asmStandard = false;
        _asmOverlapping = 0;
      };

      void SetAsmBlockPreconditioner (const PreconditionerType &preconditioner) {
        _asmBlockPreconditioner[0] = preconditioner;
        _asmBlockPreconditioner[1] = preconditioner;
      };

      void SetAsmBlockPreconditionerSolid (const PreconditionerType &preconditioner) {
        _asmBlockPreconditioner[0] = preconditioner;
      };

      void SetAsmBlockPreconditionerFluid (const PreconditionerType &preconditioner) {
        _asmBlockPreconditioner[1] = preconditioner;
      };

      void SetAsmNumeberOfSchurVariables (const unsigned &SchurVariableNumber) {
        _asmSchurVariableNumber = SchurVariableNumber;
      };
      void SetAsmOverlapping (const unsigned &overlapping) {
        _asmOverlapping = overlapping;
      }
    private:
      std::vector< std::vector < std::vector <PetscInt> > > _asmOverlappingIsIndex;
      std::vector< std::vector < std::vector <PetscInt> > > _asmLocalIsIndex;
      std::vector< std::vector <IS> > _asmOverlappingIs;
      std::vector< std::vector <IS> > _asmLocalIs;
      std::vector< std::vector <unsigned> > _asmBlockMaterialRange;
      std::vector < unsigned > _asmBlockSize;
      std::vector < PreconditionerType > _asmBlockPreconditioner;
      unsigned _asmSchurVariableNumber;

      //std::vector< PetscInt >  _nlocal;
      bool _asmStandard;
      unsigned _asmOverlapping;





  };


}


#endif
