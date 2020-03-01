/*=========================================================================

 Program: FEMUS
 Module: LinearImplicitSystem
 Authors: Simone Bnà

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_equations_LinearImplicitSystem_hpp__
#define __femus_equations_LinearImplicitSystem_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "ImplicitSystem.hpp"
#include "LinearEquationSolver.hpp"
#include "MgTypeEnum.hpp"
#include "DirichletBCTypeEnum.hpp"
#include "LinearEquationSolverEnum.hpp"
#include "FemusDefault.hpp"

#include <petscksp.h>

namespace femus {


//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------

  class LinearImplicitSystem : public ImplicitSystem {

    public:

      /** Constructor.  Optionally initializes required data structures. */
      LinearImplicitSystem (MultiLevelProblem& ml_probl, const std::string& name, const unsigned int number, const LinearEquationSolverType & smoother_type);

      /** Destructor */
      virtual ~LinearImplicitSystem();

      /** The type of the parent. */
      typedef ImplicitSystem Parent;

      /** Init the system PDE structures */
      virtual void init();

      /** @deprecated Init the system PDE structures */
      virtual void init_two();

      /** @deprecated Multigrid routine */
      void MGSolve (double Eps, int MaxIter, const uint Gamma = DEFAULT_MG_GAMMA, const uint Nc_pre = DEFAULT_NC_PRE, const uint Nc_coarse = DEFAULT_NC_COARSE, const uint Nc_post = DEFAULT_NC_POST);

      /** @deprecated Multigrid step routine */
      double MGStep (int Level, double Eps1, int MaxIter, const uint Gamma, const uint Nc_pre, const uint Nc_coarse, const uint Nc_post);


      /** Add a system level */
      void AddSystemLevel();
      /**
       * @returns \p "LinearImplicit".  Helps in identifying
       * the system type in an equation system file.
       */
      virtual std::string system_type() const {
        return "LinearImplicit";
      }

      /**
      * The \p LinearSolver defines the default interface used to
      * solve the linear_implicit system.  This class handles all the
      * details of interfacing with various linear algebra packages
      * like PETSc or LASPACK. Up to now also for the nonlinear case we use linear_solvers, in future we will add the nonlinear solver
      */
      vector < LinearEquationSolver*> _LinSolver;

      void SetNumberOfGlobalVariables (const unsigned &numberOfGlobalVariables) {
        _numberOfGlobalVariables = numberOfGlobalVariables;
      }

      /** Set the max number of linear iterationsfor solving Ax=b */
      void SetMaxNumberOfLinearIterations (unsigned int max_lin_it) {
        _n_max_linear_iterations = max_lin_it;
      };

      /** Get the final Linear Residual of the linear problem Ax=b*/
      double GetFinalLinearResidual() const {
        return _final_linear_residual;
      };

      /** Flag to print fields to file after each linear iteration */
      void SetDebugLinear(const bool my_value); 

      /** Get the absolute convergence tolerance for the linear problem Ax=b*/
      double GetAbsoluteConvergenceTolerance() const {
        return _linearAbsoluteConvergenceTolerance;
      };

      /** Set the absolute convergence tolerance for the linear problem Ax=b*/
      void SetAbsoluteLinearConvergenceTolerance (const double & absolute_convergence_tolerance) {
        _linearAbsoluteConvergenceTolerance = absolute_convergence_tolerance;
      };

      /** */
      bool IsLinearConverged (const unsigned igridn);

      /** Set the type of multigrid */
      void SetMgType (const MgType mgtype) {
        _mg_type = mgtype;
      };

      /** Set the modality of handling the BC boundary condition (penalty or elimination)*/
      void SetDirichletBCsHandling (const DirichletBCType DirichletMode);

      /** Add the variable solname to the variable set to be solved by using the Vanka smoother */
      void AddVariableToBeSolved (const char solname[]);

      /** Clear the Vanka index */
      void ClearVariablesToBeSolved();

      /** Set the multigrid smoother, gmres or Vanka (in future AMS (additive schwartz preconditioner))*/
      void SetLinearEquationSolverType (const LinearEquationSolverType LinearEquationSolverType, const CoarseLevelInclude &includeCoarseLevel = INCLUDE_COARSE_LEVEL_FALSE);      /** Set the PCFIELDSPLIT structure in linear solver */
      void SetFieldSplitTree (FieldSplitTree *fieldSplitTree);

      /** Set if the solver has to output convergence information **/
      void PrintSolverInfo (const bool & printInfo = true);

      /** Set the number of elements of a Vanka block. The formula is nelem = (2^dim)^dim_vanka_block */
      void SetElementBlockNumber (unsigned const &dim_vanka_block);

      /** Set the number of elements of a Vanka block. The formula is nelem = (2^dim)^dim_vanka_block */
      void SetElementBlockNumber (const char all[], const unsigned & overlap = 1);

      /** Set the Ksp smoother solver on the fine grids. At the coarse solver we always use the LU (Mumps) direct solver */
      void SetSolverCoarseGrid (const SolverType &solvertype);

      /** Set the preconditioner for the Ksp smoother solver on the fine grids */
      void SetPreconditionerCoarseGrid (const PreconditionerType &preconditioner_type);


      /** Set the Ksp smoother solver on the fine grids. At the coarse solver we always use the LU (Mumps) direct solver */
      void SetSolverFineGrids (const SolverType &solvertype);

      /** Set the preconditioner for the Ksp smoother solver on the fine grids */
      void SetPreconditionerFineGrids (const PreconditionerType &preconditioner_type);

      /** Set the tolerances for the ksp solver on fine grids: rtol, atol, divtol, maxits */
      void SetTolerances (const double &rtol, const double &atol,
                          const double &divtol, const unsigned &maxits,
                          const unsigned &restart = 30);

      void SetRichardsonScaleFactor (const double &richardsonScaleFactor) {
        _richardsonScaleFactor = richardsonScaleFactor;
        _richardsonScaleFactorDecrease = 0;
        _richardsonScaleFactorIsSet = true;
        for (unsigned i = 0; i < _gridn; i++) {
          _LinSolver[i]->SetRichardsonScaleFactor (_richardsonScaleFactor);
        }
      }

      void SetRichardsonScaleFactor (const double &richardsonScaleFactorMin, const double &richardsonScaleFactorMax) {
        _richardsonScaleFactor = richardsonScaleFactorMax;
        _richardsonScaleFactorDecrease = (_gridn > 1) ? (richardsonScaleFactorMin - richardsonScaleFactorMax) / (_gridn - 2) : 0;
        _richardsonScaleFactorIsSet = true;
        _LinSolver[0]->SetRichardsonScaleFactor (_richardsonScaleFactor);
        for (unsigned i = 1; i < _gridn; i++) {
          _LinSolver[i]->SetRichardsonScaleFactor (_richardsonScaleFactor + _richardsonScaleFactorDecrease * (i - 1));
        }
      }

      void ResetComputationalTime() {
        _totalAssemblyTime = 0.;
        _totalSolverTime = 0.;
      }

      void PrintComputationalTime() {
        std::cout << "Total Assembly Time = " << _totalAssemblyTime << std::endl;
        std::cout << "Total Solver Time = " << _totalSolverTime << std::endl;
        std::cout << "Total Computational Time = " << _totalAssemblyTime + _totalSolverTime << std::endl;
      }

      void SetOuterSolver (const SolverType & mgOuterSolver) {
        _mgOuterSolver = mgOuterSolver;
      };

      /** Set AMR options */
      void SetAMRSetOptions (const std::string& AMR, const unsigned &AMRlevels,
                             const std::string& AMRnorm, const double &AMRthreshold,
                             bool (* SetRefinementFlag) (const std::vector < double > &x,
                                                         const int &ElemGroupNumber, const int &level) = NULL);

      void SetAMRghborThresholdValue (const double &neighborThresholdValue) {
        _AMReighborThresholdValue = neighborThresholdValue;
      }

      /** Set the options of the Schur-Vanka smoother */
      //void SetVankaSchurOptions(bool Schur, short unsigned NSchurVar);
      void SetNumberOfSchurVariables (const unsigned short &NSchurVar);


      /** Set the number of pre-smoothing step of a Multigrid cycle */
      void SetNumberPreSmoothingStep (const unsigned &npre) {
        _npre = npre;
      };
      
      void SetNumberSmoothingStepCoarseGrid (const unsigned &npre0) {
        _npre0 = npre0;
      };

      /** Set the number of post-smoothing step of a Multigrid cycle */
      void SetNumberPostSmoothingStep (const unsigned &npost) {
        _npost = npost;
      };
      
      void GetSystemInfo();

      /** enforce sparcity pattern for setting uncoupled variables and save on memory allocation **/
      void SetSparsityPattern (vector < bool > other_sparcity_pattern);

      void SetSparsityPatternMinimumSize (const unsigned & minimumSize);

      bool GetAssembleMatrix() {
        return _assembleMatrix;
      }

      vector < SparseMatrix* > &GetProjectionMatrix() {
        return _PP;
      }

      vector < SparseMatrix* > &GetRestrictionMatrix() {
        return _RR;
      }
     
      /** Solves the system. */
      virtual void MGsolve (const MgSmootherType& mgSmootherType = MULTIPLICATIVE);
    protected:

      vector < SparseMatrix* > _PP, _RR;
      vector < SparseMatrix* > _PPamr, _RRamr;

      bool _printSolverInfo;
      bool _assembleMatrix;
      void AddAMRLevel (unsigned &AMRCounter);

      bool Vcycle (const unsigned & gridn, const MgSmootherType& mgSmootherType);


      /** Create the Prolongator matrix for the Multigrid solver */
      void Prolongator (const unsigned &gridf);

      /** Create the Restrictor matrix for the Multigrid solver */
      virtual void Restrictor (const unsigned &gridf);

      /** Prolongate the solution to a finer level */
      void ProlongatorSol (unsigned gridf);

      /** Create the Prolongator Operator in order to get the coarser matrix for the Algebraic Multigrid Solver */
      virtual void BuildProlongatorMatrix (unsigned gridf);
      virtual void BuildAmrProlongatorMatrix (unsigned level);
      void ZeroInterpolatorDirichletNodes (const unsigned &level);

      // member data
      /** The number of linear iterations required to solve the linear system Ax=b. */
      //unsigned int _n_linear_iterations;

      /** Flag for printing fields at each linear iteration */
      bool _debug_linear;

      /** The final residual for the linear system Ax=b. */
      double _final_linear_residual;

      /** The threshold residual for the linear system Ax=b. */
      double _linearAbsoluteConvergenceTolerance;

      /** The max number of linear iterations */
      unsigned int _n_max_linear_iterations;

      /** The ksp outer solver*/
      SolverType _mgOuterSolver;

      /** The type of multigrid, F-cyle, V-cycle, M-cycle */
      MgType _mg_type;

      /** To be Added */
      unsigned _npre;
      unsigned _npre0;
      unsigned _npost;

      /** To be Added */
      LinearEquationSolverType _smootherType;
      CoarseLevelInclude _includeCoarseLevelSmoother;
      bool _MGmatrixFineReuse;
      bool _MGmatrixCoarseReuse;

      /** To be Added */
      vector <unsigned> _VariablesToBeSolvedIndex;

      SolverType _finegridsolvertype;
      unsigned int _DirichletBCsHandlingMode;
      double _rtol, _atol, _divtol, _maxits, _restart;
      PreconditionerType _finegridpreconditioner;

      bool _numblock_test;
      unsigned _num_block;

      bool _numblock_all_test;
      bool _overlap;

      bool _NSchurVar_test;
      unsigned short _NSchurVar;
      bool _AMRtest;
      unsigned _maxAMRlevels;
      short _AMRnorm;
      double _AMReighborThresholdValue;
      std::vector <double> _AMRthreshold;

      vector <bool> _SparsityPattern;

      double _richardsonScaleFactor;
      double _richardsonScaleFactorDecrease;
      bool _richardsonScaleFactorIsSet;

      double _totalSolverTime;
      double _totalAssemblyTime;

      bool _bitFlipOccurred;
      unsigned _bitFlipCounter;

      unsigned _numberOfGlobalVariables;
      unsigned _sparsityPatternMinimumSize;

  };

} //end namespace femus

#endif
