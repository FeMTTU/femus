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
#include "System.hpp"
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

  class LinearImplicitSystem : public System {


//==== Constr/Destr - BEGIN ========
    public:

      /** Constructor.  Optionally initializes required data structures. */
      LinearImplicitSystem (MultiLevelProblem& ml_probl, const std::string& name, const unsigned int number, const LinearEquationSolverType & smoother_type);

      /** Destructor */
      virtual ~LinearImplicitSystem();
//==== Constr/Destr - END ========

      
//==== Basic - BEGIN ========
    public:

      /**
       * @returns \p "LinearImplicit".  Helps in identifying
       * the system type in an equation system file.
       */
      virtual std::string system_type() const {
        return "LinearImplicit";
      }
      
      /** The type of the parent. */
      typedef System Parent;

      /** Init the system PDE structures */
      virtual void init();
      
      void GetSystemInfo();
//==== Basic - END ========


//==== Matrix, Assemble - BEGIN ========
    public:

      bool GetAssembleMatrix() {
        return _assembleMatrix;
      }

      
    protected:
      
      bool _assembleMatrix;
      
//==== Matrix, Assemble - END ========


//==== Matrix, Unknown types & Sparsity Pattern - BEGIN ========
    public:

      void SetNumberOfGlobalVariables (const unsigned &numberOfGlobalVariables) {
        _numberOfGlobalVariables = numberOfGlobalVariables;
      }

      /** enforce sparcity pattern for setting uncoupled variables and save on memory allocation **/
      void SetSparsityPattern (std::vector < bool > other_sparcity_pattern);

      void SetSparsityPatternMinimumSize (const unsigned &minimumSize, const std::string variableName = "All");

    protected:
      
      unsigned _numberOfGlobalVariables;
      
      std::vector <bool> _SparsityPattern;

      std::vector< std::string > _sparsityPatternSolName;
      std::vector < unsigned > _sparsityPatternMinimumSize;

//==== Matrix, Unknown types & Sparsity Pattern - END ========


//==== Matrix, Boundary Conditions - BEGIN ========
    public:

      /** Set the modality of handling the BC boundary condition (penalty or elimination)*/
      void SetDirichletBCsHandling (const DirichletBCType DirichletMode);
      
    protected:

      unsigned int _DirichletBCsHandlingMode;   ///@todo it is not used anywhere
      
      void ZeroInterpolatorDirichletNodes (const unsigned &level);

//==== Matrix, Boundary Conditions - END ========


      
//==== Unknowns - BEGIN ========
    public:
      
      /** Add the variable solname to the variable set to be solved by using the Vanka smoother  ? Is it only for Vanka or for all?  */
      void AddVariableToBeSolved (const char solname[]);

      /** Clear the Vanka index ? Is it only for Vanka or for all? */
      void ClearVariablesToBeSolved();

    protected:
      
      /** To be Added */
      std::vector <unsigned> _VariablesToBeSolvedIndex;

      
//==== Unknowns - END ========

      

//==== Solver - BEGIN ========
    public:

      /** @deprecated Multigrid routine */
      void MGSolve (double Eps, int MaxIter, const uint Gamma = DEFAULT_MG_GAMMA, const uint Nc_pre = DEFAULT_NC_PRE, const uint Nc_coarse = DEFAULT_NC_COARSE, const uint Nc_post = DEFAULT_NC_POST);

      /** @deprecated Multigrid step routine */
      double MGStep (int Level, double Eps1, int MaxIter, const uint Gamma, const uint Nc_pre, const uint Nc_coarse, const uint Nc_post);

      /** Solves the system. */
      virtual void MGsolve (const MgSmootherType& mgSmootherType = MULTIPLICATIVE);

      /** Set the max number of linear iterationsfor solving Ax=b */
      void SetMaxNumberOfLinearIterations (unsigned int max_lin_it) {
        _n_max_linear_iterations = max_lin_it;
      };


      /** Get the absolute convergence tolerance for the linear problem Ax=b*/
      double GetAbsoluteLinearConvergenceTolerance() const {
        return _linearAbsoluteConvergenceTolerance;
      };

      /** Set the absolute convergence tolerance for the linear problem Ax=b*/
      void SetAbsoluteLinearConvergenceTolerance (const double & absolute_convergence_tolerance) {
        _linearAbsoluteConvergenceTolerance = absolute_convergence_tolerance;
      };
      

    protected:
      
      bool _bitFlipOccurred;
      unsigned _bitFlipCounter;

      bool _MGmatrixFineReuse;
      bool _MGmatrixCoarseReuse;

      /** Convergence control */
      bool HasLinearConverged (const unsigned igridn);
      
      /** The threshold residual for the linear system Ax=b. */
      double _linearAbsoluteConvergenceTolerance;

      /** The max number of linear iterations */
      unsigned int _n_max_linear_iterations;

      
//==== Solver - END ========


//==== Solver, Debug - BEGIN ========
    public:
      
      void ResetComputationalTime() {
        _totalAssemblyTime = 0.;
        _totalSolverTime = 0.;
      }

      void PrintComputationalTime() {
        std::cout << "Total Assembly Time = " << _totalAssemblyTime << std::endl;
        std::cout << "Total Solver Time = " << _totalSolverTime << std::endl;
        std::cout << "Total Computational Time = " << _totalAssemblyTime + _totalSolverTime << std::endl;
      }

      /** Set if the solver has to output convergence information **/
      void PrintSolverInfo (const bool & printInfo = true);

      /** Flag to print fields to file after each linear iteration */
      void SetDebugLinear(const bool my_value); 

    protected:

      bool _printSolverInfo;
      
      /** Flag for printing fields at each linear iteration */
      bool _debug_linear;
      
      double _totalSolverTime;
      double _totalAssemblyTime;

//==== Solver, Debug - END ========



//==== Solver, Multigrid - BEGIN ========
    public:

      /**
      * The \p LinearSolver defines the default interface used to
      * solve the linear_implicit system.  This class handles all the
      * details of interfacing with various linear algebra packages
      * like PETSc or LASPACK. Up to now also for the nonlinear case we use linear_solvers, in future we will add the nonlinear solver
      */
      std::vector < LinearEquationSolver*> _LinSolver;   /* vector of number of levels */

      /** Set the tolerances for the ksp solver on each grid: rtol, atol, divtol, maxits */
      void SetTolerances (const double &rtol, const double &atol,
                          const double &divtol, const unsigned &maxits,
                          const unsigned &restart = 30);
      
      void SetOuterSolver (const SolverType & mgOuterSolver) {
        _mgOuterSolver = mgOuterSolver;
      };

      /** Set the type of multigrid */
      void SetMgType (const MgType mgtype) {
        _mg_type = mgtype;
      };

      /** Set the multigrid smoother, gmres or Vanka (in future AMS (additive schwartz preconditioner))*/
      void SetLinearEquationSolverType (const LinearEquationSolverType LinearEquationSolverType, const CoarseLevelInclude &includeCoarseLevel = INCLUDE_COARSE_LEVEL_FALSE);


    protected:


      double _rtol, _atol, _divtol, _maxits, _restart;
      
      /** The ksp outer solver*/
      SolverType _mgOuterSolver;

      /** The type of multigrid, F-cyle, V-cycle, M-cycle */
      MgType _mg_type;
      
      
      bool Vcycle (const unsigned & gridn, const MgSmootherType& mgSmootherType);

//==== Solver, Multigrid - END ========


//==== Solver, Multigrid, Prol/Rest Operators - BEGIN ========
    public:
      
      std::vector < SparseMatrix* > &GetProjectionMatrix() {
        return _PP;
      }

      std::vector < SparseMatrix* > &GetRestrictionMatrix() {
        return _RR;
      }
     
     
    protected:
      
      /** Create the Prolongator matrix for the Multigrid solver */
      void Prolongator (const unsigned &gridf);

      /** Create the Restrictor matrix for the Multigrid solver */
      virtual void Restrictor (const unsigned &gridf);

      /** Prolongate the solution to a finer level */
      void ProlongatorSol (unsigned gridf);

      /** Create the Prolongator Operator in order to get the coarser matrix for the Algebraic Multigrid Solver */
      virtual void BuildProlongatorMatrix (unsigned gridf);
      
      /** Prolongator and Restrictor */
      std::vector < SparseMatrix* > _PP, _RR;
      
//==== Solver, Multigrid, Prol/Rest Operators - END ========


//==== Solver, Multigrid, Smoothing - BEGIN ========
    public:
      
      /** Set the Ksp smoother solver on the fine grids. At the coarse solver we always use the LU (Mumps) direct solver */
      void SetSolverCoarseGrid (const SolverType &solvertype);

      /** Set the preconditioner for the Ksp smoother solver on the fine grids */
      void SetPreconditionerCoarseGrid (const PreconditionerType &preconditioner_type);


      /** Set the Ksp smoother solver on the fine grids. At the coarse solver we always use the LU (Mumps) direct solver */
      void SetSolverFineGrids (const SolverType &solvertype);

      /** Set the preconditioner for the Ksp smoother solver on the fine grids */
      void SetPreconditionerFineGrids (const PreconditionerType &preconditioner_type);

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
      

      
    protected:
      
      SolverType _finegridsolvertype;

      PreconditionerType _finegridpreconditioner;

      /** To be Added */
      LinearEquationSolverType _smootherType;
      CoarseLevelInclude _includeCoarseLevelSmoother;
      
      /** To be Added */
      unsigned _npre;
      unsigned _npre0;
      unsigned _npost;

      
//==== Solver, Multigrid, Smoothing - END ========
      
      
//==== Solver, Multigrid, AMR - BEGIN ========
    public:
      
      /** Add a system level, for AMR */
      void AddSystemLevel();

      /** Set AMR options */
      void SetAMRSetOptions (const std::string& AMR, const unsigned &AMRlevels,
                             const std::string& AMRnorm, const double &AMRthreshold,
                             bool (* SetRefinementFlag) (const std::vector < double > &x,
                                                         const int &ElemGroupNumber, const int &level) = NULL);

      void SetAMRghborThresholdValue (const double &neighborThresholdValue) {
        _AMReighborThresholdValue = neighborThresholdValue;
      }
      
      
    protected:
      
      void AddAMRLevel (unsigned &AMRCounter);
      

      bool _AMRtest;
      unsigned _maxAMRlevels;
      short _AMRnorm;
      double _AMReighborThresholdValue;
      std::vector <double> _AMRthreshold;

//==== Solver, Multigrid, AMR - END ========

      
      
//==== Solver, Multigrid, AMR, Prol/Rest Operators - BEGIN ========
    public:
      


    protected:
      
      /** Prolongator and Restrictor for AMR */
      std::vector < SparseMatrix* > _PPamr, _RRamr;
      
      virtual void BuildAmrProlongatorMatrix (unsigned level);

//==== Solver, Multigrid, AMR, Prol/Rest Operators - END ========



      
//==== Solver, Preconditioner - BEGIN ========
    public:

      /** Set the PCFIELDSPLIT structure in linear solver */
      void SetFieldSplitTree (FieldSplitTree *fieldSplitTree);

      
      
//==== Solver, Preconditioner - END ========


      
//==== ASM - BEGIN ========
    public:
      

      /** Set the number of elements of a Vanka block. The formula is nelem = (2^dim)^dim_vanka_block */
      void SetElementBlockNumber (unsigned const &dim_vanka_block);

      /** Set the number of elements of a Vanka block. The formula is nelem = (2^dim)^dim_vanka_block */
      void SetElementBlockNumber (const char all[], const unsigned & overlap = 1);

      /** Set the options of the Schur-Vanka smoother */
      //void SetVankaSchurOptions(bool Schur, short unsigned NSchurVar);
      void SetNumberOfSchurVariables (const unsigned short &NSchurVar);
      
    protected:
      
      
      bool _numblock_test;
      unsigned _num_block;

      bool _numblock_all_test;
      bool _overlap;

      bool _NSchurVar_test;
      unsigned short _NSchurVar;
      
//==== ASM - END ========


//==== Richardson - BEGIN ========
    public:
      
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
      
    protected:

      double _richardsonScaleFactor;
      double _richardsonScaleFactorDecrease;
      bool _richardsonScaleFactorIsSet;
      
//==== Richardson - END ========
      

  };

} //end namespace femus

#endif
