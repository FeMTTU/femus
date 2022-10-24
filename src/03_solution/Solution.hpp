/*=========================================================================

 Program: FEMuS
 Module: Solution
 Authors: Simone Bnà

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/**
 * This class is useful to handle multilevel variables. It's going to be split
 * into meshvariables and solutionvariables
*/


#ifndef __femus_solution_Solution_hpp__
#define __femus_solution_Solution_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "ParallelObject.hpp"
#include "Mesh.hpp"
#include "FElemTypeEnum.hpp"
#include "Math.hpp"

#include "petscmat.h"

#include <vector>


namespace femus {

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
  class NumericVector;
  class SparseMatrix;
  class Mesh;
  class MultiLevelSolution;


  class Solution : public ParallelObject {


// === CONSTR/DESTR  - BEGIN =================
    public:
      /* Constructor */
      Solution(Mesh *other_msh);

      /** Destructor */
      ~Solution();
// === CONSTR/DESTR  - END =================

// === BASIC SOL MANIPULATION - BEGIN =================
    public:
      /** Add a new variable called 'name' */
      void AddSolution(const char name[], const FEFamily fefamily, const FEOrder order, const unsigned& tmorder = 0, const bool &Pde_type = 1);

      void AddSolution(const char name[], const FEFamily fefamily, const FEOrder order_v, const FEOrder order_b,
                             const unsigned& tmorder = 0, const bool &Pde_type = 1);

      /** Add a number of solutions */
      void AddSolution_par(const int n_sols, const char name[], const FEFamily fefamily, const FEOrder order, const unsigned& tmorder = 0, const bool &Pde_type = 1);
      
      void ResizeSolution_par(const int new_size);
  
      /** Resize the solution vector */
      void ResizeSolutionVector(const char name[]);

      /** Free the solution vectors */
      void FreeSolutionVectors();
// === BASIC SOL MANIPULATION - END =================

// === NAME & INDEX  - BEGIN =================
  public:

        /** Get the index of the variable -name- */
      unsigned GetIndex(const char name[]) const;
      
    private:
        
      /** Vector size: number of added Solutions. Solution name */
      std::vector <char*>    _SolName;
// === NAME & INDEX  - END =================

// === MESH - BEGIN =================
    public:
        
      Mesh * GetMesh() {
        return _msh;
      }
      
    private:
      /** Pointer to underlying mesh object */
      Mesh *_msh;
      
// === MESH - END =================

// === SPACE DISCRETIZATION (FE) - BEGIN =================
    public:
      
      unsigned GetSolutionType(const unsigned &index ) {
	return _SolType[index];
      }
      
     const unsigned GetSolutionType(const unsigned &index ) const {
	return _SolType[index];
      }
      
    /// compute the sequential index for the FE family 
    static int compute_fe_sol_type( const FEFamily fefamily,  const FEOrder order_v, const FEOrder order_b) {
        return   order_v - ((fefamily == LAGRANGE) ? 1 : 0) + fefamily * 3;
    }
      
    private:
                
      /** Vector size: number of added Solutions. Tells the FE index */
      std::vector <int>      _SolType;
      /** Vector size: number of added Solutions. FE family */
      std::vector <FEFamily> _family;
      /** Vector size: number of added Solutions. FE order within a family */
      std::vector <FEOrder>  _order;
// === SPACE DISCRETIZATION (FE) - END =================
      
// === Solution, NumericVector - BEGIN =================
    public:
      /** Get a const solution (Numeric Vector) by name @todo make the _Sol object private */
      const NumericVector& GetSolutionName(const char* var) const {
        return *_Sol[GetIndex(var)];
      }

      /** Get a solution (Numeric Vector) by name @todo make the _Sol object private */
      NumericVector& GetSolutionName(const char* var) {
        return *_Sol[GetIndex(var)];
      }
      
      /** Add another solution */
      void add_solution(const unsigned index_read, const unsigned index_write);
      


    public:
      /** Vector size: number of added Solutions. */
      std::vector <NumericVector*> _Sol;
// === Solution, NumericVector - END =================

// === TIME EVOLUTION (NOT DISCRETIZATION)  - BEGIN =================
// (why here? or is it only the fact of changing in time or not? or maybe it is just SOLUTION with ITERATION (so you have stages, like a time-dep algorithm or a nonlinear algoritm)
  public:
      /** Update the solution */
      void CopySolutionToOldSolution();
      
      void ResetSolutionToOldSolution();

      unsigned GetSolutionTimeOrder(unsigned i) {
        return _SolTmOrder[i];
      }
      
      /** Vector size: number of added Solutions. Used only if the Solution is time-dependent */
      std::vector <NumericVector*> _SolOld;

    private:
      /** Vector size: number of added Solutions. Time type of Solution: 0 = steady, 2 = time dependent */
      std::vector <unsigned> _SolTmOrder;
// === TIME EVOLUTION (NOT DISCRETIZATION)  - END =================

// === INITIALIZATION ANALYTICAL FUNCTION - BEGIN =================
public:
    
  void set_analytical_function(const char * name,  Math::Function< double > * func_in);

   Math::Function< double > * get_analytical_function(const char * name) const;
  
protected:
    
   std::vector< Math::Function< double > * >  _analytical_function;
// === INITIALIZATION ANALYTICAL FUNCTION - END =================
    
      
//=========== 
///==== @todo from now on, it is more stuff that shouldn't be in a basic Solution object      
//=========== 
      
// === Solution as Unknown of System - BEGIN =================
    public:
//       /** Sum to Solution vector the Epsilon vector. It is used inside the multigrid cycle */
//       void UpdateSolAndRes(const vector <unsigned> &_SolPdeIndex,  NumericVector* EPS, NumericVector* RES, const vector <vector <unsigned> > &KKoffset);

      void UpdateSol(const std::vector <unsigned> &_SolPdeIndex,  NumericVector* EPS, const std::vector <std::vector <unsigned> > &KKoffset);
      /** */
      void UpdateRes(const std::vector <unsigned> &_SolPdeIndex, NumericVector* _RES, const std::vector <std::vector <unsigned> > &KKoffset);


      /** Vector size: number of added Solutions. Used only if the Solution is an unknown to some PDE */
      std::vector <NumericVector*> _Res;
      /** Vector size: number of added Solutions. Used only if the Solution is an unknown to some PDE */
      std::vector <NumericVector*> _Eps;
      /** Vector size: number of added Solutions. Used only if the Solution is an unknown to some PDE */
      std::vector <NumericVector*> _Bdc;
      /** Vector size: number of added Solutions. Tells if the Solution is an unknown to some PDE */
      std::vector <bool> _ResEpsBdcFlag;
// === Solution as Unknown of System - END =================


// === NULL SPACE (of what here? for pressure variable, when you pinpoint it) - BEGIN =================
    public:
        
      void RemoveNullSpace(const unsigned &index) {
        _removeNullSpace[index] = true ;
      }
      bool GetIfRemoveNullSpace(const unsigned &index) {
        return _removeNullSpace[index];
      }
      
    private:
        
      /** Vector size: number of added Solutions. */
      std::vector <bool>     _removeNullSpace;
      
// === NULL SPACE (of what here? for pressure variable, when you pinpoint it) - END =================


// === FSI - BEGIN =================
    public:

    void SetIfFSI(const bool &FSI = true){
	_FSI = FSI; 
      }
      
      bool GetIfFSI(){
	return _FSI; 
      }
      
    private:
      bool _FSI;

// === FSI - END =================


// === MESH, AMR - BEGIN =================
// does this depend on a System or can it be also equation-independent?
  public:
//      /** Flag the elemets to be refined in the AMR alghorithm based on the epsilon*/
//     bool FlagAMRRegionBasedOnl2(const vector <unsigned> &_SolPdeIndex, const double &AMRthreshold);
//
//      /** Flag the elemets to be refined in the AMR alghorithm based on the solution gradient*/
//     bool FlagAMRRegionBasedOnSemiNorm(const vector <unsigned> &SolIndex,const double &AMRthreshold);

      /** Flag the elemets to be refined in the AMR alghorithm based on the solution gradient*/
      bool FlagAMRRegionBasedOnErroNorm(const std::vector <unsigned> &solIndex, std::vector <double> &AMRthreshold, const unsigned& normType);
      bool FlagAMRRegionBasedOnErroNormAdaptive(const std::vector <unsigned> &solIndex, std::vector <double> &AMRthreshold, 
						const unsigned& normType, const double &neighborThresholdValue);
      /** Init and set to zero The AMR Eps vector */
      void InitAMREps();
      
      /** Vector size: number of added Solutions. */
      std::vector <NumericVector*> _AMREps;
      bool _AMR_flag;
      
      
      
    public:
      std::vector < std::vector <NumericVector*> > _GradVec;

      /** Build Grad Matrix structure for SolType 0,1,2 */
      void BuildGradMatrixStructure(unsigned SolType);

      std::vector <SparseMatrix*> _GradMat[5];
      // bool _GradMatFlag[5];

// === MESH, AMR - END =================

  };


} //end namespace femus



#endif
