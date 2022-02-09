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
#include "Mesh.hpp"
#include "FElemTypeEnum.hpp"
#include "ParallelObject.hpp"

#include "petscmat.h"

#include <vector>


namespace femus {

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
  class elem_type;
  class NumericVector;
  class SparseMatrix;
  class Mesh;
  class MultiLevelSolution;


  class Solution : public ParallelObject {

    public:

      /* Constructor */
      Solution(Mesh *other_msh);

      /** Destructor */
      ~Solution();

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

      unsigned GetSolutionTimeOrder(unsigned i) {
        return _SolTmOrder[i];
      };
      
      
//       /** Sum to Solution vector the Epsilon vector. It is used inside the multigrid cycle */
//       void UpdateSolAndRes(const vector <unsigned> &_SolPdeIndex,  NumericVector* EPS, NumericVector* RES, const vector <vector <unsigned> > &KKoffset);

      void UpdateSol(const std::vector <unsigned> &_SolPdeIndex,  NumericVector* EPS, const std::vector <std::vector <unsigned> > &KKoffset);
      /** */
      void UpdateRes(const std::vector <unsigned> &_SolPdeIndex, NumericVector* _RES, const std::vector <std::vector <unsigned> > &KKoffset);

      /** Update the solution */
      void CopySolutionToOldSolution();
      
      void ResetSolutionToOldSolution();

      /** Get a const solution (Numeric Vector) by name @todo make the _Sol object private */
      const NumericVector& GetSolutionName(const char* var) const {
        return *_Sol[GetIndex(var)];
      };

      /** Get a solution (Numeric Vector) by name @todo make the _Sol object private */
      NumericVector& GetSolutionName(const char* var) {
        return *_Sol[GetIndex(var)];
      };

//      /** Flag the elemets to be refined in the AMR alghorithm based on the epsilon*/
//     bool FlagAMRRegionBasedOnl2(const vector <unsigned> &_SolPdeIndex, const double &AMRthreshold);
//
//      /** Flag the elemets to be refined in the AMR alghorithm based on the solution gradient*/
//     bool FlagAMRRegionBasedOnSemiNorm(const vector <unsigned> &SolIndex,const double &AMRthreshold);

      /** Flag the elemets to be refined in the AMR alghorithm based on the solution gradient*/
      bool FlagAMRRegionBasedOnErroNorm(const std::vector <unsigned> &solIndex, std::vector <double> &AMRthreshold, const unsigned& normType);
      bool FlagAMRRegionBasedOnErroNormAdaptive(const std::vector <unsigned> &solIndex, std::vector <double> &AMRthreshold, 
						const unsigned& normType, const double &neighborThresholdValue);


      /** Build Grad Matrix structure for SolType 0,1,2 */
      void BuildGradMatrixStructure(unsigned SolType);

      /** Init and set to zero The AMR Eps vector */
      void InitAMREps();
      /** Vector size: number of added Solutions. */
      std::vector <NumericVector*> _Sol;
      /** Vector size: number of added Solutions. Used only if the Solution is time-dependent */
      std::vector <NumericVector*> _SolOld;
      /** Vector size: number of added Solutions. Used only if the Solution is an unknown to some PDE */
      std::vector <NumericVector*> _Res;
      /** Vector size: number of added Solutions. Used only if the Solution is an unknown to some PDE */
      std::vector <NumericVector*> _Eps;
      /** Vector size: number of added Solutions. */
      std::vector <NumericVector*> _AMREps;
      bool _AMR_flag;
      /** Vector size: number of added Solutions. Used only if the Solution is an unknown to some PDE */
      std::vector <NumericVector*> _Bdc;
      /** Vector size: number of added Solutions. Tells if the Solution is an unknown to some PDE */
      std::vector <bool> _ResEpsBdcFlag;

      std::vector < std::vector <NumericVector*> > _GradVec;

      std::vector <SparseMatrix*> _GradMat[5];
      // bool _GradMatFlag[5];

      void RemoveNullSpace(const unsigned &index) {
        _removeNullSpace[index] = true ;
      };
      bool GetIfRemoveNullSpace(const unsigned &index) {
        return _removeNullSpace[index];
      };

      Mesh * GetMesh() {
        return _msh;
      }
      
      /** Get the index of the variable -name- */
      unsigned GetIndex(const char name[]) const;
      
      
      unsigned GetSolutionType(const unsigned &index ) {
	return _SolType[index];
      }
      
      void SetIfFSI(const bool &FSI = true){
	_FSI = FSI; 
      }
      
      bool GetIfFSI(){
	return _FSI; 
      }
    
    /// compute the sequential index for the FE family 
    static int compute_fe_sol_type( const FEFamily fefamily,  const FEOrder order_v, const FEOrder order_b) {
        return   order_v - ((fefamily == LAGRANGE) ? 1 : 0) + fefamily * 3;
    }
      
    private:
        
      //member data
        
      /** Vector size: number of added Solutions. Tells the FE index */
      std::vector <int>      _SolType;
      /** Vector size: number of added Solutions. Solution name */
      std::vector <char*>    _SolName;
      /** Vector size: number of added Solutions. Time type of Solution: 0 = steady, 2 = time dependent */
      std::vector <unsigned> _SolTmOrder;
      /** Vector size: number of added Solutions. FE family */
      std::vector <FEFamily> _family;
      /** Vector size: number of added Solutions. FE order within a family */
      std::vector <FEOrder>  _order;
      /** Vector size: number of added Solutions. */
      std::vector <bool>     _removeNullSpace;
      
      /** Pointer to underlying mesh object */
      Mesh *_msh;
      
      bool _FSI;

  };


} //end namespace femus



#endif
