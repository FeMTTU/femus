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
#include "petscmat.h"
#include "FElemTypeEnum.hpp"
#include "ParallelObject.hpp"
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

  using std::vector;

  class Solution : public ParallelObject {

    public:

      /** Constructor */
      Solution(Mesh *other_msh);

      /** Destructor */
      ~Solution();

      /** Add a new variable called 'name' */
      void AddSolution(const char name[], const FEFamily fefamily, const FEOrder order, const unsigned& tmorder = 0, const bool &Pde_type = 1);

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

      void UpdateSol(const vector <unsigned> &_SolPdeIndex,  NumericVector* EPS, const vector <vector <unsigned> > &KKoffset);
      /** */
      void UpdateRes(const vector <unsigned> &_SolPdeIndex, NumericVector* _RES, const vector <vector <unsigned> > &KKoffset);

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
      bool FlagAMRRegionBasedOnErroNorm(const vector <unsigned> &solIndex, std::vector <double> &AMRthreshold, const unsigned& normType);
      bool FlagAMRRegionBasedOnErroNormAdaptive(const vector <unsigned> &solIndex, std::vector <double> &AMRthreshold, 
						const unsigned& normType, const double &neighborThresholdValue);


      /** Build Grad Matrix structure for SolType 0,1,2 */
      void BuildGradMatrixStructure(unsigned SolType);

      /** Init and set to zero The AMR Eps vector */
      void InitAMREps();
      /** member data - one for each variable - */
      vector <NumericVector*> _Sol;
      vector <NumericVector*> _SolOld;
      vector <NumericVector*> _Res;
      vector <NumericVector*> _Eps;
      vector <NumericVector*> _AMREps;
      bool _AMR_flag;
      vector <NumericVector*> _Bdc;
      vector <bool> _ResEpsBdcFlag;

      vector < vector <NumericVector*> > _GradVec;

      vector <SparseMatrix*> _GradMat[5];
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
      
    private:
        
      //member data
        
      /** This group of vectors has the size of the number of added solutions */
      vector <int>      _SolType;
      vector <char*>    _SolName;
      vector <unsigned> _SolTmOrder;
      vector <FEFamily> _family;
      vector <FEOrder>  _order;
      vector <bool>     _removeNullSpace;
      
      Mesh *_msh;
      
      bool _FSI;

  };


} //end namespace femus



#endif
