/*=========================================================================

 Program: FEMuS
 Module: LinearEquation
 Authors: Eugenio Aulisa, Simone Bnà

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_algebra_LinearEquation_hpp__
#define __femus_algebra_LinearEquation_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Mesh.hpp"
#include "petscmat.h"
#include "ParallelObject.hpp"


namespace femus {

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class elem_type;
class NumericVector;
class SparseMatrix;
class Mesh;
class Solution;

/**
 This class is a container that holds linear operators and other structure for solving a linear equation system
*/

class LinearEquation : public ParallelObject {

public:

  /** costructor */
  LinearEquation(Solution *other_solution);

  /** destructor */
  ~LinearEquation();

  /** To be Added */
  void InitPde(const vector <unsigned> &_SolPdeIndex,const  vector <int> &SolType,
               const vector <char*> &SolName, vector <NumericVector*> *Bdc_other,
               const unsigned &other_gridn, vector < bool > &SparsityPattern_other);

  void GetSparsityPatternSize();

  /** To be Added */
  void DeletePde();

  /** To be Added */
 // unsigned GetSystemDof(const unsigned &index_sol, const unsigned &kkindex_sol,const unsigned &idof_gmt) const;

  unsigned GetSystemDof(const unsigned &index_sol, const unsigned &kkindex_sol,
			const unsigned &i, const unsigned &iel) const;

  unsigned GetSystemDof(const unsigned &index_sol, const unsigned &kkindex_sol,
	                const unsigned &ielc, const unsigned &i0,const unsigned &i1,
		        const Mesh* mshc) const;
			
  unsigned GetSystemDof(const unsigned &soltype, const unsigned &kkindex_sol,
			const unsigned &i, const unsigned &iel, const vector < vector <unsigned> > &otherKKoffset) const;
			

  /** To be Added */
  void SetResZero();

  /** To be Added */
  void SetEpsZero();

  /** To be Added */
  void SumEpsCToEps();

  /** To be Added */
  void UpdateResidual();

  /** AddLevel */
  void AddLevel();

  void SwapMatrices(){
    SparseMatrix *Temp = _KK;
    _KK = _KKamr;
    _KKamr = Temp;
  }
  
  void SetNumberOfGlobalVariables(const unsigned &numberOfGlobalVariables){
    _numberOfGlobalVariables = numberOfGlobalVariables;
  }
  
  // member data
  Mesh *_msh;
  Solution *_solution;
  NumericVector *_EPS, *_EPSC, *_RES, *_RESC;
  SparseMatrix *_KK;
  SparseMatrix *_KKamr;
  
  
  vector < vector <unsigned> > KKoffset;     // size [_SolPdeIndex.size() + 1][_nprocs]
  vector < unsigned > KKghostsize;           // size [_nprocs]
  vector < vector < int> > KKghost_nd;       // size [_nprocs][KKghostsize[i]]
  vector <int> KKIndex;                      // size [_SolPdeIndex.size() + 1]
  unsigned _gridn;  
  vector < int > d_nnz;
  vector < int > o_nnz;
  
  void SetSparsityPatternMinimumSize (const std::vector < unsigned> &minimumSize, const std::vector < unsigned > &variableIndex);

protected:

  /** To be Added */
  unsigned GetIndex(const char name[]);

  // member data
  vector <unsigned> _SolPdeIndex;
  vector <int> _SolType;
  vector <char*> _SolName;
  const vector <NumericVector*> *_Bdc;
  vector <bool> _SparsityPattern;
  std::vector < unsigned > _sparsityPatternMinimumSize; 
  std::vector <unsigned> _sparsityPatternVariableIndex;
  unsigned _numberOfGlobalVariables;

};

} //end namespace femus



#endif
