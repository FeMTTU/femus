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
#include "ParallelObject.hpp"


#include "petscmat.h"


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
  void InitPde(const std::vector <unsigned> &_SolPdeIndex,const  std::vector <int> &SolType,
               const std::vector <char*> &SolName, std::vector <NumericVector*> *Bdc_other,
               const unsigned &other_gridn, std::vector < bool > &SparsityPattern_other);

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
			const unsigned &i, const unsigned &iel, const std::vector < std::vector <unsigned> > &otherKKoffset) const;

            
  /** Sparsity pattern: print number of nonzeros per row */
  void sparsity_pattern_print_nonzeros(const std::string filename_base, const std::string on_or_off);

  /** Print numeric vector with structure */
  void print_with_structure_matlab_friendly(const unsigned iproc, const std::string filename_base, NumericVector * num_vec_in) const;

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
  
  
  /** Pointer to underlying mesh */
  Mesh *_msh;
  
  /** To be Added */
  Solution *_solution;
  
  /** Vectors for error and residual */
  NumericVector *_EPS, *_EPSC, *_RES, *_RESC;
  
  /** Matrix */
  SparseMatrix *_KK;
  
  /** AMR Matrix */
  SparseMatrix *_KKamr;
  
  
  /** size [_SolPdeIndex.size() + 1][_nprocs] */
  std::vector < std::vector <unsigned> > KKoffset;
  
  /** size [_nprocs] */
  std::vector < unsigned > KKghostsize;
  
  /** size [_nprocs][KKghostsize[i]] */
  std::vector < std::vector < int> > KKghost_nd;
  
  /** size [_SolPdeIndex.size() + 1] : number of dofs for each variable, summed over all processors (expressed in offset mode) */
  std::vector < int > KKIndex;
  
  /** number of levels */
  unsigned _gridn;
  
  /** number of non-zeros per row, on-diagonal */
  std::vector < int > d_nnz;
  
  /** number of non-zeros per row, off-diagonal */
  std::vector < int > o_nnz;
  
  void SetSparsityPatternMinimumSize (const std::vector < unsigned> &minimumSize, const std::vector < unsigned > &variableIndex);

protected:

  /** To be Added */
  unsigned GetIndex(const char name[]);

  // member data
  /** size: number of unknowns */
  std::vector <unsigned> _SolPdeIndex;
  /** size: number of unknowns */
  std::vector <int> _SolType;
  /** size: number of unknowns */
  std::vector <char*> _SolName;
  /** size: number of unknowns */
  const std::vector <NumericVector*> *_Bdc;
  /** size: number of unknowns */
  std::vector <bool> _SparsityPattern;
  /** size: number of unknowns */
  std::vector < unsigned > _sparsityPatternMinimumSize; 
  /** size: number of unknowns */
  std::vector <unsigned> _sparsityPatternVariableIndex;
  
  unsigned _numberOfGlobalVariables;

};

} //end namespace femus



#endif
