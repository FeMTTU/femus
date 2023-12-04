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
#include "ParallelObject.hpp"

#include <vector>
#include <string>



namespace femus {

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class NumericVector;
class SparseMatrix;

class Mesh;
class Solution;

/**
 This class is a container that holds linear operators and other structure for solving a linear equation system
*/

class LinearEquation : public ParallelObject {


//==== Constructors / Destructor - BEGIN ========
 public:
  /** costructor */
  LinearEquation(Solution *other_solution);

  /** destructor */
  ~LinearEquation();
//==== Constructors / Destructor - END ========

  
//==== All equation - BEGIN ========
 public:

  /** To be Added */
  void InitPde(const std::vector <unsigned> &_SolPdeIndex,
               const std::vector <int> &SolType,
               const std::vector <char*> &SolName,
               std::vector <NumericVector*> *Bdc_other,
               const unsigned &other_gridn,
               std::vector < bool > &SparsityPattern_other);

  /** To be Added */
  void DeletePde();
//==== All equation - END ========


//==== Number of levels - BEGIN ========
 public:
 
  /** AddLevel */
  void AddLevel();
  
 protected:

  /** number of levels */
  unsigned _gridn;
  
//==== Number of levels - END ========


//==== DOFMAP - BEGIN ========
 public:
   
  /** To be Added */
 // unsigned GetSystemDof(const unsigned &index_sol, const unsigned &kkindex_sol,const unsigned &idof_gmt) const;

  unsigned GetSystemDof(const unsigned &index_sol, const unsigned &kkindex_sol,
			const unsigned &i, const unsigned &iel) const;

  unsigned GetSystemDof(const unsigned &index_sol, const unsigned &kkindex_sol,
	                const unsigned &ielc, const unsigned &i0,const unsigned &i1,
		        const Mesh* mshc) const;
			
  unsigned GetSystemDof(const unsigned &soltype, const unsigned &kkindex_sol,
			const unsigned &i, const unsigned &iel, const std::vector < std::vector <unsigned> > &otherKKoffset) const;

 public:
   
   inline const Mesh *  GetMeshFromLinEq() const { return _msh; }
   
 private:
  
  /** Pointer to underlying mesh */
  const Mesh *_msh;
  
//==== DOFMAP - END ========

  
//==== Linear Algebra - BEGIN ============

//==== Indices (both for Vectors and Matrices) - BEGIN ========
 
 public:
  
  /** size [_SolPdeIndex.size() + 1][_nprocs] */
  std::vector < std::vector <unsigned> > KKoffset;
  
  /** size [_SolPdeIndex.size() + 1] : number of dofs for each variable, summed over all processors (expressed in offset mode) */
  std::vector < int > KKIndex;
  
 private:
    
  /** size [_nprocs] */
  std::vector < unsigned > KKghostsize;
  
  /** size [_nprocs][KKghostsize[i]] */
  std::vector < std::vector < int> > KKghost_nd;
  
//==== Indices (both for Vectors and Matrices) - END ========
  
  
  
//==== MATRIX - BEGIN ========
 public:
            
  void SwapMatrices() {
    SparseMatrix *Temp = _KK;
    _KK = _KKamr;
    _KKamr = Temp;
  }
  
  /** Matrix */
  SparseMatrix *_KK;
  
  /** AMR Matrix */
  SparseMatrix *_KKamr;
  
 protected:

//==== MATRIX - END ========
            
//==== MATRIX, Sparsity Pattern - BEGIN ========
 public:
   
  void GetSparsityPatternSize();

  void SetSparsityPatternMinimumSize (const std::vector < unsigned> &minimumSize, const std::vector < unsigned > &variableIndex);
  
  /** Sparsity pattern: print number of nonzeros per row */
  void sparsity_pattern_print_nonzeros(const std::string filename_base, const std::string on_or_off);
  
 protected:
   
  /** size: number of unknowns */
  std::vector <bool> _SparsityPattern;
  /** size: number of unknowns */
  std::vector < unsigned > _sparsityPatternMinimumSize; 
  /** size: number of unknowns */
  std::vector <unsigned> _sparsityPatternVariableIndex;
  
 private:
   
  /** number of non-zeros per row, on-diagonal */
  std::vector < int > d_nnz;
  
  /** number of non-zeros per row, off-diagonal */
  std::vector < int > o_nnz;
  
//==== MATRIX, Sparsity Pattern - END ========

//==== Residual - BEGIN ========
 public:

  /** Print numeric vector with structure */
  void print_residual_with_structure_matlab_friendly(const unsigned iproc, const std::string filename_base, NumericVector * num_vec_in) const;

  /** To be Added */
  void SetResZero();
   
  /** To be Added */
  void UpdateResidual();
  
  /** Vectors for  residual */
  NumericVector *_RES, *_RESC;

 protected:

//==== Residual - END ========
  

//==== Unknown, Error - BEGIN ========
 public:
   
  /** To be Added */
  void SetEpsZero();

  /** To be Added */
  void SumEpsCToEps();
  
  /** Vectors for error  */
  NumericVector *_EPS, *_EPSC;

 protected:
   
  // member data
  /** size: number of Added Unknowns */
  std::vector <unsigned> _SolPdeIndex;
  
   
//==== Unknown, Error - END ========

 
//==== Unknown, Global Variables - BEGIN ========
 public:
   
   void SetNumberOfGlobalVariables(const unsigned &numberOfGlobalVariables){
    _numberOfGlobalVariables = numberOfGlobalVariables;
  }
  
  
 protected:
   
  unsigned _numberOfGlobalVariables;
//==== Unknown, Global Variables - END ========

//==== Linear Algebra - END =============

 
//==== Solution - BEGIN ========
 protected:
   
  /** Only used in one place to retrieve RemoveNullSpace */
  Solution *_solution;

  /** To be Added */
  unsigned GetIndex(const char name[]);

  /** size: number of Solutions */
  std::vector <int> _SolType;
  /** size: number of Solutions */
  std::vector <char*> _SolName;
  
   
   
//==== Solution - END ========


//==== Solution, Boundary Conditions - BEGIN ========
 protected:
   
  /** size: number of Solutions */
  const std::vector <NumericVector*> *_Bdc;
//==== Solution, Boundary Conditions - END ========
  


};

} //end namespace femus



#endif
