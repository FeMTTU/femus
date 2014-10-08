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

#ifndef __linear_equation_hpp__
#define __linear_equation_hpp__

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
class mesh;

/**
 This class is a container that holds linear operators and other structure for solving a linear equation system
*/

class LinearEquation : public ParallelObject {

public:   

  
public:  
  
  /** costructor */
  LinearEquation(mesh *other_msh);
  
  /** destructor */
  ~LinearEquation();
  
  /** To be Added */
  void InitPde(const vector <unsigned> &_SolPdeIndex,const  vector <int> &SolType,  
               const vector <char*> &SolName, vector <NumericVector*> *Bdc_other, 
               const unsigned &other_gridr, const unsigned &other_gridn);
  
  void GetSparsityPatternSize();
  
  /** To be Added */
  void DeletePde();
  
  /** To be Added */
  unsigned GetKKDof(const unsigned &index_sol, const unsigned &kkindex_sol,const unsigned &idof_gmt) const;
  
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
  
  // member data
  mesh *_msh; 
  NumericVector *_EPS, *_EPSC, *_RES, *_RESC;
  SparseMatrix *_KK, *_PP,*_RR, *_CC; 
  vector < vector <unsigned> > KKoffset;
  vector < unsigned > KKghostsize;
  vector < vector < int> > KKghost_nd;
  vector <int> KKIndex;
  bool _CC_flag; 
  unsigned _gridr,_gridn;
  
  vector < int > d_nnz;
  vector < int > o_nnz;

protected:
  
  /** To be Added */
  unsigned GetIndex(const char name[]);
  
  // member data 
  vector <unsigned> _SolPdeIndex;
  vector <int> _SolType;  
  vector <char*> _SolName;
  const vector <NumericVector*> *_Bdc;
    
};

} //end namespace femus



#endif