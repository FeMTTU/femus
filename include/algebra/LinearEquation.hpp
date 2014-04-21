/*=========================================================================

 Program: FEMUS
 Module: LinearEquation
 Authors: Eugenio Aulisa, Simone Bnà
 
 Copyright (c) FEMTTU
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

class LinearEquation {

  //Data 
private:
  bool _is_symmetric;
  bool _stabilization;
  double _compressibility;
  
protected:
  vector <unsigned> _SolPdeIndex;
  vector <int> _SolType;  
  vector <char*> _SolName;
  const vector <NumericVector*> *_Bdc;
  unsigned int _DirichletBCsHandlingMode; //* 0 Penalty method,  1 Elimination method */
  
public:   
  mesh *_msh; 
  vector < vector <unsigned> > KKoffset;
  vector < unsigned > KKghostsize;
  vector < vector < int> > KKghost_nd;
  vector <int> KKIndex;
  
  
    
  NumericVector *_EPS, *_EPSC, *_RES, *_RESC;
  SparseMatrix *_KK, *_PP,*_RR, *_CC; //will become SparseMatrix ASAP
  bool _CC_flag; 
  unsigned _gridr,_gridn;
  
  //Functions
  
public:  
  
  /** costructor */
  LinearEquation(mesh *other_msh);
  
  /** destructor */
  ~LinearEquation();
  
  void InitPde(const vector <unsigned> &_SolPdeIndex,const  vector <int> &SolType,  
	      const vector <char*> &SolName, vector <NumericVector*> *Bdc_other, 
	      const unsigned &other_gridr, const unsigned &other_gridn);
  
  void DeletePde();
  unsigned GetKKDof(const unsigned &index_sol, const unsigned &kkindex_sol,const unsigned &idof_gmt) const;
  
  void set_dirichletBCsHandling(unsigned int DirichletBCsHandlingMode);
  void SetMatrixProperties(const bool property);
  bool GetMatrixProperties();
  void AddStabilization(const bool stab, const double compressibility);
  double GetCompressibility();
  bool GetStabilization();
  void SetResZero();
  void SetEpsZero();
  void SumEpsCToEps();
  void UpdateResidual();

protected:
  unsigned GetIndex(const char name[]);
};

} //end namespace femus



#endif