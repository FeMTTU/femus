/*=========================================================================

 Program: FEMUS
 Module: Solution
 Authors: Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/**
 * This class is useful to handle multilevel variables. It's going to be splitted
 * into meshvariables and solutionvariables 
*/


#ifndef __solution_hpp__
#define __solution_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Mesh.hpp"
#include "petscmat.h"

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class elem_type;
class NumericVector;
class SparseMatrix;


class Solution {
  
  //member data
private:
  vector <int> _SolType;  
  vector <char*> _SolName;
  vector <unsigned> _SolTmOrder;
  mesh *_msh;
public:
  vector <NumericVector*> _Sol;      //one for every variable
  vector <NumericVector*> _SolOld;   //one for every variable
  vector <NumericVector*> _Res;      //one for every variable
  vector <NumericVector*> _Eps;      //one for every variable
  vector <NumericVector*> _Bdc;      //one for every variable
  vector <bool> _ResEpsBdcFlag;      //one for every variable
    
  SparseMatrix* _ProjMat[5];  //Projection matrices (one for every type of variable)
  bool _ProjMatFlag[5];			 //One for every type of variable

  // Constructor - Destructor
  Solution(mesh *other_msh);
  ~Solution();
  
  // functions
private:
  unsigned GetIndex(const char name[]) const;
public:
  void AddSolution( const char name[], const char order[],const unsigned& tmorder, const bool &Pde_type=1);
  void ResizeSolutionVector(const char name[]);
  void FreeSolutionVectors();
  void SetCoarseCoordinates( vector < vector < double> > &vt);
  void SumEpsToSol(const vector <unsigned> &_SolPdeIndex,  NumericVector* EPS, NumericVector* RES, const vector <vector <unsigned> > &KKoffset);
  void UpdateSolution();
  void SetElementRefiniement(const unsigned &test=0);

};

#endif