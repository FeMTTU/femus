#ifndef __solution_hpp__
#define __solution_hpp__

#include "Mesh.hpp"
#include "petscmat.h"

//Forward Declaration
class elem_type;
class NumericVector;
class SparseRectangularMatrix;

//class lsysPDE: public mesh {

class Solution{
  
public:
  
  mesh *_msh;
  
  vector <int> _SolType;  
  vector <char*> _SolName;
  vector <unsigned> _SolTmOrder;
  vector <NumericVector*> _Sol;      //one for every variable
  vector <NumericVector*> _SolOld;   //one for every variable
  vector <NumericVector*> _Res;      //one for every variable
  vector <NumericVector*> _Eps;      //one for every variable
  vector <NumericVector*> _Bdc;      //one for every variable
  vector <bool> _ResEpsBdcFlag;
    
  //Projection matrices (for every kind of variable)
  SparseRectangularMatrix* _ProjMat[5];
  bool _ProjMatFlag[5];

  static const unsigned _END_IND[5];

  // Constructor - Destructor
  Solution(mesh *other_msh);
  ~Solution();

  // function
  void AddSolutionVector( const char name[], const char order[],const unsigned& tmorder, const bool &PDE_type=1);
  unsigned GetIndex(const char name[]) const;
  void ResizeSolutionVector(const char name[]);
  void FreeSolutionVectors();
  void SetCoarseCoordinates( vector < vector < double> > &vt);
  int SumEpsToSol(const vector <unsigned> &MGIndex, const Vec &EPS, const Vec &RES, const vector <vector <unsigned> > &KKoffset);
  void UpdateSolution();
  void set_elr(const unsigned &test=0);
};

#endif