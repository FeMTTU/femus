#ifndef __lsyspde_hpp__
#define __lsyspde_hpp__

#include "Mesh.hpp"
#include "petscmat.h"

//Forward Declaration
class elem_type;
class NumericVector;
class SparseRectangularMatrix;

//class lsysPDE: public mesh {

class lsysPDE{
  
private:
  bool _is_symmetric;
  bool _stabilization;
  double _compressibility;

public:
  
  mesh *_msh;
  
 
  //vector <vector <unsigned short> > Bdc;  
  vector <int> SolType;  
  vector <char*> SolName;
  //vector <unsigned> SolTmOrder;
  //vector <NumericVector*> Sol_;      //one for every variable
  //vector <NumericVector*> Sol_old_;  //one for every variable
  //vector <NumericVector*> Res_;      //one for every variable
  //vector <NumericVector*> Eps_;      //one for every variable
  vector <NumericVector*> Bdc_;      //one for every variable
  vector <bool> ResEpsBdc_flag_;
    
  vector <NumericVector*> *_Bdc;
 
  void GetBoundaryCondition(vector <NumericVector*> *Bdc_other);
  
  vector <PetscInt> DrchKKdofs;
  vector < vector <unsigned> > KKoffset;
  vector < unsigned > KKghostsize;
  vector < vector < int> > KKghost_nd;
  vector <int> KKIndex;
  Mat PP;
  Mat KK,CC;
  Vec RES,RESC;
  Vec EPS,EPSC;
  bool CC_flag;
  
  //Projection matrices (for every kind of variable)
//   SparseRectangularMatrix* Proj_mat[5];
//   bool Proj_mat_flag[5];

  static const unsigned END_IND[5];

  // Constructor - Destructor
  lsysPDE(const char infile[], vector < vector < double> > &vt,const double Lref);
  lsysPDE(const unsigned &igrid,elem *elc);
  lsysPDE(mesh *other_msh);
  ~lsysPDE();

  // function
  int InitMultigrid(const vector <unsigned> &MGIndex);

  int DeallocateMatrix(const vector <unsigned> &MGIndex);
  int AllocateMatrix(const vector <unsigned> &MGIndex);

  void AddSolutionVector( const char name[], const char order[],const unsigned& tmorder, const bool &PDE_type=1);
  void FreeSolutionVectors();
  void ResizeSolutionVector(const char name[]);

  unsigned GetIndex(const char name[]);
  unsigned GetKKDof(const unsigned &index_sol, const unsigned &kkindex_sol,const unsigned &idof_gmt) const;
  void SetMatrixProperties(const bool property);
  bool GetMatrixProperties();
  void AddStabilization(const bool stab, const double compressibility);
  double GetCompressibility();
  bool GetStabilization();

  int SetResZero(const vector <unsigned> &MGIndex);
  int SetEpsZero(const vector <unsigned> &MGIndex);
  int SumEpsCToEps(const vector <unsigned> &MGIndex);
  //int SumEpsToSol(const vector <unsigned> &MGIndex);

  int UpdateResidual();
  void UpdateSolution();

  //void GetCoordinates( vector < vector < double> > &vt);
  //void set_elr(const unsigned &test=0);
};


/*----------------------- functions ----------------------------------*/
//--------------------------------------------------------------------------------------------
inline unsigned lsysPDE::GetKKDof(const unsigned &index_sol, const unsigned &kkindex_sol, 
				  const unsigned &idof_gmt) const {
  
   //return KKIndex[kkindex_sol]+idof_gmt;
     
   unsigned soltype =  SolType[index_sol]; 
   unsigned isubdom = (soltype<3)?_msh->npart[idof_gmt]:(_msh->epart[idof_gmt % _msh->GetElementNumber()]);
   unsigned idof_metis = _msh->GetMetisDof(idof_gmt,soltype);   
   return KKoffset[kkindex_sol][isubdom] + idof_metis - _msh->MetisOffset[soltype][isubdom];
}
//--------------------------------------------------------------------------------------------

#endif