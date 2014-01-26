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
 //Data 
private:
  bool _is_symmetric;
  bool _stabilization;
  double _compressibility;
  
protected:
  vector <int> _SolType;  
  vector <char*> _SolName;
  const vector <NumericVector*> *_Bdc;
  
public:   
  mesh *_msh; 
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
    
  //Functions
  
public:  
  lsysPDE(mesh *other_msh);
  ~lsysPDE();

  void AddSolutionVector(const char name[], const char order[],const unsigned& tmorder, const bool &PDE_type=1);
  void SetBdcPointer(vector <NumericVector*> *Bdc_other);
  int InitMultigrid(const vector <unsigned> &MGIndex);
  int DeallocateMatrix();
  int AllocateMatrix();
  unsigned GetKKDof(const unsigned &index_sol, const unsigned &kkindex_sol,const unsigned &idof_gmt) const;
  void SetMatrixProperties(const bool property);
  bool GetMatrixProperties();
  void AddStabilization(const bool stab, const double compressibility);
  double GetCompressibility();
  bool GetStabilization();
  int SetResZero();
  int SetEpsZero();
  int SumEpsCToEps();
  int UpdateResidual();

protected:
  unsigned GetIndex(const char name[]);
};


/*----------------------- functions ----------------------------------*/
//--------------------------------------------------------------------------------------------
inline unsigned lsysPDE::GetKKDof(const unsigned &index_sol, const unsigned &kkindex_sol, 
				  const unsigned &idof_gmt) const {
  
   //return KKIndex[kkindex_sol]+idof_gmt;
     
   unsigned soltype =  _SolType[index_sol]; 
   unsigned isubdom = (soltype<3)?_msh->npart[idof_gmt]:(_msh->epart[idof_gmt % _msh->GetElementNumber()]);
   unsigned idof_metis = _msh->GetMetisDof(idof_gmt,soltype);   
   return KKoffset[kkindex_sol][isubdom] + idof_metis - _msh->MetisOffset[soltype][isubdom];
}
//--------------------------------------------------------------------------------------------

#endif