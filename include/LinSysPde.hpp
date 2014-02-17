#ifndef __lsyspde_hpp__
#define __lsyspde_hpp__

#include "Mesh.hpp"
#include "petscmat.h"

//Forward Declaration
class elem_type;
class NumericVector;
class SparseMatrix;

//class lsysPde: public mesh {

class lsysPde{
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
  vector <int> _KKdofs_tobe_Solved;
  vector < vector <unsigned> > KKoffset;
  vector < unsigned > KKghostsize;
  vector < vector < int> > KKghost_nd;
  vector <int> KKIndex;
  
  
    
  NumericVector *_EPS, *_EPSC, *_RES, *_RESC;
  SparseMatrix *_KK, *_PP, *_CC; //will become SparseMatrix ASAP
  bool _CC_flag;
  
  //Functions
  
public:  
  lsysPde(mesh *other_msh);
  ~lsysPde();
  
  int InitPde(const vector <unsigned> &_SolPdeIndex,const  vector <int> &SolType,  
	      const vector <char*> &SolName, vector <NumericVector*> *Bdc_other);
  void DeletePde();
  unsigned GetKKDof(const unsigned &index_sol, const unsigned &kkindex_sol,const unsigned &idof_gmt) const;
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


/*----------------------- functions ----------------------------------*/
//--------------------------------------------------------------------------------------------
inline unsigned lsysPde::GetKKDof(const unsigned &index_sol, const unsigned &kkindex_sol, 
				  const unsigned &idof_gmt) const {
  
   //return KKIndex[kkindex_sol]+idof_gmt;
     
   unsigned soltype =  _SolType[index_sol]; 
   unsigned isubdom = (soltype<3)?_msh->npart[idof_gmt]:(_msh->epart[idof_gmt % _msh->GetElementNumber()]);
   unsigned idof_metis = _msh->GetMetisDof(idof_gmt,soltype);   
   return KKoffset[kkindex_sol][isubdom] + idof_metis - _msh->MetisOffset[soltype][isubdom];
}
//--------------------------------------------------------------------------------------------

#endif