/*=========================================================================

 Program: FEMuS
 Module: ElemType
 Authors: Eugenio Aulisa
 
 Copyright (c) FEMuS
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/**
 * The ElemType class
*/

#ifndef __elem_type_hpp__
#define __elem_type_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Basis.hpp"
#include "SparseMatrix.hpp"
#include "Mesh.hpp"
#include "LinearEquation.hpp"
#include "adept.h"

namespace femus {



//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class elem;
class LinearEquation;

class elem_type {
  
public:

  /** constructor */
  elem_type(const char *solid,const char *order, const char* gauss_order);
  
  /** destructor */
  ~elem_type();

  /** To be Added */
  void BuildProlongation(const LinearEquation &lspdef,const LinearEquation &lspdec, const int& ielc, SparseMatrix* Projmat, 
                         const unsigned &index_sol, const unsigned &kkindex_sol) const;

  /** To be Added */  
  void BuildRestrictionTranspose(const LinearEquation &lspdef,const LinearEquation &lspdec, const int& ielc, SparseMatrix* Projmat, 
                                 const unsigned &index_sol, const unsigned &kkindex_sol, const bool &TestDisp) const;		    

  /** To be Added */
  void BuildProlongation(const mesh &meshf, const mesh &meshc, const int& ielc, SparseMatrix* Projmat) const;
  
  /** To be Added */
  void BuildProlongation(const mesh& mymesh, const int& iel, SparseMatrix* Projmat, const unsigned &itype) const;
  
	
  virtual void Jacobian_AD(const vector < vector < adept::adouble > > &vt,const unsigned &ig, adept::adouble &Weight, 
			   vector < double > &phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &nablaphi) const = 0;			 
			 
  virtual void Jacobian(const vector < vector < double > > &vt,const unsigned &ig, double &Weight, 
			vector < double > &other_phi, vector < double > &gradphi, vector < double > &nablaphi) const = 0;
			
  virtual void JacobianSur_AD(const vector < vector < adept::adouble > > &vt, const unsigned &ig, adept::adouble &Weight, 
			      vector < double > &other_phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &normal) const = 0;
			      
  virtual void JacobianSur(const vector < vector < double > > &vt, const unsigned &ig, double &Weight, 
			   vector < double > &other_phi, vector < double > &gradphi, vector < double > &normal) const = 0;
  /** To be Added */
  double* GetPhi(const unsigned &ig) const;
  
  /** To be Added */
  double* GetDPhiDXi(const unsigned &ig) const;
  
  /** To be Added */
  double* GetDPhiDEta(const unsigned &ig) const;
  
  /** To be Added */
  double* GetDPhiDZeta(const unsigned &ig) const;
  
  /** To be Added */
  void GetArea(const double *vt,const double *vty, const double *vtz, const unsigned &ig,
               double &Weight, double *other_phi) const;

  /** To be Added */
  double  GetGaussWeight(const unsigned ig) const {
    return GaussWeight[ig];
  };
  
  /** To be Added */
  unsigned GetGaussPointNumber() const {
    return GaussPoints;
  };
  
  // member data
  static unsigned _refindex;
				  
  void (elem_type::*Jacobian_sur_ptr)(const vector < vector < double > > &vt, const unsigned &ig,
                                      double &Weight, vector < double > &other_phi, vector < double > &gradphi, vector < double > &normal) const;

  void GetSparsityPatternSize(const LinearEquation &lspdef,const LinearEquation &lspdec, const int& ielc,  
			      NumericVector* NNZ_d, NumericVector* NNZ_o,
			      const unsigned &index_sol, const unsigned &kkindex_sol) const; 
  void GetSparsityPatternSize(const mesh &meshf,const mesh &meshc, const int& ielc, NumericVector* NNZ_d, NumericVector* NNZ_o) const;					 
  					 
  void GetSparsityPatternSize(const mesh& mesh,const int& iel, NumericVector* NNZ_d, NumericVector* NNZ_o, const unsigned &itype) const;
				      
protected:
  
  // member data
  void test_prol_and_rest();
  int nc_,nf_,ncf_[3];
  unsigned type_;
  unsigned SolType_;
  const double **X;
  const int **IND;
  const int **KVERT_IND;
  double** rest_val;
  int** rest_ind;
  double* mem_rest_val;
  int * mem_rest_ind;
  double** prol_val;
  int** prol_ind;
  double* mem_prol_val;
  int * mem_prol_ind;
  basis *pt_basis;
  
  const double *GaussWeight;
  unsigned GaussPoints;
  double **phi;
  double *phi_memory;
  double **dphidxi;
  double *dphidxi_memory;
  double **dphideta;
  double *dphideta_memory;
  double **dphidzeta;
  double *dphidzeta_memory;
  
  double **d2phidxi2;
  double *d2phidxi2_memory;
  double **d2phideta2;
  double *d2phideta2_memory;
  double **d2phidzeta2;
  double *d2phidzeta2_memory;
  
  double **d2phidxideta;
  double *d2phidxideta_memory;
  double **d2phidetadzeta;
  double *d2phidetadzeta_memory;
  double **d2phidzetadxi;
  double *d2phidzetadxi_memory;
  
  const double *weight;
};


class elem_type_1D : public elem_type {
public:
  /** constructor */
  elem_type_1D(const char *solid,const char *order, const char* gauss_order):
    elem_type(solid,order,gauss_order){};
	
  void Jacobian_AD(const vector < vector < adept::adouble > > &vt,const unsigned &ig, adept::adouble &Weight, 
		   vector < double > &phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &nablaphi) const;
		   
  void Jacobian(const vector < vector < double > > &vt,const unsigned &ig, double &Weight, 
		vector < double > &other_phi, vector < double > &gradphi, vector < double > &nablaphi) const;
		
  void JacobianSur_AD(const vector < vector < adept::adouble > > &vt, const unsigned &ig, adept::adouble &Weight, 
	              vector < double > &other_phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &normal) const;
		      
  void JacobianSur(const vector < vector < double > > &vt, const unsigned &ig, double &Weight, 
	           vector < double > &other_phi, vector < double > &gradphi, vector < double > &normal) const;
};

class elem_type_2D : public elem_type {
public:
  /** constructor */
  elem_type_2D(const char *solid,const char *order, const char* gauss_order):
	elem_type(solid,order,gauss_order){};
	
  void Jacobian_AD(const vector < vector < adept::adouble > > &vt,const unsigned &ig, adept::adouble &Weight, 
		   vector < double > &phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &nablaphi) const;	
		   
  void Jacobian(const vector < vector < double > > &vt,const unsigned &ig, double &Weight, 
		vector < double > &other_phi, vector < double > &gradphi, vector < double > &nablaphi) const;
		
  void JacobianSur_AD(const vector < vector < adept::adouble > > &vt, const unsigned &ig, adept::adouble &Weight, 
	              vector < double > &other_phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &normal) const;
		      
  void JacobianSur(const vector < vector < double > > &vt, const unsigned &ig, double &Weight, 
	           vector < double > &other_phi, vector < double > &gradphi, vector < double > &normal) const;
};

class elem_type_3D : public elem_type {
public:
  /** constructor */
  elem_type_3D(const char *solid,const char *order, const char* gauss_order):
	elem_type(solid,order,gauss_order){};
	
  void Jacobian_AD(const vector < vector < adept::adouble > > &vt,const unsigned &ig, adept::adouble &Weight, 
		   vector < double > &phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &nablaphi) const;
		   
  void Jacobian(const vector < vector < double > > &vt,const unsigned &ig, double &Weight, 
		vector < double > &other_phi, vector < double > &gradphi, vector < double > &nablaphi) const;
		
  void JacobianSur_AD(const vector < vector < adept::adouble > > &vt, const unsigned &ig, adept::adouble &Weight, 
	              vector < double > &other_phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &normal) const
	              {};
		      
  void JacobianSur(const vector < vector < double > > &vt, const unsigned &ig, double &Weight, 
	           vector < double > &other_phi, vector < double > &gradphi, vector < double > &normal) const
	           {};
};



} //end namespace femus



#endif


