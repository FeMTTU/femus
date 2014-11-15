/*=========================================================================

 Program: FEMuS
 Module: ElemType
 Authors: Eugenio Aulisa, Giorgio Bornia
 
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
  elem_type(){};
  
  /** destructor */
  virtual ~elem_type();

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
  virtual double* GetPhi(const unsigned &ig) const = 0;
  
  /** To be Added */
  virtual double* GetDPhiDXi(const unsigned &ig) const = 0;
  
  /** To be Added */
  virtual double* GetDPhiDEta(const unsigned &ig) const{
    std::cout<<"GetDPhiDEta does not apply for this element dimension\n"; 
    abort();
  }
  
  /** To be Added */
  virtual double* GetDPhiDZeta(const unsigned &ig) const{
    std::cout<<"GetDPhiDZeta does not apply for this element dimension\n"; 
    abort();
  };
  
//   /** To be Added */
//   void GetArea(const double *vt,const double *vty, const double *vtz, const unsigned &ig,
//                double &Weight, double *other_phi) const;

  /** To be Added */
  double  GetGaussWeight(const unsigned ig) const {
    return _GaussPointValues[ig];
  };
  
  /** To be Added */
  unsigned GetGaussPointNumber() const {
    return _GaussPointNumber;
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
  int _nc,_nf,_nlag[3];
  unsigned _SolType;
  const double **_X;
  const int **_IND;
  const int **_KVERT_IND;
  
  double **_prol_val;
  int **_prol_ind;
  double *_mem_prol_val;
  int *_mem_prol_ind;
  basis *_pt_basis;
  
  const double *_GaussPointValues;
  unsigned _GaussPointNumber;
};


class elem_type_1D : public elem_type {
public:
  /** constructor */
  elem_type_1D(const char *solid,const char *order, const char* gauss_order);
  
  /** destructor */
  ~elem_type_1D(){
    
    delete [] _phi;
    delete [] _phi_memory;
    delete [] _dphidxi;
    delete [] _dphidxi_memory;
      
    delete [] _d2phidxi2;
    delete [] _d2phidxi2_memory;
           
  };
   
  void Jacobian_AD(const vector < vector < adept::adouble > > &vt,const unsigned &ig, adept::adouble &Weight, 
		   vector < double > &phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &nablaphi) const;	
		   
  void Jacobian(const vector < vector < double > > &vt,const unsigned &ig, double &Weight, 
		vector < double > &other_phi, vector < double > &gradphi, vector < double > &nablaphi) const;
		
  void JacobianSur_AD(const vector < vector < adept::adouble > > &vt, const unsigned &ig, adept::adouble &Weight, 
	              vector < double > &other_phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &normal) const;
		      
  void JacobianSur(const vector < vector < double > > &vt, const unsigned &ig, double &Weight, 
	           vector < double > &other_phi, vector < double > &gradphi, vector < double > &normal) const;
		    
  inline double* GetPhi(const unsigned &ig) const { return _phi[ig]; }
  inline double* GetDPhiDXi(const unsigned &ig) const { return _dphidxi[ig]; }
		   
  double **_phi;
  double *_phi_memory;
  double **_dphidxi;
  double *_dphidxi_memory;
  
  double **_d2phidxi2;
  double *_d2phidxi2_memory;
};

class elem_type_2D : public elem_type {
public:
  /** constructor */
  elem_type_2D(const char *solid,const char *order, const char* gauss_order);
  
  /** destructor */
  ~elem_type_2D(){
    
    delete [] _phi;
    delete [] _phi_memory;
    delete [] _dphidxi;
    delete [] _dphidxi_memory;
    delete [] _dphideta;
    delete [] _dphideta_memory;
  
    delete [] _d2phidxi2;
    delete [] _d2phidxi2_memory;
    delete [] _d2phideta2;
    delete [] _d2phideta2_memory;
  
    delete [] _d2phidxideta;
    delete [] _d2phidxideta_memory;
    
  };
	
  void Jacobian_AD(const vector < vector < adept::adouble > > &vt,const unsigned &ig, adept::adouble &Weight, 
		   vector < double > &phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &nablaphi) const;	
		   
  void Jacobian(const vector < vector < double > > &vt,const unsigned &ig, double &Weight, 
		vector < double > &other_phi, vector < double > &gradphi, vector < double > &nablaphi) const;
		
  void JacobianSur_AD(const vector < vector < adept::adouble > > &vt, const unsigned &ig, adept::adouble &Weight, 
	              vector < double > &other_phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &normal) const;
		      
  void JacobianSur(const vector < vector < double > > &vt, const unsigned &ig, double &Weight, 
	           vector < double > &other_phi, vector < double > &gradphi, vector < double > &normal) const;

  inline double* GetPhi(const unsigned &ig) const { return _phi[ig]; }
  inline double* GetDPhiDXi(const unsigned &ig) const { return _dphidxi[ig]; }
  inline double* GetDPhiDEta(const unsigned &ig) const { return _dphideta[ig]; }
	   
private:    
  double **_phi;
  double *_phi_memory;
  double **_dphidxi;
  double *_dphidxi_memory;
  double **_dphideta;
  double *_dphideta_memory;
  
  double **_d2phidxi2;
  double *_d2phidxi2_memory;
  double **_d2phideta2;
  double *_d2phideta2_memory;
  
  double **_d2phidxideta;
  double *_d2phidxideta_memory;
};

class elem_type_3D : public elem_type {
public:
  /** constructor */
  elem_type_3D(const char *solid,const char *order, const char* gauss_order);
   /** destructor */
  ~elem_type_3D(){
    delete [] _phi;
    delete [] _phi_memory;
    delete [] _dphidxi;
    delete [] _dphidxi_memory;
    delete [] _dphideta;
    delete [] _dphideta_memory;
    delete [] _dphidzeta;
    delete [] _dphidzeta_memory;
  
    delete [] _d2phidxi2;
    delete [] _d2phidxi2_memory;
    delete [] _d2phideta2;
    delete [] _d2phideta2_memory;
    delete [] _d2phidzeta2;
    delete [] _d2phidzeta2_memory;
  
    delete [] _d2phidxideta;
    delete [] _d2phidxideta_memory;
    delete [] _d2phidetadzeta;
    delete [] _d2phidetadzeta_memory;
    delete [] _d2phidzetadxi;
    delete [] _d2phidzetadxi_memory;
      
  };
	
  void Jacobian_AD(const vector < vector < adept::adouble > > &vt,const unsigned &ig, adept::adouble &Weight, 
		   vector < double > &phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &nablaphi) const;
		   
  void Jacobian(const vector < vector < double > > &vt,const unsigned &ig, double &Weight, 
		vector < double > &other_phi, vector < double > &gradphi, vector < double > &nablaphi) const;
		
  void JacobianSur_AD(const vector < vector < adept::adouble > > &vt, const unsigned &ig, adept::adouble &Weight, 
	              vector < double > &other_phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &normal) const{};
		      
  void JacobianSur(const vector < vector < double > > &vt, const unsigned &ig, double &Weight, 
	           vector < double > &other_phi, vector < double > &gradphi, vector < double > &normal) const{};
		   
  
  //---------------------------------------------------------------------------------------------------------
  inline double* GetPhi(const unsigned &ig) const { return _phi[ig]; }
  inline double* GetDPhiDXi(const unsigned &ig) const { return _dphidxi[ig]; }
  inline double* GetDPhiDEta(const unsigned &ig) const { return _dphideta[ig]; }
  inline double* GetDPhiDZeta(const unsigned &ig) const { return _dphidzeta[ig];}
	
private:	   		   
  double **_phi;
  double *_phi_memory;
  double **_dphidxi;
  double *_dphidxi_memory;
  double **_dphideta;
  double *_dphideta_memory;
  double **_dphidzeta;
  double *_dphidzeta_memory;
  
  double **_d2phidxi2;
  double *_d2phidxi2_memory;
  double **_d2phideta2;
  double *_d2phideta2_memory;
  double **_d2phidzeta2;
  double *_d2phidzeta2_memory;
  
  double **_d2phidxideta;
  double *_d2phidxideta_memory;
  double **_d2phidetadzeta;
  double *_d2phidetadzeta_memory;
  double **_d2phidzetadxi;
  double *_d2phidzetadxi_memory;
};



} //end namespace femus



#endif


