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
  void prolongation(const mesh &meshf,const mesh &meshc, const int& ielc, SparseMatrix* Projmat) const;
  
  /** To be Added */
  void ProlQitoQj(const mesh& mymesh,const int& iel, SparseMatrix* Projmat, 
                  bool testnode[],const unsigned &itype) const;

  /** To be Added */ 
  void JacobianSur2D(const vector < vector < double > > &vt, const unsigned &ig,
                     double &Weight, vector < double > &other_phi, vector < double > &gradphi, vector < double > &normal)const;
   
  /** To be Added */
  void JacobianSur1D(const vector < vector < double > > &vt, const unsigned &ig,
                     double &Weight, vector < double > &other_phi, vector < double > &gradphi, vector < double > &normal) const;

  /** To be Added */
  void Jacobian3D(const vector < vector < double > > &vt,const unsigned &ig,
                  double &Weight, vector < double > &other_phi, vector < double > &gradphi)const;
  
  /** To be Added */
  void Jacobian2D(const vector < vector < double > > &vt,const unsigned &ig,
                  double &Weight, vector < double > &other_phi, vector < double > &gradphi) const;

  /** To be Added */
  void Jacobian1D(const vector < vector < double > > &vt,const unsigned &ig,
                  double &Weight, vector < double > &other_phi, vector < double > &gradphi) const;

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
 
  void (elem_type::*Jacobian_ptr)(const vector < vector < double > > &vt, const unsigned &ig,
                                  double &Weight, vector < double > &other_phi, vector < double > &gradphi) const;
  
  void (elem_type::*Jacobian_sur_ptr)(const vector < vector < double > > &vt, const unsigned &ig,
                                      double &Weight, vector < double > &other_phi, vector < double > &gradphi, vector < double > &normal) const;
  
private:
  
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
  hex0 hex_0;
  hexpwl hex_pwl;
  hex1 hex_1;
  hexth hex_th;
  hex2 hex_2;
  wedge1 wedge_1;
  wedgeth wedge_th;
  wedge2 wedge_2;
  tet1 tet_1;
  tet2 tet_2;
  quad0 quad_0;
  quadpwl quad_pwl;
  quad1 quad_1;
  quadth quad_th;
  quad2 quad_2;
  tri1 tri_1;
  tri2 tri_2;
  line1 line_1;
  line2 line_2;
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
  const double *weight;
  
};


} //end namespace femus



#endif


