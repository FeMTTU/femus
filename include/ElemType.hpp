#ifndef __elem_type_hpp__
#define __elem_type_hpp__

#include "Basis.hpp"
#include "petscksp.h"
#include "petscvec.h" 
#include "SparseRectangularMatrix.hpp"
#include "Mesh.hpp"
#include "LinSysPde.hpp"

class elem;

class elem_type {
private:
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

public:
  static unsigned _refindex;
  elem_type(const char *solid,const char *order, const char* gauss_order);
  ~elem_type();
  
  int prolongation(const elem* el,const elem* elc, const int& ielc, 
		   Mat &PP,const PetscInt& istart=0,const PetscInt& jstart=0) const ;
		   
  void prolongation(const lsysPDE &lspdef,const lsysPDE &lspdec, const int& ielc, Mat& PP, 
		    const unsigned &index_sol, const unsigned &kkindex_sol) const;
		   
  void prolongation(const mesh &meshf,const mesh &meshc, const int& ielc, SparseRectangularMatrix* Projmat) const;
  
  void ProlQitoQj(const mesh& mymesh,const int& iel, SparseRectangularMatrix* Projmat, 
		  bool testnode[],const unsigned &itype) const;

  // Jacobian over Riemaniann Manifold (2D or 3D)
  void JacobianSur2D(const double vt[][27],const unsigned &ig,
                     double &Weight, double *other_phi, double gradphi[][3], double normal[3])const;
  void JacobianSur1D(const double vt[][27],const unsigned &ig,
                     double &Weight, double *other_phi,double gradphi[][3], double normal[3]) const;
  void (elem_type::*Jacobian_sur_ptr)(const double vtx[][27],const unsigned &ig,
                                      double &Weight, double *other_phi, double gradphi[][3], double normal[3]) const;


  void Jacobian3D(const double vt[][27],const unsigned &ig,
                  double &Weight, double *other_phi, double gradphi[][3])const;
  void Jacobian2D(const double vt[][27],const unsigned &ig,
                  double &Weight, double *other_phi,double gradphi[][3]) const;
  void Jacobian1D(const double vt[][27],const unsigned &ig,
                  double &Weight, double *other_phi, double gradphi[][3]) const;
  void (elem_type::*Jacobian_ptr)(const double vtx[][27],const unsigned &ig,
                                  double &Weight, double *other_phi, double gradphi[][3]) const;

  double* GetPhi(const unsigned &ig) const;

  void GetArea(const double *vt,const double *vty, const double *vtz, const unsigned &ig,
               double &Weight, double *other_phi) const;

  double  GetGaussWeight(const unsigned ig) const {
    return GaussWeight[ig];
  };
  unsigned GetGaussPointNumber() const {
    return GaussPoints;
  };
};

#endif


