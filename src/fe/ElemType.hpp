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

#ifndef __femus_fe_ElemType_hpp__
#define __femus_fe_ElemType_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Basis.hpp"
#include "SparseMatrix.hpp"
#include "Mesh.hpp"
#include "LinearEquation.hpp"
#include "GaussPoints.hpp"
#include "adept.h"
#include "FETypeEnum.hpp"


namespace femus {

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class elem;
class LinearEquation;

class elem_type {

public:

  /** constructor that receives Geometric Element and Gauss info */
  elem_type(const char *geom_elem, const char *order_gauss);

  /** destructor */
  virtual ~elem_type();

  /** To be Added */
  void BuildProlongation(const LinearEquation &lspdef,const LinearEquation &lspdec, const int& ielc, SparseMatrix* Projmat,
                         const unsigned &index_sol, const unsigned &kkindex_sol) const;

  /** To be Added */
  void BuildRestrictionTranspose(const LinearEquation &lspdef,const LinearEquation &lspdec, const int& ielc, SparseMatrix* Projmat,
                                 const unsigned &index_sol, const unsigned &kkindex_sol,
				 const unsigned &index_pair_sol, const unsigned &kkindex_pair_sol) const;

  /** To be Added */
  void BuildProlongation(const Mesh &meshf, const Mesh &meshc, const int& ielc, SparseMatrix* Projmat) const;

  /** To be Added */
  void BuildProlongation(const Mesh& mymesh, const int& iel, SparseMatrix* Projmat, const unsigned &itype) const;


  virtual void Jacobian(const vector < vector < adept::adouble > > &vt,const unsigned &ig, adept::adouble &Weight,
 			vector < double > &phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &nablaphi) const = 0;

  virtual void Jacobian(const vector < vector < double > > &vt,const unsigned &ig, double &Weight,
			vector < double > &other_phi, vector < double > &gradphi, vector < double > &nablaphi) const = 0;

  virtual void JacobianSur(const vector < vector < adept::adouble > > &vt, const unsigned &ig, adept::adouble &Weight,
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

  /** @deprecated  Function pointer for DPhiDXEZ */
  typedef double* (elem_type::*_FunctionPointer)(const unsigned & ig) const;  //you need "elem_type::" for some reason
  std::vector<_FunctionPointer> _DPhiXiEtaZetaPtr;

  /** @deprecated Evaluate shape functions at all quadrature points */
  virtual void EvaluateShapeAtQP(const std::string geomel_id_in,const std::string fe_in);

  /**  @deprecated Get shape functions */
  inline const double GetPhi(const uint qp, const uint dof ) const {
     return _phi_mapGD[qp][dof];
    }

  /**  @deprecated Get shape function first derivatives */
  inline const double GetDPhiDxez(const uint qp, const uint dof ) const {
     return _dphidxez_mapGD[qp][dof];
    }

  /** To be Added */
  inline const Gauss GetGaussRule() const {
    return _gauss;
  };

  /** To be Added */
  inline double  GetGaussWeight(const unsigned ig) const {
    return _gauss.GetGaussWeightsPointer()[ig];
  };

  /** To be Added */
  inline unsigned GetGaussPointNumber() const {
    return _gauss.GetGaussPointsNumber();
  };

  /** Retrieve the number of dofs for this element */
  inline int  GetNDofs() const {
    return _nc;
  };

  /** Retrieve the dimension of the underlying geometric element */
  inline unsigned  GetDim() const {
    return _dim;
  };

  // member data
  static unsigned _refindex;

  void GetSparsityPatternSize(const LinearEquation &lspdef,const LinearEquation &lspdec, const int& ielc,
			      NumericVector* NNZ_d, NumericVector* NNZ_o,
			      const unsigned &index_sol, const unsigned &kkindex_sol) const;
  void GetSparsityPatternSize(const Mesh &meshf,const Mesh &meshc, const int& ielc, NumericVector* NNZ_d, NumericVector* NNZ_o) const;

  void GetSparsityPatternSize(const Mesh& Mesh,const int& iel, NumericVector* NNZ_d, NumericVector* NNZ_o, const unsigned &itype) const;

  static const unsigned _fe_old_to_new[QL];

  static const unsigned _fe_new_to_old[NFE_FAMS];

protected:

  // member data
  unsigned _dim; /*Spatial dimension of the geometric element*/
  int _nc,_nf,_nlag[3];
  unsigned _SolType;   /*Finite Element Family flag*/
  const double **_X;
  const int **_IND;
  const int **_KVERT_IND;

  double **_prol_val;
  int **_prol_ind;
  double *_mem_prol_val;
  int *_mem_prol_ind;
  basis *_pt_basis;

//  Gauss
  const Gauss _gauss;

  /**  @deprecated */
  bool isMpGDAllocated;
  double**      _phi_mapGD;
  double** _dphidxez_mapGD;

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

  template <class type>
  void Jacobian_type( const vector < vector < type > > &vt,const unsigned &ig, type &Weight,
		      vector < double > &phi, vector < type > &gradphi, vector < type > &nablaphi) const;

  void Jacobian( const vector < vector < adept::adouble > > &vt,const unsigned &ig, adept::adouble &Weight,
 		 vector < double > &phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &nablaphi) const{
		 Jacobian_type(vt, ig, Weight, phi, gradphi, nablaphi);
		}
  void Jacobian(const vector < vector < double > > &vt,const unsigned &ig, double &Weight,
		vector < double > &phi, vector < double > &gradphi, vector < double > &nablaphi) const{
		Jacobian_type(vt, ig, Weight, phi, gradphi, nablaphi);
		}

  template <class type>
  void JacobianSur_type(const vector < vector < type > > &vt, const unsigned &ig, type &Weight,
			vector < double > &phi, vector < type > &gradphi, vector < type > &normal) const;

  void JacobianSur(const vector < vector < adept::adouble > > &vt, const unsigned &ig, adept::adouble &Weight,
	           vector < double > &phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &normal) const{
		   JacobianSur_type(vt, ig, Weight, phi, gradphi, normal);
		  }

  void JacobianSur(const vector < vector < double > > &vt, const unsigned &ig, double &Weight,
	           vector < double > &phi, vector < double > &gradphi, vector < double > &normal) const{
		     JacobianSur_type(vt, ig, Weight, phi, gradphi, normal);
		  }

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

  template <class type>
  void Jacobian_type( const vector < vector < type > > &vt,const unsigned &ig, type &Weight,
		      vector < double > &phi, vector < type > &gradphi, vector < type > &nablaphi) const;

  void Jacobian( const vector < vector < adept::adouble > > &vt,const unsigned &ig, adept::adouble &Weight,
 		 vector < double > &phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &nablaphi) const{
		 Jacobian_type(vt, ig, Weight, phi, gradphi, nablaphi);
		}

  void Jacobian(const vector < vector < double > > &vt,const unsigned &ig, double &Weight,
		vector < double > &phi, vector < double > &gradphi, vector < double > &nablaphi) const{
		Jacobian_type(vt, ig, Weight, phi, gradphi, nablaphi);
		}

  template <class type>
  void JacobianSur_type(const vector < vector < type > > &vt, const unsigned &ig, type &Weight,
			vector < double > &phi, vector < type > &gradphi, vector < type > &normal) const;

  void JacobianSur(const vector < vector < adept::adouble > > &vt, const unsigned &ig, adept::adouble &Weight,
	           vector < double > &phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &normal) const{
		   JacobianSur_type(vt, ig, Weight, phi, gradphi, normal);
		  }

  void JacobianSur(const vector < vector < double > > &vt, const unsigned &ig, double &Weight,
	           vector < double > &phi, vector < double > &gradphi, vector < double > &normal) const{
		     JacobianSur_type(vt, ig, Weight, phi, gradphi, normal);
		  }

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


  template <class type>
  void Jacobian_type( const vector < vector < type > > &vt,const unsigned &ig, type &Weight,
		      vector < double > &phi, vector < type > &gradphi, vector < type > &nablaphi) const;

  void Jacobian( const vector < vector < adept::adouble > > &vt,const unsigned &ig, adept::adouble &Weight,
 		 vector < double > &phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &nablaphi) const{
		 Jacobian_type(vt, ig, Weight, phi, gradphi, nablaphi);
		}

  void Jacobian(const vector < vector < double > > &vt,const unsigned &ig, double &Weight,
		vector < double > &phi, vector < double > &gradphi, vector < double > &nablaphi) const{
		Jacobian_type(vt, ig, Weight, phi, gradphi, nablaphi);
		}

  void JacobianSur(const vector < vector < adept::adouble > > &vt, const unsigned &ig, adept::adouble &Weight,
	           vector < double > &other_phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &normal) const{
		   std::cout<<"Jacobian surface non-defined for elem_type_3D objects"<<std::endl;
		   abort();
		  }

  void JacobianSur(const vector < vector < double > > &vt, const unsigned &ig, double &Weight,
	           vector < double > &other_phi, vector < double > &gradphi, vector < double > &normal) const{
		   std::cout<<"Jacobian surface non-defined for elem_type_3D objects"<<std::endl;
		   abort();
		  }


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


