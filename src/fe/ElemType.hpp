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
#include <boost/optional.hpp>


namespace femus
{

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
  class elem;
  class LinearEquation;

  class elem_type
  {

    public:

      /** constructor that receives Geometric Element and Gauss info */
      elem_type(const char* geom_elem, const char* fe_order, const char* order_gauss);

      /** destructor */
      virtual ~elem_type();

      /** To be Added */
      void BuildProlongation(const LinearEquation& lspdef, const LinearEquation& lspdec, const int& ielc, SparseMatrix* Projmat,
                             const unsigned& index_sol, const unsigned& kkindex_sol) const;

      /** To be Added */
      void BuildRestrictionTranspose(const LinearEquation& lspdef, const LinearEquation& lspdec, const int& ielc, SparseMatrix* Projmat,
                                     const unsigned& index_sol, const unsigned& kkindex_sol,
                                     const unsigned& index_pair_sol, const unsigned& kkindex_pair_sol) const;

      /** To be Added */
      void BuildProlongation(const Mesh& meshf, const Mesh& meshc, const int& ielc, SparseMatrix* Projmat, const char el_dofs[]) const;
      /** To be Added */
      void BuildProlongation(const Mesh& mymesh, const int& iel, SparseMatrix* Projmat, NumericVector* NNZ_d, NumericVector* NNZ_o, const unsigned& itype) const;


      virtual void GetJacobian(const vector < vector < adept::adouble > >& vt, const unsigned& ig, adept::adouble& Weight,
                               vector< vector < adept::adouble > >& jacobianMatrix) const = 0;

      virtual void GetJacobian(const vector < vector < double > >& vt, const unsigned& ig, double& Weight,
                               vector< vector < double > >& jacobianMatrix) const = 0;

      /* mixed adept-double */                        
      virtual void Jacobian(const vector < vector < adept::adouble > >& vt, const unsigned& ig, adept::adouble& Weight,
                            vector < double >& phi, vector < adept::adouble >& gradphi,
                            boost::optional < vector < adept::adouble > & > nablaphi = boost::none) const = 0;

      /* all double */                        
      virtual void Jacobian(const vector < vector < double > >& vt, const unsigned& ig, double& Weight,
                            vector < double >& other_phi, vector < double >& gradphi,
                            boost::optional < vector < double > & > nablaphi = boost::none) const = 0;

      /* Gauss-coordinate based - mixed adept-double */                        
      virtual void Jacobian(const vector < vector < adept::adouble > >& vt, const vector <double >& xi, adept::adouble& Weight,
                            vector < double >& phi, vector < adept::adouble >& gradphi,
                            boost::optional < vector < adept::adouble > & > nablaphi = boost::none) const = 0;

      /* Gauss-coordinate based - all double*/                        
      virtual void Jacobian(const vector < vector < double > >& vt, const vector <double >& xi, double& Weight,
                            vector < double >& other_phi, vector < double >& gradphi,
                            boost::optional < vector < double > & > nablaphi = boost::none) const = 0;

      virtual void JacobianSur(const vector < vector < adept::adouble > >& vt, const unsigned& ig, adept::adouble& Weight,
                               vector < double >& other_phi, vector < adept::adouble >& gradphi, vector < adept::adouble >& normal) const = 0;

      virtual void JacobianSur(const vector < vector < double > >& vt, const unsigned& ig, double& Weight,
                               vector < double >& other_phi, vector < double >& gradphi, vector < double >& normal) const = 0;
      /** To be Added */
      virtual double* GetPhi(const unsigned& ig) const = 0;

      /** To be Added */
      virtual double* GetDPhiDXi(const unsigned& ig) const = 0;

      /** To be Added */
      virtual double* GetDPhiDEta(const unsigned& ig) const {
        std::cout << "GetDPhiDEta does not apply for this element dimension\n";
        abort();
      }

      /** To be Added */
      virtual double* GetDPhiDZeta(const unsigned& ig) const {
        std::cout << "GetDPhiDZeta does not apply for this element dimension\n";
        abort();
      };

//   /** To be Added */
//   void GetArea(const double *vt,const double *vty, const double *vtz, const unsigned &ig,
//                double &Weight, double *other_phi) const;

      /** @deprecated  Function pointer for DPhiDXEZ */
      typedef double* (elem_type::*_FunctionPointer)(const unsigned& ig) const;   //you need "elem_type::" for some reason
      std::vector<_FunctionPointer> _DPhiXiEtaZetaPtr;

      /** @deprecated Evaluate shape functions at all quadrature points */
      virtual void EvaluateShapeAtQP(const std::string geomel_id_in, const std::string fe_in);

      /**  @deprecated Get shape functions */
      inline const double GetPhi(const uint qp, const uint dof) const {
        return _phi_mapGD[qp][dof];
      }

      /**  @deprecated Get shape function first derivatives */
      inline const double GetDPhiDxez(const uint qp, const uint dof) const {
        return _dphidxez_mapGD[qp][dof];
      }

      /** To be Added */
      inline const Gauss GetGaussRule() const {
        return _gauss;
      };

      /** To be Added */
      inline const Gauss* GetGaussRule_bdry() const {
        return _gauss_bdry;
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

      /** Set numbers of coarse and fine dofs for 1 element */
      void set_coarse_and_fine_elem_data(const basis* pt_basis_in);
      
      void allocate_and_set_IND(const basis* pt_basis_in);

      void allocate_coordinates_and_KVERT_IND();
      
      /** Set node coordinates and fine node indices */
      void set_coordinates_and_KVERT_IND(const basis* pt_basis_in);
      
      /** Compute node coordinates in basis object */
      void set_coordinates_in_Basis_object(basis* pt_basis_in, const basis* linearElement) const;
   
      /** Compute element prolongation operator */
      void set_element_prolongation(const basis* linearElement);
      
     // member data
      static unsigned _refindex;

      void GetSparsityPatternSize(const LinearEquation& lspdef, const LinearEquation& lspdec, const int& ielc,
                                  NumericVector* NNZ_d, NumericVector* NNZ_o,
                                  const unsigned& index_sol, const unsigned& kkindex_sol) const;

      void GetSparsityPatternSize(const Mesh& meshf, const Mesh& meshc, const int& ielc, NumericVector* NNZ_d, NumericVector* NNZ_o, const char el_dofs[]) const;

      void GetSparsityPatternSize(const Mesh& Mesh, const int& iel, NumericVector* NNZ_d, NumericVector* NNZ_o, const unsigned& itype) const;

      static const unsigned _fe_old_to_new[QL];
 
      static const int _fe_new_to_old[NFE_FAMS];
  
      virtual void VolumeShapeAtBoundary(const vector < vector < double > > &vt, const vector < vector < double> > & vt_bdry,  const unsigned& jface, const unsigned &ig, vector < double > &phi, vector < double > &gradphi) const {
           std::cout << "Implemented only for quad4 now" << std::endl; abort(); 
      }

      
      basis* GetBasis() const {
        return _pt_basis;
      }


   protected:

      // member data
      unsigned _dim; /* Spatial dimension of the geometric element */
      int _nc, _nf, _nlag[4];  /* _nc: number of dofs of 1 element;  _nf: number of dofs in that element after refinement; 
                                  _nlag[0] = number of linear dofs in 1 element;
                                  _nlag[1] = number of serendipity dofs in 1 element; 
                                  _nlag[2] = number of tensor-product quadratic dofs in 1 element; 
                                  _nlag[3] = number of tensor-product quadratic dofs in that element after 1 refinement; 
                                  */
      unsigned _SolType;       /* Finite Element Family flag */
      const double** _X;       /* [_nf][_dim] coordinates of the _nf nodes in the refined elements */ 
      const int** _IND;        /* [_nc][_dim] */
      const int** _KVERT_IND;  /* [_nf][2] For each _nf: 0 = id of the subdivision of the fine element, 1 = local id node on the subdivision of the fine element*/

      double** _prol_val;
      int** _prol_ind;
      double* _mem_prol_val;
      int* _mem_prol_ind;
      
      basis* _pt_basis;  /* FE basis functions*/

//  Gauss
      const Gauss _gauss;
            Gauss* _gauss_bdry;

      /**  @deprecated */
      bool isMpGDAllocated;
      double**      _phi_mapGD;
      double** _dphidxez_mapGD;

  };


  class elem_type_1D : public elem_type
  {

    public:
      /** constructor */
      elem_type_1D(const char* solid, const char* order, const char* gauss_order);

      /** destructor */
      ~elem_type_1D() {

        delete [] _phi;
        delete [] _phi_memory;
        delete [] _dphidxi;
        delete [] _dphidxi_memory;

        delete [] _d2phidxi2;
        delete [] _d2phidxi2_memory;

      };

      template <class type>
      void GetJacobian_type(const vector < vector < type > >& vt, const unsigned& ig, type& Weight,
                            vector < vector < type > >& jacobianMatrix) const;

      void GetJacobian(const vector < vector < adept::adouble > >& vt, const unsigned& ig, adept::adouble& Weight,
                       vector < vector < adept::adouble > >& jacobianMatrix) const {
        GetJacobian_type(vt, ig, Weight, jacobianMatrix);
      }
      void GetJacobian(const vector < vector < double > >& vt, const unsigned& ig, double& Weight,
                       vector < vector < double > >& jacobianMatrix) const {
        GetJacobian_type(vt, ig, Weight, jacobianMatrix);
      }


      /* all type minus a double */
      template <class type>
      void Jacobian_type(const vector < vector < type > >& vt, const unsigned& ig, type& Weight,
                         vector < double >& phi, vector < type >& gradphi,
                         boost::optional < vector < type > & > nablaphi) const;

      /* mixed adept - double */                        
      void Jacobian(const vector < vector < adept::adouble > >& vt, const unsigned& ig, adept::adouble& Weight,
                    vector < double >& phi, vector < adept::adouble >& gradphi,
                    boost::optional < vector < adept::adouble > & > nablaphi = boost::none) const {
        Jacobian_type(vt, ig, Weight, phi, gradphi, nablaphi);
      }
      
      /* all double */                        
      void Jacobian(const vector < vector < double > >& vt, const unsigned& ig, double& Weight,
                    vector < double >& phi, vector < double >& gradphi,
                    boost::optional < vector < double > & > nablaphi = boost::none) const {
        Jacobian_type(vt, ig, Weight, phi, gradphi, nablaphi);
      }

      /* Gauss-coordinate based */                        
      template <class type>
      void Jacobian_type(const vector < vector < type > >& vt, const vector < double >& xi, type& Weight,
                         vector < double >& phi, vector < type >& gradphi,
                         boost::optional < vector < type > & > nablaphi) const;

      void Jacobian(const vector < vector < adept::adouble > >& vt, const vector < double >& xi, adept::adouble& Weight,
                    vector < double >& phi, vector < adept::adouble >& gradphi,
                    boost::optional < vector < adept::adouble > & > nablaphi = boost::none) const {
        Jacobian_type(vt, xi, Weight, phi, gradphi, nablaphi);
      }
      void Jacobian(const vector < vector < double > >& vt, const vector < double >& xi, double& Weight,
                    vector < double >& phi, vector < double >& gradphi,
                    boost::optional < vector < double > & > nablaphi = boost::none) const {
        Jacobian_type(vt, xi, Weight, phi, gradphi, nablaphi);
      }


      template <class type>
      void JacobianSur_type(const vector < vector < type > >& vt, const unsigned& ig, type& Weight,
                            vector < double >& phi, vector < type >& gradphi, vector < type >& normal) const;

      void JacobianSur(const vector < vector < adept::adouble > >& vt, const unsigned& ig, adept::adouble& Weight,
                       vector < double >& phi, vector < adept::adouble >& gradphi, vector < adept::adouble >& normal) const {
        JacobianSur_type(vt, ig, Weight, phi, gradphi, normal);
      }

      void JacobianSur(const vector < vector < double > >& vt, const unsigned& ig, double& Weight,
                       vector < double >& phi, vector < double >& gradphi, vector < double >& normal) const {
        JacobianSur_type(vt, ig, Weight, phi, gradphi, normal);
      }

      inline double* GetPhi(const unsigned& ig) const {
        return _phi[ig];
      }
      inline double* GetDPhiDXi(const unsigned& ig) const {
        return _dphidxi[ig];
      }

      double** _phi;
      double* _phi_memory;
      double** _dphidxi;
      double* _dphidxi_memory;

      double** _d2phidxi2;
      double* _d2phidxi2_memory;

      std::vector < std::vector < std::vector < double > > > _phiFace;
      std::vector < std::vector < std::vector < std::vector < double > > > > _gradPhiFace;
      std::vector < std::vector < std::vector < std::vector < std::vector < double > > > > > _hessianPhiFace;

  };


  class elem_type_2D : public elem_type
  {
    public:
      /** constructor */
      elem_type_2D(const char* solid, const char* order, const char* gauss_order);

      /** destructor */
      ~elem_type_2D() {

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
         
        delete [] _phi_bdry;
        delete [] _phi_memory_bdry;
        delete [] _dphidxi_bdry;
        delete [] _dphidxi_memory_bdry;
        delete [] _dphideta_bdry;
        delete [] _dphideta_memory_bdry;

      };

      template <class type>
      void GetJacobian_type(const vector < vector < type > >& vt, const unsigned& ig, type& Weight,
                            vector < vector < type > >& jacobianMatrix) const;

      void GetJacobian(const vector < vector < adept::adouble > >& vt, const unsigned& ig, adept::adouble& Weight,
                       vector < vector < adept::adouble > >& jacobianMatrix) const {
        GetJacobian_type(vt, ig, Weight, jacobianMatrix);
      }
      void GetJacobian(const vector < vector < double > >& vt, const unsigned& ig, double& Weight,
                       vector < vector < double > >& jacobianMatrix) const {
        GetJacobian_type(vt, ig, Weight, jacobianMatrix);
      }


      /* Gauss index-based : all type minus a double */
      template <class type>
      void Jacobian_type(const vector < vector < type > >& vt, const unsigned& ig, type& Weight,
                         vector < double >& phi, vector < type >& gradphi,
                         boost::optional< vector < type > & > nablaphi) const;
                         
      /* mixed adept-double */
      void Jacobian(const vector < vector < adept::adouble > >& vt, const unsigned& ig, adept::adouble& Weight,
                    vector < double >& phi, vector < adept::adouble >& gradphi,
                    boost::optional< vector < adept::adouble > & > nablaphi = boost::none) const {
        Jacobian_type(vt, ig, Weight, phi, gradphi, nablaphi);
      }

      /* all double */
      void Jacobian(const vector < vector < double > >& vt, const unsigned& ig, double& Weight,
                    vector < double >& phi, vector < double >& gradphi,
                    boost::optional< vector < double > & > nablaphi = boost::none) const {
        Jacobian_type(vt, ig, Weight, phi, gradphi, nablaphi);
      }

      /* Gauss coordinate-based */
      template <class type>
      void Jacobian_type(const vector < vector < type > >& vt, const vector < double >& xi, type& Weight,
                         vector < double >& phi, vector < type >& gradphi,
                         boost::optional < vector < type > & > nablaphi) const;

      /* mixed adept-double */
      void Jacobian(const vector < vector < adept::adouble > >& vt, const vector < double >& xi, adept::adouble& Weight,
                    vector < double >& phi, vector < adept::adouble >& gradphi,
                    boost::optional < vector < adept::adouble > & > nablaphi = boost::none) const {
        Jacobian_type(vt, xi, Weight, phi, gradphi, nablaphi);
      }
      
      /* all double */
      void Jacobian(const vector < vector < double > >& vt, const vector < double >& xi, double& Weight,
                    vector < double >& phi, vector < double >& gradphi,
                    boost::optional < vector < double > & > nablaphi = boost::none) const {
        Jacobian_type(vt, xi, Weight, phi, gradphi, nablaphi);
      }

      template <class type>
      void JacobianSur_type(const vector < vector < type > >& vt, const unsigned& ig, type& Weight,
                            vector < double >& phi, vector < type >& gradphi, vector < type >& normal) const;

      void JacobianSur(const vector < vector < adept::adouble > >& vt, const unsigned& ig, adept::adouble& Weight,
                       vector < double >& phi, vector < adept::adouble >& gradphi, vector < adept::adouble >& normal) const {
        JacobianSur_type(vt, ig, Weight, phi, gradphi, normal);
      }

      void JacobianSur(const vector < vector < double > >& vt, const unsigned& ig, double& Weight,
                       vector < double >& phi, vector < double >& gradphi, vector < double >& normal) const {
        JacobianSur_type(vt, ig, Weight, phi, gradphi, normal);
      }

      inline double* GetPhi(const unsigned& ig) const {
        return _phi[ig];
      }
      inline double* GetDPhiDXi(const unsigned& ig) const {
        return _dphidxi[ig];
      }
      inline double* GetDPhiDEta(const unsigned& ig) const {
        return _dphideta[ig];
      }

     void VolumeShapeAtBoundary(const vector < vector < double > >& vt_vol, const vector < vector < double> > & vt_bdry,  const unsigned& jface, const unsigned& ig, vector < double >& phi, vector < double >& gradphi) const;

     template <class type>
     void Jacobian_at_point(const vector < vector < double > >& vt, const vector < double >& pos_in, vector < double >& pos_out) const;

     template <class type>
     void Jacobian_type_non_isoparametric(const elem_type_2D * fe_elem_coords,
                                                const vector < vector < type > > & vt,
                                                const unsigned & ig,
                                                type & Weight,
                                                vector < type > & phi, 
                                                vector < type >   & gradphi,
                                                boost::optional< vector < type > & > nablaphi) const {
                                                    
// geometry part ==============

    type Jac[2][2] = {{0, 0}, {0, 0}};
    type JacI[2][2];
    const double* dxi_coords  = fe_elem_coords->_dphidxi[ig];
    const double* deta_coords = fe_elem_coords->_dphideta[ig];

    for(int inode = 0; inode < fe_elem_coords->_nc; inode++, dxi_coords++, deta_coords++) {
      Jac[0][0] += (*dxi_coords) * vt[0][inode];
      Jac[0][1] += (*dxi_coords) * vt[1][inode];
      Jac[1][0] += (*deta_coords) * vt[0][inode];
      Jac[1][1] += (*deta_coords) * vt[1][inode];
    }

    type det = (Jac[0][0] * Jac[1][1] - Jac[0][1] * Jac[1][0]);

    JacI[0][0] =  Jac[1][1] / det;
    JacI[0][1] = -Jac[0][1] / det;
    JacI[1][0] = -Jac[1][0] / det;
    JacI[1][1] =  Jac[0][0] / det;

    Weight = det * _gauss.GetGaussWeightsPointer()[ig];

    
// function part ================
    
    const double* dxi  = _dphidxi[ig];
    const double* deta = _dphideta[ig];

    const double* dxi2 = _d2phidxi2[ig];
    const double* deta2 = _d2phideta2[ig];
    const double* dxideta = _d2phidxideta[ig];


    phi.resize(_nc);
    gradphi.resize(_nc * 2);
    if(nablaphi) nablaphi->resize(_nc * 3);

    
    for(int inode = 0; inode < _nc; inode++, dxi++, deta++, dxi2++, deta2++, dxideta++) {

      phi[inode] = _phi[ig][inode];

      gradphi[2 * inode + 0] = (*dxi) * JacI[0][0] + (*deta) * JacI[0][1];
      gradphi[2 * inode + 1] = (*dxi) * JacI[1][0] + (*deta) * JacI[1][1];

      if(nablaphi) {
        (*nablaphi)[3 * inode + 0] =
          ((*dxi2)   * JacI[0][0] + (*dxideta) * JacI[0][1]) * JacI[0][0] +
          ((*dxideta) * JacI[0][0] + (*deta2)  * JacI[0][1]) * JacI[0][1];
        (*nablaphi)[3 * inode + 1] =
          ((*dxi2)   * JacI[1][0] + (*dxideta) * JacI[1][1]) * JacI[1][0] +
          ((*dxideta) * JacI[1][0] + (*deta2)  * JacI[1][1]) * JacI[1][1];
        (*nablaphi)[3 * inode + 2] =
          ((*dxi2)   * JacI[0][0] + (*dxideta) * JacI[0][1]) * JacI[1][0] +
          ((*dxideta) * JacI[0][0] + (*deta2)  * JacI[0][1]) * JacI[1][1];
      }
    }



};


  private:
      double** _phi;
      double* _phi_memory;
      double** _dphidxi;
      double* _dphidxi_memory;
      double** _dphideta;
      double* _dphideta_memory;

      double** _d2phidxi2;
      double* _d2phidxi2_memory;
      double** _d2phideta2;
      double* _d2phideta2_memory;

      double** _d2phidxideta;
      double* _d2phidxideta_memory;

      std::vector < std::vector < std::vector < double > > > _phiFace;
      std::vector < std::vector < std::vector < std::vector < double > > > > _gradPhiFace;
      std::vector < std::vector < std::vector < std::vector < std::vector < double > > > > > _hessianPhiFace;
      
        // values at boundary gauss points
      double **_phi_bdry;
      double *_phi_memory_bdry;
      double **_dphidxi_bdry;
      double *_dphidxi_memory_bdry;
      double **_dphideta_bdry;
      double *_dphideta_memory_bdry;


  };

  class elem_type_3D : public elem_type
  {
    public:
      /** constructor */
      elem_type_3D(const char* solid, const char* order, const char* gauss_order);
      /** destructor */
      ~elem_type_3D() {
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
      void GetJacobian_type(const vector < vector < type > >& vt, const unsigned& ig, type& Weight,
                            vector < vector < type > >& jacobianMatrix) const;

      void GetJacobian(const vector < vector < adept::adouble > >& vt, const unsigned& ig, adept::adouble& Weight,
                       vector < vector < adept::adouble > >& jacobianMatrix) const {
        GetJacobian_type(vt, ig, Weight, jacobianMatrix);
      }
      void GetJacobian(const vector < vector < double > >& vt, const unsigned& ig, double& Weight,
                       vector < vector < double > >& jacobianMatrix) const {
        GetJacobian_type(vt, ig, Weight, jacobianMatrix);
      }

      /* templated: mixed type - double */
      template <class type>
      void Jacobian_type(const vector < vector < type > >& vt, const unsigned& ig, type& Weight,
                         vector < double >& phi, vector < type >& gradphi,
                         boost::optional< vector < type > & > nablaphi) const;

      /* mixed adept-double */
      void Jacobian(const vector < vector < adept::adouble > >& vt, const unsigned& ig, adept::adouble& Weight,
                    vector < double >& phi, vector < adept::adouble >& gradphi,
                    boost::optional< vector < adept::adouble > & > nablaphi = boost::none) const {
        Jacobian_type(vt, ig, Weight, phi, gradphi, nablaphi);
      }
      
      /* all double */
      void Jacobian(const vector < vector < double > >& vt, const unsigned& ig, double& Weight,
                    vector < double >& phi, vector < double >& gradphi,
                    boost::optional< vector < double > & > nablaphi = boost::none) const {
        Jacobian_type(vt, ig, Weight, phi, gradphi, nablaphi);
      }

      template <class type>
      void Jacobian_type(const vector < vector < type > >& vt, const vector < double >& xi, type& Weight,
                         vector < double >& phi, vector < type >& gradphi,
                         boost::optional < vector < type > & > nablaphi) const;

      void Jacobian(const vector < vector < adept::adouble > >& vt, const vector < double >& xi, adept::adouble& Weight,
                    vector < double >& phi, vector < adept::adouble >& gradphi,
                    boost::optional < vector < adept::adouble > & > nablaphi = boost::none) const {
        Jacobian_type(vt, xi, Weight, phi, gradphi, nablaphi);
      }
      void Jacobian(const vector < vector < double > >& vt, const vector < double >& xi, double& Weight,
                    vector < double >& phi, vector < double >& gradphi,
                    boost::optional < vector < double > & > nablaphi = boost::none) const {
        Jacobian_type(vt, xi, Weight, phi, gradphi, nablaphi);
      }

      void JacobianSur(const vector < vector < adept::adouble > >& vt, const unsigned& ig, adept::adouble& Weight,
                       vector < double >& other_phi, vector < adept::adouble >& gradphi, vector < adept::adouble >& normal) const {
        std::cout << "Jacobian surface non-defined for elem_type_3D objects" << std::endl;
        abort();
      }

      void JacobianSur(const vector < vector < double > >& vt, const unsigned& ig, double& Weight,
                       vector < double >& other_phi, vector < double >& gradphi, vector < double >& normal) const {
        std::cout << "Jacobian surface non-defined for elem_type_3D objects" << std::endl;
        abort();
      }


      //---------------------------------------------------------------------------------------------------------
      inline double* GetPhi(const unsigned& ig) const {
        return _phi[ig];
      }
      inline double* GetDPhiDXi(const unsigned& ig) const {
        return _dphidxi[ig];
      }
      inline double* GetDPhiDEta(const unsigned& ig) const {
        return _dphideta[ig];
      }
      inline double* GetDPhiDZeta(const unsigned& ig) const {
        return _dphidzeta[ig];
      }

    private:
      double** _phi;
      double* _phi_memory;
      double** _dphidxi;
      double* _dphidxi_memory;
      double** _dphideta;
      double* _dphideta_memory;
      double** _dphidzeta;
      double* _dphidzeta_memory;

      double** _d2phidxi2;
      double* _d2phidxi2_memory;
      double** _d2phideta2;
      double* _d2phideta2_memory;
      double** _d2phidzeta2;
      double* _d2phidzeta2_memory;

      double** _d2phidxideta;
      double* _d2phidxideta_memory;
      double** _d2phidetadzeta;
      double* _d2phidetadzeta_memory;
      double** _d2phidzetadxi;
      double* _d2phidzetadxi_memory;

      std::vector < std::vector < std::vector < double > > > _phiFace;
      std::vector < std::vector < std::vector < std::vector < double > > > > _gradPhiFace;
      std::vector < std::vector < std::vector < std::vector < std::vector < double > > > > > _hessianPhiFace;

  };




} //end namespace femus



#endif


