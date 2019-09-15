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


      virtual void compute_normal(const std::vector< std::vector< double > > & Jac, std::vector< double > & normal) const = 0;
      
      virtual void compute_normal(const std::vector< std::vector< adept::adouble > > & Jac, std::vector< adept::adouble > & normal) const = 0;

      
//==============================
     /* adept */                        
     virtual void Jacobian_geometry(const vector < vector < adept::adouble > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector < adept::adouble > > & Jac,
                            std::vector < std::vector < adept::adouble > > & JacI,
                            adept::adouble & detJac,
                            const unsigned dim,
                            const unsigned space_dim) const = 0;
                            
     /* double */                        
     virtual void Jacobian_geometry(const vector < vector < double > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector < double > > & Jac,
                            std::vector < std::vector < double > > & JacI,
                            double & detJac,
                            const unsigned dim,
                            const unsigned space_dim) const = 0;
                            
                            
     /* adept-adept */                        
  virtual void JacobianSur_non_isoparametric(const elem_type * fe_elem_coords_in,
                                             const vector < vector < adept::adouble > > & vt,
                                             const unsigned & ig,
                                             adept::adouble & Weight,
                                             vector < double > & phi,
                                             vector < adept::adouble >   & gradphi,
                                             vector < adept::adouble >& normal,
                                             const unsigned dim,
                                             const unsigned space_dim) const = 0;
     
     /* adept-double */                        
  virtual void JacobianSur_non_isoparametric(const elem_type * fe_elem_coords_in,
                                             const vector < vector < double > > & vt,
                                             const unsigned & ig,
                                             double & Weight,
                                             vector < double > & phi, 
                                             vector < adept::adouble >   & gradphi,
                                              vector < double >& normal,
                                             const unsigned dim,
                                             const unsigned space_dim) const = 0;
         
     /* all double */                        
  virtual void JacobianSur_non_isoparametric(const elem_type * fe_elem_coords_in,
                                             const vector < vector < double > > & vt,
                                             const unsigned & ig,
                                             double & Weight,
                                             vector < double > & phi, 
                                             vector < double >   & gradphi,
                                             vector < double >& normal,
                                             const unsigned dim,
                                             const unsigned space_dim) const = 0;

//==============================
                                             
      /* adept-adept */                        
     virtual void Jacobian_non_isoparametric(const elem_type * fe_elem_coords_in,
                                             const vector < vector < adept::adouble > > & vt,
                                             const unsigned & ig,
                                             adept::adouble & Weight,
                                             vector < double > & phi,
                                             vector < adept::adouble >   & gradphi,
                                             boost::optional< vector < adept::adouble > & > nablaphi,
                                             const unsigned dim,
                                             const unsigned space_dim) const = 0;
     
     /* adept-double */                        
     virtual void Jacobian_non_isoparametric(const elem_type * fe_elem_coords_in,
                                             const vector < vector < double > > & vt,
                                             const unsigned & ig,
                                             double & Weight,
                                             vector < double > & phi, 
                                             vector < adept::adouble >   & gradphi,
                                             boost::optional< vector < adept::adouble > & > nablaphi,
                                             const unsigned dim,
                                             const unsigned space_dim) const = 0;
         
     /* all double */                        
     virtual void Jacobian_non_isoparametric(const elem_type * fe_elem_coords_in,
                                             const vector < vector < double > > & vt,
                                             const unsigned & ig,
                                             double & Weight,
                                             vector < double > & phi, 
                                             vector < double >   & gradphi,
                                             boost::optional< vector < double > & > nablaphi,
                                             const unsigned dim,
                                             const unsigned space_dim) const = 0;
     
//==============================
                                             
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
        std::cout << "GetDPhiDEta does not apply to this element dimension\n";
        abort();
      }

      /** To be Added */
      virtual double* GetDPhiDZeta(const unsigned& ig) const {
        std::cout << "GetDPhiDZeta does not apply to this element dimension\n";
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
      
     // member data
      static unsigned _refindex;

      void GetSparsityPatternSize(const LinearEquation& lspdef, const LinearEquation& lspdec, const int& ielc,
                                  NumericVector* NNZ_d, NumericVector* NNZ_o,
                                  const unsigned& index_sol, const unsigned& kkindex_sol) const;

      void GetSparsityPatternSize(const Mesh& meshf, const Mesh& meshc, const int& ielc, NumericVector* NNZ_d, NumericVector* NNZ_o, const char el_dofs[]) const;

      void GetSparsityPatternSize(const Mesh& Mesh, const int& iel, NumericVector* NNZ_d, NumericVector* NNZ_o, const unsigned& itype) const;

      static const unsigned _fe_old_to_new[QL];
 
      static const int _fe_new_to_old[NFE_FAMS];
  
      virtual void VolumeShapeAtBoundary(const vector < vector < double > > &vt, const vector < vector < double> > & vt_bdry,  const unsigned& jface, const unsigned &ig, vector < double > &phi, vector < double > &gradphi) const = 0;

      
      basis* GetBasis() const {
        return _pt_basis;
      }


   protected:

      // ====================================
      // member functions
      // ====================================
      
      virtual const basis* set_FE_family_and_linear_element(const char* geom_elem, unsigned int FEType_in) = 0;

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
      
      virtual void allocate_and_fill_shape_at_quadrature_points() = 0;
      
      virtual void deallocate_shape_at_quadrature_points() = 0; /*they say you shouldn't call virtual function from constr/destr, so for now I am not using them in there*/
      
      virtual void allocate_and_fill_shape_at_quadrature_points_on_faces(const char* order_gauss) = 0;

      virtual void allocate_volume_shape_at_boundary_quadrature_points() = 0;

      virtual void deallocate_volume_shape_at_boundary_quadrature_points() = 0; /*they say you shouldn't call virtual function from constr/destr, so for now I am not using them in there*/
      
      // ====================================
      // member data
      // ====================================
      
      unsigned _dim; /* Spatial dimension of the geometric element */
      int _nc, _nf, _nlag[4];  /* _nc: number of dofs of 1 element;  _nf: number of dofs in that element after refinement; 
                                  _nlag[0] = number of linear dofs in 1 element;
                                  _nlag[1] = number of serendipity dofs in 1 element; 
                                  _nlag[2] = number of tensor-product quadratic dofs in 1 element; 
                                  _nlag[3] = number of tensor-product quadratic dofs in that element after 1 refinement; 
                                  */
      unsigned _SolType;       /* Finite Element Family flag */
      const double** _X;       /* [_nf][_dim] coordinates of the _nf nodes in the refined elements */ 
      const int** _IND;        /* [_nc][_dim] */ /*///@todo This is only used to evaluate the phi and derivatives */
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
          
          deallocate_shape_at_quadrature_points();
            
          deallocate_volume_shape_at_boundary_quadrature_points();
       }
      
      
     void compute_normal(const std::vector< std::vector< double > > & Jac, std::vector< double > & normal) const {
                 compute_normal< double >(Jac, normal);
     }
     
     void compute_normal(const std::vector< std::vector< adept::adouble > > & Jac, std::vector< adept::adouble > & normal) const {
                 compute_normal< adept::adouble >(Jac, normal);
     }

     void Jacobian_geometry(const vector < vector < adept::adouble > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector < adept::adouble > > & Jac,
                            std::vector < std::vector < adept::adouble > > & JacI,
                            adept::adouble & detJac,
                            const unsigned dim,
                            const unsigned space_dim) const {
         
              Jacobian_type_geometry(vt, ig, Jac, JacI, detJac, dim, space_dim);
     }
                            
     /* double */                        
     void Jacobian_geometry(const vector < vector < double > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector < double > > & Jac,
                            std::vector < std::vector < double > > & JacI,
                            double & detJac,
                            const unsigned dim,
                            const unsigned space_dim) const {
                                
              Jacobian_type_geometry(vt, ig, Jac, JacI, detJac, dim, space_dim);
     }
                            
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

  private:
      
      // ====================================
      // member data
      // ====================================
      
      double** _phi;
      double* _phi_memory;
      double** _dphidxi;
      double* _dphidxi_memory;

      double** _d2phidxi2;
      double* _d2phidxi2_memory;

      std::vector < std::vector < std::vector < double > > > _phiFace;
      std::vector < std::vector < std::vector < std::vector < double > > > > _gradPhiFace;
      std::vector < std::vector < std::vector < std::vector < std::vector < double > > > > > _hessianPhiFace;


      // ====================================
      // member functions
      // ====================================
      
      const basis* set_FE_family_and_linear_element(const char* geom_elem, unsigned int FEType_in);
      
      
      void allocate_and_fill_shape_at_quadrature_points();
      
      void deallocate_shape_at_quadrature_points();
      
      void allocate_and_fill_shape_at_quadrature_points_on_faces(const char* order_gauss);
      
      void allocate_volume_shape_at_boundary_quadrature_points();

      void deallocate_volume_shape_at_boundary_quadrature_points();

       virtual void VolumeShapeAtBoundary(const vector < vector < double > > &vt, const vector < vector < double> > & vt_bdry,  const unsigned& jface, const unsigned &ig, vector < double > &phi, vector < double > &gradphi) const { std::cout << "Not implemented"; abort(); };
      

   template <class type, class type_mov>
      void JacobianSur_type_non_isoparametric(const elem_type * fe_elem_coords_in,
                                              const vector < vector < type_mov > >& vt,
                                              const unsigned & ig, 
                                              type_mov & Weight,
                                              vector < double >& phi, 
                                              vector < type >& gradphi, 
                                              vector < type_mov >& normal,
                                              const unsigned dim,
                                              const unsigned space_dim) const;
                                              

    template <class type_mov>
      void compute_normal(const std::vector< std::vector< type_mov > > & Jac, std::vector< type_mov > & normal) const {
    
    // if you want to compute the normal to a 1d element, you need to know to what plane the boundary element belongs...                                   

    constexpr unsigned int space_dim = 3;
    
    //JacI ====================
    type_mov JacJacT[1][1]; JacJacT[0][0] = 0.; //1x1
    for (unsigned d = 0; d < space_dim; d++) JacJacT[0][0] += Jac[0][d]*Jac[0][d];
    type_mov detJac = sqrt(JacJacT[0][0]);
                                
   //   normal module, also equal to the transformation area....
   normal.resize(space_dim); ///@todo this must change based on how my domain is oriented
    
    normal[0] =  Jac[0][1] / detJac;
    normal[1] = -Jac[0][0] / detJac;
    normal[2] = 0.;

    //The derivative of x with respect to eta (dx/deta) has the opposite sign with respect to the normal
    //obtained as cross product between (dx/deta , dy/deta, 0) x (0,0,1)
//     (dx/deta , dy/deta, 0)  is the tangent vector (not normalized)
//     (0,0,1) is the unit vector going out of the plane
//     their cross product gives the (non-normalized) normal vector:    
//       i       , j      , k    
// det   dx/deta , dy/deta, 0    =   i (dy/deta)  -j (dx/deta) 
//        0      , 0      , 1    

// More: if you take the SCALAR TRIPLE PRODUCT of the (non-normalized) tangent, the unit normal and (0,0,1),
// that has the meaning of VOLUME, but since two vectors out of three have length 1,
// that is the same as taking the LENGTH of the segment.

//       n_x    , n_y     , 0 
// det   dx/deta , dy/deta, 0   =     n_x (dy/deta)  - n_y (dx/deta)   
//        0      , 0      , 1 
    
      }
      
      
    template <class type_mov>
     void JacobianSur_type_geometry(const vector < vector < type_mov > > & vt,
                               const unsigned & ig,
                               std::vector < std::vector <type_mov> > & Jac,
                               std::vector < std::vector <type_mov> > & JacI,
                               type_mov & detJac,
                               const unsigned dim,
                               const unsigned space_dim) const;

                               
    template <class type_mov>
     void Jacobian_type_geometry(const vector < vector < type_mov > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector <type_mov> > & Jac,
                            std::vector < std::vector <type_mov> > & JacI,
                            type_mov & detJac,
                            const unsigned dim,
                            const unsigned space_dim) const;
                            

    template <class type, class type_mov>
     void Jacobian_type_non_isoparametric(const elem_type * fe_elem_coords_in,
                                                const vector < vector < type_mov > > & vt,
                                                const unsigned & ig,
                                                type_mov & Weight,
                                                vector < double > & phi, 
                                                vector < type >   & gradphi,
                                                boost::optional< vector < type > & > nablaphi,
                                                const unsigned dim,
                                                const unsigned space_dim) const;
                                                

    /* adept-adept */                        
  virtual void JacobianSur_non_isoparametric(const elem_type * fe_elem_coords_in,
                                             const vector < vector < adept::adouble > > & vt,
                                             const unsigned & ig,
                                             adept::adouble & Weight,
                                             vector < double > & phi,
                                             vector < adept::adouble >   & gradphi,
                                             vector < adept::adouble >& normal,
                                             const unsigned dim,
                                             const unsigned space_dim) const {
                                              
      JacobianSur_type_non_isoparametric< adept::adouble, adept::adouble >(fe_elem_coords_in, vt, ig, Weight, phi, gradphi, normal, dim, space_dim);
                                                 
   }
     
     /* adept-double */                        
  virtual void JacobianSur_non_isoparametric(const elem_type * fe_elem_coords_in,
                                             const vector < vector < double > > & vt,
                                             const unsigned & ig,
                                             double & Weight,
                                             vector < double > & phi, 
                                             vector < adept::adouble >   & gradphi,
                                              vector < double >& normal,
                                             const unsigned dim,
                                             const unsigned space_dim) const  {
                                              
      JacobianSur_type_non_isoparametric< adept::adouble, double >(fe_elem_coords_in, vt, ig, Weight, phi, gradphi, normal, dim, space_dim);
                                                 
   }
         
     /* all double */                        
  virtual void JacobianSur_non_isoparametric(const elem_type * fe_elem_coords_in,
                                             const vector < vector < double > > & vt,
                                             const unsigned & ig,
                                             double & Weight,
                                             vector < double > & phi, 
                                             vector < double >   & gradphi,
                                             vector < double >& normal,
                                             const unsigned dim,
                                             const unsigned space_dim) const  {
                                              
      JacobianSur_type_non_isoparametric< double, double >(fe_elem_coords_in, vt, ig, Weight, phi, gradphi, normal, dim, space_dim);
                                                 
   }
   
                                             
 //================
                                             
/* adept-adept */                        
     void Jacobian_non_isoparametric(const elem_type * fe_elem_coords_in,
                                     const vector < vector < adept::adouble > > & vt,
                                     const unsigned & ig,
                                     adept::adouble & Weight,
                                     vector < double > & phi,
                                     vector < adept::adouble >   & gradphi,
                                     boost::optional< vector < adept::adouble > & > nablaphi,
                                     const unsigned dim,
                                     const unsigned space_dim) const {
                                              
        Jacobian_type_non_isoparametric< adept::adouble, adept::adouble >(fe_elem_coords_in, vt, ig, Weight, phi, gradphi, nablaphi, dim, space_dim);
                                              
      }
     
     /* adept-double */                        
     void Jacobian_non_isoparametric(const elem_type * fe_elem_coords_in,
                                     const vector < vector < double > > & vt,
                                     const unsigned & ig,
                                     double & Weight,
                                     vector < double > & phi, 
                                     vector < adept::adouble >   & gradphi,
                                     boost::optional< vector < adept::adouble > & > nablaphi,
                                     const unsigned dim,
                                     const unsigned space_dim) const {
                                                    
        Jacobian_type_non_isoparametric< adept::adouble, double >(fe_elem_coords_in, vt, ig, Weight, phi, gradphi, nablaphi, dim, space_dim);
                                        
       }
         
     /* all double */                        
     void Jacobian_non_isoparametric(const elem_type * fe_elem_coords_in,
                                     const vector < vector < double > > & vt,
                                     const unsigned & ig,
                                     double & Weight,
                                     vector < double > & phi, 
                                     vector < double >   & gradphi,
                                     boost::optional< vector < double > & > nablaphi,
                                     const unsigned dim,
                                     const unsigned space_dim) const {
                                                    
        Jacobian_type_non_isoparametric< double, double >(fe_elem_coords_in, vt, ig, Weight, phi, gradphi, nablaphi, dim, space_dim);
                                        
       }
      
  };


  class elem_type_2D : public elem_type
  {
      
    public:
        
      /** constructor */
      elem_type_2D(const char* solid, const char* order, const char* gauss_order);

      /** destructor */
      ~elem_type_2D() {
          
          deallocate_shape_at_quadrature_points();
            
          deallocate_volume_shape_at_boundary_quadrature_points();
       }

      
     void compute_normal(const std::vector< std::vector< double > > & Jac, std::vector< double > & normal) const {
                 compute_normal< double >(Jac, normal);
     }
     
     void compute_normal(const std::vector< std::vector< adept::adouble > > & Jac, std::vector< adept::adouble > & normal) const {
                 compute_normal< adept::adouble >(Jac, normal);
     }

     /* adept */                        
     void Jacobian_geometry(const vector < vector < adept::adouble > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector < adept::adouble > > & Jac,
                            std::vector < std::vector < adept::adouble > > & JacI,
                            adept::adouble & detJac,
                            const unsigned dim,
                            const unsigned space_dim) const {
         
              Jacobian_type_geometry(vt, ig, Jac, JacI, detJac, dim, space_dim);
     }
                            
     /* double */                        
     void Jacobian_geometry(const vector < vector < double > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector < double > > & Jac,
                            std::vector < std::vector < double > > & JacI,
                            double & detJac,
                            const unsigned dim,
                            const unsigned space_dim) const {
                                
              Jacobian_type_geometry(vt, ig, Jac, JacI, detJac, dim, space_dim);
     }
     
     
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

     
  private:
      
     const basis* set_FE_family_and_linear_element(const char* geom_elem, unsigned int FEType_in);
     
     void allocate_and_fill_shape_at_quadrature_points();
      
     void deallocate_shape_at_quadrature_points();
      
     void allocate_and_fill_shape_at_quadrature_points_on_faces(const char* order_gauss);
     
     void allocate_volume_shape_at_boundary_quadrature_points();

     void deallocate_volume_shape_at_boundary_quadrature_points();

      
  void VolumeShapeAtBoundary(const vector < vector < double > >& vt_vol, const vector < vector < double> > & vt_bdry,  const unsigned& jface, const unsigned& ig, vector < double >& phi, vector < double >& gradphi) const;

     
    template <class type_mov>
       void compute_normal(const std::vector< std::vector< type_mov > > & Jac, std::vector< type_mov > & normal) const {
           
//   normal  ===================
//     Cross product
//       i         , j        , k    
// det   d x/d xi  , d y/d xi , d z/d xi    =   i (dy/dxi dz/deta - dz/dxi dy/deta)  -j (d x/d xi  d z/d eta - d z/d xi  d x/d eta) + k (d x/d xi d y/d eta - d y/d xi d x/d eta)
//       d x/d eta , d y/d eta, d z/d eta    
    
    
    ///@todo How are we guaranteed that this normal is OUTWARDS??? For instance in 2d anticlockwise order of the edges guarantees outward normal. In 3d it must be the anticlockwise order of the edges of the boundary face taken from the volume... (if you look at the surface from outside, you must have dx/dxi and then dx/deta in anticlockwise order)
    const type_mov nx = Jac[0][1] * Jac[1][2] - Jac[1][1] * Jac[0][2];
    const type_mov ny = Jac[1][0] * Jac[0][2] - Jac[1][2] * Jac[0][0];
    const type_mov nz = Jac[0][0] * Jac[1][1] - Jac[1][0] * Jac[0][1];
    const type_mov invModn = 1. / sqrt(nx * nx + ny * ny + nz * nz);

    normal.resize(3);
    normal[0] = (nx) * invModn;
    normal[1] = (ny) * invModn;
    normal[2] = (nz) * invModn;

    
// ======== COMPUTATION of ELEMENT AREA as TRIPLE PRODUCT of two tangent vectors with UNIT normal vector 
//     Jac[2][0] = normal[0];
//     Jac[2][1] = normal[1];
//     Jac[2][2] = normal[2];
// 
//     //the determinant of the matrix is the area  ///@todo This is the triple scalar product of three vectors following the right-hand rule, so it is for sure positive
//     detJac = (Jac[0][0] * (Jac[1][1] * Jac[2][2] - Jac[1][2] * Jac[2][1]) +
//               Jac[0][1] * (Jac[1][2] * Jac[2][0] - Jac[1][0] * Jac[2][2]) +
//               Jac[0][2] * (Jac[1][0] * Jac[2][1] - Jac[1][1] * Jac[2][0]));


      }
      
      
      template <class type, class type_mov>
      void JacobianSur_type_non_isoparametric(const elem_type * fe_elem_coords_in,
                                              const vector < vector < type_mov > >& vt,
                                              const unsigned & ig, 
                                              type_mov & Weight,
                                              vector < double >& phi, 
                                              vector < type >& gradphi, 
                                              vector < type_mov >& normal,
                                                const unsigned dim,
                                                const unsigned space_dim) const;
                                                
                                                
     template <class type_mov>
     void JacobianSur_type_geometry(const vector < vector < type_mov > > & vt,
                               const unsigned & ig,
                               std::vector < std::vector <type_mov> > & Jac,
                               std::vector < std::vector <type_mov> > & JacI,
                               type_mov & detJac,
                               const unsigned dim,
                               const unsigned space_dim) const;
                               
                               
     template <class type_mov>
     void Jacobian_type_geometry(const vector < vector < type_mov > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector <type_mov> > & Jac,
                            std::vector < std::vector <type_mov> > & JacI,
                            type_mov & detJac,
                            const unsigned dim,
                            const unsigned space_dim) const;
                            
                            
     template <class type, class type_mov>
     void Jacobian_type_non_isoparametric(const elem_type * fe_elem_coords_in,
                                                const vector < vector < type_mov > > & vt,
                                                const unsigned & ig,
                                                type_mov & Weight,
                                                vector < double > & phi, 
                                                vector < type >   & gradphi,
                                                boost::optional< vector < type > & > nablaphi,
                                                const unsigned dim,
                                                const unsigned space_dim) const;                            
                                                  
     /* adept-adept */                        
  virtual void JacobianSur_non_isoparametric(const elem_type * fe_elem_coords_in,
                                             const vector < vector < adept::adouble > > & vt,
                                             const unsigned & ig,
                                             adept::adouble & Weight,
                                             vector < double > & phi,
                                             vector < adept::adouble >   & gradphi,
                                             vector < adept::adouble >& normal,
                                             const unsigned dim,
                                             const unsigned space_dim) const {
                                              
      JacobianSur_type_non_isoparametric< adept::adouble, adept::adouble >(fe_elem_coords_in, vt, ig, Weight, phi, gradphi, normal, dim, space_dim);
                                                 
   }
     
     /* adept-double */                        
  virtual void JacobianSur_non_isoparametric(const elem_type * fe_elem_coords_in,
                                             const vector < vector < double > > & vt,
                                             const unsigned & ig,
                                             double & Weight,
                                             vector < double > & phi, 
                                             vector < adept::adouble >   & gradphi,
                                              vector < double >& normal,
                                             const unsigned dim,
                                             const unsigned space_dim) const  {
                                              
      JacobianSur_type_non_isoparametric< adept::adouble, double >(fe_elem_coords_in, vt, ig, Weight, phi, gradphi, normal, dim, space_dim);
                                                 
   }
         
     /* all double */                        
  virtual void JacobianSur_non_isoparametric(const elem_type * fe_elem_coords_in,
                                             const vector < vector < double > > & vt,
                                             const unsigned & ig,
                                             double & Weight,
                                             vector < double > & phi, 
                                             vector < double >   & gradphi,
                                             vector < double >& normal,
                                             const unsigned dim,
                                             const unsigned space_dim) const  {
                                              
      JacobianSur_type_non_isoparametric< double, double >(fe_elem_coords_in, vt, ig, Weight, phi, gradphi, normal, dim, space_dim);
                                                 
   }
   
                                             
 //================
                                             
     /* adept-adept */     //adept, moving domain
     void Jacobian_non_isoparametric(const elem_type * fe_elem_coords_in,
                                          const vector < vector < adept::adouble > > & vt,
                                          const unsigned & ig,
                                          adept::adouble & Weight,
                                          vector < double > & phi,
                                          vector < adept::adouble >   & gradphi,
                                          boost::optional< vector < adept::adouble > & > nablaphi,
                                                const unsigned dim,
                                                const unsigned space_dim) const {
         
         Jacobian_type_non_isoparametric<adept::adouble, adept::adouble>(fe_elem_coords_in,vt,ig,Weight, phi, gradphi, nablaphi, dim, space_dim);
         
     }
     
     /* adept-double */    //adept, fixed domain
     void Jacobian_non_isoparametric(const elem_type * fe_elem_coords_in,
                                                const vector < vector < double > > & vt,
                                                const unsigned & ig,
                                                double & Weight,
                                                vector < double > & phi, 
                                                vector < adept::adouble >   & gradphi,
                                                boost::optional< vector < adept::adouble > & > nablaphi,
                                                const unsigned dim,
                                                const unsigned space_dim) const {
         
         Jacobian_type_non_isoparametric<adept::adouble, double>(fe_elem_coords_in,vt,ig,Weight, phi, gradphi, nablaphi, dim, space_dim);
         
     }
     
     /* all double */   //no adept, fixed domain 
     void Jacobian_non_isoparametric(const elem_type * fe_elem_coords_in,
                                                const vector < vector < double > > & vt,
                                                const unsigned & ig,
                                                double & Weight,
                                                vector < double > & phi, 
                                                vector < double >   & gradphi,
                                                boost::optional< vector < double > & > nablaphi,
                                                const unsigned dim,
                                                const unsigned space_dim) const {
         
         Jacobian_type_non_isoparametric<double, double>(fe_elem_coords_in,vt,ig,Weight, phi, gradphi, nablaphi, dim, space_dim);
         
     }



      
  void volume_shape_functions_at_reference_boundary_quadrature_points(const vector < vector < double> > & vt_bdry,  
                                           const unsigned jface) const;
                                           
      // ====================================
      // member data
      // ====================================
      
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
          
          deallocate_shape_at_quadrature_points();
            
          deallocate_volume_shape_at_boundary_quadrature_points();
      }
      
      
   
     void compute_normal(const std::vector< std::vector< double > > & Jac, std::vector< double > & normal) const {
           std::cout << "Normal non-defined for 3D objects" << std::endl; abort();
     }
     
     void compute_normal(const std::vector< std::vector< adept::adouble > > & Jac, std::vector< adept::adouble > & normal) const {
           std::cout << "Normal non-defined for 3D objects" << std::endl; abort();
     }
      
     /* adept */                        
     void Jacobian_geometry(const vector < vector < adept::adouble > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector < adept::adouble > > & Jac,
                            std::vector < std::vector < adept::adouble > > & JacI,
                            adept::adouble & detJac,
                            const unsigned dim,
                            const unsigned space_dim) const {
         
              Jacobian_type_geometry(vt, ig, Jac, JacI, detJac, dim, space_dim);
     }
                            
     /* double */                        
     void Jacobian_geometry(const vector < vector < double > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector < double > > & Jac,
                            std::vector < std::vector < double > > & JacI,
                            double & detJac,
                            const unsigned dim,
                            const unsigned space_dim) const {
                                
              Jacobian_type_geometry(vt, ig, Jac, JacI, detJac, dim, space_dim);
     }
     
     
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

     const basis* set_FE_family_and_linear_element(const char* geom_elem, unsigned int FEType_in);
     
     void allocate_and_fill_shape_at_quadrature_points();

     void deallocate_shape_at_quadrature_points();
     
     void allocate_and_fill_shape_at_quadrature_points_on_faces(const char* order_gauss);
     
     void allocate_volume_shape_at_boundary_quadrature_points();

     void deallocate_volume_shape_at_boundary_quadrature_points();

     void VolumeShapeAtBoundary(const vector < vector < double > >& vt_vol, const vector < vector < double> > & vt_bdry,  const unsigned& jface, const unsigned& ig, vector < double >& phi, vector < double >& gradphi) const;

     
     template <class type_mov>
     void Jacobian_type_geometry(const vector < vector < type_mov > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector <type_mov> > & Jac,
                            std::vector < std::vector <type_mov> > & JacI,
                            type_mov & detJac,
                            const unsigned dim,
                            const unsigned space_dim) const;
                            
                            
     template <class type, class type_mov>
     void Jacobian_type_non_isoparametric(const elem_type * fe_elem_coords_in,
                                          const vector < vector < type_mov > > & vt,
                                          const unsigned & ig,
                                          type_mov & Weight,
                                          vector < double > & phi,
                                          vector < type >   & gradphi,
                                          boost::optional< vector < type > & > nablaphi,
                                                const unsigned dim,
                                                const unsigned space_dim) const;
        
        
     /* adept-adept */                        
  virtual void JacobianSur_non_isoparametric(const elem_type * fe_elem_coords_in,
                                             const vector < vector < adept::adouble > > & vt,
                                             const unsigned & ig,
                                             adept::adouble & Weight,
                                             vector < double > & phi,
                                             vector < adept::adouble >   & gradphi,
                                             vector < adept::adouble >& normal,
                                             const unsigned dim,
                                             const unsigned space_dim) const {
                                              
        std::cout << "Jacobian surface non-defined for elem_type_3D objects" << std::endl;
        abort();
                                                 
   }
     
     /* adept-double */                        
  virtual void JacobianSur_non_isoparametric(const elem_type * fe_elem_coords_in,
                                             const vector < vector < double > > & vt,
                                             const unsigned & ig,
                                             double & Weight,
                                             vector < double > & phi, 
                                             vector < adept::adouble >   & gradphi,
                                              vector < double >& normal,
                                             const unsigned dim,
                                             const unsigned space_dim) const  {
                                              
        std::cout << "Jacobian surface non-defined for elem_type_3D objects" << std::endl;
        abort();
                                                 
   }
         
     /* all double */                        
  virtual void JacobianSur_non_isoparametric(const elem_type * fe_elem_coords_in,
                                             const vector < vector < double > > & vt,
                                             const unsigned & ig,
                                             double & Weight,
                                             vector < double > & phi, 
                                             vector < double >   & gradphi,
                                             vector < double >& normal,
                                             const unsigned dim,
                                             const unsigned space_dim) const  {
                                              
        std::cout << "Jacobian surface non-defined for elem_type_3D objects" << std::endl;
        abort();
                                                 
   }
   
                                             
 //================
                                             
     /* adept-adept */                        
     void Jacobian_non_isoparametric(const elem_type * fe_elem_coords_in,
                                          const vector < vector < adept::adouble > > & vt,
                                          const unsigned & ig,
                                          adept::adouble & Weight,
                                          vector < double > & phi,
                                          vector < adept::adouble >   & gradphi,
                                          boost::optional< vector < adept::adouble > & > nablaphi,
                                                const unsigned dim,
                                                const unsigned space_dim) const {
                                              
         Jacobian_type_non_isoparametric< adept::adouble, adept::adouble >(fe_elem_coords_in, vt, ig, Weight, phi, gradphi, nablaphi, dim, space_dim);
   
      }
     
     /* adept-double */                        
     void Jacobian_non_isoparametric(const elem_type * fe_elem_coords_in,
                                                const vector < vector < double > > & vt,
                                                const unsigned & ig,
                                                double & Weight,
                                                vector < double > & phi, 
                                                vector < adept::adouble >   & gradphi,
                                                boost::optional< vector < adept::adouble > & > nablaphi,
                                                const unsigned dim,
                                                const unsigned space_dim) const {
                                              
         Jacobian_type_non_isoparametric< adept::adouble, double >(fe_elem_coords_in, vt, ig, Weight, phi, gradphi, nablaphi, dim, space_dim);
   
      }
         
     /* all double */                        
     void Jacobian_non_isoparametric(const elem_type * fe_elem_coords_in,
                                                const vector < vector < double > > & vt,
                                                const unsigned & ig,
                                                double & Weight,
                                                vector < double > & phi, 
                                                vector < double >   & gradphi,
                                                boost::optional< vector < double > & > nablaphi,
                                                const unsigned dim,
                                                const unsigned space_dim) const {
                                              
         Jacobian_type_non_isoparametric< double, double >(fe_elem_coords_in, vt, ig, Weight, phi, gradphi, nablaphi, dim, space_dim);
   
      }
        
        
      // ====================================
      // member data
      // ====================================
        
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
      
        // values at boundary gauss points
      double **_phi_bdry;
      double *_phi_memory_bdry;
      double **_dphidxi_bdry;
      double *_dphidxi_memory_bdry;
      double **_dphideta_bdry;
      double *_dphideta_memory_bdry;
      double **_dphidzeta_bdry;
      double *_dphidzeta_memory_bdry;
      ///@todo do the same also in elem_type_1D
 
      
      
      
  };


//===========================================================
//============== elem_type_1D ===============================
//===========================================================
  
       template <class type, class type_mov>
      void elem_type_1D::JacobianSur_type_non_isoparametric(const elem_type * fe_elem_coords_in,
                                              const vector < vector < type_mov > >& vt,
                                              const unsigned & ig, 
                                              type_mov & weight,
                                              vector < double >& phi, 
                                              vector < type >& gradphi, 
                                              vector < type_mov >& normal,
                                              const unsigned dim,
                                              const unsigned space_dim) const  {
   
     std::vector < std::vector <type_mov> >  JacI;
     std::vector < std::vector <type_mov> >  Jac;
     type_mov detJac;

     const elem_type_1D *   fe_elem_coords_cast =  static_cast<const elem_type_1D*> (fe_elem_coords_in);
     
     fe_elem_coords_cast->JacobianSur_type_geometry<type_mov>(vt, ig, Jac, JacI, detJac, dim, space_dim);
     
     compute_normal<type_mov>(Jac, normal);

    // function part ====================
    weight = detJac * _gauss.GetGaussWeightsPointer()[ig];
    
    phi.resize(_nc);
    
    gradphi.resize(_nc * space_dim);  std::fill(gradphi.begin(),gradphi.end(),0.);
    const double* dxi = _dphidxi[ig];

    for(int inode = 0; inode < _nc; inode++, dxi++) {
        
      phi[inode] = _phi[ig][inode];

      for (unsigned d = 0; d < space_dim; d++) gradphi[ inode * space_dim + d] = (*dxi) * JacI[d][0];

    }
    
}



//return det without quadrature weight, pure geometric information
//     std::vector < std::vector <type_mov> > Jac(dim);
//     
//     for (unsigned d = 0; d < dim; d++) { Jac[d].resize(space_dim);	std::fill(Jac[d].begin(), Jac[d].end(), 0.); }
// 
//     for (unsigned d = 0; d < space_dim; d++) {
//     const double* dxi_coords  = _dphidxi[ig];
//        for(int inode = 0; inode < _nc; inode++, dxi_coords++) {
//           Jac[0][d] += (*dxi_coords) * vt[d][inode];
//         }
//      }
//    
//     type_mov JacJacT = 0.; //1x1
//     for (unsigned d = 0; d < space_dim; d++) JacJacT += Jac[0][d]*Jac[0][d];
//     
//     for (unsigned d = 0; d < space_dim; d++) JacI[d][0] = Jac[0][d] * 1. / JacJacT;
// 
//     detJac = sqrt(JacJacT)/*Jac[0][0]*/;  ///@todo in the old implementation shouldn't we take the absolute value??? I'd say we don't because it goes both on the lhs and on the rhs...

    template <class type_mov>
     void elem_type_1D::JacobianSur_type_geometry(const vector < vector < type_mov > > & vt,
                               const unsigned & ig,
                               std::vector < std::vector <type_mov> > & Jac,
                               std::vector < std::vector <type_mov> > & JacI,
                               type_mov & detJac,
                               const unsigned dim,
                               const unsigned space_dim) const {
                                   
    //Jac =====================
    Jac.resize(dim);
    for (unsigned d = 0; d < dim; d++) { Jac[d].resize(space_dim);	std::fill(Jac[d].begin(), Jac[d].end(), 0.); }

      for (unsigned d = 0; d < space_dim; d++) {
    const double* dxi_coords = _dphidxi[ig];
    for(int inode = 0; inode < _nc; inode++, dxi_coords++) {
      Jac[0][d] += (*dxi_coords) * vt[d][inode];
         }
      }
      
    
    //JacI ====================
    type_mov JacJacT[1][1]; JacJacT[0][0] = 0.; //1x1
    for (unsigned d = 0; d < space_dim; d++) JacJacT[0][0] += Jac[0][d]*Jac[0][d];
    detJac = sqrt(JacJacT[0][0]);

    JacI.resize(space_dim);
    for (unsigned d = 0; d < space_dim; d++) JacI[d].resize(dim);

    for (unsigned d = 0; d < space_dim; d++) JacI[d][0] = Jac[0][d] * 1. / JacJacT[0][0];

                               
    }
    
    
     template <class type_mov>
     void elem_type_1D::Jacobian_type_geometry(const vector < vector < type_mov > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector <type_mov> > & Jac,
                            std::vector < std::vector <type_mov> > & JacI,
                            type_mov & detJac,
                            const unsigned dim,
                            const unsigned space_dim) const {
//here the convention for the Jacobian is that the real coordinates are put along a COLUMN, so you have
 //    J = [ d x_1/ d xi  |  d x_2/ d xi  | d x_3 / d xi  ]     (1x3 matrix)
 // Then, if we denote differentials with D, we have
 //  [ D x_1 |  D x_2 | D x_3 ] =  D \xi [ d x_1/ d xi  |  d x_2/ d xi  | d x_3 / d xi  ]                                

 //  [ D x_1 |  D x_2 | D x_3 ] =   D \xi  J                               
                                
 //  [ D x_1 |  D x_2 | D x_3 ] J^T               =  D \xi  J  J^T                              

 //  [ D x_1 |  D x_2 | D x_3 ] J^T (J  J^T)^{-1} =  D \xi                              
                                
 //  [ D x_1 |  D x_2 | D x_3 ] | d xi / dx_1 | =  D \xi                              
 //                             | d xi / dx_2 |
 //                             | d xi / dx_3 |   
// 
// | d xi / dx_1 | 
// | d xi / dx_2 | =  J^T (J J^T)^{-1}
// | d xi / dx_3 |                              

                                
    //Jac =================
    Jac.resize(dim);
    
    for (unsigned d = 0; d < dim; d++) { 
        Jac[d].resize(space_dim);	std::fill(Jac[d].begin(), Jac[d].end(), 0.); }

    for (unsigned d = 0; d < space_dim; d++) {
    const double* dxi_coords  = _dphidxi[ig];
       for(int inode = 0; inode < _nc; inode++, dxi_coords++) {
          Jac[0][d] += (*dxi_coords) * vt[d][inode];
        }
     }
     
    //JacI  =================
    type_mov JacJacT[1][1];  JacJacT[0][0] = 0.; //1x1
    for (unsigned d = 0; d < space_dim; d++) JacJacT[0][0] += Jac[0][d]*Jac[0][d];
    detJac = sqrt(JacJacT[0][0]);
    
    JacI.resize(space_dim);
    for (unsigned d = 0; d < space_dim; d++) JacI[d].resize(dim);

    for (unsigned d = 0; d < space_dim; d++) JacI[d][0] = Jac[0][d] * 1. / JacJacT[0][0];

  ///@todo in the old implementation shouldn't we take the absolute value??? I'd say we don't because it goes both on the lhs and on the rhs...
             
     }
     
     
   template <class type, class type_mov>
     void elem_type_1D::Jacobian_type_non_isoparametric(const elem_type * fe_elem_coords_in,
                                                const vector < vector < type_mov > > & vt,
                                                const unsigned & ig,
                                                type_mov & Weight,
                                                vector < double > & phi, 
                                                vector < type >   & gradphi,
                                                boost::optional< vector < type > & > nablaphi,
                                                const unsigned dim,
                                                const unsigned space_dim) const {
                                                    
                                                    
     
// geometry part ================
     const elem_type_1D *   fe_elem_coords_cast =  static_cast<const elem_type_1D*> (fe_elem_coords_in);
     
     std::vector < std::vector <type_mov> >  JacI;
     std::vector < std::vector <type_mov> >  Jac;
     type_mov detJac;
   
     fe_elem_coords_cast->Jacobian_type_geometry<type_mov>(vt, ig, Jac, JacI, detJac, dim, space_dim);

// function part ================
    Weight = detJac * _gauss.GetGaussWeightsPointer()[ig];
    
    const double* dxi  = _dphidxi[ig];
    const double* dxi2 = _d2phidxi2[ig];

    phi.resize(_nc);
    gradphi.resize(_nc * space_dim);  std::fill(gradphi.begin(),gradphi.end(),0.);
    if(nablaphi) nablaphi->resize(_nc * space_dim);   ///@todo fix this: once space_dim was only 1

    
    for(int inode = 0; inode < _nc; inode++, dxi++, dxi2++) {

      phi[inode] = _phi[ig][inode];
      
      for (unsigned d = 0; d < space_dim; d++) gradphi[ inode * space_dim + d] = (*dxi) * JacI[d][0];

      if(nablaphi)(*nablaphi)[inode] = (*dxi2) * JacI[0][0] * JacI[0][0]; ///@todo fix this

    }

// function part - end ================


}

  
  
  
//===========================================================
//============== elem_type_2D ===============================
//===========================================================
  
  
     template <class type, class type_mov>
      void elem_type_2D::JacobianSur_type_non_isoparametric(const elem_type * fe_elem_coords_in,
                                              const vector < vector < type_mov > >& vt,
                                              const unsigned & ig, 
                                              type_mov & weight,
                                              vector < double >& phi, 
                                              vector < type >& gradphi, 
                                              vector < type_mov >& normal,
                                                const unsigned dim,
                                                const unsigned space_dim) const  {

     std::vector < std::vector <type_mov> >  JacI;
     std::vector < std::vector <type_mov> >  Jac;
     type_mov detJac;
     
     const elem_type_2D *   fe_elem_coords_cast =  static_cast<const elem_type_2D*> (fe_elem_coords_in);
     
     fe_elem_coords_cast->JacobianSur_type_geometry<type_mov>(vt, ig, Jac, JacI, detJac, dim, space_dim);
    
     compute_normal<type_mov>(Jac, normal);
     
    // function part ============
    weight = detJac * _gauss.GetGaussWeightsPointer()[ig];
    
    phi.resize(_nc);
    
    gradphi.resize(_nc * space_dim);  std::fill(gradphi.begin(),gradphi.end(),0.);
    const double* dxi  = _dphidxi[ig];
    const double* deta = _dphideta[ig];
    
    for(int inode = 0; inode < _nc; inode++, dxi++, deta++) {
        
      phi[inode] = _phi[ig][inode];
      
      for (unsigned d = 0; d < space_dim; d++) gradphi[ inode * space_dim + d] = (*dxi) * JacI[d][0] + (*deta) * JacI[d][1];
 
       }                                                  
                                                  
     }

     
     
//return det without quadrature weight, pure geometric information
     template <class type_mov>
     void elem_type_2D::JacobianSur_type_geometry(const vector < vector < type_mov > > & vt,
                               const unsigned & ig,
                               std::vector < std::vector <type_mov> > & Jac,
                               std::vector < std::vector <type_mov> > & JacI,
                               type_mov & detJac,
                               const unsigned dim,
                               const unsigned space_dim) const {
                
    //Jac ===================
    Jac.resize(2/*dim*/);
    for (unsigned d = 0; d < 2/*dim*/; d++) { Jac[d].resize(3/*space_dim*/);	std::fill(Jac[d].begin(), Jac[d].end(), 0.); }

    
  for (unsigned d = 0; d < 3/*space_dim*/; d++) {
      
    const double* dxi_coords = _dphidxi[ig];
    const double* deta_coords = _dphideta[ig];

    for(int inode = 0; inode < _nc; inode++, dxi_coords++, deta_coords++) {
      Jac[0][d] += (*dxi_coords)  * vt[d][inode];
      Jac[1][d] += (*deta_coords) * vt[d][inode];
    }
    
  }
  
  
    //JacI ===================
    type_mov JacJacT[2/*dim*/][2/*dim*/] = {{0., 0.}, {0., 0.}};
    type_mov JacJacT_inv[2/*dim*/][2/*dim*/] = {{0., 0.}, {0., 0.}};
    
    for (unsigned i = 0; i < 2/*dim*/; i++)
        for (unsigned j = 0; j < 2/*dim*/; j++)
            for (unsigned k = 0; k < space_dim; k++) JacJacT[i][j] += Jac[i][k]*Jac[j][k];
            
    type_mov detJacJacT = (JacJacT[0][0] * JacJacT[1][1] - JacJacT[0][1] * JacJacT[1][0]); ///@todo How are sure that this is POSITIVE?
            
    detJac = sqrt(abs(detJacJacT));
    
    JacJacT_inv[0][0] =  JacJacT[1][1] / detJacJacT;
    JacJacT_inv[0][1] = -JacJacT[0][1] / detJacJacT;
    JacJacT_inv[1][0] = -JacJacT[1][0] / detJacJacT;
    JacJacT_inv[1][1] =  JacJacT[0][0] / detJacJacT;
        
    JacI.resize(space_dim);
    for (unsigned d = 0; d < space_dim; d++) { JacI[d].resize(dim);  std::fill(JacI[d].begin(),JacI[d].end(),0.); }
    
    for (unsigned i = 0; i < space_dim; i++)
        for (unsigned j = 0; j < dim; j++)
            for (unsigned k = 0; k < dim; k++) JacI[i][j] += Jac[k][i]*JacJacT_inv[k][j];

    }
                            
                            
//return det without quadrature weight, pure geometric information
     template <class type_mov>
     void elem_type_2D::Jacobian_type_geometry(const vector < vector < type_mov > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector <type_mov> > & Jac,
                            std::vector < std::vector <type_mov> > & JacI,
                            type_mov & detJac,
                            const unsigned dim,
                            const unsigned space_dim) const {
                                                      
     //Jac ===============
    Jac.resize(dim);
    
    for (unsigned d = 0; d < dim; d++) { 
        Jac[d].resize(space_dim);	std::fill(Jac[d].begin(), Jac[d].end(), 0.); }

    for (unsigned d = 0; d < space_dim; d++) {
    const double* dxi_coords  = _dphidxi[ig];
    const double* deta_coords = _dphideta[ig];
    
      for(int inode = 0; inode < _nc; inode++, dxi_coords++, deta_coords++) {
          Jac[0][d] += (*dxi_coords)  * vt[d][inode];
          Jac[1][d] += (*deta_coords) * vt[d][inode];
        }
     }
     

     //JacI ===============
    type_mov JacJacT[2/*dim*/][2/*dim*/] = {{0., 0.}, {0., 0.}};
    type_mov JacJacT_inv[2/*dim*/][2/*dim*/] = {{0., 0.}, {0., 0.}};
    
    for (unsigned i = 0; i < 2/*dim*/; i++)
        for (unsigned j = 0; j < 2/*dim*/; j++)
            for (unsigned k = 0; k < space_dim; k++) JacJacT[i][j] += Jac[i][k]*Jac[j][k];
            
    type_mov detJacJacT = (JacJacT[0][0] * JacJacT[1][1] - JacJacT[0][1] * JacJacT[1][0]);
            
    detJac = sqrt(abs(detJacJacT));
    
    JacJacT_inv[0][0] =  JacJacT[1][1] / detJacJacT;
    JacJacT_inv[0][1] = -JacJacT[0][1] / detJacJacT;
    JacJacT_inv[1][0] = -JacJacT[1][0] / detJacJacT;
    JacJacT_inv[1][1] =  JacJacT[0][0] / detJacJacT;
        
    JacI.resize(space_dim);
    for (unsigned d = 0; d < space_dim; d++) { JacI[d].resize(dim);  std::fill(JacI[d].begin(),JacI[d].end(),0.); }
    
    for (unsigned i = 0; i < space_dim; i++)
        for (unsigned j = 0; j < dim; j++)
            for (unsigned k = 0; k < dim; k++) JacI[i][j] += Jac[k][i]*JacJacT_inv[k][j];

    
     }
     
                                                      
     template <class type, class type_mov>
     void elem_type_2D::Jacobian_type_non_isoparametric(const elem_type * fe_elem_coords_in,
                                                const vector < vector < type_mov > > & vt,
                                                const unsigned & ig,
                                                type_mov & Weight,
                                                vector < double > & phi, 
                                                vector < type >   & gradphi,
                                                boost::optional< vector < type > & > nablaphi,
                                                const unsigned dim,
                                                const unsigned space_dim) const {
                                                    

// geometry part ================
   const elem_type_2D *   fe_elem_coords_cast =  static_cast<const elem_type_2D*> (fe_elem_coords_in);                                                  

   std::vector < std::vector <type_mov> >  JacI;
   std::vector < std::vector <type_mov> >  Jac;
   type_mov detJac;
   
   fe_elem_coords_cast->Jacobian_type_geometry<type_mov>(vt, ig, Jac, JacI, detJac, dim, space_dim);

// function part ================
    Weight = detJac * _gauss.GetGaussWeightsPointer()[ig];
    
    phi.resize(_nc);
    
    gradphi.resize(_nc * space_dim);  std::fill(gradphi.begin(),gradphi.end(),0.);
    const double* dxi  = _dphidxi[ig];
    const double* deta = _dphideta[ig];
    
    if(nablaphi) nablaphi->resize(_nc * 3);
    const double* dxi2 = _d2phidxi2[ig];
    const double* deta2 = _d2phideta2[ig];
    const double* dxideta = _d2phidxideta[ig];

    
    for(int inode = 0; inode < _nc; inode++, dxi++, deta++, dxi2++, deta2++, dxideta++) {

      phi[inode] = _phi[ig][inode];

      for (unsigned d = 0; d < space_dim; d++) gradphi[ inode * space_dim + d] = (*dxi) * JacI[d][0] + (*deta) * JacI[d][1];

//       gradphi[inode * 2 + 0] = (*dxi) * JacI[0][0] + (*deta) * JacI[0][1];
//       gradphi[inode * 2 + 1] = (*dxi) * JacI[1][0] + (*deta) * JacI[1][1];

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



}


  
//===========================================================
//============== elem_type_3D ===============================
//===========================================================


//return det without quadrature weight, pure geometric information
     template <class type_mov>
     void elem_type_3D::Jacobian_type_geometry(const vector < vector < type_mov > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector <type_mov> > & Jac,
                            std::vector < std::vector <type_mov> > & JacI,
                            type_mov & detJac,
                            const unsigned dim,
                            const unsigned space_dim) const {
                                
     //Jac ===============
    Jac.resize(3/*dim*/);
    for (unsigned d = 0; d < 3/*dim*/; d++) {
        Jac[d].resize(3/*space_dim*/);   std::fill(Jac[d].begin(), Jac[d].end(), 0.); }
    
    for (unsigned d = 0; d < 3/*dim*/; d++) {
    const double * dxi_coords   = _dphidxi[ig];
    const double * deta_coords  = _dphideta[ig];
    const double * dzeta_coords = _dphidzeta[ig];

    for(int inode = 0; inode < _nc; inode++, dxi_coords++, deta_coords++, dzeta_coords++) {
      Jac[0][d] += (*dxi_coords)   * vt[d][inode];
      Jac[1][d] += (*deta_coords)  * vt[d][inode];
      Jac[2][d] += (*dzeta_coords) * vt[d][inode];
       }
    }
    
     //JacI ===============
    JacI.resize(space_dim);
    for (unsigned d = 0; d < space_dim; d++) JacI[d].resize(dim);
                                
    detJac = (Jac[0][0] * (Jac[1][1] * Jac[2][2] - Jac[1][2] * Jac[2][1]) +
                    Jac[0][1] * (Jac[1][2] * Jac[2][0] - Jac[1][0] * Jac[2][2]) +
                    Jac[0][2] * (Jac[1][0] * Jac[2][1] - Jac[1][1] * Jac[2][0]));

    JacI[0][0] = (-Jac[1][2] * Jac[2][1] + Jac[1][1] * Jac[2][2]) / detJac;
    JacI[0][1] = (Jac[0][2] * Jac[2][1] - Jac[0][1] * Jac[2][2]) / detJac;
    JacI[0][2] = (-Jac[0][2] * Jac[1][1] + Jac[0][1] * Jac[1][2]) / detJac;
    JacI[1][0] = (Jac[1][2] * Jac[2][0] - Jac[1][0] * Jac[2][2]) / detJac;
    JacI[1][1] = (-Jac[0][2] * Jac[2][0] + Jac[0][0] * Jac[2][2]) / detJac;
    JacI[1][2] = (Jac[0][2] * Jac[1][0] - Jac[0][0] * Jac[1][2]) / detJac;
    JacI[2][0] = (-Jac[1][1] * Jac[2][0] + Jac[1][0] * Jac[2][1]) / detJac;
    JacI[2][1] = (Jac[0][1] * Jac[2][0] - Jac[0][0] * Jac[2][1]) / detJac;
    JacI[2][2] = (-Jac[0][1] * Jac[1][0] + Jac[0][0] * Jac[1][1]) / detJac;

    }
     
     
     template <class type, class type_mov>
     void elem_type_3D::Jacobian_type_non_isoparametric(const elem_type * fe_elem_coords_in,
                                          const vector < vector < type_mov > > & vt,
                                          const unsigned & ig,
                                          type_mov & Weight,
                                          vector < double > & phi,
                                          vector < type >   & gradphi,
                                          boost::optional< vector < type > & > nablaphi,
                                                const unsigned dim,
                                                const unsigned space_dim) const {
                                                    
// geometry part ==============
     
   const elem_type_3D *   fe_elem_coords_cast =  static_cast<const elem_type_3D*> (fe_elem_coords_in);                                                  

   std::vector < std::vector <type_mov> >  JacI;
   std::vector < std::vector <type_mov> >  Jac;
   type_mov detJac;
   
   fe_elem_coords_cast->Jacobian_type_geometry<type_mov>(vt, ig, Jac, JacI, detJac, dim, space_dim);
// geometry part - end ==============

    
// function part ================
    Weight = detJac * _gauss.GetGaussWeightsPointer()[ig];

    const double* dxi = _dphidxi[ig];
    const double* deta = _dphideta[ig];
    const double* dzeta = _dphidzeta[ig];

    const double* dxi2 = _d2phidxi2[ig];
    const double* deta2 = _d2phideta2[ig];
    const double* dzeta2 = _d2phidzeta2[ig];
    const double* dxideta = _d2phidxideta[ig];
    const double* detadzeta = _d2phidetadzeta[ig];
    const double* dzetadxi = _d2phidzetadxi[ig];

    phi.resize(_nc);
    gradphi.resize(_nc * 3);  std::fill(gradphi.begin(),gradphi.end(),0.);
    if(nablaphi) nablaphi->resize(_nc * 6);
    
    
    for(int inode = 0; inode < _nc; inode++, dxi++, deta++, dzeta++, dxi2++, deta2++, dzeta2++, dxideta++, detadzeta++, dzetadxi++) {

      phi[inode] = _phi[ig][inode];

      gradphi[3 * inode + 0] = (*dxi) * JacI[0][0] + (*deta) * JacI[0][1] + (*dzeta) * JacI[0][2];
      gradphi[3 * inode + 1] = (*dxi) * JacI[1][0] + (*deta) * JacI[1][1] + (*dzeta) * JacI[1][2];
      gradphi[3 * inode + 2] = (*dxi) * JacI[2][0] + (*deta) * JacI[2][1] + (*dzeta) * JacI[2][2];

      if(nablaphi) {
        (*nablaphi)[6 * inode + 0] =
          ((*dxi2)    * JacI[0][0] + (*dxideta)  * JacI[0][1] + (*dzetadxi) * JacI[0][2]) * JacI[0][0] +
          ((*dxideta) * JacI[0][0] + (*deta2)    * JacI[0][1] + (*detadzeta) * JacI[0][2]) * JacI[0][1] +
          ((*dzetadxi) * JacI[0][0] + (*detadzeta) * JacI[0][1] + (*dzeta2)   * JacI[0][2]) * JacI[0][2];
        (*nablaphi)[6 * inode + 1] =
          ((*dxi2)    * JacI[1][0] + (*dxideta)  * JacI[1][1] + (*dzetadxi) * JacI[1][2]) * JacI[1][0] +
          ((*dxideta) * JacI[1][0] + (*deta2)    * JacI[1][1] + (*detadzeta) * JacI[1][2]) * JacI[1][1] +
          ((*dzetadxi) * JacI[1][0] + (*detadzeta) * JacI[1][1] + (*dzeta2)   * JacI[1][2]) * JacI[1][2];
        (*nablaphi)[6 * inode + 2] =
          ((*dxi2)    * JacI[2][0] + (*dxideta)  * JacI[2][1] + (*dzetadxi) * JacI[2][2]) * JacI[2][0] +
          ((*dxideta) * JacI[2][0] + (*deta2)    * JacI[2][1] + (*detadzeta) * JacI[2][2]) * JacI[2][1] +
          ((*dzetadxi) * JacI[2][0] + (*detadzeta) * JacI[2][1] + (*dzeta2)   * JacI[2][2]) * JacI[2][2];
        (*nablaphi)[6 * inode + 3] =
          ((*dxi2)    * JacI[0][0] + (*dxideta)  * JacI[0][1] + (*dzetadxi) * JacI[0][2]) * JacI[1][0] +
          ((*dxideta) * JacI[0][0] + (*deta2)    * JacI[0][1] + (*detadzeta) * JacI[0][2]) * JacI[1][1] +
          ((*dzetadxi) * JacI[0][0] + (*detadzeta) * JacI[0][1] + (*dzeta2)   * JacI[0][2]) * JacI[1][2];
        (*nablaphi)[6 * inode + 4] =
          ((*dxi2)    * JacI[1][0] + (*dxideta)  * JacI[1][1] + (*dzetadxi) * JacI[1][2]) * JacI[2][0] +
          ((*dxideta) * JacI[1][0] + (*deta2)    * JacI[1][1] + (*detadzeta) * JacI[1][2]) * JacI[2][1] +
          ((*dzetadxi) * JacI[1][0] + (*detadzeta) * JacI[1][1] + (*dzeta2)   * JacI[1][2]) * JacI[2][2];
        (*nablaphi)[6 * inode + 5] =
          ((*dxi2)    * JacI[2][0] + (*dxideta)  * JacI[2][1] + (*dzetadxi) * JacI[2][2]) * JacI[0][0] +
          ((*dxideta) * JacI[2][0] + (*deta2)    * JacI[2][1] + (*detadzeta) * JacI[2][2]) * JacI[0][1] +
          ((*dzetadxi) * JacI[2][0] + (*detadzeta) * JacI[2][1] + (*dzeta2)   * JacI[2][2]) * JacI[0][2];
      }
    }

// function part - end ================


}




} //end namespace femus



#endif


