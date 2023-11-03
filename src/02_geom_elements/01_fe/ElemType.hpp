/* =========================================================================

 Program: FEMuS
 Module: ElemType
 Authors: Eugenio Aulisa, Giorgio Bornia

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

========================================================================= */


#ifndef __femus_fe_ElemType_hpp__
#define __femus_fe_ElemType_hpp__

#include "Basis.hpp"
#include "SparseMatrix.hpp"
#include "GaussPoints.hpp"
#include "FElemTypeEnum_list.hpp"
#include "GeomElTypeEnum.hpp"

#include "adept.h"
#include <boost/optional.hpp>


namespace femus
{


//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class elem;

class LinearEquation;
class Mesh;

  
      /** Class for Finite Element on 1 single Geometric Element 
       @todo must factorize this class a lot 
       It should have only abstract stuff about FE on 1 elem */
  class elem_type
  {

// =========================================
// ===  Constructors / Destructor - BEGIN =================
// =========================================
    public:

      /** constructor that receives Geometric Element, Finite Element and Gauss info */
      elem_type(const char* geom_elem, const char* fe_order, const char* order_gauss);

      /** constructor that receives Geometric Element, Finite Element */
      elem_type(const char* geom_elem, const char* fe_order);
  
      /** destructor */
      virtual ~elem_type();
// =========================================
// ===  Constructors / Destructor - END =================
// =========================================

// =========================================
// ===  Geometry - related - BEGIN =================
// =========================================
    public:
        
      /** Retrieve the dimension of the underlying geometric element */
      inline unsigned  GetDim() const {
        return _dim;
      }

   protected:
       
      void initialize_geom_elem(const char* geom_elem);
      
      unsigned _dim; /* Spatial dimension of the geometric element */
      
      GeomElType _GeomElemType;  /* Geometric Element flag */
// =========================================
// ===  Geometry - related - END =================
// =========================================

      
// =========================================
// ===   FE (without Quadrature evaluations) - BEGIN =================
// =========================================
    public:
        
      basis* GetBasis() const {
        return _pt_basis;
      }
      
      
   protected:
       
      void initialize_fe_and_multigrid_parts(const char* geom_elem);
      
      void deallocate_fe_and_multigrid_parts();
      
      void initialize_fe_soltype(const char* fe_order);
      
      /** Finite Element Family flag */
      unsigned _SolType;
      
      /** FE basis functions*/
      basis* _pt_basis;
      
      virtual const basis* set_current_FE_family_and_underlying_linear_FE_family(const char* geom_elem, unsigned int FEType_in) = 0;
      
      void allocate_and_set_coarse_node_indices(const basis* pt_basis_in);
      
      /** [_nc][_dim] */ /*///@todo This is only used to evaluate the phi and derivatives */
      const int** _IND;
// =========================================
// ===   FE (without Quadrature evaluations) - END =================
// =========================================


      
// =========================================
// ===  FE with MG - BEGIN =================
// =========================================
    public:
        
      /** Retrieve the number of dofs for this element */
      inline int  GetNDofs() const {
        return _nc;
      }
      
   protected:
       
      /** Compute element prolongation operator */
      void set_element_prolongation(const basis* linearElement);
      
      void allocate_fine_coordinates_and_KVERT_IND();
      
      /** Compute node coordinates in basis object */
      void set_fine_coordinates_in_Basis_object(basis* pt_basis_in, const basis* linearElement) const;
      
      /** Set fine node coordinates and fine node indices */
      void set_fine_coordinates_and_KVERT_IND(const basis* pt_basis_in);
      
      /** Set numbers of coarse and fine dofs for 1 element */
      void set_coarse_and_fine_num_dofs(const basis* pt_basis_in);
      
      /** _nc: number of dofs of 1 element;  */
      int _nc;
      
      /** _nf: number of dofs in the element after refinement; */
      int _nf;
      
      /**  _nlag[0] = number of linear dofs in 1 element;
           _nlag[1] = number of serendipity dofs in 1 element; 
           _nlag[2] = number of tensor-product quadratic dofs in 1 element; 
           _nlag[3] = number of tensor-product quadratic dofs in that element after 1 refinement; 
       */
      int _nlag[4];
      
      /** [_nf][_dim] coordinates of the _nf nodes in the refined elements ... @todo in what order? */ 
      const double** _X;
      
      /** [_nf][2] For each _nf: 0 = id of the subdivision of the fine element, 1 = local id node on the subdivision of the fine element */
      const int** _KVERT_IND;

      /** Prolongator value */
      double** _prol_val;
      
      /** Prolongator value, array for contiguous memory */
      double* _mem_prol_val;
      
      /** Prolongator index */
      int** _prol_ind;
 
      /** Prolongator index, array for contiguous memory */
      int* _mem_prol_ind;

// =========================================
// ===  FE with MG - END =================
// =========================================



// =========================================
// ===  Quadrature, without FE evaluations - BEGIN =================
// =========================================
  public:
      
      /** To be Added */
      inline const Gauss* GetGaussRule() const {
        return _gauss;
      }

      /** To be Added */
      inline double  GetGaussWeight(const unsigned ig) const {
        return GetGaussRule()->GetGaussWeightsPointer()[ig];
      }

      /** To be Added */
      inline unsigned GetGaussPointNumber() const {
        return GetGaussRule()->GetGaussPointsNumber();
      }
      
  protected:
      
      void initialize_quadrature_all(const char* geom_elem, const char* order_gauss);
     
      void deallocate_quadrature_all();
      
      void initialize_quadrature(const char* geom_elem, const char* order_gauss);
      
      void deallocate_quadrature();
      
      void initialize_quadrature_boundary(const char* geom_elem, const char* order_gauss);
      
      void deallocate_quadrature_boundary();
      
      Gauss * _gauss;
      ///@todo this must become a std::vector because for a Wedge there are 2 boundary quadrature rules, since there are 2 types of geom elems
      Gauss * _gauss_bdry; 
// =========================================
// ===  Quadrature, without FE evaluations - END =================
// =========================================

      
// =========================================
// ===  Quadrature, with FE evaluations - BEGIN =================
// =========================================
  public:
        
      virtual void GetJacobian(const std::vector < std::vector < adept::adouble > >& vt,
                               const unsigned& ig, 
                               adept::adouble& Weight,
                               std::vector< std::vector < adept::adouble > >& jacobianMatrix) const = 0;

      virtual void GetJacobian(const std::vector < std::vector < double > >& vt,
                               const unsigned& ig, double& Weight,
                               std::vector< std::vector < double > >& jacobianMatrix) const = 0;


                                             
      /* mixed adept-double */                        
      virtual void Jacobian(const std::vector < std::vector < adept::adouble > >& vt,
                            const unsigned& ig, 
                            adept::adouble& Weight,
                            std::vector < double >& phi,
                            std::vector < adept::adouble >& gradphi,
                            boost::optional < std::vector < adept::adouble > & > nablaphi = boost::none) const = 0;

      /* all double */                        
      virtual void Jacobian(const std::vector < std::vector < double > >& vt,
                            const unsigned& ig,
                            double& Weight,
                            std::vector < double >& other_phi,
                            std::vector < double >& gradphi,
                            boost::optional < std::vector < double > & > nablaphi = boost::none) const = 0;

      /* Gauss-coordinate based - mixed adept-double */                        
      virtual void Jacobian(const std::vector < std::vector < adept::adouble > >& vt,
                            const std::vector <double >& xi,
                            adept::adouble& Weight,
                            std::vector < double >& phi, 
                            std::vector < adept::adouble >& gradphi,
                            boost::optional < std::vector < adept::adouble > & > nablaphi = boost::none) const = 0;

      /* Gauss-coordinate based - all double*/                        
      virtual void Jacobian(const std::vector < std::vector < double > >& vt,
                            const std::vector <double >& xi, double& Weight,
                            std::vector < double >& other_phi, std::vector < double >& gradphi,
                            boost::optional < std::vector < double > & > nablaphi = boost::none) const = 0;

      virtual void JacobianSur(const std::vector < std::vector < adept::adouble > >& vt,
                               const unsigned& ig, 
                               adept::adouble& Weight,
                               std::vector < double >& other_phi,
                               std::vector < adept::adouble >& gradphi,
                               std::vector < adept::adouble >& normal) const = 0;

      virtual void JacobianSur(const std::vector < std::vector < double > >& vt,
                               const unsigned& ig,
                               double& Weight,
                               std::vector < double >& other_phi,
                               std::vector < double >& gradphi,
                               std::vector < double >& normal) const = 0;
                               
      /** To be Added */
      virtual double* GetPhi(const unsigned& ig) const = 0;
      
      virtual void GetPhi(std::vector<double> &phi, const std::vector < double >& xi ) const = 0;

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
      }

      
      virtual void fill_volume_shape_funcs_at_boundary_quadrature_points_on_current_elem(const std::vector < std::vector < double > > &vt,
                                                                                         const std::vector < std::vector < double> > & vt_bdry, 
                                                                                         const unsigned & jface, 
                                                                                         const unsigned & ig,
                                                                                         std::vector < double > & phi,
                                                                                         std::vector < double > & gradphi) const = 0;

      
     virtual  void initialize_quadrature_with_fe_evals_from_child(const char* geom_elem, const char* order_gauss) = 0;
      
     
   protected:
       
       
      void initialize_fe_quadrature_evaluations(const char* order_gauss);
      
      void initialize_to_null_fe_quadrature_evaluations_all();
      
      virtual void initialize_to_null_fe_quadrature_evaluations_vol_at_vol() = 0;
      
      virtual void initialize_to_null_fe_quadrature_evaluations_vol_at_bdry() = 0;
      
      virtual void allocate_and_fill_shape_at_quadrature_points() = 0;
      
      virtual void deallocate_shape_at_quadrature_points() = 0; /* you shouldn't call virtual function from constr/destr, so I am not using them in there*/
      
      virtual void allocate_and_fill_volume_shape_at_reference_boundary_quadrature_points_on_faces(const char* order_gauss) = 0;

      virtual void allocate_volume_shape_at_reference_boundary_quadrature_points_per_current_face() = 0;

      virtual void fill_volume_shape_at_reference_boundary_quadrature_points_per_face(/*const std::vector < std::vector < double> > & vt_bdry,  */const unsigned jface) const = 0;

      virtual void deallocate_volume_shape_at_reference_boundary_quadrature_points() = 0;
      
      ///  Quadrature - FE values at quadrature points on faces - for each face, for each Gauss point on the face, for each shape function of the volume
      std::vector < std::vector < std::vector < double > > > _phiFace;
      
      /// for each face, for each Gauss point on the face, for each shape function of the volume, for each reference direction (xi, eta, zeta)      
      std::vector < std::vector < std::vector < std::vector < double > > > > _gradPhiFace;
      
      /// for each face, for each Gauss point on the face, for each shape function of the volume, for each reference direction (xi, eta, zeta), for  each reference direction (xi, eta, zeta) again
      std::vector < std::vector < std::vector < std::vector < std::vector < double > > > > > _hessianPhiFace;

      // // //         std::vector < std::vector <  std::vector < double > > > _dphidxi_vol_at_bdry_templ;
      /// for every Direction, for every Quadrature Point, for every Dof  
                                                                                ///@todo unfortunately I cannot put this only once in the templated father, because it is not found...
                                                                                /// But, if I put it in this templated father, it works... Why?
      std::vector < std::vector <  std::vector < double > > > _dphidxi_templ;     
// =========================================
// ===  Quadrature, with FE evaluations - END =================
// =========================================
      
// =========================================
// ===  Equation, Sparsity pattern and Multigrid - BEGIN =================
// =========================================
    public:

      /** @todo move away from here */
      void GetSparsityPatternSize(const LinearEquation& lspdef, 
                                  const LinearEquation& lspdec,
                                  const int& ielc,
                                  NumericVector* NNZ_d,
                                  NumericVector* NNZ_o,
                                  const unsigned& index_sol, 
                                  const unsigned& kkindex_sol) const;

      /** @todo move away from here */
      void GetSparsityPatternSize(const Mesh& meshf,
                                  const Mesh& meshc,
                                  const int& ielc,
                                  NumericVector* NNZ_d,
                                  NumericVector* NNZ_o,
                                  const char el_dofs[]) const;

      /** for solution printing @todo move away from here */
      void GetSparsityPatternSize(const Mesh& Mesh,
                                  const int& iel, 
                                  NumericVector* NNZ_d,
                                  NumericVector* NNZ_o,
                                  const unsigned& itype) const;

        
      /** @todo move away from here */
      void BuildProlongation(const LinearEquation& lspdef,
                             const LinearEquation& lspdec,
                             const int& ielc, 
                             SparseMatrix* Projmat,
                             const unsigned& index_sol, 
                             const unsigned& kkindex_sol) const;

      /** @todo move away from here */
      void BuildProlongation(const Mesh& meshf,
                             const Mesh& meshc,
                             const int& ielc,
                             SparseMatrix* Projmat, 
                             const char el_dofs[]) const;
                             
      /** for solution printing @todo move away from here */
      void BuildProlongation(const Mesh& mesh,
                             const int& iel,
                             SparseMatrix* Projmat,
                             NumericVector* NNZ_d,
                             NumericVector* NNZ_o,
                             const unsigned& itype) const;

      /** @todo move away from here */
      void BuildRestrictionTranspose(const LinearEquation& lspdef,
                                     const LinearEquation& lspdec,
                                     const int& ielc,
                                     SparseMatrix* Projmat,
                                     const unsigned& index_sol,
                                     const unsigned& kkindex_sol,
                                     const unsigned& index_pair_sol,
                                     const unsigned& kkindex_pair_sol) const;
                                     
// =========================================
// ===  Equation, Sparsity pattern and Multigrid - END =================
// =========================================


  };
  
  


  class elem_type_1D : public elem_type
  {

    public:
        
// ===  Constructors / Destructor - BEGIN =================
      /** constructor */
      elem_type_1D(const char* geom_elem, const char* fe_order, const char* gauss_order);

      /** constructor */
      elem_type_1D(const char* geom_elem, const char* fe_order);

      /** destructor */
      ~elem_type_1D();      
// ===  Constructors / Destructor - END =================

// =========================================
// ===   FE (without Quadrature evaluations) - BEGIN =================
// =========================================
  protected:
      
      const basis* set_current_FE_family_and_underlying_linear_FE_family(const char* geom_elem, unsigned int FEType_in);
// ===   FE (without Quadrature evaluations) - END =================
      
                            
// =========================================
// ===  Quadrature, with FE evaluations - BEGIN =================
// =========================================
    public:
        
        
      template <class type>
      void GetJacobian_type(const std::vector < std::vector < type > >& vt, const unsigned& ig, type& Weight,
                            std::vector < std::vector < type > >& jacobianMatrix) const;

      void GetJacobian(const std::vector < std::vector < adept::adouble > >& vt,
                       const unsigned& ig,
                       adept::adouble& Weight,
                       std::vector < std::vector < adept::adouble > >& jacobianMatrix) const {
        GetJacobian_type(vt, ig, Weight, jacobianMatrix);
      }
      void GetJacobian(const std::vector < std::vector < double > >& vt,
                       const unsigned& ig,
                       double& Weight,
                       std::vector < std::vector < double > >& jacobianMatrix) const {
        GetJacobian_type(vt, ig, Weight, jacobianMatrix);
      }
      

      /* all type minus a double */
      template <class type>
      void Jacobian_type(const std::vector < std::vector < type > >& vt,
                         const unsigned& ig, 
                         type& Weight,
                         std::vector < double >& phi, 
                         std::vector < type >& gradphi,
                         boost::optional < std::vector < type > & > nablaphi) const;

      /* mixed adept - double */                        
      void Jacobian(const std::vector < std::vector < adept::adouble > >& vt,
                    const unsigned& ig, 
                    adept::adouble& Weight,
                    std::vector < double >& phi,
                    std::vector < adept::adouble >& gradphi,
                    boost::optional < std::vector < adept::adouble > & > nablaphi = boost::none) const {
        Jacobian_type(vt, ig, Weight, phi, gradphi, nablaphi);
      }
      
      /* all double */                        
      void Jacobian(const std::vector < std::vector < double > >& vt, 
                    const unsigned& ig, 
                    double& Weight,
                    std::vector < double >& phi,
                    std::vector < double >& gradphi,
                    boost::optional < std::vector < double > & > nablaphi = boost::none) const {
        Jacobian_type(vt, ig, Weight, phi, gradphi, nablaphi);
      }

      /* Gauss-coordinate based */                        
      template <class type>
      void Jacobian_type(const std::vector < std::vector < type > >& vt, 
                         const std::vector < double >& xi,
                         type& Weight,
                         std::vector < double >& phi, 
                         std::vector < type >& gradphi,
                         boost::optional < std::vector < type > & > nablaphi) const;

      void Jacobian(const std::vector < std::vector < adept::adouble > >& vt, 
                    const std::vector < double >& xi, 
                    adept::adouble& Weight,
                    std::vector < double >& phi, 
                    std::vector < adept::adouble >& gradphi,
                    boost::optional < std::vector < adept::adouble > & > nablaphi = boost::none) const {
        Jacobian_type(vt, xi, Weight, phi, gradphi, nablaphi);
      }
      void Jacobian(const std::vector < std::vector < double > >& vt, 
                    const std::vector < double >& xi, 
                    double& Weight,
                    std::vector < double >& phi, 
                    std::vector < double >& gradphi,
                    boost::optional < std::vector < double > & > nablaphi = boost::none) const {
        Jacobian_type(vt, xi, Weight, phi, gradphi, nablaphi);
      }


      template <class type>
      void JacobianSur_type(const std::vector < std::vector < type > >& vt,
                            const unsigned& ig,
                            type& Weight,
                            std::vector < double >& phi, 
                            std::vector < type >& gradphi, 
                            std::vector < type >& normal) const;

      void JacobianSur(const std::vector < std::vector < adept::adouble > >& vt, 
                       const unsigned& ig,
                       adept::adouble& Weight,
                       std::vector < double >& phi,
                       std::vector < adept::adouble >& gradphi,
                       std::vector < adept::adouble >& normal) const {
        JacobianSur_type(vt, ig, Weight, phi, gradphi, normal);
      }

      void JacobianSur(const std::vector < std::vector < double > >& vt,
                       const unsigned& ig,
                       double& Weight,
                       std::vector < double >& phi, 
                       std::vector < double >& gradphi, 
                       std::vector < double >& normal) const {
        JacobianSur_type(vt, ig, Weight, phi, gradphi, normal);
      }

      inline double* GetPhi(const unsigned& ig) const {
        return _phi[ig];
      }
      
      void GetPhi(std::vector<double> &phi, const std::vector < double >& xi ) const {
        phi.resize(_nc);
        for(unsigned i = 0; i < _nc; i++) {
          phi[i] = _pt_basis->eval_phi(_IND[i], &xi[0]);
        }
      }
      
      inline double* GetDPhiDXi(const unsigned& ig) const {
        return _dphidxi[ig];
      }

      void initialize_quadrature_with_fe_evals_from_child(const char* geom_elem, const char* order_gauss);
      
  protected:
      
       
      void fill_volume_shape_at_reference_boundary_quadrature_points_per_face(const unsigned  jface) const;
      
      void initialize_to_null_fe_quadrature_evaluations_vol_at_vol();
      
      void initialize_to_null_fe_quadrature_evaluations_vol_at_bdry();
      
      void allocate_and_fill_shape_at_quadrature_points();
      
      void deallocate_shape_at_quadrature_points();
      
      void allocate_and_fill_volume_shape_at_reference_boundary_quadrature_points_on_faces(const char* order_gauss);
      
      void allocate_volume_shape_at_reference_boundary_quadrature_points_per_current_face();

      void deallocate_volume_shape_at_reference_boundary_quadrature_points();

      void fill_volume_shape_funcs_at_boundary_quadrature_points_on_current_elem(const std::vector < std::vector < double > > &vt,
                                                                                 const std::vector < std::vector < double> > & vt_bdry,  
                                                                                 const unsigned& jface, 
                                                                                 const unsigned &ig,
                                                                                 std::vector < double > &phi,
                                                                                 std::vector < double > &gradphi) const { std::cout << "Not implemented"; abort(); }
            
      
      double** _phi;
      double*  _phi_memory;
      double** _dphidxi;
      double*  _dphidxi_memory;

      double** _d2phidxi2;
      double*  _d2phidxi2_memory;


        // values at boundary gauss points
      double ** _phi_vol_at_bdry;
      double *  _phi_vol_at_bdry_memory;
      double ** _dphidxi_vol_at_bdry;
      double *  _dphidxi_vol_at_bdry_memory;
// ===  Quadrature, with FE evaluations - END =================
      
  };


  class elem_type_2D : public elem_type
  {
      
// ===  Constructors / Destructor - BEGIN =================
    public:
        
      /** constructor */
      elem_type_2D(const char* geom_elem, const char* fe_order, const char* gauss_order);

      /** constructor */
      elem_type_2D(const char* geom_elem, const char* fe_order);

      /** destructor */
      ~elem_type_2D();
// ===  Constructors / Destructor - END =================

      
// =========================================
// ===   FE (without Quadrature evaluations) - BEGIN =================
// =========================================
    protected:
      
    const basis* set_current_FE_family_and_underlying_linear_FE_family(const char* geom_elem, unsigned int FEType_in);
// ===   FE (without Quadrature evaluations) - END =================

// =========================================
// ===  Quadrature, with FE evaluations - BEGIN =================
// =========================================
    public:
        
      template <class type>
      void GetJacobian_type(const std::vector < std::vector < type > >& vt, const unsigned& ig, type& Weight,
                            std::vector < std::vector < type > >& jacobianMatrix) const;

      void GetJacobian(const std::vector < std::vector < adept::adouble > >& vt, const unsigned& ig, adept::adouble& Weight,
                       std::vector < std::vector < adept::adouble > >& jacobianMatrix) const {
        GetJacobian_type(vt, ig, Weight, jacobianMatrix);
      }
      void GetJacobian(const std::vector < std::vector < double > >& vt, const unsigned& ig, double& Weight,
                       std::vector < std::vector < double > >& jacobianMatrix) const {
        GetJacobian_type(vt, ig, Weight, jacobianMatrix);
      }


      /* Gauss index-based : all type minus a double */
      template <class type>
      void Jacobian_type(const std::vector < std::vector < type > >& vt, const unsigned& ig, type& Weight,
                         std::vector < double >& phi, std::vector < type >& gradphi,
                         boost::optional< std::vector < type > & > nablaphi) const;
                         
      /* mixed adept-double */
      void Jacobian(const std::vector < std::vector < adept::adouble > >& vt, const unsigned& ig, adept::adouble& Weight,
                    std::vector < double >& phi, std::vector < adept::adouble >& gradphi,
                    boost::optional< std::vector < adept::adouble > & > nablaphi = boost::none) const {
        Jacobian_type(vt, ig, Weight, phi, gradphi, nablaphi);
      }

      /* all double */
      void Jacobian(const std::vector < std::vector < double > >& vt, const unsigned& ig, double& Weight,
                    std::vector < double >& phi, std::vector < double >& gradphi,
                    boost::optional< std::vector < double > & > nablaphi = boost::none) const {
        Jacobian_type(vt, ig, Weight, phi, gradphi, nablaphi);
      }

      /* Gauss coordinate-based */
      template <class type>
      void Jacobian_type(const std::vector < std::vector < type > >& vt, const std::vector < double >& xi, type& Weight,
                         std::vector < double >& phi, std::vector < type >& gradphi,
                         boost::optional < std::vector < type > & > nablaphi) const;

      /* mixed adept-double */
      void Jacobian(const std::vector < std::vector < adept::adouble > >& vt, const std::vector < double >& xi, adept::adouble& Weight,
                    std::vector < double >& phi, std::vector < adept::adouble >& gradphi,
                    boost::optional < std::vector < adept::adouble > & > nablaphi = boost::none) const {
        Jacobian_type(vt, xi, Weight, phi, gradphi, nablaphi);
      }
      
      /* all double */
      void Jacobian(const std::vector < std::vector < double > >& vt, const std::vector < double >& xi, double& Weight,
                    std::vector < double >& phi, std::vector < double >& gradphi,
                    boost::optional < std::vector < double > & > nablaphi = boost::none) const {
        Jacobian_type(vt, xi, Weight, phi, gradphi, nablaphi);
      }

      template <class type>
      void JacobianSur_type(const std::vector < std::vector < type > >& vt, const unsigned& ig, type& Weight,
                            std::vector < double >& phi, std::vector < type >& gradphi, std::vector < type >& normal) const;

      void JacobianSur(const std::vector < std::vector < adept::adouble > >& vt, const unsigned& ig, adept::adouble& Weight,
                       std::vector < double >& phi, std::vector < adept::adouble >& gradphi, std::vector < adept::adouble >& normal) const {
        JacobianSur_type(vt, ig, Weight, phi, gradphi, normal);
      }

      void JacobianSur(const std::vector < std::vector < double > >& vt, const unsigned& ig, double& Weight,
                       std::vector < double >& phi, std::vector < double >& gradphi, std::vector < double >& normal) const {
        JacobianSur_type(vt, ig, Weight, phi, gradphi, normal);
      }

      inline double* GetPhi(const unsigned& ig) const {
        return _phi[ig];
      }
      
      void GetPhi(std::vector<double> &phi, const std::vector < double >& xi ) const{
        phi.resize(_nc);
        for(unsigned i = 0; i < _nc; i++) {
          phi[i] = _pt_basis->eval_phi(_IND[i], &xi[0]);
        }
      }
      
      inline double* GetDPhiDXi(const unsigned& ig) const {
        return _dphidxi[ig];
      }
      inline double* GetDPhiDEta(const unsigned& ig) const {
        return _dphideta[ig];
      }

     
      void initialize_quadrature_with_fe_evals_from_child(const char* geom_elem, const char* order_gauss);
      
  protected:

     void fill_volume_shape_at_reference_boundary_quadrature_points_per_face(/*const std::vector < std::vector < double> > & vt_bdry,  */const unsigned jface) const;
                                           
      void initialize_to_null_fe_quadrature_evaluations_vol_at_vol();
      
      void initialize_to_null_fe_quadrature_evaluations_vol_at_bdry();
      
     void allocate_and_fill_shape_at_quadrature_points();
      
     void deallocate_shape_at_quadrature_points();
      
     void allocate_and_fill_volume_shape_at_reference_boundary_quadrature_points_on_faces(const char* order_gauss);
     
     void allocate_volume_shape_at_reference_boundary_quadrature_points_per_current_face();

     void deallocate_volume_shape_at_reference_boundary_quadrature_points();

      
     void fill_volume_shape_funcs_at_boundary_quadrature_points_on_current_elem(const std::vector < std::vector < double > >& vt_vol, const std::vector < std::vector < double> > & vt_bdry,  const unsigned& jface, const unsigned& ig, std::vector < double >& phi, std::vector < double >& gradphi) const;

      
      double** _phi;
      double*  _phi_memory;
      double** _dphidxi;
      double*  _dphidxi_memory;
      double** _dphideta;
      double*  _dphideta_memory;

      double** _d2phidxi2;
      double*  _d2phidxi2_memory;
      double** _d2phideta2;
      double*  _d2phideta2_memory;

      double** _d2phidxideta;
      double*  _d2phidxideta_memory;

      // values at boundary gauss points ///@todo probably remove
      double ** _phi_vol_at_bdry;
      double *  _phi_vol_at_bdry_memory;
      double ** _dphidxi_vol_at_bdry;
      double *  _dphidxi_vol_at_bdry_memory;
      double ** _dphideta_vol_at_bdry;
      double *  _dphideta_vol_at_bdry_memory;
// ===  Quadrature, with FE evaluations - END =================
      

  };
  
  
  

  class elem_type_3D : public elem_type
  {
      
// ===  Constructors / Destructor - BEGIN =================
    public:
        
      /** constructor */
      elem_type_3D(const char* geom_elem, const char* fe_order, const char* gauss_order);
      
      /** constructor */
      elem_type_3D(const char* geom_elem, const char* fe_order);

      /** destructor */
      ~elem_type_3D();
// ===  Constructors / Destructor - END =================
      
// =========================================
// ===   FE (without Quadrature evaluations) - BEGIN =================
// =========================================
    protected:
     
     const basis* set_current_FE_family_and_underlying_linear_FE_family(const char* geom_elem, unsigned int FEType_in);
// ===   FE (without Quadrature evaluations) - END =================
     
      
// =========================================
// ===  Quadrature, with FE evaluations - BEGIN =================
// =========================================
    public:
     
      template <class type>
      void GetJacobian_type(const std::vector < std::vector < type > >& vt, const unsigned& ig, type& Weight,
                            std::vector < std::vector < type > >& jacobianMatrix) const;

      void GetJacobian(const std::vector < std::vector < adept::adouble > >& vt, const unsigned& ig, adept::adouble& Weight,
                       std::vector < std::vector < adept::adouble > >& jacobianMatrix) const {
        GetJacobian_type(vt, ig, Weight, jacobianMatrix);
      }
      void GetJacobian(const std::vector < std::vector < double > >& vt, const unsigned& ig, double& Weight,
                       std::vector < std::vector < double > >& jacobianMatrix) const {
        GetJacobian_type(vt, ig, Weight, jacobianMatrix);
      }
     
                                                
      /* templated: mixed type - double */
      template <class type>
      void Jacobian_type(const std::vector < std::vector < type > >& vt, const unsigned& ig, type& Weight,
                         std::vector < double >& phi, std::vector < type >& gradphi,
                         boost::optional< std::vector < type > & > nablaphi) const;

      /* mixed adept-double */
      void Jacobian(const std::vector < std::vector < adept::adouble > >& vt, const unsigned& ig, adept::adouble& Weight,
                    std::vector < double >& phi, std::vector < adept::adouble >& gradphi,
                    boost::optional< std::vector < adept::adouble > & > nablaphi = boost::none) const {
        Jacobian_type(vt, ig, Weight, phi, gradphi, nablaphi);
      }
      
      /* all double */
      void Jacobian(const std::vector < std::vector < double > >& vt, const unsigned& ig, double& Weight,
                    std::vector < double >& phi, std::vector < double >& gradphi,
                    boost::optional< std::vector < double > & > nablaphi = boost::none) const {
        Jacobian_type(vt, ig, Weight, phi, gradphi, nablaphi);
      }

      template <class type>
      void Jacobian_type(const std::vector < std::vector < type > >& vt, const std::vector < double >& xi, type& Weight,
                         std::vector < double >& phi, std::vector < type >& gradphi,
                         boost::optional < std::vector < type > & > nablaphi) const;

      void Jacobian(const std::vector < std::vector < adept::adouble > >& vt, const std::vector < double >& xi, adept::adouble& Weight,
                    std::vector < double >& phi, std::vector < adept::adouble >& gradphi,
                    boost::optional < std::vector < adept::adouble > & > nablaphi = boost::none) const {
        Jacobian_type(vt, xi, Weight, phi, gradphi, nablaphi);
      }
      void Jacobian(const std::vector < std::vector < double > >& vt, const std::vector < double >& xi, double& Weight,
                    std::vector < double >& phi, std::vector < double >& gradphi,
                    boost::optional < std::vector < double > & > nablaphi = boost::none) const {
        Jacobian_type(vt, xi, Weight, phi, gradphi, nablaphi);
      }

      void JacobianSur(const std::vector < std::vector < adept::adouble > >& vt, const unsigned& ig, adept::adouble& Weight,
                       std::vector < double >& other_phi, std::vector < adept::adouble >& gradphi, std::vector < adept::adouble >& normal) const {
        std::cout << "Jacobian surface non-defined for elem_type_3D objects" << std::endl;
        abort();
      }

      void JacobianSur(const std::vector < std::vector < double > >& vt, const unsigned& ig, double& Weight,
                       std::vector < double >& other_phi, std::vector < double >& gradphi, std::vector < double >& normal) const {
        std::cout << "Jacobian surface non-defined for elem_type_3D objects" << std::endl;
        abort();
      }


      //---------------------------------------------------------------------------------------------------------
      inline double* GetPhi(const unsigned& ig) const {
        return _phi[ig];
      }
      
      void GetPhi(std::vector<double> &phi, const std::vector < double >& xi ) const {
        phi.resize(_nc);
        for(unsigned i = 0; i < _nc; i++) {
          phi[i] = _pt_basis->eval_phi(_IND[i], &xi[0]);
        }
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

      void initialize_quadrature_with_fe_evals_from_child(const char* geom_elem, const char* order_gauss);
      
    protected:

     void fill_volume_shape_at_reference_boundary_quadrature_points_per_face(const unsigned  jface) const;
     
      void initialize_to_null_fe_quadrature_evaluations_vol_at_vol();
      
      void initialize_to_null_fe_quadrature_evaluations_vol_at_bdry();
      
     void allocate_and_fill_shape_at_quadrature_points();

     void deallocate_shape_at_quadrature_points();
     
     void allocate_and_fill_volume_shape_at_reference_boundary_quadrature_points_on_faces(const char* order_gauss);
     
     void allocate_volume_shape_at_reference_boundary_quadrature_points_per_current_face();

     void deallocate_volume_shape_at_reference_boundary_quadrature_points();

     void fill_volume_shape_funcs_at_boundary_quadrature_points_on_current_elem(const std::vector < std::vector < double > >& vt_vol, const std::vector < std::vector < double> > & vt_bdry,  const unsigned& jface, const unsigned& ig, std::vector < double >& phi, std::vector < double >& gradphi) const;

     
      double** _phi;
      double*  _phi_memory;
      double** _dphidxi;
      double*  _dphidxi_memory;
      double** _dphideta;
      double*  _dphideta_memory;
      double** _dphidzeta;
      double*  _dphidzeta_memory;

      double** _d2phidxi2;
      double*  _d2phidxi2_memory;
      double** _d2phideta2;
      double*  _d2phideta2_memory;
      double** _d2phidzeta2;
      double*  _d2phidzeta2_memory;

      double** _d2phidxideta;
      double*  _d2phidxideta_memory;
      double** _d2phidetadzeta;
      double*  _d2phidetadzeta_memory;
      double** _d2phidzetadxi;
      double*  _d2phidzetadxi_memory;
      
        // values at boundary gauss points
      double ** _phi_vol_at_bdry;
      double *  _phi_vol_at_bdry_memory;
      double ** _dphidxi_vol_at_bdry;
      double *  _dphidxi_vol_at_bdry_memory;
      double ** _dphideta_vol_at_bdry;
      double *  _dphideta_vol_at_bdry_memory;
      double ** _dphidzeta_vol_at_bdry;
      double *  _dphidzeta_vol_at_bdry_memory;
// ===  Quadrature, with FE evaluations - END =================
 
  };

  
  
  
  
  //-----------elem_type_1D - BEGIN ----------------------------------------------------------------------------------------------

  template <class type>
  void elem_type_1D::GetJacobian_type(const std::vector < std::vector < type > >& vt, const unsigned& ig, type& Weight,
                                      std::vector < std::vector < type > >& jacobianMatrix) const
  {

    jacobianMatrix.resize(1);
    jacobianMatrix[1].resize(1);


    type Jac = 0.;

    const double* dxi = _dphidxi[ig];

    for(int inode = 0; inode < _nc; inode++, dxi++) {
      Jac += (*dxi) * vt[0][inode];
    }

    jacobianMatrix[0][0] = 1 / Jac;

    Weight = Jac * _gauss->GetGaussWeightsPointer()[ig];

  }


//---------------------------------------------------------------------------------------------------------

  template <class type>
  void elem_type_1D::Jacobian_type(const std::vector < std::vector < type > >& vt, const unsigned& ig, type& Weight,
                                   std::vector < double >& phi, std::vector < type >& gradphi,
                                   boost::optional< std::vector < type > & > nablaphi) const
  {

//    bool hermitianMatrix = true;
//     if(&nablaphi == NULL) {
//       hermitianMatrix = false;
//     }

    phi.resize(_nc);
    gradphi.resize(_nc * 1);
    if(nablaphi) nablaphi->resize(_nc * 1);

    type Jac = 0.;
    type JacI;

    const double* dxi = _dphidxi[ig];

    for(int inode = 0; inode < _nc; inode++, dxi++) {
      Jac += (*dxi) * vt[0][inode];
    }

    Weight = Jac * _gauss->GetGaussWeightsPointer()[ig];

    JacI = 1 / Jac;

    dxi = _dphidxi[ig];
    const double* dxi2 = _d2phidxi2[ig];

    for(int inode = 0; inode < _nc; inode++, dxi++, dxi2++) {
      phi[inode] = _phi[ig][inode];
      gradphi[inode] = (*dxi) * JacI;
      if(nablaphi)(*nablaphi)[inode] = (*dxi2) * JacI * JacI;
    }

  }


  //---------------------------------------------------------------------------------------------------------

  template <class type>
  void elem_type_1D::Jacobian_type(const std::vector < std::vector < type > >& vt, const std::vector<double>& xi, type& Weight,
                                   std::vector < double >& phi, std::vector < type >& gradphi,
                                   boost::optional < std::vector < type > & > nablaphi) const
  {

    phi.resize(_nc);
    gradphi.resize(_nc * 1);

    if(nablaphi) nablaphi->resize(_nc * 1);

    std::vector <double> dphidxi(_nc);
    std::vector <double> d2phidxi2(_nc);

    for(int j = 0; j < _nc; j++) {
      phi[j] = _pt_basis->eval_phi(_IND[j], &xi[0]);
      dphidxi[j] = _pt_basis->eval_dphidx(_IND[j], &xi[0]);
      d2phidxi2[j] = _pt_basis->eval_d2phidx2(_IND[j], &xi[0]);
    }

    type Jac = 0.;
    type JacI;

    const double* dxi = &dphidxi[0];

    for(int inode = 0; inode < _nc; inode++, dxi++) {
      Jac += (*dxi) * vt[0][inode];
    }

    Weight = Jac;

    JacI = 1 / Jac;

    dxi = &dphidxi[0];
    const double* dxi2 = &d2phidxi2[0];

    for(int inode = 0; inode < _nc; inode++, dxi++, dxi2++) {
      gradphi[inode] = (*dxi) * JacI;
      if(nablaphi) {
        (*nablaphi)[inode] = (*dxi2) * JacI * JacI;
      }
    }
  }

//---------------------------------------------------------------------------------------------------------



  template <class type>
  void elem_type_1D::JacobianSur_type(const std::vector < std::vector < type > >& vt, const unsigned& ig, type& Weight,
                                      std::vector < double >& phi, std::vector < type >& gradphi, std::vector < type >& normal) const
  {

    phi.resize(_nc);
    normal.resize(2);

    type Jac[2][2] = {{0., 0.}, {0., 0.}};
    type JacI[2][2];

    const double* dfeta = _dphidxi[ig];

    for(int inode = 0; inode < _nc; inode++, dfeta++) {
      Jac[0][0] += (*dfeta) * vt[0][inode];
      Jac[1][0] += (*dfeta) * vt[1][inode];
    }

//   normal module
    type modn = sqrt(Jac[0][0] * Jac[0][0] + Jac[1][0] * Jac[1][0]);

    normal[0] =  Jac[1][0] / modn;
    normal[1] = -Jac[0][0] / modn;

    //The derivative of x with respect to eta (dx/deta) has the opposite sign with respect to the normal
    //obtained as cross product between (dx/deta , dy/deta, 0) x (0,0,1)
    //The Jacobian has the structure
    // |dx/deta  -nx|
    // |dy/deta  -ny|
    Jac[0][1] = -normal[0];
    Jac[1][1] = -normal[1];

    //The determinant of that matrix is the area
    type det = (Jac[0][0] * Jac[1][1] - Jac[0][1] * Jac[1][0]);

    JacI[0][0] =  Jac[1][1] / det;
    JacI[0][1] = -Jac[0][1] / det;
    JacI[1][0] = -Jac[1][0] / det;
    JacI[1][1] =  Jac[0][0] / det;

    Weight = det * _gauss->GetGaussWeightsPointer()[ig];

    for(int inode = 0; inode < _nc; inode++) {
      phi[inode] = _phi[ig][inode];
    }

    ///@todo warning the surface gradient is missing!!!!!!!!!!!!!!!

}



  //-----------elem_type_1D - END  ----------------------------------------------------------------------------------------------


  //-----------elem_type_2D - BEGIN  ----------------------------------------------------------------------------------------------

  template <class type>
  void elem_type_2D::GetJacobian_type(const std::vector < std::vector < type > >& vt, const unsigned& ig, type& Weight,
                                      std::vector < std::vector < type > >& jacobianMatrix) const
  {

    jacobianMatrix.resize(2);
    jacobianMatrix[0].resize(2);
    jacobianMatrix[1].resize(2);

    type Jac[2][2] = {{0, 0}, {0, 0}};
    const double* dxi = _dphidxi[ig];
    const double* deta = _dphideta[ig];

    for(int inode = 0; inode < _nc; inode++, dxi++, deta++) {
      Jac[0][0] += (*dxi) * vt[0][inode];
      Jac[0][1] += (*dxi) * vt[1][inode];
      Jac[1][0] += (*deta) * vt[0][inode];
      Jac[1][1] += (*deta) * vt[1][inode];
    }

    type det = (Jac[0][0] * Jac[1][1] - Jac[0][1] * Jac[1][0]);

    jacobianMatrix[0][0] = Jac[1][1] / det;
    jacobianMatrix[0][1] = -Jac[0][1] / det;
    jacobianMatrix[1][0] = -Jac[1][0] / det;
    jacobianMatrix[1][1] = Jac[0][0] / det;

    Weight = det * _gauss->GetGaussWeightsPointer()[ig];

  }


//---------------------------------------------------------------------------------------------------------

  template <class type>
  void elem_type_2D::Jacobian_type(const std::vector < std::vector < type > >& vt, const unsigned& ig, type& Weight,
                                   std::vector < double >& phi, std::vector < type >& gradphi,
                                   boost::optional< std::vector < type > & > nablaphi) const
  {

//     bool hermitianMatrix = true;phi_x
//     if( &nablaphi == NULL ) {
//       hermitianMatrix = false;
//     }



    phi.resize(_nc);
    gradphi.resize(_nc * 2);
    if(nablaphi) nablaphi->resize(_nc * 3);

    type Jac[2][2] = {{0, 0}, {0, 0}};
    type JacI[2][2];
    const double* dxi = _dphidxi[ig];
    const double* deta = _dphideta[ig];

    for(int inode = 0; inode < _nc; inode++, dxi++, deta++) {
      Jac[0][0] += (*dxi) * vt[0][inode];
      Jac[0][1] += (*dxi) * vt[1][inode];
      Jac[1][0] += (*deta) * vt[0][inode];
      Jac[1][1] += (*deta) * vt[1][inode];
    }

    type det = (Jac[0][0] * Jac[1][1] - Jac[0][1] * Jac[1][0]);

    JacI[0][0] = Jac[1][1] / det;
    JacI[0][1] = -Jac[0][1] / det;
    JacI[1][0] = -Jac[1][0] / det;
    JacI[1][1] = Jac[0][0] / det;

    Weight = det * _gauss->GetGaussWeightsPointer()[ig];

    dxi = _dphidxi[ig];
    deta = _dphideta[ig];

    const double* dxi2 = _d2phidxi2[ig];
    const double* deta2 = _d2phideta2[ig];
    const double* dxideta = _d2phidxideta[ig];

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
  }

//---------------------------------------------------------------------------------------------------------

  template <class type>
  void elem_type_2D::Jacobian_type(const std::vector < std::vector < type > >& vt, const std::vector <double>& xi, type& Weight,
                                   std::vector < double >& phi, std::vector < type >& gradphi,
                                   boost::optional < std::vector < type > & > nablaphi) const
  {

    phi.resize(_nc);
    gradphi.resize(_nc * 2);
    if(nablaphi) nablaphi->resize(_nc * 3);

    std::vector <double> dphidxi(_nc);
    std::vector <double> dphideta(_nc);
    std::vector <double> d2phidxi2(_nc);
    std::vector <double> d2phideta2(_nc);
    std::vector <double> d2phidxideta(_nc);

    for(int j = 0; j < _nc; j++) {
      phi[j] = _pt_basis->eval_phi(_IND[j], &xi[0]);
      dphidxi[j] = _pt_basis->eval_dphidx(_IND[j], &xi[0]);
      dphideta[j] = _pt_basis->eval_dphidy(_IND[j], &xi[0]);
      d2phidxi2[j] = _pt_basis->eval_d2phidx2(_IND[j], &xi[0]);
      d2phideta2[j] = _pt_basis->eval_d2phidy2(_IND[j], &xi[0]);
      d2phidxideta[j] = _pt_basis->eval_d2phidxdy(_IND[j], &xi[0]);
    }

    type Jac[2][2] = {{0, 0}, {0, 0}};
    type JacI[2][2];
    const double* dxi = &dphidxi[0];
    const double* deta = &dphideta[0];

    for(int inode = 0; inode < _nc; inode++, dxi++, deta++) {
      Jac[0][0] += (*dxi) * vt[0][inode];
      Jac[0][1] += (*dxi) * vt[1][inode];
      Jac[1][0] += (*deta) * vt[0][inode];
      Jac[1][1] += (*deta) * vt[1][inode];
    }

    type det = (Jac[0][0] * Jac[1][1] - Jac[0][1] * Jac[1][0]);

    JacI[0][0] = Jac[1][1] / det;
    JacI[0][1] = -Jac[0][1] / det;
    JacI[1][0] = -Jac[1][0] / det;
    JacI[1][1] = Jac[0][0] / det;

    Weight = det;

    dxi = &dphidxi[0];
    deta = &dphideta[0];

    const double* dxi2 = &d2phidxi2[0];
    const double* deta2 = &d2phideta2[0];
    const double* dxideta = &d2phidxideta[0];

    for(int inode = 0; inode < _nc; inode++, dxi++, deta++, dxi2++, deta2++, dxideta++) {

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
  }

//---------------------------------------------------------------------------------------------------------



  template <class type>
  void elem_type_2D::JacobianSur_type(const std::vector < std::vector < type > >& vt, const unsigned& ig, type& Weight,
                                      std::vector < double >& phi, std::vector < type >& gradphi, std::vector < type >& normal) const
  {
    phi.resize(_nc);
    normal.resize(3);

    type Jac[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

    const double* dfx = _dphidxi[ig];
    const double* dfy = _dphideta[ig];

    for(int inode = 0; inode < _nc; inode++, dfx++, dfy++) {
      Jac[0][0] += (*dfx) * vt[0][inode];
      Jac[1][0] += (*dfx) * vt[1][inode];
      Jac[2][0] += (*dfx) * vt[2][inode];

      Jac[0][1] += (*dfy) * vt[0][inode];
      Jac[1][1] += (*dfy) * vt[1][inode];
      Jac[2][1] += (*dfy) * vt[2][inode];
    }

    //   normal module
    type nx = Jac[1][0] * Jac[2][1] - Jac[1][1] * Jac[2][0];
    type ny = Jac[0][1] * Jac[2][0] - Jac[2][1] * Jac[0][0];
    type nz = Jac[0][0] * Jac[1][1] - Jac[0][1] * Jac[1][0];
    type invModn = 1. / sqrt(nx * nx + ny * ny + nz * nz);

    normal[0] = (nx) * invModn;
    normal[1] = (ny) * invModn;
    normal[2] = (nz) * invModn;

    Jac[0][2] = normal[0];
    Jac[1][2] = normal[1];
    Jac[2][2] = normal[2];

    //the determinant of the matrix is the area
    type det = (Jac[0][0] * (Jac[1][1] * Jac[2][2] - Jac[1][2] * Jac[2][1]) +
                Jac[0][1] * (Jac[1][2] * Jac[2][0] - Jac[1][0] * Jac[2][2]) +
                Jac[0][2] * (Jac[1][0] * Jac[2][1] - Jac[1][1] * Jac[2][0]));

    Weight = det * _gauss->GetGaussWeightsPointer()[ig];

    for(int inode = 0; inode < _nc; inode++) {
      phi[inode] = _phi[ig][inode];
    }

    ///@todo warning the surface gradient is missing!!!!!!!!!!!!!!!
  }
  

  //-----------elem_type_2D - END  ----------------------------------------------------------------------------------------------
  
  
  //-----------elem_type_3D - BEGIN  ----------------------------------------------------------------------------------------------
   
 
  template <class type>
  void elem_type_3D::GetJacobian_type(const std::vector < std::vector < type > >& vt, const unsigned& ig, type& Weight,
                                      std::vector< std::vector < type > >& jacobianMatrix) const
  {

    jacobianMatrix.resize(3);
    jacobianMatrix[0].resize(3);
    jacobianMatrix[1].resize(3);
    jacobianMatrix[2].resize(3);

    type Jac[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

    const double* dxi = _dphidxi[ig];
    const double* deta = _dphideta[ig];
    const double* dzeta = _dphidzeta[ig];

    for(int inode = 0; inode < _nc; inode++, dxi++, deta++, dzeta++) {
      Jac[0][0] += (*dxi) * vt[0][inode];
      Jac[0][1] += (*dxi) * vt[1][inode];
      Jac[0][2] += (*dxi) * vt[2][inode];
      Jac[1][0] += (*deta) * vt[0][inode];
      Jac[1][1] += (*deta) * vt[1][inode];
      Jac[1][2] += (*deta) * vt[2][inode];
      Jac[2][0] += (*dzeta) * vt[0][inode];
      Jac[2][1] += (*dzeta) * vt[1][inode];
      Jac[2][2] += (*dzeta) * vt[2][inode];
    }

    type det = (Jac[0][0] * (Jac[1][1] * Jac[2][2] - Jac[1][2] * Jac[2][1]) +
                Jac[0][1] * (Jac[1][2] * Jac[2][0] - Jac[1][0] * Jac[2][2]) +
                Jac[0][2] * (Jac[1][0] * Jac[2][1] - Jac[1][1] * Jac[2][0]));

    jacobianMatrix[0][0] = (-Jac[1][2] * Jac[2][1] + Jac[1][1] * Jac[2][2]) / det;
    jacobianMatrix[0][1] = (Jac[0][2] * Jac[2][1] - Jac[0][1] * Jac[2][2]) / det;
    jacobianMatrix[0][2] = (-Jac[0][2] * Jac[1][1] + Jac[0][1] * Jac[1][2]) / det;
    jacobianMatrix[1][0] = (Jac[1][2] * Jac[2][0] - Jac[1][0] * Jac[2][2]) / det;
    jacobianMatrix[1][1] = (-Jac[0][2] * Jac[2][0] + Jac[0][0] * Jac[2][2]) / det;
    jacobianMatrix[1][2] = (Jac[0][2] * Jac[1][0] - Jac[0][0] * Jac[1][2]) / det;
    jacobianMatrix[2][0] = (-Jac[1][1] * Jac[2][0] + Jac[1][0] * Jac[2][1]) / det;
    jacobianMatrix[2][1] = (Jac[0][1] * Jac[2][0] - Jac[0][0] * Jac[2][1]) / det;
    jacobianMatrix[2][2] = (-Jac[0][1] * Jac[1][0] + Jac[0][0] * Jac[1][1]) / det;

    Weight = det * _gauss->GetGaussWeightsPointer()[ig];

  }

  
//---------------------------------------------------------------------------------------------------------
  template <class type>
  void elem_type_3D::Jacobian_type(const std::vector < std::vector < type > >& vt, const unsigned& ig, type& Weight,
                                   std::vector < double >& phi, std::vector < type >& gradphi,
                                   boost::optional< std::vector < type > & > nablaphi) const
  {

//     bool hermitianMatrix = true;
//     if(&nablaphi == NULL) {
//       hermitianMatrix = false;
//     }

    phi.resize(_nc);
    gradphi.resize(_nc * 3);
    if(nablaphi) nablaphi->resize(_nc * 6);


    type Jac[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    type JacI[3][3];

    const double* dxi = _dphidxi[ig];
    const double* deta = _dphideta[ig];
    const double* dzeta = _dphidzeta[ig];

    for(int inode = 0; inode < _nc; inode++, dxi++, deta++, dzeta++) {
      Jac[0][0] += (*dxi) * vt[0][inode];
      Jac[0][1] += (*dxi) * vt[1][inode];
      Jac[0][2] += (*dxi) * vt[2][inode];
      Jac[1][0] += (*deta) * vt[0][inode];
      Jac[1][1] += (*deta) * vt[1][inode];
      Jac[1][2] += (*deta) * vt[2][inode];
      Jac[2][0] += (*dzeta) * vt[0][inode];
      Jac[2][1] += (*dzeta) * vt[1][inode];
      Jac[2][2] += (*dzeta) * vt[2][inode];
    }

    type det = (Jac[0][0] * (Jac[1][1] * Jac[2][2] - Jac[1][2] * Jac[2][1]) +
                Jac[0][1] * (Jac[1][2] * Jac[2][0] - Jac[1][0] * Jac[2][2]) +
                Jac[0][2] * (Jac[1][0] * Jac[2][1] - Jac[1][1] * Jac[2][0]));

    JacI[0][0] = (-Jac[1][2] * Jac[2][1] + Jac[1][1] * Jac[2][2]) / det;
    JacI[0][1] = (Jac[0][2] * Jac[2][1] - Jac[0][1] * Jac[2][2]) / det;
    JacI[0][2] = (-Jac[0][2] * Jac[1][1] + Jac[0][1] * Jac[1][2]) / det;
    JacI[1][0] = (Jac[1][2] * Jac[2][0] - Jac[1][0] * Jac[2][2]) / det;
    JacI[1][1] = (-Jac[0][2] * Jac[2][0] + Jac[0][0] * Jac[2][2]) / det;
    JacI[1][2] = (Jac[0][2] * Jac[1][0] - Jac[0][0] * Jac[1][2]) / det;
    JacI[2][0] = (-Jac[1][1] * Jac[2][0] + Jac[1][0] * Jac[2][1]) / det;
    JacI[2][1] = (Jac[0][1] * Jac[2][0] - Jac[0][0] * Jac[2][1]) / det;
    JacI[2][2] = (-Jac[0][1] * Jac[1][0] + Jac[0][0] * Jac[1][1]) / det;

    Weight = det * _gauss->GetGaussWeightsPointer()[ig];

    dxi = _dphidxi[ig];
    deta = _dphideta[ig];
    dzeta = _dphidzeta[ig];

    const double* dxi2 = _d2phidxi2[ig];
    const double* deta2 = _d2phideta2[ig];
    const double* dzeta2 = _d2phidzeta2[ig];
    const double* dxideta = _d2phidxideta[ig];
    const double* detadzeta = _d2phidetadzeta[ig];
    const double* dzetadxi = _d2phidzetadxi[ig];

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

  }


  
//---------------------------------------------------------------------------------------------------------
  template <class type>
  void elem_type_3D::Jacobian_type(const std::vector < std::vector < type > >& vt, const std::vector <double>& xi, type& Weight,
                                   std::vector < double >& phi, std::vector < type >& gradphi,
                                   boost::optional < std::vector < type > & > nablaphi) const
  {

    phi.resize(_nc);
    gradphi.resize(_nc * 3);
    if(nablaphi) nablaphi->resize(_nc * 6);

    std::vector < double > dphidxi(_nc);
    std::vector < double > dphideta(_nc);
    std::vector < double > dphidzeta(_nc);

    std::vector < double > d2phidxi2(_nc);
    std::vector < double > d2phideta2(_nc);
    std::vector < double > d2phidzeta2(_nc);

    std::vector < double > d2phidxideta(_nc);
    std::vector < double > d2phidetadzeta(_nc);
    std::vector < double > d2phidzetadxi(_nc);

    for(int j = 0; j < _nc; j++) {
      phi[j] = _pt_basis->eval_phi(_IND[j], &xi[0]);
      dphidxi[j] = _pt_basis->eval_dphidx(_IND[j], &xi[0]);
      dphideta[j] = _pt_basis->eval_dphidy(_IND[j], &xi[0]);
      dphidzeta[j] = _pt_basis->eval_dphidz(_IND[j], &xi[0]);

      d2phidxi2[j] = _pt_basis->eval_d2phidx2(_IND[j], &xi[0]);
      d2phideta2[j] = _pt_basis->eval_d2phidy2(_IND[j], &xi[0]);
      d2phidzeta2[j] = _pt_basis->eval_d2phidz2(_IND[j], &xi[0]);

      d2phidxideta[j] = _pt_basis->eval_d2phidxdy(_IND[j], &xi[0]);
      d2phidetadzeta[j] = _pt_basis->eval_d2phidydz(_IND[j], &xi[0]);
      d2phidzetadxi[j] = _pt_basis->eval_d2phidzdx(_IND[j], &xi[0]);
    }


    type Jac[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    type JacI[3][3];

    const double* dxi = &dphidxi[0];
    const double* deta = &dphideta[0];
    const double* dzeta = &dphidzeta[0];

    for(int inode = 0; inode < _nc; inode++, dxi++, deta++, dzeta++) {
      Jac[0][0] += (*dxi) * vt[0][inode];
      Jac[0][1] += (*dxi) * vt[1][inode];
      Jac[0][2] += (*dxi) * vt[2][inode];
      Jac[1][0] += (*deta) * vt[0][inode];
      Jac[1][1] += (*deta) * vt[1][inode];
      Jac[1][2] += (*deta) * vt[2][inode];
      Jac[2][0] += (*dzeta) * vt[0][inode];
      Jac[2][1] += (*dzeta) * vt[1][inode];
      Jac[2][2] += (*dzeta) * vt[2][inode];
    }

    type det = (Jac[0][0] * (Jac[1][1] * Jac[2][2] - Jac[1][2] * Jac[2][1]) +
                Jac[0][1] * (Jac[1][2] * Jac[2][0] - Jac[1][0] * Jac[2][2]) +
                Jac[0][2] * (Jac[1][0] * Jac[2][1] - Jac[1][1] * Jac[2][0]));

    JacI[0][0] = (-Jac[1][2] * Jac[2][1] + Jac[1][1] * Jac[2][2]) / det;
    JacI[0][1] = (Jac[0][2] * Jac[2][1] - Jac[0][1] * Jac[2][2]) / det;
    JacI[0][2] = (-Jac[0][2] * Jac[1][1] + Jac[0][1] * Jac[1][2]) / det;
    JacI[1][0] = (Jac[1][2] * Jac[2][0] - Jac[1][0] * Jac[2][2]) / det;
    JacI[1][1] = (-Jac[0][2] * Jac[2][0] + Jac[0][0] * Jac[2][2]) / det;
    JacI[1][2] = (Jac[0][2] * Jac[1][0] - Jac[0][0] * Jac[1][2]) / det;
    JacI[2][0] = (-Jac[1][1] * Jac[2][0] + Jac[1][0] * Jac[2][1]) / det;
    JacI[2][1] = (Jac[0][1] * Jac[2][0] - Jac[0][0] * Jac[2][1]) / det;
    JacI[2][2] = (-Jac[0][1] * Jac[1][0] + Jac[0][0] * Jac[1][1]) / det;

    Weight = det;

    dxi = &dphidxi[0];
    deta = &dphideta[0];
    dzeta = &dphidzeta[0];

    const double* dxi2 = &d2phidxi2[0];
    const double* deta2 = &d2phideta2[0];
    const double* dzeta2 = &d2phidzeta2[0];
    const double* dxideta = &d2phidxideta[0];
    const double* detadzeta = &d2phidetadzeta[0];
    const double* dzetadxi = &d2phidzetadxi[0];

    for(int inode = 0; inode < _nc; inode++, dxi++, deta++, dzeta++, dxi2++, deta2++, dzeta2++, dxideta++, detadzeta++, dzetadxi++) {

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
  }

  //-----------elem_type_3D - END  ----------------------------------------------------------------------------------------------
   
  
  
  
  

} //end namespace femus



#endif


