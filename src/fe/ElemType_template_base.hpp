#ifndef __femus_fe_ElemType_templ_base_hpp__
#define __femus_fe_ElemType_templ_base_hpp__

#include <vector>
#include <string>
#include <boost/optional.hpp>


namespace femus

{
    
    
    
 template <class type, class type_mov>
    class elem_type_templ_base {
          
      public:

     elem_type_templ_base(){};
         
     virtual const double get_dphidxi_ref(const unsigned idim, const unsigned qp, const unsigned dof) = 0;

     virtual void fill_dphidxi_at_quadrature_points() = 0;     
     virtual void fill_dphidxi_at_quadrature_points_vol_at_bdry() = 0;     
                                   
     virtual /*inline*/ void JacJacInv_vol_at_bdry_new(const std::vector < std::vector < type_mov > > & vt,
                                                       const unsigned & ig,
                                                       const unsigned jface,
                                                       std::vector < std::vector <type_mov> > & Jac,
                                                       std::vector < std::vector <type_mov> > & JacI,
                                                       type_mov & detJac,
                                                       const unsigned space_dimension) /*const*/ = 0;

     virtual /*inline*/ void JacJacInv(const std::vector < std::vector < type_mov > > & vt,
                                       const unsigned & ig,
                                       std::vector < std::vector <type_mov> > & Jac,
                                       std::vector < std::vector <type_mov> > & JacI,
                                       type_mov & detJac,
                                       const unsigned space_dimension) const = 0;

     virtual /*inline*/ void compute_normal(const std::vector< std::vector< type_mov > > & Jac, std::vector< type_mov > & normal) const = 0;

     virtual /*inline*/ void shape_funcs_current_elem_flexible(const unsigned & ig,
                                                               const std::vector < std::vector <type_mov> > & JacI,
                                                               double **     phi_ref,
                                                               const std::vector < double ** > & dphidxi_ref,
                                                               const std::vector < double ** > & d2phidxi2_ref,
                                                               std::vector < double > & phi,
                                                               std::vector < type >   & gradphi,
                                                               boost::optional< std::vector < type > & > nablaphi,
                                                               const unsigned space_dimension) const = 0;
                                             
     virtual /*inline*/ void shape_funcs_current_elem(const unsigned & ig,
                                                      const std::vector < std::vector <type_mov> > & JacI,
                                                      std::vector < double > & phi,
                                                      std::vector < type >   & gradphi,
                                                      boost::optional< std::vector < type > & > nablaphi,
                                                      const unsigned space_dimension) const = 0;
                                             
    virtual /*inline*/ void shape_funcs_vol_at_bdry_current_elem(const unsigned ig,
                                                                 const unsigned jface,
                                                                 const std::vector < std::vector <type_mov> > & JacI_qp,
                                                                 std::vector < double > & phi_vol_at_bdry,
                                                                 std::vector < type >   & phi_x_vol_at_bdry,
                                                                 boost::optional< std::vector < type > & > nablaphi_vol_at_bdry,
                                                                 const unsigned space_dimension) const = 0;
          
// run-time selection
      static elem_type_templ_base<type, type_mov> * build(const std::string geom_elem, /*dimension is contained in the Geometric Element*/
                                                          const std::string fe_fam,
                                                          const std::string order_gauss,
                                                          const unsigned space_dimension); 

     virtual /*inline*/ void jac_jacT(const std::vector < std::vector <type_mov> > & Jac,
                                      std::vector < std::vector <type_mov> > & JacJacT,
                                      const unsigned space_dimension) const = 0;

     virtual /*inline*/ void jac_jacT_inv(const std::vector < std::vector <type_mov> > & JacJacT,
                                          std::vector < std::vector <type_mov> > & JacJacT_inv,
                                          const unsigned space_dimension) const = 0;
                                         
    protected:
        
     virtual /*inline*/ void jacobian_flexible(const std::vector < std::vector < type_mov > > & vt,
                                               const unsigned & ig,
                                               const std::vector < std::vector < std::vector < double > > > & dphidxi,
                                               std::vector < std::vector <type_mov> > & Jac,
                                               const unsigned space_dimension) const = 0;

     virtual /*inline*/ void measure_transf(const std::vector < std::vector <type_mov> > & JacJacT,
                                            type_mov & area,
                                            const unsigned space_dimension) const = 0;

     virtual /*inline*/ void jacobian_inv(const std::vector < std::vector <type_mov> > & Jac,
                                          const std::vector < std::vector <type_mov> > & JacJacT_inv,
                                          std::vector < std::vector <type_mov> > & Jac_inv,
                                          const unsigned space_dimension) const = 0;
                                          
    
      };    



}

#endif
