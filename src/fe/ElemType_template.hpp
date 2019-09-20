#ifndef __femus_fe_ElemType_templ_hpp__
#define __femus_fe_ElemType_templ_hpp__



#include "ElemType.hpp"
#include "ElemType_template_base.hpp"



namespace femus {


//these classes are meant to be used only for Evaluation at Quadrature Points of jacobians/shape functions/derivatives
// The goal is to put together, for a given Geometric Element, the FE family and the Quadrature Rule.
// We should try to distinguish well the part of the information that is Reference (Abstract) from the part that is Real.

//The old ElemType is meant to be Abstract
// My problem is that ElemType is non-templated, while I would like to make it templated
// There are two kinds of templates that I need:
    // a template with rispect to Dim and Space_dim
    // a template with respect to Type and Type_mov
// The additional templating with respect to Type and Type_Mov should only refer to the Jacobian (and old JacobianSur) routines
// The templating wrt. Type and Type_Mov should be implemented without template specialization, but it should all be abstract
// The templating wrt. Dim and Space_dim should be implemented with template specialization

    
// The original ElemType has two main parts: a part related to Multigrid, a part related to Quadrature. We should split them.

    
    
    template <class type, class type_mov, unsigned int dim, unsigned int space_dim>
     class elem_type_templ  : public elem_type_templ_base<type, type_mov>  {
      
  public: 
      

      ~elem_type_templ(){ }
     
     elem_type_templ(const std::string geom_elem, const std::string fe_elem, const std::string order_gauss) 
//      : elem_type_templ_base<type, type_mov, dim, space_dim>(geom_elem, order_gauss)
     { }
      

     inline void JacJacInv(const std::vector < std::vector < type_mov > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector <type_mov> > & Jac,
                            std::vector < std::vector <type_mov> > & JacI,
                            type_mov & detJac,
                            const unsigned space_dimension) const;

     inline void compute_normal(const std::vector< std::vector< type_mov > > & Jac, std::vector< type_mov > & normal) const;

     inline void shape_funcs_current_elem(const unsigned & ig,
                                             const std::vector < std::vector <type_mov> > & JacI,
                                             vector < double > & phi, 
                                             vector < type >   & gradphi,
                                             boost::optional< vector < type > & > nablaphi,
                                             const unsigned space_dimension) const;
                                             
     inline void shape_funcs_volume_at_bdry_current_elem(const unsigned ig, 
                                                         const unsigned jface, 
                                                         const std::vector < std::vector <type_mov> > & JacI_qp, 
                                                         std::vector < double > & phi_vol_at_bdry,
                                                         std::vector < type >   & phi_x_vol_at_bdry, 
                                                         boost::optional< std::vector < type > & > nablaphi_vol_at_bdry,
                                                         const unsigned space_dimension) const;                                        

     private:

     inline void jacobian(const std::vector < std::vector < type_mov > > & vt,
                           const unsigned & ig,
                           const std::vector < double ** > & dphidxi,
                           std::vector < std::vector <type_mov> > & Jac,
                           const unsigned space_dimension) const;

     inline void jac_jacT(const unsigned & ig,
                          const std::vector < std::vector <type_mov> > & Jac,
                          std::vector < std::vector <type_mov> > & JacJacT,
                          const unsigned space_dimension) const;

     inline void jac_jacT_inv(const unsigned & ig,
                          const std::vector < std::vector <type_mov> > & Jac,
                          std::vector < std::vector <type_mov> > & JacJacT,
                          const unsigned space_dimension) const;

     inline void area_transf(const unsigned & ig,
                          const std::vector < std::vector <type_mov> > & JacJacT,
                          type_mov & area,
                          const unsigned space_dimension) const;

                           
  };

  

// Each class was transformed into a templated one, and the virtuality of the dimensions  is implemented as template specialization
// PARTIAL *CLASS* SPECIALIZATION! (function specialization can only be full)
  template <class type, class type_mov>
   class  elem_type_templ<type, type_mov, 1, 3>  : public elem_type_1D,  public elem_type_templ_base<type, type_mov/*, 1, 3>*/>  {
       
   private: 
       
   void area_transf(const unsigned & ig,
                          const std::vector < std::vector <type_mov> > & JacJacT,
                          type_mov & area,
                          const unsigned space_dimension) const {

      area = sqrt(JacJacT[0][0]);

  }
                          
    void jac_jacT_inv(const unsigned & ig,
                          const std::vector < std::vector <type_mov> > & JacJacT,
                          std::vector < std::vector <type_mov> > & JacJacT_inv,
                          const unsigned space_dimension) const {

    /*const type_mov*/ JacJacT_inv[0][0] = 1. / JacJacT[0][0];
                           
    }
                          
     void jac_jacT(const unsigned & ig,
                   const std::vector < std::vector <type_mov> > & Jac,
                   std::vector < std::vector <type_mov> > & JacJacT,
                   const unsigned space_dimension) const {
                              
//     type_mov JacJacT[1][1];
    JacJacT[0][0] = 0.; //1x1
    for (unsigned d = 0; d < space_dimension; d++) JacJacT[0][0] += Jac[0][d]*Jac[0][d];
    
   }
                            
                            
     void jacobian(const std::vector < std::vector < type_mov > > & vt,
                           const unsigned & ig,
                           const std::vector < double ** > & dphidxi,
                           std::vector < std::vector <type_mov> > & Jac,
                           const unsigned space_dimension) const {
                               
    constexpr unsigned int dim = 1;

    //Jac =================
    Jac.resize(dim);
    
    for (unsigned d = 0; d < dim; d++) { 
        Jac[d].resize(space_dimension);	
        std::fill(Jac[d].begin(), Jac[d].end(), 0.); 
    }

    for (unsigned d = 0; d < space_dimension; d++) {
    const double* dxi_coords  = /*_dphidxi*/dphidxi[0][ig];
       for(int inode = 0; inode < _nc; inode++, dxi_coords++) {
          Jac[0][d] += (*dxi_coords) * vt[d][inode];
        }
     }
                                    
                               
         }

         
   public: 

             
             
             
       elem_type_templ(const std::string geom_elem, const std::string fe_elem, const std::string order_gauss) 
       : elem_type_1D(geom_elem.c_str(), fe_elem.c_str(), order_gauss.c_str() )
       {   }
     
          ~elem_type_templ(){ }

                           
                           
     void JacJacInv(const std::vector < std::vector < type_mov > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector <type_mov> > & Jac,
                            std::vector < std::vector <type_mov> > & JacI,
                            type_mov & detJac,
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

    constexpr unsigned int dim = 1;
    
    //Jac =================
    Jac.resize(dim);
    
    for (unsigned d = 0; d < dim; d++) { 
        Jac[d].resize(space_dim);	
        std::fill(Jac[d].begin(), Jac[d].end(), 0.); 
    }

    for (unsigned d = 0; d < space_dim; d++) {
    const double* dxi_coords  = _dphidxi[ig];
       for(int inode = 0; inode < _nc; inode++, dxi_coords++) {
          Jac[0][d] += (*dxi_coords) * vt[d][inode];
        }
     }
     
    //JacJacT  =================
    type_mov JacJacT[1][1];  JacJacT[0][0] = 0.; //1x1
    for (unsigned d = 0; d < space_dim; d++) JacJacT[0][0] += Jac[0][d]*Jac[0][d];
    
    //area  =================
    detJac = sqrt(JacJacT[0][0]);
    
    //JacJacT_inv  =================
    const type_mov JacJacT_inv = 1. / JacJacT[0][0];
    
    //JacI  =================
    JacI.resize(space_dim);
    for (unsigned d = 0; d < space_dim; d++) JacI[d].resize(dim);

    for (unsigned d = 0; d < space_dim; d++) JacI[d][0] = Jac[0][d] * JacJacT_inv;

  ///@todo in the old implementation shouldn't we take the absolute value??? I'd say we don't because it goes both on the lhs and on the rhs...
             
     }
     
     
     
      void compute_normal(const std::vector< std::vector< type_mov > > & Jac, std::vector< type_mov > & normal) const {
    
    // if you want to compute the normal to a 1d element, you need to know to what plane the boundary element belongs...                                   

    constexpr unsigned int space_dim = 3;
    
    //JacJacT  =================
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
      
      
     void shape_funcs_current_elem(const unsigned & ig,
                                             const std::vector < std::vector <type_mov> > & JacI,
                                             vector < double > & phi, 
                                             vector < type >   & gradphi,
                                             boost::optional< vector < type > & > nablaphi,
                                             const unsigned space_dimension) const {
                                                 
    const double* dxi  = _dphidxi[ig];
    const double* dxi2 = _d2phidxi2[ig];

    phi.resize(_nc);
    gradphi.resize(_nc * space_dimension);  std::fill(gradphi.begin(),gradphi.end(),0.);
    if(nablaphi) nablaphi->resize(_nc * space_dimension);   ///@todo fix this: once space_dim was only 1

    
    for(int inode = 0; inode < _nc; inode++, dxi++, dxi2++) {

      phi[inode] = _phi[ig][inode];
      
      for (unsigned d = 0; d < space_dimension; d++) gradphi[ inode * space_dimension + d] = (*dxi) * JacI[d][0];

      if(nablaphi)(*nablaphi)[inode] = (*dxi2) * JacI[0][0] * JacI[0][0]; ///@todo fix this

    }                                           
                                                 
                                                 
  }
 
 
 
     void shape_funcs_volume_at_bdry_current_elem(const unsigned ig, 
                                                         const unsigned jface, 
                                                         const std::vector < std::vector <type_mov> > & JacI_qp, 
                                                         std::vector < double > & phi_vol_at_bdry,
                                                         std::vector < type >   & phi_x_vol_at_bdry, 
                                                         boost::optional< std::vector < type > & > nablaphi_vol_at_bdry,
                                                         const unsigned space_dimension) const {
 //1d                                                            
         abort();                                                    
                                                             
         }
 
 
};




  template <class type, class type_mov>
   class  elem_type_templ<type, type_mov, 2, 3>  : public elem_type_2D,  public elem_type_templ_base<type, type_mov>   {
       
   private: 
       
  void area_transf(const unsigned & ig,
                          const std::vector < std::vector <type_mov> > & JacJacT,
                          type_mov & area,
                          const unsigned space_dimension) const {
                              
    type_mov detJacJacT = (JacJacT[0][0] * JacJacT[1][1] - JacJacT[0][1] * JacJacT[1][0]);
            
    area = sqrt(abs(detJacJacT));
                              
    }
                          
                          
  void jacobian(const std::vector < std::vector < type_mov > > & vt,
                           const unsigned & ig,
                           const std::vector < double ** > & dphidxi,
                           std::vector < std::vector <type_mov> > & Jac,
                           const unsigned space_dimension) const {

     constexpr unsigned int dim = 2;                           
     //Jac ===============
    Jac.resize(dim);
    
    for (unsigned d = 0; d < dim; d++) { 
        Jac[d].resize(space_dimension);	
        std::fill(Jac[d].begin(), Jac[d].end(), 0.); 
    }

    for (unsigned d = 0; d < space_dimension; d++) {
    const double* dxi_coords  = /*_dphidxi*/ dphidxi[0][ig];
    const double* deta_coords = /*_dphideta*/dphidxi[1][ig];
    
      for(int inode = 0; inode < _nc; inode++, dxi_coords++, deta_coords++) {
          Jac[0][d] += (*dxi_coords)  * vt[d][inode];
          Jac[1][d] += (*deta_coords) * vt[d][inode];
        }
     }                               

  }
       
       
     void jac_jacT(const unsigned & ig,
                          const std::vector < std::vector <type_mov> > & Jac,
                          std::vector < std::vector <type_mov> > & JacJacT,
                          const unsigned space_dimension) const {

//     type_mov JacJacT[2/*dim*/][2/*dim*/] = {{0., 0.}, {0., 0.}};
    
    for (unsigned i = 0; i < 2/*dim*/; i++) std::fill( JacJacT[i].begin(), JacJacT[i].end(), 0.);
    
    for (unsigned i = 0; i < 2/*dim*/; i++)
        for (unsigned j = 0; j < 2/*dim*/; j++)
            for (unsigned k = 0; k < space_dimension; k++) JacJacT[i][j] += Jac[i][k]*Jac[j][k];
                              
   }
   

     void jac_jacT_inv(const unsigned & ig,
                          const std::vector < std::vector <type_mov> > & JacJacT,
                          std::vector < std::vector <type_mov> > & JacJacT_inv,
                          const unsigned space_dimension) const {
                                                            
    const type_mov detJacJacT = (JacJacT[0][0] * JacJacT[1][1] - JacJacT[0][1] * JacJacT[1][0]);
            
    JacJacT_inv[0][0] =  JacJacT[1][1] / detJacJacT;
    JacJacT_inv[0][1] = -JacJacT[0][1] / detJacJacT;
    JacJacT_inv[1][0] = -JacJacT[1][0] / detJacJacT;
    JacJacT_inv[1][1] =  JacJacT[0][0] / detJacJacT;
        
   }
   
   
   public: 

       elem_type_templ(const std::string geom_elem, const std::string fe_elem, const std::string order_gauss) 
       : elem_type_2D(geom_elem.c_str(), fe_elem.c_str(), order_gauss.c_str() )
     { }
     
          ~elem_type_templ(){ }
          
          
                           
     void JacJacInv(const std::vector < std::vector < type_mov > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector <type_mov> > & Jac,
                            std::vector < std::vector <type_mov> > & JacI,
                            type_mov & detJac,
                            const unsigned space_dim) const {

     constexpr unsigned int dim = 2;                           
     //Jac ===============
    Jac.resize(dim);
    
    for (unsigned d = 0; d < dim; d++) { 
        Jac[d].resize(space_dim);	
        std::fill(Jac[d].begin(), Jac[d].end(), 0.); 
    }

    for (unsigned d = 0; d < space_dim; d++) {
    const double* dxi_coords  = _dphidxi[ig];
    const double* deta_coords = _dphideta[ig];
    
      for(int inode = 0; inode < _nc; inode++, dxi_coords++, deta_coords++) {
          Jac[0][d] += (*dxi_coords)  * vt[d][inode];
          Jac[1][d] += (*deta_coords) * vt[d][inode];
        }
     }
     

    //JacJacT  =================
    type_mov JacJacT[2/*dim*/][2/*dim*/] = {{0., 0.}, {0., 0.}};
    
    for (unsigned i = 0; i < 2/*dim*/; i++)
        for (unsigned j = 0; j < 2/*dim*/; j++)
            for (unsigned k = 0; k < space_dim; k++) JacJacT[i][j] += Jac[i][k]*Jac[j][k];
            
    //area  =================
    type_mov detJacJacT = (JacJacT[0][0] * JacJacT[1][1] - JacJacT[0][1] * JacJacT[1][0]);
            
    detJac = sqrt(abs(detJacJacT));
    
    //JacJacT_inv  =================
    type_mov JacJacT_inv[2/*dim*/][2/*dim*/];/* = {{0., 0.}, {0., 0.}};*/
    
    JacJacT_inv[0][0] =  JacJacT[1][1] / detJacJacT;
    JacJacT_inv[0][1] = -JacJacT[0][1] / detJacJacT;
    JacJacT_inv[1][0] = -JacJacT[1][0] / detJacJacT;
    JacJacT_inv[1][1] =  JacJacT[0][0] / detJacJacT;
        
     //JacI ===============
    JacI.resize(space_dim);
    for (unsigned d = 0; d < space_dim; d++) { 
        JacI[d].resize(dim);  
        std::fill(JacI[d].begin(),JacI[d].end(),0.); 
    }
    
    for (unsigned i = 0; i < space_dim; i++)
        for (unsigned j = 0; j < dim; j++)
            for (unsigned k = 0; k < dim; k++) JacI[i][j] += Jac[k][i]*JacJacT_inv[k][j];

    
     }
    
          
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

 
 
      void shape_funcs_current_elem(const unsigned & ig,
                                             const std::vector < std::vector <type_mov> > & JacI,
                                             vector < double > & phi, 
                                             vector < type >   & gradphi,
                                             boost::optional< vector < type > & > nablaphi,
                                             const unsigned space_dimension) const {
                                                 

    phi.resize(_nc);
    
    gradphi.resize(_nc * space_dimension);  std::fill(gradphi.begin(),gradphi.end(),0.);
    const double* dxi  = _dphidxi[ig];
    const double* deta = _dphideta[ig];
    
    if(nablaphi) nablaphi->resize(_nc * 3);
    const double* dxi2 = _d2phidxi2[ig];
    const double* deta2 = _d2phideta2[ig];
    const double* dxideta = _d2phidxideta[ig];

    
    for(int inode = 0; inode < _nc; inode++, dxi++, deta++, dxi2++, deta2++, dxideta++) {

      phi[inode] = _phi[ig][inode];

      for (unsigned d = 0; d < space_dimension; d++) gradphi[ inode * space_dimension + d] = (*dxi) * JacI[d][0] + (*deta) * JacI[d][1];

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
 
   
   
    void shape_funcs_volume_at_bdry_current_elem(const unsigned ig, 
                                                         const unsigned jface, 
                                                         const std::vector < std::vector <type_mov> > & JacI_qp, 
                                                         std::vector < double > & phi_vol_at_bdry,
                                                         std::vector < type >   & phi_x_vol_at_bdry, 
                                                         boost::optional< std::vector < type > & > nablaphi_vol_at_bdry,
                                                         const unsigned space_dimension) const {
 //2d                                                            
         abort();                                                    
                                                             
         }
   
   
};



  template <class type, class type_mov>
   class  elem_type_templ<type, type_mov, 3, 3>   : public elem_type_3D,  public elem_type_templ_base<type, type_mov>   {

   private:
       
  void area_transf(const unsigned & ig,
                          const std::vector < std::vector <type_mov> > & Jac,
                          type_mov & area,
                          const unsigned space_dimension) const {
                              
  //here you don't pass JacJacT, but more simply just Jac                         
                              
    area = (Jac[0][0] * (Jac[1][1] * Jac[2][2] - Jac[1][2] * Jac[2][1]) +
            Jac[0][1] * (Jac[1][2] * Jac[2][0] - Jac[1][0] * Jac[2][2]) +
            Jac[0][2] * (Jac[1][0] * Jac[2][1] - Jac[1][1] * Jac[2][0]));
    
  }
                          
                          
   void jac_jacT_inv(const unsigned & ig,
                          const std::vector < std::vector <type_mov> > & JacJacT,
                          std::vector < std::vector <type_mov> > & JacJacT_inv,
                          const unsigned space_dimension) const {
                              
       std::cout << "Not needed with 3d in 3d space"; abort();                              

       }
                                                            
      void jac_jacT(const unsigned & ig,
                          const std::vector < std::vector <type_mov> > & Jac,
                          std::vector < std::vector <type_mov> > & JacJacT,
                          const unsigned space_dimension) const {
                              
       std::cout << "Not needed with 3d in 3d space"; abort();                              

       }
                          
       
     void jacobian(const std::vector < std::vector < type_mov > > & vt,
                           const unsigned & ig,
                           const std::vector < double ** > & dphidxi,
                           std::vector < std::vector <type_mov> > & Jac,
                           const unsigned space_dimension) const {
                               
     constexpr unsigned int dim = 3;
     
     //Jac ===============
//     Jac.resize(3/*dim*/);
//     for (unsigned d = 0; d < 3/*dim*/; d++) {
// //         Jac[d].resize(3/*space_dim*/);   
//         std::fill(Jac[d].begin(), Jac[d].end(), 0.); 
//     }
        std::fill(Jac[0].begin(), Jac[0].end(), 0.); 
        std::fill(Jac[1].begin(), Jac[1].end(), 0.); 
        std::fill(Jac[2].begin(), Jac[2].end(), 0.); 
    

    
//     for (unsigned d = 0; d < 3/*dim*/; d++) {
    const double * dxi_coords   = /*_dphidxi*/  dphidxi[0][ig];
    const double * deta_coords  = /*_dphideta*/ dphidxi[1][ig];
    const double * dzeta_coords = /*_dphidzeta*/dphidxi[2][ig];

    for(int inode = 0; inode < _nc; inode++, dxi_coords++, deta_coords++, dzeta_coords++) {
      Jac[0][0] += (*dxi_coords) * vt[0][inode];
      Jac[0][1] += (*dxi_coords) * vt[1][inode];
      Jac[0][2] += (*dxi_coords) * vt[2][inode];
      Jac[1][0] += (*deta_coords) * vt[0][inode];
      Jac[1][1] += (*deta_coords) * vt[1][inode];
      Jac[1][2] += (*deta_coords) * vt[2][inode];
      Jac[2][0] += (*dzeta_coords) * vt[0][inode];
      Jac[2][1] += (*dzeta_coords) * vt[1][inode];
      Jac[2][2] += (*dzeta_coords) * vt[2][inode];
      
//       Jac[0][d] += (*dxi_coords)   * vt[d][inode];
//       Jac[1][d] += (*deta_coords)  * vt[d][inode];
//       Jac[2][d] += (*dzeta_coords) * vt[d][inode];
       }
//    }                               
                               
                               
    }
    
    
         public: 

       elem_type_templ(const std::string geom_elem, const std::string fe_elem, const std::string order_gauss)
       : elem_type_3D(geom_elem.c_str(), fe_elem.c_str(), order_gauss.c_str() )
     {}
     
          ~elem_type_templ(){}
          
                           
                               
  void JacJacInv(const std::vector < std::vector < type_mov > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector <type_mov> > & Jac,
                            std::vector < std::vector <type_mov> > & JacI,
                            type_mov & detJac,
                            const unsigned space_dim) const {
                                                                
     constexpr unsigned int dim = 3;
     
     //Jac ===============
//     Jac.resize(3/*dim*/);
//     for (unsigned d = 0; d < 3/*dim*/; d++) {
// //         Jac[d].resize(3/*space_dim*/);   
//         std::fill(Jac[d].begin(), Jac[d].end(), 0.); 
//     }
        std::fill(Jac[0].begin(), Jac[0].end(), 0.); 
        std::fill(Jac[1].begin(), Jac[1].end(), 0.); 
        std::fill(Jac[2].begin(), Jac[2].end(), 0.); 
    
// Jac[0][0] = 0.;
// Jac[0][1] = 0.; 
// Jac[0][2] = 0.; 
// Jac[1][0] = 0.; 
// Jac[1][1] = 0.; 
// Jac[1][2] = 0.; 
// Jac[2][0] = 0.; 
// Jac[2][1] = 0.; 
// Jac[2][2] = 0.; 
    
    
//     for (unsigned d = 0; d < 3/*dim*/; d++) {
    const double * dxi_coords   = _dphidxi[ig];
    const double * deta_coords  = _dphideta[ig];
    const double * dzeta_coords = _dphidzeta[ig];

    for(int inode = 0; inode < _nc; inode++, dxi_coords++, deta_coords++, dzeta_coords++) {
      Jac[0][0] += (*dxi_coords) * vt[0][inode];
      Jac[0][1] += (*dxi_coords) * vt[1][inode];
      Jac[0][2] += (*dxi_coords) * vt[2][inode];
      Jac[1][0] += (*deta_coords) * vt[0][inode];
      Jac[1][1] += (*deta_coords) * vt[1][inode];
      Jac[1][2] += (*deta_coords) * vt[2][inode];
      Jac[2][0] += (*dzeta_coords) * vt[0][inode];
      Jac[2][1] += (*dzeta_coords) * vt[1][inode];
      Jac[2][2] += (*dzeta_coords) * vt[2][inode];
      
//       Jac[0][d] += (*dxi_coords)   * vt[d][inode];
//       Jac[1][d] += (*deta_coords)  * vt[d][inode];
//       Jac[2][d] += (*dzeta_coords) * vt[d][inode];
       }
//    }
    
     //JacI ===============
//     JacI.resize(space_dim);
//     for (unsigned d = 0; d < space_dim; d++) JacI[d].resize(dim);
                                
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
     
     
   
    void compute_normal(const std::vector< std::vector< type_mov > > & Jac, std::vector< type_mov > & normal) const {
           std::cout << "Normal non-defined for 3D objects" << std::endl; abort();
     }
     
     

      void shape_funcs_current_elem(const unsigned & ig,
                                             const std::vector < std::vector <type_mov> > & JacI,
                                             vector < double > & phi, 
                                             vector < type >   & gradphi,
                                             boost::optional< vector < type > & > nablaphi,
                                             const unsigned space_dimension) const {

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
    gradphi.resize(_nc * 3);  /*std::fill(gradphi.begin(),gradphi.end(),0.);*/
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


                                                 
}


     void shape_funcs_volume_at_bdry_current_elem(const unsigned ig, 
                                                         const unsigned jface, 
                                                         const std::vector < std::vector <type_mov> > & JacI_qp, 
                                                         std::vector < double > & phi_vol_at_bdry,
                                                         std::vector < type >   & phi_x_vol_at_bdry, 
                                                         boost::optional< std::vector < type > & > nablaphi_vol_at_bdry,
                                                         const unsigned space_dimension) const {
 //3d                                                            
         abort();                                                    
                                                             
         }
         
};


    
} //end namespace femus

#endif
