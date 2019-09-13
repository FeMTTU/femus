#ifndef __femus_fe_ElemType_Jac_templ_hpp__
#define __femus_fe_ElemType_Jac_templ_hpp__



#include "GaussPoints.hpp"

namespace femus
{

    
  template <class type, class type_mov, unsigned int dim, unsigned int space_dim>
      class elem_type_jac_templ_base {
          
          
  protected:
      
     elem_type_jac_templ_base(const std::string geom_elem, const std::string order_gauss)
      
     : _gauss(geom_elem.c_str(), order_gauss.c_str())  { }
         
      
      const Gauss _gauss;
      
      };
 
      
      
      
    
    template <class type, class type_mov, unsigned int dim, unsigned int space_dim>
     class elem_type_jac_templ  : public elem_type_jac_templ_base<type, type_mov,dim, space_dim>    {
      
  public: 
      
       
      ~elem_type_jac_templ(){ }
     
     elem_type_jac_templ(const std::string geom_elem, const std::string order_gauss) :
     elem_type_jac_templ_base<type, type_mov, dim, space_dim>(geom_elem, order_gauss) { }
      
  };

  


// PARTIAL *CLASS* SPECIALIZATION! (function specialization can only be full)
  template <class type, class type_mov>
   class  elem_type_jac_templ<type, type_mov, 1, 3>  : public elem_type_jac_templ_base<type, type_mov, 1, 3>  {
       
   private: 
       
       
       
         public: 

             
             
             
       elem_type_jac_templ(const std::string geom_elem, const std::string order_gauss) :
     elem_type_jac_templ_base<type, type_mov, 1, 3>(geom_elem, order_gauss)
     
     { 
         
         
         
         
         
    }
     
          ~elem_type_jac_templ(){ }
          
          
     void Jacobian_type_geometry(const std::vector < std::vector < type_mov > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector <type_mov> > & Jac,
                            std::vector < std::vector <type_mov> > & JacI,
                            type_mov & detJac,
                            const unsigned dim,
                            const unsigned space_dim) const {
// // // //here the convention for the Jacobian is that the real coordinates are put along a COLUMN, so you have
// // //  //    J = [ d x_1/ d xi  |  d x_2/ d xi  | d x_3 / d xi  ]     (1x3 matrix)
// // //  // Then, if we denote differentials with D, we have
// // //  //  [ D x_1 |  D x_2 | D x_3 ] =  D \xi [ d x_1/ d xi  |  d x_2/ d xi  | d x_3 / d xi  ]                                
// // // 
// // //  //  [ D x_1 |  D x_2 | D x_3 ] =   D \xi  J                               
// // //                                 
// // //  //  [ D x_1 |  D x_2 | D x_3 ] J^T               =  D \xi  J  J^T                              
// // // 
// // //  //  [ D x_1 |  D x_2 | D x_3 ] J^T (J  J^T)^{-1} =  D \xi                              
// // //                                 
// // //  //  [ D x_1 |  D x_2 | D x_3 ] | d xi / dx_1 | =  D \xi                              
// // //  //                             | d xi / dx_2 |
// // //  //                             | d xi / dx_3 |   
// // // // 
// // // // | d xi / dx_1 | 
// // // // | d xi / dx_2 | =  J^T (J J^T)^{-1}
// // // // | d xi / dx_3 |                              
// // // 
// // //                                 
// // //     //Jac =================
// // //     Jac.resize(dim);
// // //     
// // //     for (unsigned d = 0; d < dim; d++) { 
// // //         Jac[d].resize(space_dim);	std::fill(Jac[d].begin(), Jac[d].end(), 0.); }
// // // 
// // //     for (unsigned d = 0; d < space_dim; d++) {
// // //     const double* dxi_coords  = _dphidxi[ig];
// // //        for(int inode = 0; inode < _nc; inode++, dxi_coords++) {
// // //           Jac[0][d] += (*dxi_coords) * vt[d][inode];
// // //         }
// // //      }
// // //      
// // //     //JacI  =================
// // //     type_mov JacJacT[1][1];  JacJacT[0][0] = 0.; //1x1
// // //     for (unsigned d = 0; d < space_dim; d++) JacJacT[0][0] += Jac[0][d]*Jac[0][d];
// // //     detJac = sqrt(JacJacT[0][0]);
// // //     
// // //     JacI.resize(space_dim);
// // //     for (unsigned d = 0; d < space_dim; d++) JacI[d].resize(dim);
// // // 
// // //     for (unsigned d = 0; d < space_dim; d++) JacI[d][0] = Jac[0][d] * 1. / JacJacT[0][0];
// // // 
// // //   ///@todo in the old implementation shouldn't we take the absolute value??? I'd say we don't because it goes both on the lhs and on the rhs...
             
     }
     
     
   
};




  template <class type, class type_mov>
   class  elem_type_jac_templ<type, type_mov, 2, 3>  : public elem_type_jac_templ_base<type, type_mov, 2, 3>  {
       
         public: 

       elem_type_jac_templ(const std::string geom_elem, const std::string order_gauss) :
     elem_type_jac_templ_base<type, type_mov, 2, 3>(geom_elem, order_gauss)
     { }
     
          ~elem_type_jac_templ(){ }
   
};



  template <class type, class type_mov>
   class  elem_type_jac_templ<type, type_mov, 3, 3>  : public elem_type_jac_templ_base<type, type_mov, 3, 3>  {
       
         public: 

       elem_type_jac_templ(const std::string geom_elem, const std::string order_gauss) :
     elem_type_jac_templ_base<type, type_mov, 3, 3>(geom_elem, order_gauss)
     {}
     
          ~elem_type_jac_templ(){}
   
   
        void Jacobian_type_geometry(const std::vector < std::vector < type_mov > > & vt,
                            const unsigned & ig,
                            std::vector < std::vector <type_mov> > & Jac,
                            std::vector < std::vector <type_mov> > & JacI,
                            type_mov & detJac,
                            const unsigned dim,
                            const unsigned space_dim) const {
                                                                
//      //Jac ===============
//     Jac.resize(3/*dim*/);
//     for (unsigned d = 0; d < 3/*dim*/; d++) {
//         Jac[d].resize(3/*space_dim*/);   std::fill(Jac[d].begin(), Jac[d].end(), 0.); }
//     
//     for (unsigned d = 0; d < 3/*dim*/; d++) {
//     const double * dxi_coords   = _dphidxi[ig];
//     const double * deta_coords  = _dphideta[ig];
//     const double * dzeta_coords = _dphidzeta[ig];
// 
//     for(int inode = 0; inode < _nc; inode++, dxi_coords++, deta_coords++, dzeta_coords++) {
//       Jac[0][d] += (*dxi_coords)   * vt[d][inode];
//       Jac[1][d] += (*deta_coords)  * vt[d][inode];
//       Jac[2][d] += (*dzeta_coords) * vt[d][inode];
//        }
//     }
//     
//      //JacI ===============
//     JacI.resize(space_dim);
//     for (unsigned d = 0; d < space_dim; d++) JacI[d].resize(dim);
//                                 
//     detJac = (Jac[0][0] * (Jac[1][1] * Jac[2][2] - Jac[1][2] * Jac[2][1]) +
//                     Jac[0][1] * (Jac[1][2] * Jac[2][0] - Jac[1][0] * Jac[2][2]) +
//                     Jac[0][2] * (Jac[1][0] * Jac[2][1] - Jac[1][1] * Jac[2][0]));
// 
//     JacI[0][0] = (-Jac[1][2] * Jac[2][1] + Jac[1][1] * Jac[2][2]) / detJac;
//     JacI[0][1] = (Jac[0][2] * Jac[2][1] - Jac[0][1] * Jac[2][2]) / detJac;
//     JacI[0][2] = (-Jac[0][2] * Jac[1][1] + Jac[0][1] * Jac[1][2]) / detJac;
//     JacI[1][0] = (Jac[1][2] * Jac[2][0] - Jac[1][0] * Jac[2][2]) / detJac;
//     JacI[1][1] = (-Jac[0][2] * Jac[2][0] + Jac[0][0] * Jac[2][2]) / detJac;
//     JacI[1][2] = (Jac[0][2] * Jac[1][0] - Jac[0][0] * Jac[1][2]) / detJac;
//     JacI[2][0] = (-Jac[1][1] * Jac[2][0] + Jac[1][0] * Jac[2][1]) / detJac;
//     JacI[2][1] = (Jac[0][1] * Jac[2][0] - Jac[0][0] * Jac[2][1]) / detJac;
//     JacI[2][2] = (-Jac[0][1] * Jac[1][0] + Jac[0][0] * Jac[1][1]) / detJac;                               

             
     }
   
   
};


    
} //end namespace femus

#endif
