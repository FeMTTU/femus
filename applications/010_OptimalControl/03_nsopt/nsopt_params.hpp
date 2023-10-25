#ifndef NSOPT_PARAMETERS
#define NSOPT_PARAMETERS


#include   "../param.hpp"

namespace femus {

//*********************** Physics *******************************************************
#define FLUID_DENSITY  1

 double force[3] = {0., 0., 0.};

//******************************** switch between stokes and navier stokes *********************************************
 
 const int advection_flag = 0;
 const int advection_Picard = 0;
 
//  Stokes: advection_flag = 0; advection_Picard = 0;
//  Newton: advection_flag = 1; advection_Picard = 0;
//  Picard: advection_flag = 1; advection_Picard = 1;


 
} //end namespace

 
 
 //****** Convergence ********************************
#define exact_sol_flag 0 // 1 = if we want to use manufactured solution; 0 = if we use regular convention
#define compute_conv_flag 0 // 1 = if we want to compute the convergence and error ; 0 =  no error computation
 //****** Convergence ********************************

 
 namespace ns_state_only  {
    
const unsigned no_of_norms = 5;   // for L2 norm of U,V,P and H1 norm of U,V
    
}


namespace pure_boundary_norms { 
 
 const unsigned no_of_l2_norms = 9;   //U,V,P,UADJ,VADJ,PADJ,ctrl_0,ctrl_1,THETA
 const unsigned no_of_h1_norms = 6;    //U,V,UADJ,VADJ,ctrl_0, ctrl_1

}
 
 
 
  namespace lifting_internal_norms { 

 const unsigned no_of_l2_norms = 11;   //U,V,P,adj_0,adj_1,PADJ,ctrl_0,ctrl_1,PCTRL,U+U0,V+V0
 const unsigned no_of_h1_norms = 8;    //U,V,adj_0,adj_1,ctrl_0,ctrl_1,U+U0,V+V0

  }

  
 
#endif
