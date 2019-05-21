#ifndef SOLIDOPT_PARAMETERS
#define SOLIDOPT_PARAMETERS


//*********************** Sets Number of subdivisions in X and Y direction *****************************************

#define NSUB  2

//*********************** Model *****************************************
#define MODEL "Linear_elastic"
// #define MODEL "Mooney-Rivlin" 
// #define MODEL "Neo-Hookean"

//*********************** Parameters for Solid Model *****************************************
#define YOUNGS_MODULUS 1.5e+6
#define SOLID_DENSITY  1.
#define NI 0.5
//******************************************* Desired Target  and RHS function*******************************************************
    
 double _gravity[3] = {0., 0., 0.};  // gravity
 double TargetDisp[3] = {0., 0., 1.e-1};

//*********************** Sets the regularization parameters *******************************************************

 double alpha_val = 1.;
 double beta_val  = 1.e-3;
 double gamma_val = 1.e-3;
 
 
 
//*********************** Find volume elements that contain a  Target domain element ********************************
template < class var_type >
int ElementTargetFlag(const std::vector<var_type> & elem_center) {

 //***** set target domain flag ********************************** 
  int target_flag = 0;
  
    if ( sqrt(elem_center[0] * elem_center[0] + elem_center[1] * elem_center[1]) < 1.0 + 1.e-5  &&
	  elem_center[2] > 0.75 - 1.e-5  &&  elem_center[2] < 1.0  + 1.e-5
  ) //target for cylinder

//    if ( elem_center[0] > 0.   - 1.e-5  &&  elem_center[0] < 0.25  + 1.e-5  && 
//         elem_center[1] > 0.25 - 1.e-5  &&  elem_center[1] < 0.75  + 1.e-5
//   ) //target on left 
   
//     if (  elem_center[0] > 0.25 - 1.e-5  &&  elem_center[0] < 0.75  + 1.e-5  && 
// 	  elem_center[1] > 0.75  - 1.e-5  &&  elem_center[1] < 1.0   + 1.e-5
//   ) //target on top

//    if ( elem_center[0] > 0.75  - 1.e-5  &&  elem_center[0] < 1.0   + 1.e-5  && 
//         elem_center[1] > 0.25 - 1.e-5  &&  elem_center[1] < 0.75  + 1.e-5
//   ) //target on right 
   
//     if (  elem_center[0] > 0.25 - 1.e-5  &&  elem_center[0] < 0.75  + 1.e-5  && 
// 	  elem_center[1] > 0.   - 1.e-5  &&  elem_center[1] < 0.25  + 1.e-5
//   ) //target on bottom
  
//     if (  elem_center[0] > 0. - 1.e-5  &&  elem_center[0] < 0.25  + 1.e-5  && 
// 	  elem_center[1] > 0.75   - 1.e-5  &&  elem_center[1] < 1.0  + 1.e-5
//   ) //target box  on NW

//     if (  elem_center[0] > 0.75 - 1.e-5  &&  elem_center[0] < 1.0  + 1.e-5  && 
// 	  elem_center[1] > 0.   - 1.e-5  &&  elem_center[1] < 0.25  + 1.e-5
//   ) //target box  on SE


  {
     
     target_flag = 1;
     
  }
  
     return target_flag;

}


//*********************** Find volume elements that contain a Control Face element *********************************

int ControlDomainFlag(const std::vector<double> & elem_center) {

 //***** set control domain flag ***** 
  double mesh_size = 1./NSUB/*_Y*/;
  int control_el_flag = 0;
   if ( elem_center[1] >  1. - mesh_size ) { control_el_flag = 1; }

     return control_el_flag;
}




#endif
