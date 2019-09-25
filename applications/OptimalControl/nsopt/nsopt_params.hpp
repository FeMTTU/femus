#ifndef NSOPT_PARAMETERS
#define NSOPT_PARAMETERS


#include "Mesh.hpp"


//*********************** Sets Number of subdivisions in X and Y direction *****************************************

#define NSUB_X  2
#define NSUB_Y  2
#define NSUB_Z  2

#define FLUID_DENSITY  1

//*********************** Control boundary extremes *******************************************************

#define GAMMA_CONTROL_LOWER 0.25
#define GAMMA_CONTROL_UPPER 0.75

//******************************************* Desired Target  and RHS function*******************************************************

 double force[3] = {0.,0.,0.};
 double Vel_desired[3] = {0.,1.,0.};

//*********************** Sets the regularization parameters *******************************************************

 double alpha_val = 1.;
 double beta_val  = 1.e-3;
 double gamma_val = 1.e-3;
 
 
//******************************** switch between stokes and navier stokes *********************************************
 
 int advection_flag = 1;
 int advection_Picard = 0;
 
//  Newton: advection_flag = 1; advection_Picard = 0;
//  Picard: advection_flag = 1; advection_Picard = 1;
 
 
const unsigned int axis_direction_Gamma_control(const unsigned int face_index) {
    
    int axis_dir;
    
        if (face_index == 1 || face_index == 2) { axis_dir = 1; }
   else if (face_index == 3 || face_index == 4) { axis_dir = 0; }
   else if (face_index == 5 || face_index == 6) { axis_dir = 0; }

    return axis_dir;
    
}


const unsigned int axis_direction_target_reg(const unsigned int face_index) {
    
    int axis_dir;
    
        if (face_index == 1 || face_index == 2) { axis_dir = 0; }
   else if (face_index == 3 || face_index == 4) { axis_dir = 1; }
   else if (face_index == 5 || face_index == 6) { axis_dir = 2; }

    return axis_dir;
    
}



   
const  int target_line_sign_func(const unsigned int face_index) {
    
   int  target_line_sign;
  
        if (face_index == 1 || face_index == 3 || face_index == 5) { target_line_sign = 1;  }
   else if (face_index == 2 || face_index == 4 || face_index == 6) { target_line_sign = -1; }
   
   return target_line_sign;
   
}



const double extreme_position(const unsigned int face_index) {
    
  double extreme_pos;
  
        if (face_index == 1 || face_index == 3 || face_index == 5) {  extreme_pos = 0.; }
   else if (face_index == 2 || face_index == 4 || face_index == 6) {  extreme_pos = 1.; }
   
   return extreme_pos;
   
}


//*********************** Find volume elements that contain a  Target domain element **************************************

int ElementTargetFlag(const std::vector<double> & elem_center) {

 //***** set target domain flag ******
  int target_flag = 0; //set 0 to 1 to get the entire domain
  
  const double offset_to_include_line = 1.e-5;
   
  double target_region_width = 0.25;
  
  const int  target_line_sign = target_line_sign_func(FACE_FOR_CONTROL);
  
  const double extreme_pos = extreme_position(FACE_FOR_CONTROL);
   
  const unsigned int axis_dir = axis_direction_Gamma_control(FACE_FOR_CONTROL);
 
  if (   ( target_line_sign * elem_center[1 - axis_dir] <   target_line_sign * ( extreme_pos + target_line_sign * target_region_width + target_line_sign * offset_to_include_line ) )
      && ( target_line_sign * elem_center[1 - axis_dir] > - 0.5 + target_line_sign * (0.5 - target_line_sign * offset_to_include_line))
      && ( elem_center[axis_dir] > GAMMA_CONTROL_LOWER - offset_to_include_line ) 
      && ( elem_center[axis_dir] < GAMMA_CONTROL_UPPER + offset_to_include_line ) 
     )
   {  target_flag = 1;  }
  
     return target_flag;

}


// // // //*********************** Find volume elements that contain a  Target domain element ********************************
// // // 
// // // int ElementTargetFlag(const std::vector<double> & elem_center) {
// // // 
// // //  //***** set target domain flag ********************************** 
// // //   int target_flag = 0;
// // // 
// // // //     if ( sqrt(elem_center[0] * elem_center[0] + elem_center[1] * elem_center[1]) < 0.5 + 1.e-5  &&
// // // // 	  elem_center[2] > 0.75 - 1.e-5  &&  elem_center[2] < 1.0  + 1.e-5
// // // //   ) //target for cylinder
// // // 
// // // //     if (  elem_center[0] > 0.25 - 1.e-5  &&  elem_center[0] < 0.75  + 1.e-5  && 
// // // // 	  elem_center[1] > 0.75  - 1.e-5  &&  elem_center[1] < 1.0   + 1.e-5 /*&&
// // // // 	  elem_center[2] > 0.0 - 1.e-5  &&  elem_center[2] < 1.0  + 1.e-5*/
// // // //   ) //target on top for cube
// // //  
// // // //    if ( elem_center[0] > 0.   - 1.e-5  &&  elem_center[0] < 0.25  + 1.e-5  && 
// // // //         elem_center[1] > 0.25 - 1.e-5  &&  elem_center[1] < 0.75  + 1.e-5
// // // //   ) //target on left 
// // //    
// // // //     if (  elem_center[0] > 0.25 - 1.e-5  &&  elem_center[0] < 0.75  + 1.e-5  && 
// // // // 	      elem_center[1] > 0.75  - 1.e-5  &&  elem_center[1] < 1.0   + 1.e-5 
// // // //   ) //target on top
// // // 
// // // //    if ( elem_center[0] > 0.75  - 1.e-5  &&  elem_center[0] < 1.0   + 1.e-5  && 
// // // //         elem_center[1] > 0.25 - 1.e-5  &&  elem_center[1] < 0.75  + 1.e-5
// // // //   ) //target on right 
// // //    
// // //     if (  elem_center[0] > 0.25 - 1.e-5  &&  elem_center[0] < 0.75  + 1.e-5  && 
// // // 	  elem_center[1] > 0.   - 1.e-5  &&  elem_center[1] < 0.25  + 1.e-5
// // //   ) //target on bottom
// // //   
// // // //     if (  elem_center[0] > 0. - 1.e-5  &&  elem_center[0] < 0.25  + 1.e-5  && 
// // // // 	  elem_center[1] > 0.75   - 1.e-5  &&  elem_center[1] < 1.0  + 1.e-5
// // // //   ) //target box  on NW
// // // 
// // // //     if (  elem_center[0] > 0.75 - 1.e-5  &&  elem_center[0] < 1.0  + 1.e-5  && 
// // // // 	  elem_center[1] > 0.   - 1.e-5  &&  elem_center[1] < 0.25  + 1.e-5
// // // //   ) //target box  on SE
// // // 
// // // 
// // //   {
// // //      
// // //      target_flag = 1;
// // //      
// // //   }
// // //   
// // //      return target_flag;
// // // 
// // // }


//*********************** Find volume elements that contain a Control Face element *********************************

int ControlDomainFlag_bdry(const std::vector<double> & elem_center) {

  const double mesh_size = 1./*/NSUB_X*/;  //this picks a lot more elements, but then the if on the faces only gets the control boundary
   
  int control_el_flag = 0;
  
  const double offset_to_include_line = 1.e-5;

     
   const int  target_line_sign = target_line_sign_func(FACE_FOR_CONTROL);
 
   const unsigned int axis_dir = axis_direction_target_reg(FACE_FOR_CONTROL);

   const double target_line = 0.5 + target_line_sign * offset_to_include_line; 

  
     if (     (  target_line_sign * elem_center[axis_dir] < target_line_sign * target_line ) 
           && (  target_line_sign * elem_center[axis_dir] > - 0.5 + target_line_sign * (0.5 - target_line_sign * offset_to_include_line)))
            { control_el_flag = 1; }

     return control_el_flag;
}




// // // //*********************** Find volume elements that contain a Control Face element *********************************
// // // 
// // // int ControlDomainFlag(const std::vector<double> & elem_center) {
// // // 
// // //  //***** set control domain flag ***** 
// // //   double mesh_size = 1./NSUB_Y;
// // //   int control_el_flag = 0;
// // //    if ( elem_center[1] >  1. - mesh_size ) { control_el_flag = 1; }
// // // 
// // //      return control_el_flag;
// // // }




#endif
