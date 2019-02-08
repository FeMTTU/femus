#ifndef ELLIPTIC_PARAMETERS
#define ELLIPTIC_PARAMETERS


//*********************** Sets Number of subdivisions in X and Y direction *****************************************

#define NSUB_X  128
#define NSUB_Y  128


//*********************** Sets the regularization parameters *******************************************************
#define ALPHA_CTRL_BDRY 1.e-3
#define BETA_CTRL_BDRY 0 //1.e-2


#define ALPHA_CTRL_VOL 1.e-3
#define BETA_CTRL_VOL 1.e-2


//*********************** Control box constraints *******************************************************
#define  INEQ_FLAG 0
#define  C_COMPL 1.


 double InequalityConstraint(const std::vector<double> & dof_obj_coord, const bool upper) {

     double constr_value = 0.;
     double constr_value_upper = 0.5; //0.2 + dof_obj_coord[0]*(1. - dof_obj_coord[0]);
     double constr_value_lower = -1000; //-3.e-13;
     assert(constr_value_lower < constr_value_upper); 
     
    if (upper)   constr_value = constr_value_upper;
    else         constr_value = constr_value_lower; 
    
    
  return constr_value;
     
}
   


//*********************** Find volume elements that contain a  Target domain element **************************************

int ElementTargetFlag(const std::vector<double> & elem_center) {

 //***** set target domain flag ******
  int target_flag = 0; //set 0 to 1 to get the entire domain
  
  double target_line_sign;
  
        if (FACE_FOR_CONTROL == 3 || FACE_FOR_CONTROL == 2) { target_line_sign = -1; }
   else if (FACE_FOR_CONTROL == 1 || FACE_FOR_CONTROL == 4) { target_line_sign = 1;  }
   
   const double offset_to_include_line = 1.e-5;
   const double target_line = 0.5 + target_line_sign * offset_to_include_line; 
   
      if (  target_line_sign * elem_center[1-AXIS_DIRECTION_CONTROL_SIDE] < target_line_sign * target_line ) {  target_flag = 1;  }
  
     return target_flag;

}


//******************************************* Desired Target *******************************************************

double DesiredTarget()
{
   return 1.;
}




//*********************** Find volume elements that contain a Control Face element *********************************

int ControlDomainFlag_bdry(const std::vector<double> & elem_center) {

  const double mesh_size = 1./NSUB_X;
   
  int control_el_flag = 0;

  double target_line_sign;
  double extreme_pos;
  
        if (FACE_FOR_CONTROL == 3 || FACE_FOR_CONTROL == 2) { target_line_sign = -1; extreme_pos = 1.; }
   else if (FACE_FOR_CONTROL == 1 || FACE_FOR_CONTROL == 4) { target_line_sign = 1;  extreme_pos = 0.; }

  
   if (  target_line_sign * elem_center[1-AXIS_DIRECTION_CONTROL_SIDE] <   target_line_sign * (  extreme_pos  + target_line_sign * mesh_size) ) { control_el_flag = 1; }

     return control_el_flag;
}



//*********************** Find volume elements that contain a Control domain element *********************************

int ControlDomainFlag_internal_restriction(const std::vector<double> & elem_center) {

 //***** set target domain flag ******
 // flag = 1: we are in the lifting nonzero domain
  int control_el_flag = 1.;
   if ( elem_center[0] >  0.7) { control_el_flag = 1; }

     return control_el_flag;

}


//*********************** Find volume elements that contain a Control domain element *********************************

int ControlDomainFlag_external_restriction(const std::vector<double> & elem_center) {

 //***** set target domain flag ******
 // flag = 1: we are in the lifting nonzero domain
  int exterior_el_flag = 0.;
   if ( elem_center[0] >  1. -  1.e-5) { exterior_el_flag = 1; }

     return exterior_el_flag;

}




#endif
