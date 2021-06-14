#ifndef ELLIPTIC_NONLIN_PARAMETERS
#define ELLIPTIC_NONLIN_PARAMETERS


//*********************** Sets Number of subdivisions in X and Y direction *****************************************

#define NSUB_X  64
#define NSUB_Y  64


//*********************** Sets the regularization parameters *******************************************************
#define ALPHA_CTRL_VOL 1.e-5
#define BETA_CTRL_VOL 1.e-5


//*********************** Control box constraints *******************************************************
#define  INEQ_FLAG 1
#define  C_COMPL 1.


 double InequalityConstraint(const std::vector<double> & dof_obj_coord, const bool upper) {

     double constr_value = 0.;
     double constr_value_upper = 0.2 + dof_obj_coord[0]*(1. - dof_obj_coord[0]);
     double constr_value_lower = -1000.;
     assert(constr_value_lower < constr_value_upper); 
     
    if (upper)   constr_value = constr_value_upper;
    else         constr_value = constr_value_lower; 
    
    
  return constr_value;
     
}


//*********************** Find volume elements that contain a  Target domain element **************************************

int ElementTargetFlag(const std::vector<double> & elem_center) {

 //***** set target domain flag ******
  int target_flag = 0; //set 0 to 1 to get the entire domain
  
   if (   elem_center[0] < 0.75 + 1.e-5    &&  elem_center[0] > 0.25  - 1.e-5  &&  
          elem_center[1] < 0.9  + 1.e-5   &&  elem_center[1] > 0.75 -  1.e-5  /*(1./16. + 1./64.)*/
  ) {
     
     target_flag = 1;
     
  }
  
     return target_flag;

}


//******************************************* Desired Target *******************************************************

double DesiredTarget(const std::vector<double> & xyz)
{
   return /*sin(M_PI * xyz[0]) * sin(M_PI * xyz[1]);*/ 2.;
}




//*********************** Find volume elements that contain a Control Face element *********************************

int ControlDomainFlag_bdry(const std::vector<double> & elem_center) {

 //***** set control domain flag ***** 
  double mesh_size = 1./NSUB_Y;
  int control_el_flag = 0;
   if ( elem_center[1] >  1. - mesh_size ) { control_el_flag = 1; }

     return control_el_flag;
}



//*********************** Find volume elements that contain a Control domain element *********************************

int ControlDomainFlag_internal_restriction(const std::vector<double> & elem_center) {

 //***** set target domain flag ******
 // flag = 1: we are in the lifting nonzero domain
  int control_el_flag = 1.;
   if ( elem_center[1] >  0.7) { control_el_flag = 1; }

     return control_el_flag;

}





#endif
