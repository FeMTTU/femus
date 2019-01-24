#ifndef ELLIPTIC_PARAMETERS
#define ELLIPTIC_PARAMETERS

//#include </elliptic_lift_restr_bdd_ctrl_ext/input/ext_box.neu>

//*********************** Sets Number of subdivisions in X and Y direction *****************************************

#define NSUB_X  32
#define NSUB_Y  32


//*********************** Sets the regularization parameters *******************************************************
#define ALPHA_CTRL_BDRY 1.e-3
#define BETA_CTRL_BDRY 1.e-2


#define ALPHA_CTRL_VOL 1.e-3
#define BETA_CTRL_VOL 1.e-2


//*********************** Control box constraints *******************************************************
#define  INEQ_FLAG 1
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
  
   if (   /*elem_center[0] < 0.75 + 1.e-5    && */      elem_center[0] < 0.5  + 1.e-5 /* && */ 
        /*elem_center[1] <  0.5  + 1.e-5     && */    /*elem_center[1] > 0.5 -  1.e-5*/  /*(1./16. + 1./64.)*/
        //elem_center[0] > 0.5 -  1.e-5
  ) {
     
     target_flag = 1;
     
  }
  
     return target_flag;

}


//******************************************* Desired Target *******************************************************

double DesiredTarget()
{
   return 1.;
}




//*********************** Find volume elements that contain a Control Face element *********************************

int ControlDomainFlag_bdry(const std::vector<double> & elem_center) {

 //***** set control domain flag ***** 
  double mesh_size = 1./NSUB_Y;
  int control_el_flag = 0;
   if ( elem_center[0] <  0. + mesh_size ) { control_el_flag = 1; }

     return control_el_flag;
}



//*********************** Find volume elements that contain a Control domain element *********************************

int ControlDomainFlag_internal_restriction(const std::vector<double> & elem_center) {

 //***** set target domain flag ******
 // flag = 1: we are in the lifting nonzero domain
  int control_el_flag = 0.;
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
