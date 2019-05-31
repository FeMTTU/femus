#ifndef ELLIPTIC_PARAMETERS
#define ELLIPTIC_PARAMETERS



#include "Mesh.hpp"


//*********************** Sets Number of subdivisions in X and Y direction *****************************************

#define NSUB_X  32
#define NSUB_Y  32


//*********************** Sets the regularization parameters *******************************************************
#define ALPHA_CTRL_BDRY 1.e-3
#define BETA_CTRL_BDRY 1.e-2


#define ALPHA_CTRL_VOL 1.e-3
#define BETA_CTRL_VOL 1.e-2


//*********************** Control box constraints *******************************************************
#define  INEQ_FLAG 1.
#define  C_COMPL 1.


 double InequalityConstraint(const std::vector<double> & dof_obj_coord, const bool upper) {

     double constr_value = 0.;
     double constr_value_upper =  .2;//0.5 * dof_obj_coord[AXIS_DIRECTION_CONTROL_SIDE]; //dof_obj_coord[1]*(1. - dof_obj_coord[1]);
     double constr_value_lower = -1000.; //-3.e-13;
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
   
      if ((  target_line_sign * elem_center[1-AXIS_DIRECTION_CONTROL_SIDE] < target_line_sign * target_line ) && 
          (  target_line_sign * elem_center[1-AXIS_DIRECTION_CONTROL_SIDE] > - 0.5 + target_line_sign * (0.5 - target_line_sign * offset_to_include_line)))
          {  target_flag = 1;  }
  
     return target_flag;

}


//******************************************* Desired Target *******************************************************

double DesiredTarget()
{
   return 1.;
}




//*********************** Find volume elements that contain a Control Face element *********************************

int ControlDomainFlag_bdry(const std::vector<double> & elem_center) {

  const double mesh_size = 1./*/NSUB_X*/;  //this picks a lot more elements, but then the if on the faces only gets the control boundary
   
  int control_el_flag = 0;
  
  const double offset_to_include_line = 1.e-5;

  double target_line_sign;
  double extreme_pos;
  
        if (FACE_FOR_CONTROL == 3 || FACE_FOR_CONTROL == 2) { target_line_sign = -1; extreme_pos = 1.; }
   else if (FACE_FOR_CONTROL == 1 || FACE_FOR_CONTROL == 4) { target_line_sign = 1;  extreme_pos = 0.; }

  
   if ( ( target_line_sign * elem_center[1-AXIS_DIRECTION_CONTROL_SIDE] <   target_line_sign * (  extreme_pos  + target_line_sign * mesh_size) )
       && ( elem_center[AXIS_DIRECTION_CONTROL_SIDE] > 0.25 - offset_to_include_line ) 
       && ( elem_center[AXIS_DIRECTION_CONTROL_SIDE] < 0.75 + offset_to_include_line ) )
      { control_el_flag = 1; }

     return control_el_flag;
}



//*********************** Find volume elements that contain a Control domain element *********************************

int ControlDomainFlag_internal_restriction(const std::vector<double> & elem_center) {

 //***** set target domain flag ******
 // flag = 1: we are in the lifting nonzero domain
  int control_el_flag = 0.;
  
  const double offset_to_include_line = 1.e-5;
  
  double control_domain_width = 0.25;
  
  double target_line_sign;
  double extreme_pos;
  
        if (FACE_FOR_CONTROL == 3 || FACE_FOR_CONTROL == 2) { target_line_sign = -1; extreme_pos = 1.;}
   else if (FACE_FOR_CONTROL == 1 || FACE_FOR_CONTROL == 4) { target_line_sign = 1;  extreme_pos = 0.;}
   
   if ( ( target_line_sign * elem_center[1-AXIS_DIRECTION_CONTROL_SIDE] <   target_line_sign * ( extreme_pos + target_line_sign * control_domain_width ) )
       && ( elem_center[AXIS_DIRECTION_CONTROL_SIDE] > 0.25 - offset_to_include_line ) 
       && ( elem_center[AXIS_DIRECTION_CONTROL_SIDE] < 0.75 + offset_to_include_line ) )
      { control_el_flag = 1; }
   
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



 void update_active_set_flag_for_current_nonlinear_iteration(const femus::Mesh* msh,
                                                             const femus::Solution* sol,
                                                             const unsigned int iel,
                                                             const std::vector < std::vector < double > > & coords_at_dofs,
                                                             const std::vector < std::vector < double > > sol_eldofs,
                                                             const std::vector < unsigned int > & Sol_n_el_dofs,
                                                             const unsigned int pos_mu,
                                                             const unsigned int pos_ctrl,
                                                             const unsigned int c_compl,
                                                                   std::vector < double > & ctrl_lower,
                                                                   std::vector < double > & ctrl_upper,
                                                                   std::vector < double > & sol_actflag,
                                                             const unsigned int solFEType_act_flag,
                                                             const unsigned int solIndex_act_flag) {
     
        const unsigned int dim = coords_at_dofs.size();
        
// 0: inactive; 1: active_a; 2: active_b
        assert(Sol_n_el_dofs[pos_mu] == Sol_n_el_dofs[pos_ctrl]);
        sol_actflag.resize(Sol_n_el_dofs[pos_mu]);
        ctrl_lower.resize(Sol_n_el_dofs[pos_mu]);
        ctrl_upper.resize(Sol_n_el_dofs[pos_mu]);
        std::fill(sol_actflag.begin(), sol_actflag.end(), 0);
        std::fill(ctrl_lower.begin(), ctrl_lower.end(), 0.);
        std::fill(ctrl_upper.begin(), ctrl_upper.end(), 0.);

        for (unsigned i = 0; i < sol_actflag.size(); i++) {
                std::cout << sol_eldofs[pos_mu][i] << " ";
        }
        
        std::cout << std::endl;
        
        for (unsigned i = 0; i < sol_actflag.size(); i++) {
                std::cout << sol_eldofs[pos_ctrl][i] << " ";
        }
        
        for (unsigned i = 0; i < sol_actflag.size(); i++) {
            std::vector<double> node_coords_i(dim,0.);
            for (unsigned d = 0; d < dim; d++) node_coords_i[d] = coords_at_dofs[d][i];
            
            ctrl_lower[i] = InequalityConstraint(node_coords_i,false);
            ctrl_upper[i] = InequalityConstraint(node_coords_i,true);
            
            const double lower_test_value = sol_eldofs[pos_mu][i] + c_compl * ( sol_eldofs[pos_ctrl][i] - ctrl_lower[i] );
            const double upper_test_value = sol_eldofs[pos_mu][i] + c_compl * ( sol_eldofs[pos_ctrl][i] - ctrl_upper[i] );

            if      ( lower_test_value < 0 )  {
                std::cout << "Found active node below" << std::endl;
                std::cout << "The current value of mu is " <<  sol_eldofs[pos_mu][i] << std::endl;
                   sol_actflag[i] = 1;
            }
            else if ( upper_test_value > 0 )  {
                std::cout << "Found active node above" << std::endl;
                std::cout << "The current value of mu is " <<  sol_eldofs[pos_mu][i] << std::endl;
                sol_actflag[i] = 2;
            }
        }

//************** local to global act flag ***************************

        for (unsigned i = 0; i < sol_actflag.size(); i++) {
            unsigned solDof_mu = msh->GetSolutionDof(i, iel, solFEType_act_flag);
            (sol->_Sol[solIndex_act_flag])->set(solDof_mu,sol_actflag[i]);
        }

}





 void update_active_set_flag_for_current_nonlinear_iteration_bdry(const femus::Mesh* msh,
                                                             const femus::Solution* sol,
                                                             const unsigned int iel,
                                                             const unsigned int jface,
                                                             const unsigned int solType_coords,
                                                             const std::vector < std::vector < double > > & coords_at_dofs,
                                                             const std::vector < std::vector < double > > sol_eldofs,
                                                             const std::vector < unsigned int > & Sol_n_el_dofs,
                                                             const unsigned int pos_mu,
                                                             const unsigned int pos_ctrl,
                                                             const unsigned int c_compl,
                                                                   std::vector < double > & ctrl_lower,
                                                                   std::vector < double > & ctrl_upper,
                                                                   std::vector < double > & sol_actflag,
                                                             const unsigned int solFEType_act_flag,
                                                             const unsigned int solIndex_act_flag) {
     
		const unsigned nve_bdry = msh->GetElementFaceDofNumber(iel,jface,solType_coords);
        
        const unsigned dim = coords_at_dofs.size();
        
         // 0: inactive; 1: active_a; 2: active_b
        assert(Sol_n_el_dofs[pos_mu] == Sol_n_el_dofs[pos_ctrl]);
            sol_actflag.resize(nve_bdry/*nDof_mu*/);
            ctrl_lower.resize(nve_bdry/*nDof_mu*/);
            ctrl_upper.resize(nve_bdry/*nDof_mu*/);
           std::fill(sol_actflag.begin(), sol_actflag.end(), 0);
           std::fill(ctrl_lower.begin(), ctrl_lower.end(), 0.);
           std::fill(ctrl_upper.begin(), ctrl_upper.end(), 0.);

      for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
        std::vector<double> node_coords_i(dim,0.);
        for (unsigned d = 0; d < dim; d++) node_coords_i[d] = coords_at_dofs[d][i_bdry];
        
        ctrl_lower[i_bdry] = InequalityConstraint(node_coords_i,false);
        ctrl_upper[i_bdry] = InequalityConstraint(node_coords_i,true);

        const double lower_test_value = sol_eldofs[pos_mu][i_vol] + c_compl * ( sol_eldofs[pos_ctrl][i_vol] - ctrl_lower[i_bdry] );
        const double upper_test_value = sol_eldofs[pos_mu][i_vol] + c_compl * ( sol_eldofs[pos_ctrl][i_vol] - ctrl_upper[i_bdry] );
        
        if      ( lower_test_value < 0 )  sol_actflag[i_bdry] = 1;
        else if ( upper_test_value > 0 )  sol_actflag[i_bdry] = 2;
            }
            
        //************** act flag **************************** 
      for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
	    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
      unsigned solDof_actflag = msh->GetSolutionDof(i_vol, iel, solFEType_act_flag); 
      (sol->_Sol[solIndex_act_flag])->set(solDof_actflag,sol_actflag[i_bdry]);     
    }
    
    
    }
    
    

 std::vector< double > face_elem_center(const std::vector< std::vector< double > > & coords_at_dofs_bdry) {
     
     const unsigned int dim = coords_at_dofs_bdry.size();
     
            std::vector < double > elem_center_bdry(dim);
            
            for (unsigned j = 0; j < dim; j++) {  elem_center_bdry[j] = 0.;  }
            
            
            for (unsigned j = 0; j < dim; j++) {
                for (unsigned i = 0; i < coords_at_dofs_bdry[j].size(); i++) {
                    elem_center_bdry[j] += coords_at_dofs_bdry[j][i];
                }
            }
            for (unsigned j = 0; j < dim; j++) {
                elem_center_bdry[j] = elem_center_bdry[j]/coords_at_dofs_bdry[j].size();
            }
            
            return elem_center_bdry;
            
 }
 
 
#endif
