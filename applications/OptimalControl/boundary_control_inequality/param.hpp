#ifndef ELLIPTIC_PARAMETERS
#define ELLIPTIC_PARAMETERS



#include "Mesh.hpp"


//*********************** Sets Number of subdivisions in X and Y direction *****************************************

#define NSUB_X  2
#define NSUB_Y  2
#define NSUB_Z  2


//*********************** Sets the regularization parameters *******************************************************
#define ALPHA_CTRL_BDRY 1.e-3
#define BETA_CTRL_BDRY 1.e-2


#define ALPHA_CTRL_VOL 1.e-3
#define BETA_CTRL_VOL 1.e-2


//*********************** Control boundary extremes *******************************************************

#define GAMMA_CONTROL_LOWER 0.25
#define GAMMA_CONTROL_UPPER 0.75


//*********************** Control box constraints *******************************************************
#define  INEQ_FLAG 1.
#define  C_COMPL 1.


namespace femus {


 double InequalityConstraint(const std::vector<double> & dof_obj_coord, const bool upper) {

     double constr_value = 0.;
     double constr_value_upper =  .3;// dof_obj_coord[1]*(1. - dof_obj_coord[1]);
     double constr_value_lower = -1000.; //-3.e-13;
     assert(constr_value_lower < constr_value_upper); 
     
    if (upper)   constr_value = constr_value_upper;
    else         constr_value = constr_value_lower; 
    
    
  return constr_value;
     
}
   

   
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
   
  const int  target_line_sign = target_line_sign_func(FACE_FOR_CONTROL);
  
  const unsigned int axis_dir = axis_direction_target_reg(FACE_FOR_CONTROL);
   
   const double target_line = 0.5 + target_line_sign * offset_to_include_line; 
   
   
   
      if ((  target_line_sign * elem_center[axis_dir] < target_line_sign * target_line ) && 
          (  target_line_sign * elem_center[axis_dir] > - 0.5 + target_line_sign * (0.5 - target_line_sign * offset_to_include_line)))
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

     
   const int  target_line_sign = target_line_sign_func(FACE_FOR_CONTROL);

   const double extreme_pos = extreme_position(FACE_FOR_CONTROL);
   
   const unsigned int axis_dir = axis_direction_Gamma_control(FACE_FOR_CONTROL);

  
   if ( ( target_line_sign * elem_center[1 - axis_dir] <   target_line_sign * (  extreme_pos  + target_line_sign * mesh_size) )
       && ( elem_center[axis_dir] > GAMMA_CONTROL_LOWER - offset_to_include_line ) 
       && ( elem_center[axis_dir] < GAMMA_CONTROL_UPPER + offset_to_include_line ) )
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
  
   
   const int  target_line_sign = target_line_sign_func(FACE_FOR_CONTROL);

   const double extreme_pos = extreme_position(FACE_FOR_CONTROL);

   const unsigned int axis_dir = axis_direction_Gamma_control(FACE_FOR_CONTROL);

   
   if ( ( target_line_sign * elem_center[1 - axis_dir] <   target_line_sign * ( extreme_pos + target_line_sign * control_domain_width ) )
       && ( elem_center[axis_dir] > GAMMA_CONTROL_LOWER - offset_to_include_line ) 
       && ( elem_center[axis_dir] < GAMMA_CONTROL_UPPER + offset_to_include_line ) )
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
                                                             const std::vector < std::vector < double > > sol_eldofs,  ///@todo why not a reference?
                                                             const std::vector < unsigned int > & Sol_n_el_dofs,
                                                             const unsigned int pos_mu,
                                                             const unsigned int pos_ctrl,
                                                             const double c_compl,
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
                                                             const std::vector < std::vector < double > > & coords_at_dofs,
                                                             const std::vector < std::vector < double > > sol_eldofs,  ///@todo why not a reference?
                                                             const std::vector < unsigned int > & Sol_n_el_dofs,
                                                             const unsigned int pos_mu,
                                                             const unsigned int pos_ctrl,
                                                             const double c_compl,
                                                                   std::vector < double > & ctrl_lower,
                                                                   std::vector < double > & ctrl_upper,
                                                                   std::vector < double > & sol_actflag,
                                                             const unsigned int solFEType_act_flag,
                                                             const unsigned int solIndex_act_flag) {
     
		const unsigned nve_bdry = msh->GetElementFaceDofNumber(iel,jface, solFEType_act_flag);
        
        const unsigned dim = coords_at_dofs.size();
        
         // 0: inactive; 1: active_a; 2: active_b
        assert(Sol_n_el_dofs[pos_mu] == Sol_n_el_dofs[pos_ctrl]);///@todo More appropriately, 
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
    
 
 
 void node_insertion_bdry(const unsigned int iel,
                         const unsigned int jface,
                         const   Mesh* msh,
                         const     vector < vector < int > > & L2G_dofmap,
                         const double c_compl,
                         const double ineq_flag,
                         const unsigned int pos_mu,
                         const unsigned int pos_ctrl,
                         const std::vector < std::vector < double > > & sol_eldofs,
                         const std::vector < unsigned int > & Sol_n_el_dofs,
                         std::vector < double > & sol_actflag,
                         const unsigned int solFEType_act_flag,
                         const std::vector < double > & ctrl_lower,
                         const std::vector < double > & ctrl_upper,
                         SparseMatrix*             KK,
                         NumericVector* RES
                        ) {

 //============= delta_mu row ===============================
      std::vector<double> Res_mu (Sol_n_el_dofs[pos_mu]); std::fill(Res_mu.begin(),Res_mu.end(), 0.);
      
      for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
	    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
//     for (unsigned i = 0; i < sol_actflag.size(); i++) {
      if (sol_actflag[i_bdry] == 0){  //inactive
         Res_mu [i_vol] = - ineq_flag * ( 1. * sol_eldofs[pos_mu][i_vol] - 0. ); 
// 	 Res_mu [i] = Res[nDof_u + nDof_ctrl + nDof_adj + i]; 
      }
      else if (sol_actflag[i_bdry] == 1){  //active_a 
	 Res_mu [i_vol] = - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl][i_vol] - c_compl * ctrl_lower[i_bdry]);
      }
      else if (sol_actflag[i_bdry] == 2){  //active_b 
	Res_mu [i_vol]  =  - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl][i_vol] - c_compl * ctrl_upper[i_bdry]);
      }
    }

    
    RES->insert(Res_mu,  L2G_dofmap[pos_mu]);    
 //============= delta_mu row - end ===============================
    
 //============= delta_mu-delta_ctrl row ===============================
 //auxiliary volume vector for act flag
 unsigned nDof_actflag_vol  = msh->GetElementDofNumber(iel, solFEType_act_flag);
 std::vector<double> sol_actflag_vol(nDof_actflag_vol); 


 for (unsigned i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++) if (sol_actflag[i_bdry] != 0 ) sol_actflag[i_bdry] = ineq_flag * c_compl;    
 
 std::fill(sol_actflag_vol.begin(), sol_actflag_vol.end(), 0.);
    for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
       unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
       sol_actflag_vol[i_vol] = sol_actflag[i_bdry];
    }
 
 KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap[pos_mu], L2G_dofmap[pos_ctrl], sol_actflag_vol);
 //============= delta_mu-delta_ctrl row - end ===============================

 //============= delta_mu-delta_mu row ===============================
  for (unsigned i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++) sol_actflag[i_bdry] =  ineq_flag * (1 - sol_actflag[i_bdry]/c_compl)  + (1-ineq_flag) * 1.;  //can do better to avoid division, maybe use modulo operator 

 std::fill(sol_actflag_vol.begin(), sol_actflag_vol.end(), 0.);
    for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
       unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
       sol_actflag_vol[i_vol] = sol_actflag[i_bdry];
    }
  
  KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap[pos_mu], L2G_dofmap[pos_mu], sol_actflag_vol );
 //============= delta_mu-delta_mu row - end ===============================
  

}


} //end namespace
 
#endif
