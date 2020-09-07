#ifndef ELLIPTIC_PARAMETERS
#define ELLIPTIC_PARAMETERS



#include "Mesh.hpp"
#include "CurrentElem.hpp"
#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"
#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"


#include "../fractional_functions.hpp"


#include <functional>



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
     double constr_value_upper =  .03;// dof_obj_coord[1]*(1. - dof_obj_coord[1]);
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
                                                             const unsigned int iface,
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
     
		const unsigned nve_bdry = msh->GetElementFaceDofNumber(iel,iface, solFEType_act_flag);
        
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
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
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
	    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
      unsigned solDof_actflag = msh->GetSolutionDof(i_vol, iel, solFEType_act_flag); 
      (sol->_Sol[solIndex_act_flag])->set(solDof_actflag,sol_actflag[i_bdry]);     
    }
    
    
    }
    
 
 
 void node_insertion_bdry(const unsigned int iel,
                         const unsigned int iface,
                         const   Mesh* msh,
                         const     vector < vector < int > > & L2G_dofmap,
                         const unsigned int pos_mu,
                         const unsigned int pos_ctrl,
                         const std::vector < std::vector < double > > & sol_eldofs,
                         const std::vector < unsigned int > & Sol_n_el_dofs,
                         std::vector < double > & sol_actflag,
                         const unsigned int solFEType_act_flag,
                         const double ineq_flag,
                         const double c_compl,
                         const std::vector < double > & ctrl_lower,
                         const std::vector < double > & ctrl_upper,
                         SparseMatrix*             KK,
                         NumericVector* RES,
                          const bool assembleMatrix
                        ) {

// Create the L2G boundary maps from the volume ones
  std::vector < int > L2G_dofmap_mu_bdry(sol_actflag.size());
  std::vector < int > L2G_dofmap_ctrl_bdry(sol_actflag.size());

      for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
	    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
  L2G_dofmap_mu_bdry[i_bdry]   = L2G_dofmap[pos_mu][i_vol];
  L2G_dofmap_ctrl_bdry[i_bdry] = L2G_dofmap[pos_ctrl][i_vol];
      }
      
 //============= delta_mu row ===============================
      std::vector<double> Res_mu_bdry (sol_actflag.size());     std::fill(Res_mu_bdry.begin(),Res_mu_bdry.end(), 0.);
      std::vector<double> Res_mu (Sol_n_el_dofs[pos_mu]);       std::fill(Res_mu.begin(),Res_mu.end(), 0.);
      
      for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
	    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
        
      if (sol_actflag[i_bdry] == 0) {  //inactive
         Res_mu [i_vol]      = - ineq_flag * ( 1. * sol_eldofs[pos_mu][i_vol] - 0. ); 
         Res_mu_bdry[i_bdry] = - ineq_flag * ( 1. * sol_eldofs[pos_mu][i_vol] - 0. ); 
      }
      else if (sol_actflag[i_bdry] == 1) {  //active_a 
	 Res_mu [i_vol]      = - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl][i_vol] - c_compl * ctrl_lower[i_bdry]);
     Res_mu_bdry[i_bdry] = - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl][i_vol] - c_compl * ctrl_lower[i_bdry]);
          
    }
      else if (sol_actflag[i_bdry] == 2) {  //active_b 
	Res_mu [i_vol]      =  - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl][i_vol] - c_compl * ctrl_upper[i_bdry]);
    Res_mu_bdry[i_bdry] =  - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl][i_vol] - c_compl * ctrl_upper[i_bdry]);
      }
    }

    
//     RES->insert(Res_mu,  L2G_dofmap[pos_mu]);    
    RES->insert(Res_mu_bdry,  L2G_dofmap_mu_bdry);    
 //============= delta_mu row - end ===============================
    
 //============= delta_mu-delta_ctrl row ===============================
 //auxiliary volume vector for act flag
 unsigned nDof_actflag_vol  = msh->GetElementDofNumber(iel, solFEType_act_flag);
 std::vector<double> sol_actflag_vol(nDof_actflag_vol); 


 for (unsigned i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++) if (sol_actflag[i_bdry] != 0 ) sol_actflag[i_bdry] = ineq_flag * c_compl;    
 
 std::fill(sol_actflag_vol.begin(), sol_actflag_vol.end(), 0.);
    for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
       unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
       sol_actflag_vol[i_vol] = sol_actflag[i_bdry];
    }
 
//  KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap[pos_mu], L2G_dofmap[pos_ctrl], sol_actflag_vol);
 if (assembleMatrix) { KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap_mu_bdry, L2G_dofmap_ctrl_bdry, sol_actflag); }
 //============= delta_mu-delta_ctrl row - end ===============================

 //============= delta_mu-delta_mu row ===============================
 // Attention: this equation goes in contrast with \mu = 0 on \Omega \setminus \Gamma_c
 // In fact, here we shouldn't insert all VOLUME values, but only the BOUNDARY ones
 // The best way is to then do a L2G_map ON THE BOUNDARY only
 
  for (unsigned i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++) sol_actflag[i_bdry] =  ineq_flag * (1 - sol_actflag[i_bdry]/c_compl)  + (1-ineq_flag) * 1.;  //can do better to avoid division, maybe use modulo operator 

 std::fill(sol_actflag_vol.begin(), sol_actflag_vol.end(), 0.);
    for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
       unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
       sol_actflag_vol[i_vol] = sol_actflag[i_bdry];
    }
  
//   KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap[pos_mu], L2G_dofmap[pos_mu], sol_actflag_vol );
  if (assembleMatrix) { KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap_mu_bdry, L2G_dofmap_mu_bdry, sol_actflag);  }
 //============= delta_mu-delta_mu row - end ===============================
  

}


///@todo This is being added to a weak form?
 void add_one_times_mu_res_ctrl_bdry(const unsigned iproc,
                         const double ineq_flag,
                         const unsigned int pos_ctrl,
                         const unsigned int pos_mu,
                         const vector < unsigned > & SolIndex,
                         const Solution*                sol,
                         const NonLinearImplicitSystemWithPrimalDualActiveSetMethod * mlPdeSys,
                         const  LinearEquationSolver* pdeSys,
                         NumericVector* RES) {
     
 const unsigned int ctrl_index = pos_ctrl;
 const unsigned int mu_index = pos_mu;

  unsigned int ctrl_size_iproc = pdeSys->KKoffset[ctrl_index + 1][iproc] - pdeSys->KKoffset[ctrl_index][iproc];
  unsigned int mu_size_iproc = (*sol->_Sol[ SolIndex[pos_mu] ]).last_local_index() - (*sol->_Sol[ SolIndex[pos_mu] ]).first_local_index(); // pdeSys->KKoffset[mu_index + 1][iproc] - pdeSys->KKoffset[mu_index][iproc];

  assert(ctrl_size_iproc == mu_size_iproc);

  std::vector<double>  one_times_mu(ctrl_size_iproc, 0.);
  std::vector<int>    positions_ctrl_in_Res(ctrl_size_iproc);
  std::vector<int>    positions_mu_in_Sol(mu_size_iproc);      

  for (unsigned i = 0; i < positions_ctrl_in_Res.size(); i++) {
    positions_ctrl_in_Res[i] = pdeSys->KKoffset[ctrl_index][iproc] + i;
    positions_mu_in_Sol[i] = (*sol->_Sol[ SolIndex[pos_mu] ]).first_local_index()/*pdeSys->KKoffset[mu_index][iproc]*/ + i;
    //this should not come from pdeSys but from Sol ///@todo put the Dof range for Sol //actually I can take it from the Numeric Vector!
//                 unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);  //this is only if I am on an ELEMENT loop, but here I am in a NODE loop
    
    one_times_mu[i] = ineq_flag * 1. * (*sol->_Sol[ SolIndex[pos_mu] ])(positions_mu_in_Sol[i]/*i*//*position_mu_i*/) ;
    }
    RES->add_vector_blocked(one_times_mu, positions_ctrl_in_Res);
    
 }
 
 

 
void el_dofs_unknowns(const Solution*                sol,
                      const Mesh * msh,
                      const  LinearEquationSolver* pdeSys,
                      const unsigned int iel,
                      const    vector < unsigned > & SolFEType,
                      const vector < unsigned > & SolIndex,
                      const vector < unsigned > & SolPdeIndex,
                      vector < unsigned > & Sol_n_el_dofs, 
                      vector < vector < double > > & sol_eldofs,  
                      vector < vector < int > > & L2G_dofmap ) {
    
    assert(Sol_n_el_dofs.size() == sol_eldofs.size());
    
        //all vars###################################################################
        for (unsigned  k = 0; k < Sol_n_el_dofs.size(); k++) {
            unsigned  ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
            Sol_n_el_dofs[k] = ndofs_unk;
            sol_eldofs[k].resize(ndofs_unk);
            L2G_dofmap[k].resize(ndofs_unk);
            for (unsigned i = 0; i < ndofs_unk; i++) {
                unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);
                sol_eldofs[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);
                L2G_dofmap[k][i] = pdeSys->GetSystemDof(SolIndex[k], SolPdeIndex[k], i, iel);
            }
        }
        //all vars###################################################################

}


  void store_act_flag_in_old(  const NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys,
                               const MultiLevelSolution *    ml_sol,
                              Solution *                sol,
                             unsigned int & solIndex_act_flag,
                             unsigned int & solFEType_act_flag) {
      
      
  const std::string act_flag_name = mlPdeSys->GetActiveSetFlagName();
  solIndex_act_flag = ml_sol->GetIndex(act_flag_name.c_str());
  solFEType_act_flag = ml_sol->GetSolutionType(solIndex_act_flag); 
     if(sol->GetSolutionTimeOrder(solIndex_act_flag) == 2) {
       *(sol->_SolOld[solIndex_act_flag]) = *(sol->_Sol[solIndex_act_flag]);
     }
     
  }

  bool check_if_boundary_control_face(const int bdry_index_i, const unsigned int face_in_rectangle_domain_i) {
      
	    // look for boundary faces && look for control faces
		
   return ( bdry_index_i < 0 && face_in_rectangle_domain_i == FACE_FOR_CONTROL );
      
  }
  
  
  bool check_if_same_elem(const unsigned iel, const unsigned jel) {
      
   return (iel == jel);
      
  }
  
  bool check_if_same_elem_bdry(const unsigned iel, const unsigned jel, const unsigned iface, const unsigned jface) {
      
   return (iel == jel && iface == jface);
      
  }

  
  void mixed_integral(const unsigned UNBOUNDED,
                      const unsigned dim,
                      const unsigned dim_bdry,
                      const double EX_1,
                      const double EX_2,
                      const double weight_iqp_bdry,
                      const std::vector < double > & x_iqp_bdry,
                      const std::vector < double > & phi_ctrl_iel_bdry_iqp_bdry,
                      const double sol_ctrl_iqp_bdry,
                      const double s_frac,
                      const double check_limits,
                      const double C_ns,
                      const unsigned int OP_Hhalf,
                      const unsigned int nDof_iel,
                      std::vector < double > & KK_local_iel,
                      std::vector < double > & Res_local_iel
                     ) {
      
  if(UNBOUNDED == 1) {
      
      //============ Mixed Integral 1D - Analytical ==================      
      if (dim_bdry == 1) {
      
              double ex_1 = EX_1;
              double ex_2 = EX_2;
              std::vector < double > ex_1_vec(dim);
              std::vector < double > ex_2_vec(dim);
              ex_1_vec[0] = EX_1;
              ex_1_vec[1] = 0.;
              ex_2_vec[0] = EX_2;
              ex_2_vec[1] = 0.;
//               ex_1_vec[0] = 0.;
//               ex_1_vec[1] = EX_1;
//               ex_2_vec[0] = 0.;
//               ex_2_vec[1] = EX_2;
              
              double dist_1 = 0.;  //distance from node to extreme 1
              double dist_2 = 0.;  //distance from node to extreme 2
              
              for(int d = 0; d < dim; d++) {
                dist_1 += sqrt((x_iqp_bdry[d] - ex_1_vec[d]) * (x_iqp_bdry[d] - ex_1_vec[d]));
                dist_2 += sqrt((x_iqp_bdry[d] - ex_2_vec[d]) * (x_iqp_bdry[d] - ex_2_vec[d]));
              }
              
              double mixed_term = pow(dist_1, -2. * s_frac) + pow(dist_2, - 2. * s_frac);

              for(unsigned i = 0; i < nDof_iel; i++) {
                for(unsigned j = 0; j < nDof_iel /* @todo is this correct? */; j++) {
                  KK_local_iel[ i * nDof_iel + j ] += (0.5 * C_ns) * check_limits * (1. / s_frac) * OP_Hhalf * phi_ctrl_iel_bdry_iqp_bdry[i] * phi_ctrl_iel_bdry_iqp_bdry[j] * weight_iqp_bdry * mixed_term;
                }
                Res_local_iel[ i ] += (0.5 * C_ns) * check_limits * (1. / s_frac) * OP_Hhalf * phi_ctrl_iel_bdry_iqp_bdry[i] * sol_ctrl_iqp_bdry * weight_iqp_bdry * mixed_term;
              }
   
      }
      
    //============ Mixed Integral 2D - Numerical ==================      
      else if (dim_bdry == 2) {           
          std::cout << "check dim vs dim_bdry ";   abort();
// // //             double mixed_term1 = 0.;
// // // //     for(int kel = msh->_elementOffset[iproc]; kel < msh->_elementOffset[iproc + 1]; kel++) {
// // //             // *** Face Gauss point loop (boundary Integral) ***
// // //             for(unsigned jj = 0; jj < bd_face.size(); jj++) {
// // // 
// // //               int jface = bd_face[jj];
// // //               // look for boundary faces
// // // 
// // //               unsigned faceDofs = el->GetNFACENODES(ielGeom2, jface, solType);
// // // 
// // //               vector  < vector  <  double> > faceCoordinates(dim);    // A matrix holding the face coordinates rowwise.
// // //               for(int k = 0; k < dim; k++) {
// // //                 faceCoordinates[k].resize(faceDofs);
// // //               }
// // //               for(unsigned i = 0; i < faceDofs; i++) {
// // //                 unsigned inode = el->GetIG(ielGeom2, jface, i);  // face-to-element local node mapping.
// // //                 for(unsigned k = 0; k < dim; k++) {
// // //                   faceCoordinates[k][i] =  x2[k][inode] - x_iqp_bdry[k]; // We extract the local coordinates on the face from local coordinates on the element.
// // //                 }
// // //               }
// // //               const unsigned div = 10;
// // //               vector  < vector  <  double> > interpCoordinates(dim);
// // //               for(int k = 0; k < dim; k++) {
// // //                 interpCoordinates[k].resize(div + 1); // set "4" as a parameter
// // //               }
// // //               for(unsigned n = 0; n <= div; n++) {
// // //                 for(int k = 0; k < dim; k++) {
// // //                   interpCoordinates[k][n] = faceCoordinates[k][0] + n * (faceCoordinates[k][1] - faceCoordinates[k][0]) /  div ;
// // //                 }
// // //               }
// // //               for(unsigned n = 0; n < div; n++) {
// // //                 double teta2 = atan2(interpCoordinates[1][n + 1], interpCoordinates[0][n + 1]);
// // //                 double teta1 = atan2(interpCoordinates[1][n], interpCoordinates[0][n]);
// // // 
// // //                 if(teta2 < teta1) teta2 += 2. * M_PI;
// // // 
// // //                 double delta_teta = teta2 - teta1;
// // // 
// // // 
// // //                 vector <double> mid_point;
// // //                 mid_point.resize(dim);
// // //                 for(unsigned k = 0; k < dim; k++) {
// // //                   mid_point[k] = (interpCoordinates[k][n + 1] + interpCoordinates[k][n]) * 0.5;
// // //                 }
// // //                 double dist2 = 0;
// // //                 for(int k = 0; k < dim; k++) {
// // //                   dist2 += mid_point[k] * mid_point[k];
// // //                 }
// // //                 double dist = sqrt(dist2);
// // //                 mixed_term1 += 2. * pow(dist, -  2. * s_frac) * (1. / (2. * s_frac)) * delta_teta;
// // //               }
// // //             }
// // // 
// // //             for(unsigned i = 0; i < nDof1; i++) {
// // //               for(unsigned j = 0; j < nDof1; j++) {
// // //                 KK_local_iel_mixed_num[ i * nDof1 + j ] += (C_ns / 2.) * check_limits * OP_Hhalf * phi1[i] * phi1[j] * weight1 * mixed_term1;
// // //               }
// // //               Res_local_iel_mixed_num[ i ] += (C_ns / 2.) * check_limits * OP_Hhalf * weight1 * phi1[i] * solX * mixed_term1;
// // //             }
          
      }          
      
  }
  
  
 }
  
  
  
  
    //********** BEGIN FRAC CONTROL *****************************************

  void control_eqn_bdry_fractional(const unsigned iproc,
                                   const unsigned nprocs,
                        MultiLevelProblem &    ml_prob,
                        MultiLevelSolution*    ml_sol,
                        const Solution*        sol,
                        const Mesh * msh,
                        const  LinearEquationSolver* pdeSys,
                        //-----------
                        CurrentElem < double > & geom_element_iel,
                        CurrentElem < double > & geom_element_jel,
                        const unsigned int solType_coords,
                        const unsigned int dim,
                        const unsigned int space_dim,
                        const unsigned int dim_bdry,
                        //-----------
                        const unsigned int n_unknowns,
                        const    vector < std::string > & Solname_Mat,
                        const    vector < unsigned > & SolFEType_Mat,
                        const vector < unsigned > & SolIndex_Mat,
                        const vector < unsigned > & SolPdeIndex,
                        vector < unsigned > & Sol_n_el_dofs_Mat, 
                        vector < vector < double > > & sol_eldofs_Mat,  
                        vector < vector < int > > & L2G_dofmap_Mat,
                        std::vector< int >    &   L2G_dofmap_Mat_AllVars,
                        const unsigned int maxSize,
                        //-----------
                        std::vector < double > & Res,
                        std::vector < double > & Jac,
                        //-----------
                        const unsigned int n_quantities,
                        vector < unsigned > SolFEType_quantities,
                        vector < unsigned > Sol_n_el_dofs_quantities,
                        //-----------
                        std::vector < std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > > elem_all,
                        //-----------
                        std::vector < std::vector < double > >  Jac_iel_bdry_iqp_bdry,
                        std::vector < std::vector < double > >  JacI_iel_bdry_iqp_bdry,
                        double detJac_iel_bdry_iqp_bdry,
                        double weight_iqp_bdry,
                        vector <double> phi_ctrl_iel_bdry_iqp_bdry,
                        vector <double> phi_ctrl_x_iel_bdry_iqp_bdry, 
                        vector <double> phi_coords_iel_bdry_iqp_bdry,
                        vector <double> phi_coords_x_iel_bdry_iqp_bdry, 
                        //-----------
                        std::vector < std::vector < double > >  Jac_jel_bdry_jqp_bdry,
                        std::vector < std::vector < double > >  JacI_jel_bdry_jqp_bdry,
                        double detJac_jel_bdry_jqp_bdry,
//                         double weight_jqp_bdry,
//                         vector <double> phi_ctrl_jel_bdry_jqp_bdry,    //it is a stored vector below
//                         vector <double> phi_ctrl_x_jel_bdry_jqp_bdry,    //it is a stored vector below 
                        //-----------
                        const unsigned int pos_mat_ctrl,
                        const unsigned int pos_sol_ctrl,
                        const unsigned int is_block_dctrl_ctrl_inside_bdry,
                        //-----------
                        SparseMatrix*  KK,
                        NumericVector* RES,
                        const bool assembleMatrix,
                        //-----------
                        const double alpha,
                        const double beta,
                        const unsigned int Nsplit,
                        const double s_frac,
                        const double check_limits,
                        const double C_ns,
                        const unsigned int OP_Hhalf,
                        const unsigned int OP_L2,
                        const unsigned int RHS_ONE,
                        const unsigned int UNBOUNDED,
                        const double EX_1,
                        const double EX_2,
                        const unsigned qrule_i,
                        const unsigned qrule_j
                       ) {
      
      
      
  unsigned solType =   SolFEType_Mat[pos_mat_ctrl];
  unsigned soluIndex = SolIndex_Mat[pos_mat_ctrl];
  unsigned soluPdeIndex = SolPdeIndex[pos_mat_ctrl];
  
  std::vector < double > sol_ctrl_iel;
  std::vector < double > sol_ctrl_jel;
    
  //------- geometry ---------------
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

//   vector < vector < double > > x1(dim);
  vector < vector < double > > x2(dim);
  for(unsigned k = 0; k < dim; k++) {
//     x1[k].reserve(maxSize);
    x2[k].reserve(maxSize);
  }
 
  //------- geometry ---------------
 
 
 
  //-------- local to global mappings --------------
  vector< int > l2gMap_iel;  l2gMap_iel.reserve(maxSize);
  vector< int > l2gMap_jel;  l2gMap_jel.reserve(maxSize);


  //-------- Local matrices and rhs --------------
  vector < double > Res_local_iel; Res_local_iel.reserve(maxSize);
  vector < double > KK_local_iel;  KK_local_iel.reserve(maxSize * maxSize);

//   Local matrices and rhs for adaptive quadrature (iel == jel)
  vector < double > Res_local_iel_refined; Res_local_iel_refined.reserve(maxSize);
  vector < double > KK_local_iel_refined;   KK_local_iel_refined.reserve(maxSize * maxSize);

//   Local matrices and rhs for the mixed internal-external integral term (both adaptive and non-adaptive)
  vector < double > Res_local_iel_mixed_num;  Res_local_iel_mixed_num.reserve(maxSize);
  vector < double > KK_local_iel_mixed_num;   KK_local_iel_mixed_num.reserve(maxSize * maxSize);

//   Non local matrices and vectors for H^s laplacian operator
  vector< double >         Res_nonlocal_iel;  Res_nonlocal_iel.reserve(maxSize);
  vector< double >         Res_nonlocal_jel;  Res_nonlocal_jel.reserve(maxSize);

  vector < double > KK_nonlocal_iel_iel;  KK_nonlocal_iel_iel.reserve(maxSize * maxSize);
  vector < double > KK_nonlocal_iel_jel;  KK_nonlocal_iel_jel.reserve(maxSize * maxSize);
  vector < double > KK_nonlocal_jel_iel;  KK_nonlocal_jel_iel.reserve(maxSize * maxSize);
  vector < double > KK_nonlocal_jel_jel;  KK_nonlocal_jel_jel.reserve(maxSize * maxSize); 
 
 //----------------------
  KK->zero();
  RES->zero(); 
 //----------------------
  
  
   for(int kproc = 0; kproc < nprocs; kproc++) {
       
    for(int jel = msh->_elementOffset[kproc]; jel < msh->_elementOffset[kproc + 1]; jel++) {

// ***************************************
// ******* only control volume elements *************
// ***************************************
// --- geometry        
        geom_element_jel.set_coords_at_dofs_and_geom_type(jel, solType_coords);

        short unsigned jelGeom = geom_element_jel.geom_type();
              
        geom_element_jel.set_elem_center(jel, solType_coords);
// --- geometry        
        
    //************ set control flag *********************
  int control_flag_jel = 0;
        control_flag_jel = ControlDomainFlag_bdry(geom_element_jel.get_elem_center());
  std::vector<int> control_node_flag_j(Sol_n_el_dofs_Mat[pos_mat_ctrl], 0);
 //*************************************************** 
      
	// Perform face loop over elements that contain some control face
	if (control_flag_jel == 1) {
// ***************************************
// ******* only control volume elements *************
// ***************************************
         
         

      
// ***************************************
// ******* jel-related stuff *************
// ***************************************
        
// --- 1 - geometry -----------------
      
// --- geom_el type
      short unsigned jelGeom;
      unsigned nDof_jel_coords;
      if(kproc == iproc) {
          jelGeom = msh->GetElementType(jel); 
        nDof_jel_coords = msh->GetElementDofNumber(jel, solType_coords);    // number of coordinate element dofs
    }
      MPI_Bcast(&jelGeom, 1, MPI_UNSIGNED_SHORT, kproc, MPI_COMM_WORLD);
      MPI_Bcast(&nDof_jel_coords, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
// --- geom_el type
      
// --- coords - one way
      for(int k = 0; k < dim; k++) {  x2[k].resize(nDof_jel_coords);  }

      
      if(kproc == iproc) {
        for(unsigned j = 0; j < nDof_jel_coords; j++) {
          unsigned xDof  = msh->GetSolutionDof(j, jel, solType_coords);  // global to global mapping between coordinates node and coordinate dof
          for(unsigned k = 0; k < dim; k++) {
            x2[k][j] = (*msh->_topology->_Sol[k])(xDof);  // global extraction and local storage for the element coordinates
          }
        }
       }
       
      for(unsigned k = 0; k < dim; k++) {
        MPI_Bcast(& x2[k][0], nDof_jel_coords, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }
// --- coords - one way
      

// --- coords - other way
      if(kproc == iproc) {
        geom_element_jel.set_coords_at_dofs_and_geom_type(jel, solType_coords);
      }
      for(unsigned k = 0; k < dim; k++) {
        MPI_Bcast(& geom_element_jel.get_coords_at_dofs()[k][0], nDof_jel_coords, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }
      for(unsigned k = 0; k < space_dim; k++) {
        MPI_Bcast(& geom_element_jel.get_coords_at_dofs_3d()[k][0], nDof_jel_coords, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }
// --- coords - other way

// --- 1 - geometry -----------------


// --- 2 - solution -----------------
      unsigned nDof_jel;

      if(kproc == iproc) {
        nDof_jel  = msh->GetElementDofNumber(jel, solType);    // number of solution element dofs
      }

      MPI_Bcast(&nDof_jel, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
      
      
      sol_ctrl_jel.resize(nDof_jel);
      
      if(kproc == iproc) {
        for(unsigned j = 0; j < nDof_jel; j++) {
          unsigned jDof  = msh->GetSolutionDof(j, jel, solType);  // global to global mapping between coordinates node and coordinate dof
          sol_ctrl_jel[j] = (*sol->_Sol[soluIndex])(jDof);  // global extraction and local storage for the element coordinates
        }
      }
      
      MPI_Bcast(& sol_ctrl_jel[0], nDof_jel, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
// --- 2 - solution -----------------

      
// --- 3 - l2GMap -----------------
        l2gMap_jel.resize(nDof_jel);

      // local storage of global mapping and solution ********************
      if(kproc == iproc) {
        for(unsigned j = 0; j < nDof_jel; j++) {
          l2gMap_jel[j] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, j, jel);  // global to global mapping between solution node and pdeSys dof
        }
      }
      MPI_Bcast(&l2gMap_jel[0], nDof_jel, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
      // ******************************************************************
// --- 3 - l2GMap -----------------

// ***************************************
// ******* jel-related stuff - end *************
// ***************************************


      
//------------------------------------        
//------------ jface opening ---------        
//------------------------------------        
	  // loop on faces of the current element
	  for(unsigned jface = 0; jface < msh->GetElementFaceNumber(jel); jface++) {
          
       const unsigned jelGeom_bdry = msh->GetElementFaceType(jel, jface);    
       
       geom_element_jel.set_coords_at_dofs_bdry_3d(jel, jface, solType_coords);
 
       geom_element_jel.set_elem_center_bdry_3d();

	    // look for boundary faces
            const int bdry_index_j = msh->el->GetFaceElementIndex(jel, jface);
	    // look for control faces
	      const unsigned int face_in_rectangle_domain_j = -( msh->el->GetFaceElementIndex(jel,jface) + 1);
          
         
	    if( check_if_boundary_control_face(bdry_index_j, face_in_rectangle_domain_j) ) {
//------------------------------------        
//------------ jface opening ---------    
//------------------------------------        
            
            

// // // //---- Quadrature in jqp_bdry, preparation right before iel ------- 
// // // //---- Quadrature in jqp_bdry, preparation right before iel ------- 
// // // //---- Quadrature in jqp_bdry, preparation right before iel ------- 

// Wait a second... If I want to prepare the jqp_bdry loop, I must be inside jface as well...
// So now it seems to me that I have to do jel jface iel iface instead...
// Previously it was jel iel iqp jqp
// Now, it has to be jel jface  - iel iface - iqp_bdry jqp_bdry
// The two quadrature loops must be the innermost. In this way you exclude all non-needed volume elements and all non-needed faces, so that you minimize the number of inner ifs. You keep them outside as much as possible
// There will be a storage of jqp_bdry
            
      const unsigned n_jqp_bdry = ml_prob.GetQuadratureRuleMultiple(qrule_j, jelGeom_bdry).GetGaussPointsNumber();

      vector < vector < double > > x_jqp_bdry(n_jqp_bdry);
      vector < double > weight_jqp_bdry(n_jqp_bdry);
      vector < vector < double > > phi_ctrl_jel_bdry_jqp_bdry(n_jqp_bdry);
      vector < vector < double > > phi_ctrl_x_jel_bdry_jqp_bdry(n_jqp_bdry);
      vector < vector < double > > phi_coords_jel_bdry_jqp_bdry(n_jqp_bdry);
      vector < vector < double > > phi_coords_x_jel_bdry_jqp_bdry(n_jqp_bdry);
      std::vector< double > sol_ctrl_jqp_bdry(n_jqp_bdry, 0.);

         for(unsigned jqp_bdry = 0; jqp_bdry < n_jqp_bdry; jqp_bdry++) {

    elem_all[qrule_j][jelGeom_bdry][solType_coords]->JacJacInv(geom_element_jel.get_coords_at_dofs_bdry_3d(), jqp_bdry, Jac_jel_bdry_jqp_bdry, JacI_jel_bdry_jqp_bdry, detJac_jel_bdry_jqp_bdry, space_dim);
    
    weight_jqp_bdry[jqp_bdry] = detJac_jel_bdry_jqp_bdry * ml_prob.GetQuadratureRuleMultiple(qrule_j, jelGeom_bdry).GetGaussWeightsPointer()[jqp_bdry];

    elem_all[qrule_j][jelGeom_bdry][SolFEType_quantities[pos_sol_ctrl]] ->shape_funcs_current_elem(jqp_bdry, JacI_jel_bdry_jqp_bdry, phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry], phi_ctrl_x_jel_bdry_jqp_bdry[jqp_bdry], boost::none, space_dim);
            
    elem_all[qrule_j][jelGeom_bdry][solType_coords] ->shape_funcs_current_elem(jqp_bdry, JacI_jel_bdry_jqp_bdry, phi_coords_jel_bdry_jqp_bdry[jqp_bdry], phi_coords_x_jel_bdry_jqp_bdry[jqp_bdry], boost::none, space_dim);

//========== compute gauss quantities on the boundary ===============================================
//--- geom
           x_jqp_bdry[jqp_bdry].assign(dim, 0.);
          
            for(unsigned d = 0; d < dim; d++) {
	      for (int j_bdry = 0; j_bdry < geom_element_jel.get_coords_at_dofs_bdry_3d()[d].size(); j_bdry++)  {
			
              x_jqp_bdry[jqp_bdry][d] += geom_element_jel.get_coords_at_dofs_bdry_3d()[d][j_bdry] * phi_coords_jel_bdry_jqp_bdry[jqp_bdry][j_bdry];

		      }
            }
//--- geom
    
//--- solution
    sol_ctrl_jqp_bdry[jqp_bdry] = 0.;
	      for (int j_bdry = 0; j_bdry < Sol_n_el_dofs_quantities[pos_sol_ctrl]; j_bdry++)  {
		    unsigned int j_vol = msh->GetLocalFaceVertexIndex(jel, jface, j_bdry);
			
			sol_ctrl_jqp_bdry[jqp_bdry] +=  /*sol_eldofs_Mat[pos_mat_ctrl]*/sol_ctrl_jel[j_vol] * phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][j_bdry];

		      }
//--- solution
//========== compute gauss quantities on the boundary ================================================


        }  //jqp_bdry

// // // //---- Quadrature in jqp_bdry, preparation right before iel ------- 
// // // //---- Quadrature in jqp_bdry, preparation right before iel ------- 
// // // //---- Quadrature in jqp_bdry, preparation right before iel ------- 

            
// ---- boundary faces in jface: compute and broadcast - BEGIN ----
// This is only needed for when the boundary is a 2D face. We'll look at it later on
// look for what face of jface are on the boundary of the domain
// I believe we have to see deeply how we can extend this to the boundary case

// ---- boundary faces in jface: compute and broadcast - END ----    


              
       for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
           
           
// ***************************************
// ******* only control volume elements *************
// ***************************************
// --- geometry        
        geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);

        short unsigned ielGeom = geom_element_iel.geom_type();
              
        geom_element_iel.set_elem_center(iel, solType_coords);
// --- geometry        

   //************ set control flag *********************
  int control_flag_iel = 0;
        control_flag_iel = ControlDomainFlag_bdry(geom_element_iel.get_elem_center());
  std::vector<int> control_node_flag_i(Sol_n_el_dofs_Mat[pos_mat_ctrl], 0);
 //*************************************************** 
      
	// Perform face loop over elements that contain some control face
	if (control_flag_iel == 1) {
// ***************************************
// ******* only control volume elements *************
// ***************************************


        
// ***************************************
// ******* iel-related stuff *************
// ***************************************
        
      
// --- l2GMap
        unsigned nDof_iel  = msh->GetElementDofNumber(iel, solType);
        
        l2gMap_iel.resize(nDof_iel);
        for(unsigned i = 0; i < nDof_iel; i++) {
          l2gMap_iel[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);
        }
// --- l2GMap

// --- element matrix and vector resizing
        KK_nonlocal_iel_iel.assign(nDof_iel * nDof_iel, 0.);   //resize
        KK_nonlocal_iel_jel.assign(nDof_iel * nDof_jel, 0.);   //resize
        KK_nonlocal_jel_iel.assign(nDof_jel * nDof_iel, 0.);   //resize
        KK_nonlocal_jel_jel.assign(nDof_jel * nDof_jel, 0.);   //resize
//         Res_nonlocal.assign(nDof_iel, 0);    //resize
        Res_nonlocal_iel.assign(nDof_iel, 0.);    //resize
        Res_nonlocal_jel.assign(nDof_jel, 0.);    //resize

        Res_local_iel_mixed_num.assign(nDof_iel, 0.);    //resize
        KK_local_iel_mixed_num.assign(nDof_iel * nDof_iel, 0.);

        if( check_if_same_elem(iel, jel) ) {
          Res_local_iel.assign(nDof_iel, 0.);    //resize
          KK_local_iel.assign(nDof_iel * nDof_iel, 0.);
          if(Nsplit != 0) {
//             Vectors and matrices for adaptive quadrature
            Res_local_iel_refined.assign(nDof_iel, 0.);    //resize
            KK_local_iel_refined.assign(nDof_iel * nDof_iel, 0.);
          }
        }
// --- element matrix and vector resizing 


// --- solution        

       sol_ctrl_iel.resize(nDof_iel);
        
        for(unsigned i = 0; i < nDof_iel; i++) {
          unsigned iDof  = msh->GetSolutionDof(i, iel, solType);  // global to global mapping between coordinates node and coordinate dof
          sol_ctrl_iel[i] = (*sol->_Sol[soluIndex])(iDof);  // global extraction and local storage for the element coordinates
        }
// --- solution        

// ***************************************
// ******* iel-related stuff - end *************
// ***************************************
        


        
//------------ iface opening ---------        
	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
          
       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, iface);    
       
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();

	    // look for boundary faces
            const int bdry_index_i = msh->el->GetFaceElementIndex(iel, iface);
	    // look for control faces
	      const unsigned int face_in_rectangle_domain_i = -( msh->el->GetFaceElementIndex(iel,iface)+1);
          
	    if( check_if_boundary_control_face( bdry_index_i, face_in_rectangle_domain_i) ) {
//------------ iface opening ---------        
		


             
     //**** Adaptive preparation - BEGIN ******** 
     std::vector < std::vector < double > >  Jac_kel_bdry_kqp_bdry;
     std::vector < std::vector < double > >  JacI_kel_bdry_kqp_bdry;
     double detJac_kel_bdry_kqp_bdry;
      vector < double >  x_kqp_bdry;
      double  weight_kqp_bdry;
      vector < double >  phi_ctrl_kel_bdry_kqp_bdry;
      vector < double >  phi_ctrl_x_kel_bdry_kqp_bdry;
      vector < double >  phi_coords_kel_bdry_kqp_bdry;
      vector < double >  phi_coords_x_kel_bdry_kqp_bdry;
      double sol_ctrl_kqp_bdry = 0.;
//         double weight3;
//         vector < double > phi3;

//        Evaluating coarse FE functions on Quadrature Points of the "sub-elements"
       std::vector < std::vector < std::vector <double > > > aP(3);  //[NFE_FAMS][DIM==3][N_DOFS]
        if(Nsplit != 0) {
          for(unsigned fe_type = 0; fe_type < solType + 1; fe_type++) { //loop up to the FE type + 1 of the unknown
            ProjectNodalToPolynomialCoefficients(aP[fe_type], geom_element_iel.get_coords_at_dofs_bdry_3d(), ielGeom_bdry, fe_type) ;         ///@todo check this!!!!!!!!!!!!
          }
        }                      
     //**** Adaptive preparation - END ********  
                      
              //Quadrature loop initialization - BEGIN
        const unsigned n_iqp_bdry = ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussPointsNumber();
         double sol_ctrl_iqp_bdry = 0.;
              //Quadrature loop initialization - END
    
         
		for(unsigned iqp_bdry = 0; iqp_bdry < n_iqp_bdry; iqp_bdry++) {
            
    elem_all[qrule_i][ielGeom_bdry][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_bdry_3d(), iqp_bdry, Jac_iel_bdry_iqp_bdry, JacI_iel_bdry_iqp_bdry, detJac_iel_bdry_iqp_bdry, space_dim);
    
    weight_iqp_bdry = detJac_iel_bdry_iqp_bdry * ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussWeightsPointer()[iqp_bdry];

    elem_all[qrule_i][ielGeom_bdry][SolFEType_quantities[pos_sol_ctrl]] ->shape_funcs_current_elem(iqp_bdry, JacI_iel_bdry_iqp_bdry, phi_ctrl_iel_bdry_iqp_bdry, phi_ctrl_x_iel_bdry_iqp_bdry, boost::none, space_dim);
            
    elem_all[qrule_i][ielGeom_bdry][solType_coords] ->shape_funcs_current_elem(iqp_bdry, JacI_iel_bdry_iqp_bdry, phi_coords_iel_bdry_iqp_bdry, phi_coords_x_iel_bdry_iqp_bdry, boost::none, space_dim);

//========== compute gauss quantities on the boundary ===============================================
//--- geom
          vector < double > x_iqp_bdry(dim, 0.);  ///@todo is this dim or dim_bdry?

            for(unsigned d = 0; d < x_iqp_bdry.size(); d++) {
	      for (int i_bdry = 0; i_bdry < geom_element_iel.get_coords_at_dofs_bdry_3d()[d].size(); i_bdry++)  {
			
              x_iqp_bdry[d] += geom_element_iel.get_coords_at_dofs_bdry_3d()[d][i_bdry] * phi_coords_iel_bdry_iqp_bdry[i_bdry];

		      }
            }
//--- geom
    
//--- solution
    sol_ctrl_iqp_bdry = 0.;
	      for (int i_bdry = 0; i_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
			
			sol_ctrl_iqp_bdry +=  /*sol_eldofs_Mat[pos_mat_ctrl]*/sol_ctrl_iel[i_vol] * phi_ctrl_iel_bdry_iqp_bdry[i_bdry];

		      }
//--- solution
//========== compute gauss quantities on the boundary ================================================


  
      //============  Non-fractional assembly - BEGIN ==================
          
            if( check_if_same_elem_bdry(iel, jel, iface, jface) ) {
              
       //============  Mass assembly - BEGIN ==================
           for(unsigned l_bdry = 0; l_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); l_bdry++) {
               		    unsigned int l_vol = msh->GetLocalFaceVertexIndex(iel, iface, l_bdry);

              for(unsigned m_bdry = 0; m_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); m_bdry++) {
               		    unsigned int m_vol = msh->GetLocalFaceVertexIndex(iel, iface, m_bdry);
                KK_local_iel[ l_vol * nDof_iel + m_vol ] += OP_L2 * phi_ctrl_iel_bdry_iqp_bdry[l_bdry] * phi_ctrl_iel_bdry_iqp_bdry[m_bdry] * weight_iqp_bdry;
              }
              double mass_res_i = phi_ctrl_iel_bdry_iqp_bdry[l_bdry] * sol_ctrl_iqp_bdry ;
              Res_local_iel[ l_vol ] += OP_L2 * weight_iqp_bdry * mass_res_i ;
              Res_local_iel[ l_vol ] += - RHS_ONE * weight_iqp_bdry * (phi_ctrl_iel_bdry_iqp_bdry[l_bdry] * (-1.) /** ( sin(2 * acos(0.0) * x1[0][l_bdry])) * ( sin(2 * acos(0.0) * x1[1][l_bdry]))*/);
            }
        //============  Mass assembly - END ==================
         
      //============  Laplacian assembly - BEGIN ==================
      //============  Laplacian assembly - END ==================
      
    
            } 
      //============  Non-fractional assembly - END ==================
             
             
      //============  Fractional assembly - BEGIN ==================
        if(OP_Hhalf != 0) {
                 
      //============ Same elem, && Adaptive quadrature - BEGIN ==================
       if( check_if_same_elem_bdry(iel, jel, iface, jface) && Nsplit != 0) {
                
// // // // //                       /*const*/ short unsigned kelGeom_bdry = ielGeom_bdry;
// // // // //            
// // // // //                       const unsigned n_kqp_bdry = ml_prob.GetQuadratureRuleMultiple(qrule_k, kelGeom_bdry).GetGaussPointsNumber();
// // // // //                 
// // // // //                 
// // // // //              std::cout.precision(14);
// // // // //               std::vector< std::vector< std::vector<double> > > x3;
// // // // // 
// // // // //               for(unsigned split = 0; split <= Nsplit; split++) {
// // // // // 
// // // // // 
// // // // //                 if (dim_bdry/*dim*/ == 1) GetElementPartition1D(x_iqp_bdry, geom_element_iel.get_coords_at_dofs_bdry_3d(), split, Nsplit, x3, space_dim);
// // // // //                 else if (dim_bdry/*dim*/ == 2) {
// // // // //                   //GetElementPartition2D(x_iqp_bdry, geom_element_iel.get_coords_at_dofs_bdry_3d(), split, Nsplit, x3);
// // // // //                   GetElementPartitionQuad(x_iqp_bdry, geom_element_iel.get_coords_at_dofs_bdry_3d(), split, Nsplit, x3);
// // // // //                 }
// // // // // 
// // // // //                 //for(unsigned r = 0; r < size_part; r++) {
// // // // //                 for(unsigned r = 0; r < x3.size(); r++) {
// // // // // 
// // // // // 
// // // // //                   for(unsigned k_qp_bdry = 0; k_qp_bdry < n_kqp_bdry; k_qp_bdry++) {
// // // // // 
// // // // // // ********* PREPARATION PART - BEGIN ***************
// // // // //     elem_all[kelGeom_bdry][solType_coords]->JacJacInv(x3[r]/*geom_element_iel.get_coords_at_dofs_bdry_3d()*/, k_qp_bdry, Jac_kel_bdry_kqp_bdry, JacI_kel_bdry_kqp_bdry, detJac_kel_bdry_kqp_bdry, space_dim);
// // // // //     
// // // // //     weight_kqp_bdry = detJac_kel_bdry_kqp_bdry * ml_prob.GetQuadratureRuleMultiple(qrule_k, kelGeom_bdry).GetGaussWeightsPointer()[k_qp_bdry];
// // // // // 
// // // // //     elem_all[kelGeom_bdry][SolFEType_quantities[pos_sol_ctrl]] ->shape_funcs_current_elem(k_qp_bdry, JacI_kel_bdry_kqp_bdry, phi_ctrl_kel_bdry_kqp_bdry, phi_ctrl_x_kel_bdry_kqp_bdry, boost::none, space_dim);
// // // // //             
// // // // //     elem_all[kelGeom_bdry][solType_coords] ->shape_funcs_current_elem(k_qp_bdry, JacI_kel_bdry_kqp_bdry, phi_coords_kel_bdry_kqp_bdry, phi_coords_x_kel_bdry_kqp_bdry, boost::none, space_dim);
// // // // // 
// // // // // //--- geom
// // // // //           vector < double > x_kqp_bdry(dim, 0.);  ///@todo is this dim or dim_bdry?
// // // // // 
// // // // //             for(unsigned d = 0; d < x_kqp_bdry.size(); d++) {
// // // // // 	      for (int k_bdry = 0; k_bdry < /*geom_element_iel.get_coords_at_dofs_bdry_3d()[d].size()*/phi_coords_kel_bdry_kqp_bdry.size(); k_bdry++)  {
// // // // // 			
// // // // //               x_kqp_bdry[d] += x3[r][d][k_bdry]/*geom_element_iel.get_coords_at_dofs_bdry_3d()[d][i_vol]*/ * phi_coords_kel_bdry_kqp_bdry[k_bdry];
// // // // // 
// // // // // 		      }
// // // // //             }
// // // // // //--- geom
// // // // //     
// // // // //                     std::vector<double> xi3(dim, 0.);
// // // // // 
// // // // //                     GetClosestPointInReferenceElement(geom_element_iel.get_coords_at_dofs_bdry_3d(), x_kqp_bdry, kelGeom_bdry, xi3);
// // // // //                     GetInverseMapping(solType_coords, kelGeom_bdry, aP, x_kqp_bdry, xi3, 1000);  ///@todo generalize to rectangular Jacobian
// // // // // 
// // // // //                     msh->_finiteElement[kelGeom_bdry][solType]->GetPhi(phi_ctrl_kel_bdry_kqp_bdry, xi3);
// // // // // 
// // // // //                     double solY3 = 0.;
// // // // //                     for(unsigned i_bdry = 0; i_bdry < phi_ctrl_kel_bdry_kqp_bdry.size()/*nDof_iel*/; i_bdry++) {
// // // // // 		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
// // // // //                       solY3 += sol_ctrl_iel[i_vol] * phi_ctrl_kel_bdry_kqp_bdry[i_bdry];
// // // // //                     }
    
// // // 
// // // 
// // //                     msh->_finiteElement[ielGeom][solType]->Jacobian(x3[r], k_qp_bdry, weight_kqp_bdry, phi3, phi_x);
// // // 
// // //                     vector < double > xg3(dim, 0.);
// // // 
// // //                     for(unsigned d = 0; d < dim; d++) {
// // //                       for(unsigned i = 0; i < nDof_iel; i++) {
// // //                         xg3[d] += x3[r][d][i] * phi3[i];
// // //                       }
// // //                     }
// // // 
// // //                     std::vector<double> xi3(dim, 0.);
// // // 
// // //                     GetClosestPointInReferenceElement(geom_element_iel.get_coords_at_dofs_bdry_3d(), xg3, ielGeom, xi3);
// // //                     GetInverseMapping(solType, ielGeom, aP, xg3, xi3, 1000);
// // // 
// // //                     msh->_finiteElement[ielGeom][solType]->GetPhi(phi3, xi3);
// // // 
// // //                     double solY3 = 0.;
// // //                     for(unsigned i = 0; i < nDof_iel; i++) {
// // //                       solY3 += sol_ctrl_iel[i] * phi3[i];
// // //                     }
// ********* PREPARATION PART - END ***************

// ********* BOUNDED PART - BEGIN ***************
// // // 
// // //                     double dist_xyz3 = 0;
// // //                     for(unsigned k = 0; k < dim; k++) {
// // //                       dist_xyz3 += (x_iqp_bdry[k] - xg3[k]) * (x_iqp_bdry[k] - xg3[k]);
// // //                     }
// // // 
// // //                     const double denom3 = pow(dist_xyz3, (double)((dim / 2.) + s_frac));
// // // 
// // //                     for(unsigned i = 0; i < nDof_iel; i++) {
// // // 
// // //                       Res_local_iel_refined[ i ]    +=      - (0.5 * C_ns) * OP_Hhalf * check_limits *
// // //                                                         ((sol_ctrl_iqp_bdry - solY3) * (phi_ctrl_iel_bdry_iqp_bdry[i] - phi3[i]) * weight_kqp_bdry / denom3
// // //                                                         ) * weight_iqp_bdry ;
// // // 
// // //                       for(unsigned j = 0; j < nDof_jel; j++) {
// // //                         KK_local_iel_refined[ i * nDof_jel + j ] += (0.5 * C_ns) * OP_Hhalf * check_limits *
// // //                                                             ((phi_ctrl_iel_bdry_iqp_bdry[j] - phi3[j]) * (phi_ctrl_iel_bdry_iqp_bdry[i] - phi3[i]) * weight_kqp_bdry / denom3
// // //                                                             ) * weight_iqp_bdry ;
// // // 
// // //                       }
// // //                     }
// ********* BOUNDED PART - END ***************
                      
// ********* UNBOUNDED PART - BEGIN ***************
// ********* UNBOUNDED PART - END ***************
                                         
// // // // //                   } //end k_qp_bdry
// // // // //                 } //end r
// // // // //               }  //end split
              
              
                
                
                
            }  //end iel == jel && Nsplit != 0
      //============ Same elem, && Adaptive quadrature - END ==================
              
             
      //============ Either different elements, or lack of adaptivity (so all elements) - BEGIN ==================
        else {  //  if(iel != jel || Nsplit == 0) 
            
// ********* UNBOUNDED PART - BEGIN ***************
               mixed_integral(UNBOUNDED,
                              dim,
                              dim_bdry,
                              EX_1,
                              EX_2,
                              weight_iqp_bdry,
                              x_iqp_bdry,
                              phi_ctrl_iel_bdry_iqp_bdry,
                              sol_ctrl_iqp_bdry,
                              s_frac,
                              check_limits,
                              C_ns,
                              OP_Hhalf,
                              nDof_iel,
                              KK_local_iel,
                              Res_local_iel
                             ); 
              
// ********* UNBOUNDED PART - END ***************
            
// ********* BOUNDED PART - BEGIN ***************
            for(unsigned jqp_bdry = 0; jqp_bdry < n_jqp_bdry; jqp_bdry++) {

              double dist_xyz = 0.;
              for(unsigned d = 0; d < x_jqp_bdry.size(); d++) {
                dist_xyz += (x_iqp_bdry[d] - x_jqp_bdry[jqp_bdry][d]) * (x_iqp_bdry[d] - x_jqp_bdry[jqp_bdry][d]);
              }

              const double denom = pow(dist_xyz, (double)(  0.5 * /*dim*/dim_bdry + s_frac));
              
              const double common_weight = (0.5 * C_ns) * OP_Hhalf * check_limits * weight_iqp_bdry * weight_jqp_bdry[jqp_bdry]  / denom;

//               for(unsigned i = 0; i < nDof_iel; i++) {
           for(unsigned l_bdry = 0; l_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); l_bdry++) { //dofs of test function
               		    unsigned int l_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, l_bdry);
               		    unsigned int l_vol_jel = msh->GetLocalFaceVertexIndex(jel, jface, l_bdry);

                Res_nonlocal_iel[ l_vol_iel ]      +=      - common_weight * (sol_ctrl_iqp_bdry - sol_ctrl_jqp_bdry[jqp_bdry]) * (phi_ctrl_iel_bdry_iqp_bdry[l_bdry]);

                Res_nonlocal_jel[ l_vol_jel ]      +=      - common_weight * (sol_ctrl_iqp_bdry - sol_ctrl_jqp_bdry[jqp_bdry]) * (- phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][l_bdry]);

               
//                 for(unsigned j = 0; j < nDof_jel; j++) {
           for(unsigned m_bdry = 0; m_bdry < phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry].size(); m_bdry++) { //dofs of unknown function
               		    unsigned int m_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, m_bdry);
               		    unsigned int m_vol_jel = msh->GetLocalFaceVertexIndex(jel, jface, m_bdry);

                  KK_nonlocal_iel_iel[ l_vol_iel * nDof_jel + m_vol_iel ] += common_weight *          phi_ctrl_iel_bdry_iqp_bdry[m_bdry]            *    phi_ctrl_iel_bdry_iqp_bdry[l_bdry];

                  KK_nonlocal_iel_jel[ l_vol_iel * nDof_jel + m_vol_jel ] += common_weight * (- 1.) * phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][m_bdry]  *    phi_ctrl_iel_bdry_iqp_bdry[l_bdry];

                  KK_nonlocal_jel_iel[ l_vol_jel * nDof_jel + m_vol_iel ] += common_weight * (- 1.) * phi_ctrl_iel_bdry_iqp_bdry[m_bdry]            *   phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][l_bdry];

                  KK_nonlocal_jel_jel[ l_vol_jel * nDof_jel + m_vol_jel ] += common_weight *          phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][m_bdry]  *    phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][l_bdry];


                  }
                }
                
              } //endl jqp_bdry loop
// ********* BOUNDED PART - END ***************
            
            
         } //end if(iel != jel || Nsplit == 0)
      //============ Either different elements, or lack of adaptivity (so all elements) - END ==================

        
         } //OP_Hhalf != 0              
      //============  Fractional assembly - END ==================
            
      }   //end iqp_bdry
              
//----- iface ---        
        } //end if(bdry_index_i < 0) //end if(face_in_rectangle_domain_i == FACE_FOR_CONTROL)
      } //end iface
//----- iface ---        




                
//============ add to global - BEGIN ==================
// // multiply everything by -1.? Don't think so
// std::transform(KK_local_iel.begin(), KK_local_iel.end(), KK_local_iel.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(Res_local_iel.begin(), Res_local_iel.end(), Res_local_iel.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// 
// std::transform(KK_local_iel_refined.begin(), KK_local_iel_refined.end(), KK_local_iel_refined.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(Res_local_iel_refined.begin(), Res_local_iel_refined.end(), Res_local_iel_refined.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// 
// std::transform(KK_local_iel_mixed_num.begin(), KK_local_iel_mixed_num.end(), KK_local_iel_mixed_num.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(Res_local_iel_mixed_num.begin(), Res_local_iel_mixed_num.end(), Res_local_iel_mixed_num.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// 
// std::transform(KK_nonlocal_iel_iel.begin(), KK_nonlocal_iel_iel.end(), KK_nonlocal_iel_iel.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(KK_nonlocal_iel_jel.begin(), KK_nonlocal_iel_jel.end(), KK_nonlocal_iel_jel.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(KK_nonlocal_jel_iel.begin(), KK_nonlocal_jel_iel.end(), KK_nonlocal_jel_iel.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(KK_nonlocal_jel_jel.begin(), KK_nonlocal_jel_jel.end(), KK_nonlocal_jel_jel.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// 
// std::transform(Res_nonlocal_iel.begin(), Res_nonlocal_iel.end(), Res_nonlocal_iel.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(Res_nonlocal_jel.begin(), Res_nonlocal_jel.end(), Res_nonlocal_jel.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// // multiply everything by -1.

if( check_if_same_elem(iel, jel) ) {
          KK->add_matrix_blocked(KK_local_iel, l2gMap_iel, l2gMap_iel);
          RES->add_vector_blocked(Res_local_iel, l2gMap_iel);

          if(Nsplit != 0) {
            KK->add_matrix_blocked(KK_local_iel_refined, l2gMap_iel, l2gMap_iel);
            RES->add_vector_blocked(Res_local_iel_refined, l2gMap_iel);
          }
        }

        KK->add_matrix_blocked(KK_local_iel_mixed_num, l2gMap_iel, l2gMap_iel);
        RES->add_vector_blocked(Res_local_iel_mixed_num, l2gMap_iel);

        KK->add_matrix_blocked(KK_nonlocal_iel_iel, l2gMap_iel, l2gMap_iel);
        KK->add_matrix_blocked(KK_nonlocal_jel_iel, l2gMap_jel, l2gMap_iel);
        KK->add_matrix_blocked(KK_nonlocal_iel_jel, l2gMap_iel, l2gMap_jel);
        KK->add_matrix_blocked(KK_nonlocal_jel_jel, l2gMap_jel, l2gMap_jel);
// Since 1 is dense and 3 are sparse, and the dense dofs are 30, we should have at most 3x9 + 30 = 57, but in the sparsity print it shows 30. That's the problem


        RES->add_vector_blocked(Res_nonlocal_iel, l2gMap_iel);
        RES->add_vector_blocked(Res_nonlocal_jel, l2gMap_jel);
//============ add to global - END ==================
        
        
//----- iel ---        
    } //end control elem flag i (control_flag_iel == 1)
  } //end iel
//----- iel ---        


//----- jface ---        
        } //end if(bdry_index_j < 0)//end if(face_in_rectangle_domain_j == FACE_FOR_CONTROL)
      } //end jface
 //----- jface ---   
 
     
//----- jel ---        
    }  //end control elem flag j (control_flag_iel == 1
   } //end jel
//----- jel ---        
   
 } //end kproc
    
    
  
  }
    
  //********** END FRAC CONTROL *****************************************
  
  
    
  void control_eqn_bdry(const unsigned iproc,
                        MultiLevelProblem &    ml_prob,
                        MultiLevelSolution*    ml_sol,
                        const Solution*        sol,
                        const Mesh * msh,
                        const  LinearEquationSolver* pdeSys,
                        //-----------
                        CurrentElem < double > & geom_element_iel,
                        const unsigned int solType_coords,
                        const unsigned int space_dim,
                        //-----------
                        const unsigned int n_unknowns,
                        const    vector < std::string > & Solname_Mat,
                        const    vector < unsigned > & SolFEType_Mat,
                        const vector < unsigned > & SolIndex_Mat,
                        const vector < unsigned > & SolPdeIndex,
                        vector < unsigned > & Sol_n_el_dofs_Mat, 
                        vector < vector < double > > & sol_eldofs_Mat,  
                        vector < vector < int > > & L2G_dofmap_Mat,
                        std::vector< int >    &   L2G_dofmap_Mat_AllVars,
                        //-----------
                        std::vector < double > & Res,
                        std::vector < double > & Jac,
                        //-----------
                        const unsigned int n_quantities,
                        vector < unsigned > SolFEType_quantities,
                        vector < unsigned > Sol_n_el_dofs_quantities,
                        //-----------
                        std::vector < std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > > elem_all,
                        std::vector < std::vector < double > >  Jac_qp_bdry,
                        std::vector < std::vector < double > >  JacI_qp_bdry,
                        double detJac_qp_bdry,
                        double weight_bdry,
                        vector <double> phi_ctrl_iel_bdry_iqp_bdry,
                        vector <double> phi_ctrl_x_iel_bdry_iqp_bdry, 
                        //-----------
                        const unsigned int pos_mat_ctrl,
                        const unsigned int pos_sol_ctrl,
                        const unsigned int is_block_dctrl_ctrl_inside_bdry,
                        //-----------
                        SparseMatrix*             KK,
                        NumericVector* RES,
                        const bool assembleMatrix,
                        //-----------
                        const double alpha,
                        const double beta,
                        const unsigned qrule_i
                       ) {

 
    for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
        
            geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);

            const short unsigned ielGeom = geom_element_iel.geom_type();

   el_dofs_unknowns(sol, msh, pdeSys, iel,
                        SolFEType_Mat,
                        SolIndex_Mat,
                        SolPdeIndex,
                        Sol_n_el_dofs_Mat, 
                        sol_eldofs_Mat,  
                        L2G_dofmap_Mat);
        
 //***************************************************
    unsigned int nDof_max          = ElementJacRes<double>::compute_max_n_dofs(Sol_n_el_dofs_Mat);
    
    unsigned int sum_Sol_n_el_dofs = ElementJacRes<double>::compute_sum_n_dofs(Sol_n_el_dofs_Mat);

    Res.resize(sum_Sol_n_el_dofs);                        std::fill(Res.begin(), Res.end(), 0.);

    Jac.resize(sum_Sol_n_el_dofs * sum_Sol_n_el_dofs);    std::fill(Jac.begin(), Jac.end(), 0.);
    
    L2G_dofmap_Mat_AllVars.resize(0);
      for (unsigned  k = 0; k < n_unknowns; k++)     L2G_dofmap_Mat_AllVars.insert(L2G_dofmap_Mat_AllVars.end(), L2G_dofmap_Mat[k].begin(), L2G_dofmap_Mat[k].end());
 //***************************************************

   geom_element_iel.set_elem_center(iel, solType_coords);


   //************ set control flag *********************
  int control_el_flag = 0;
        control_el_flag = ControlDomainFlag_bdry(geom_element_iel.get_elem_center());
  std::vector<int> control_node_flag(Sol_n_el_dofs_Mat[pos_mat_ctrl],0);
 //*************************************************** 
      
	// Perform face loop over elements that contain some control face
	if (control_el_flag == 1) {
	  
	  double tau = 0.;
	       
	  // loop on faces of the current element

	  for(unsigned iface=0; iface < msh->GetElementFaceNumber(iel); iface++) {
          
       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, iface);    
       
       std::vector<unsigned int> Sol_el_n_dofs_current_face(n_quantities); ///@todo the active flag is not an unknown!

       for (unsigned  k = 0; k < Sol_el_n_dofs_current_face.size(); k++) {
                 if (SolFEType_quantities[k] < 3) Sol_el_n_dofs_current_face[k] = msh->GetElementFaceDofNumber(iel, iface, SolFEType_quantities[k]);  ///@todo fix this absence
       }
       
       const unsigned nDof_max_bdry = ElementJacRes<double>::compute_max_n_dofs(Sol_el_n_dofs_current_face);
       
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();

	    // look for boundary faces
            const int bdry_index = msh->el->GetFaceElementIndex(iel, iface);
            
	    if( bdry_index < 0) {
	      const unsigned int face_in_rectangle_domain = -( msh->el->GetFaceElementIndex(iel,iface)+1);
		
// 	      if( !ml_sol->_SetBoundaryConditionFunction(xx,"U",tau,face,0.) && tau!=0.){
	      if(  face_in_rectangle_domain == FACE_FOR_CONTROL) { //control face
              
 //=================================================== 
		//we use the dirichlet flag to say: if dirichlet = true, we set 1 on the diagonal. if dirichlet = false, we put the boundary equation
	      bool  dir_bool = ml_sol->GetBdcFunction()(geom_element_iel.get_elem_center_bdry(), Solname_Mat[pos_mat_ctrl].c_str(), tau, face_in_rectangle_domain, 0.);

 //=================================================== 
        
 
//========= initialize gauss quantities on the boundary ============================================
                double sol_ctrl_iqp_bdry = 0.;
                std::vector<double> sol_ctrl_x_iqp_bdry(space_dim);   std::fill(sol_ctrl_x_iqp_bdry.begin(), sol_ctrl_x_iqp_bdry.end(), 0.);

//========= initialize gauss quantities on the boundary ============================================
		
        const unsigned n_qp_bdry = ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussPointsNumber();
        
    
		for(unsigned iqp_bdry = 0; iqp_bdry < n_qp_bdry; iqp_bdry++) {
    
    elem_all[qrule_i][ielGeom_bdry][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_bdry_3d(), iqp_bdry, Jac_qp_bdry, JacI_qp_bdry, detJac_qp_bdry, space_dim);
// 	elem_all[qrule_i][ielGeom_bdry][solType_coords]->compute_normal(Jac_qp_bdry, normal);
    
    weight_bdry = detJac_qp_bdry * ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussWeightsPointer()[iqp_bdry];

    elem_all[qrule_i][ielGeom_bdry][SolFEType_quantities[pos_sol_ctrl]] ->shape_funcs_current_elem(iqp_bdry, JacI_qp_bdry, phi_ctrl_iel_bdry_iqp_bdry, phi_ctrl_x_iel_bdry_iqp_bdry, boost::none, space_dim);
  
//========== compute gauss quantities on the boundary ===============================================
		  sol_ctrl_iqp_bdry = 0.;
                  std::fill(sol_ctrl_x_iqp_bdry.begin(), sol_ctrl_x_iqp_bdry.end(), 0.);
		      for (int i_bdry = 0; i_bdry < Sol_n_el_dofs_quantities[pos_sol_ctrl]; i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
			
			sol_ctrl_iqp_bdry +=  sol_eldofs_Mat[pos_mat_ctrl][i_vol] * phi_ctrl_iel_bdry_iqp_bdry[i_bdry];
                            for (int d = 0; d < space_dim; d++) {
			      sol_ctrl_x_iqp_bdry[d] += sol_eldofs_Mat[pos_mat_ctrl][i_vol] * phi_ctrl_x_iel_bdry_iqp_bdry[i_bdry * space_dim + d];
			    }
		      }
		      
	      
//========== compute gauss quantities on the boundary ================================================

		  // *** phi_i loop ***
		  for(unsigned i_bdry=0; i_bdry < nDof_max_bdry; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);

                 double lap_rhs_dctrl_ctrl_bdry_gss_i = 0.;
                 for (unsigned d = 0; d < space_dim; d++) {
                       if ( i_vol < Sol_n_el_dofs_Mat[pos_mat_ctrl] )  lap_rhs_dctrl_ctrl_bdry_gss_i +=  phi_ctrl_x_iel_bdry_iqp_bdry[i_bdry * space_dim + d] * sol_ctrl_x_iqp_bdry[d];
                 }
                 
//=============== construct control node flag field on the go  =========================================    
	      /* (control_node_flag[i])       picks nodes on \Gamma_c
	         (1 - control_node_flag[i])   picks nodes on \Omega \setminus \Gamma_c
	       */
	      if (dir_bool == false) { 
// 		std::cout << " found boundary control nodes ==== " << std::endl;
			for(unsigned k=0; k<control_node_flag.size(); k++) {
				  control_node_flag[i_vol] = 1;
			}
              }
//=============== construct control node flag field on the go  =========================================    

		 
//============ Bdry Residuals ==================	
                Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat,pos_mat_ctrl,i_vol) ]  +=  - control_node_flag[i_vol] *  weight_bdry *
                                                                                (     ( 1 - is_block_dctrl_ctrl_inside_bdry ) * alpha * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * sol_ctrl_iqp_bdry
							                           +  ( 1 - is_block_dctrl_ctrl_inside_bdry ) * beta * lap_rhs_dctrl_ctrl_bdry_gss_i 
							                         );  //boundary optimality condition
//============ Bdry Residuals ==================    
		    
		    for(unsigned j_bdry=0; j_bdry < nDof_max_bdry; j_bdry ++) {
		         unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, iface, j_bdry);

//============ Bdry Jacobians ==================	
//============ Bdry Jacobians ==================	



// SECOND BLOCK ROW
//=========== boundary control eqn =============	    

//========block delta_control / control ========
              double  lap_mat_dctrl_ctrl_bdry_gss = 0.;
		      for (unsigned d = 0; d < space_dim; d++) {  lap_mat_dctrl_ctrl_bdry_gss += phi_ctrl_x_iel_bdry_iqp_bdry[i_bdry * space_dim + d] * phi_ctrl_x_iel_bdry_iqp_bdry[j_bdry * space_dim + d];    }

          
              Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat, sum_Sol_n_el_dofs, pos_mat_ctrl, pos_mat_ctrl, i_vol, j_vol) ] 
			+=  control_node_flag[i_vol] *  weight_bdry * ( ( 1 - is_block_dctrl_ctrl_inside_bdry ) * alpha * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * phi_ctrl_iel_bdry_iqp_bdry[j_bdry] 
			                                              + ( 1 - is_block_dctrl_ctrl_inside_bdry ) * beta *  lap_mat_dctrl_ctrl_bdry_gss);   
    
		   
//============ End Bdry Jacobians ==================	
//============ End Bdry Jacobians ==================	
				
	      }  //end j loop
	      
		  }  //end i loop
		}  //end iqp_bdry loop
	  }    //end if control face
	      
	    }  //end if boundary faces
	  }    //end loop over faces
	  
	} //end if control element flag
	
	
       RES->add_vector_blocked(Res, L2G_dofmap_Mat_AllVars);

    if (assembleMatrix) {
      KK->add_matrix_blocked(Jac, L2G_dofmap_Mat_AllVars, L2G_dofmap_Mat_AllVars);
    }   
      
    }

    
}


  
  
  
  
  
} //end namespace
 
#endif
