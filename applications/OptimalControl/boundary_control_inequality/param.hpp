#ifndef ELLIPTIC_PARAMETERS
#define ELLIPTIC_PARAMETERS



#include "Mesh.hpp"
#include "CurrentElem.hpp"
#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"
#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"


#include "../fractional_functions.hpp"


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
	    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
  L2G_dofmap_mu_bdry[i_bdry]   = L2G_dofmap[pos_mu][i_vol];
  L2G_dofmap_ctrl_bdry[i_bdry] = L2G_dofmap[pos_ctrl][i_vol];
      }
      
 //============= delta_mu row ===============================
      std::vector<double> Res_mu_bdry (sol_actflag.size());     std::fill(Res_mu_bdry.begin(),Res_mu_bdry.end(), 0.);
      std::vector<double> Res_mu (Sol_n_el_dofs[pos_mu]);       std::fill(Res_mu.begin(),Res_mu.end(), 0.);
      
      for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
	    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
        
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
       unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
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
       unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
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

  
  
  
    //********** BEGIN FRAC CONTROL *****************************************

  void control_eqn_bdry_fractional(const unsigned iproc,
                                   const unsigned nprocs,
                        MultiLevelProblem &    ml_prob,
                        MultiLevelSolution*    ml_sol,
                        const Solution*        sol,
                        const Mesh * msh,
                        const  LinearEquationSolver* pdeSys,
                        //-----------
                        CurrentElem < double > & geom_element,
                        CurrentElem < double > & geom_element_jel,
                        const unsigned int solType_coords,
                        const unsigned int dim,
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
                        const unsigned int maxSize,
                        //-----------
                        std::vector < double > & Res,
                        std::vector < double > & Jac,
                        //-----------
                        const unsigned int n_quantities,
                        vector < unsigned > SolFEType_quantities,
                        vector < unsigned > Sol_n_el_dofs_quantities,
                        //-----------
                        std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > elem_all,
                        std::vector < std::vector < double > >  Jac_qp/*_bdry*/,
                        std::vector < std::vector < double > >  JacI_qp/*_bdry*/,
                        double detJac_qp/*_bdry*/,
                        double weight/*_bdry*/,
                        vector <double> phi_ctrl/*_bdry*/,
                        vector <double> phi_ctrl_x/*_bdry*/, 
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
                        const unsigned int OP_Hhalf
                       ) {
      
//  //***************************************************
//   //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
// //***************************************************
//   std::vector < std::vector < double > >  JacI_qp(space_dim);
//   std::vector < std::vector < double > >  Jac_qp(dim);
//   for(unsigned d = 0; d < dim; d++) {
//     Jac_qp[d].resize(space_dim);
//   }
//   for(unsigned d = 0; d < space_dim; d++) {
//     JacI_qp[d].resize(dim);
//   }
// 
//   double detJac_qp;
//   std::vector < std::vector < /*const*/ elem_type_templ_base< double, double > *  > > elem_all;
//   ml_prob.get_all_abstract_fe(elem_all);
// //***************************************************
     
   vector < double > phi_x;
     
      
  unsigned solType =   SolFEType_Mat[pos_mat_ctrl];
  unsigned soluIndex = SolIndex_Mat[pos_mat_ctrl];
  unsigned soluPdeIndex = SolPdeIndex[pos_mat_ctrl];
  
  std::vector < double > solu1;
  std::vector < double > solu2;
    
  //------- geometry ---------------
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector < vector < double > > x1(dim);
  vector < vector < double > > x2(dim);
  for(unsigned k = 0; k < dim; k++) {
    x1[k].reserve(maxSize);
    x2[k].reserve(maxSize);
  }
 
  //------- geometry ---------------
 
 
 
  //-------- local to global mappings --------------
  vector< int > l2GMap1;  l2GMap1.reserve(maxSize);
  vector< int > l2GMap2;  l2GMap2.reserve(maxSize);


  //-------- Local matrices and rhs --------------
  vector < double > KK_local;  KK_local.reserve(maxSize * maxSize);
  vector < double > Res_local; Res_local.reserve(maxSize);

//   Local matrices and rhs for adaptive quadrature
  vector < double > Res_local_refined; Res_local_refined.reserve(maxSize);
  vector < double > CClocal_refined;   CClocal_refined.reserve(maxSize * maxSize);

  vector < double > KK_mixed;   KK_mixed.reserve(maxSize * maxSize);
  vector < double > Res_mixed;  Res_mixed.reserve(maxSize);

//   Non local matrices and vectors for H^s laplacian operator
//   vector< double >         Res_nonlocal;
//   Res_nonlocal.reserve(maxSize);  // local residual vector
  vector< double >         Res_nonlocalI;  Res_nonlocalI.reserve(maxSize);
  vector< double >         Res_nonlocalJ;  Res_nonlocalJ.reserve(maxSize);
//   vector < double > CClocal;
//   CClocal.reserve(maxSize * maxSize);
  vector < double > CC_nonlocal_II;  CC_nonlocal_II.reserve(maxSize * maxSize);
  vector < double > CC_nonlocal_IJ;  CC_nonlocal_IJ.reserve(maxSize * maxSize);
  vector < double > CC_nonlocal_JI;  CC_nonlocal_JI.reserve(maxSize * maxSize);
  vector < double > CC_nonlocal_JJ;  CC_nonlocal_JJ.reserve(maxSize * maxSize); 
 
 //----------------------
  KK->zero();
  RES->zero(); 
 //----------------------
  
  
   for(int kproc = 0; kproc < nprocs; kproc++) {
    for(int jel = msh->_elementOffset[kproc]; jel < msh->_elementOffset[kproc + 1]; jel++) {
 
       short unsigned ielGeom2;
      unsigned nDof2;
      unsigned nDofx2;
      unsigned nDofu2;
      //unsigned n_face;

// --- vector sizes
      if(kproc == iproc) {
        nDof2  = msh->GetElementDofNumber(jel, solType);    // number of solution element dofs
        nDofx2 = msh->GetElementDofNumber(jel, solType_coords);    // number of coordinate element dofs

      }

      MPI_Bcast(&nDof2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
      MPI_Bcast(&nDofx2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
      //MPI_Bcast(&n_face, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
// --- vector sizes

      
// --- l2GMap
        l2GMap2.resize(nDof2);

      // local storage of global mapping and solution ********************
      if(kproc == iproc) {
        for(unsigned j = 0; j < nDof2; j++) {
          l2GMap2[j] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, j, jel);  // global to global mapping between solution node and pdeSys dof
        }
      }
      MPI_Bcast(&l2GMap2[0], nDof2, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
      // ******************************************************************
// --- l2GMap

      
// --- geometry and solution
      
      
      
      
      for(int k = 0; k < dim; k++) {
        x2[k].resize(nDofx2);
      }

      solu2.resize(nDof2);
      
      if(kproc == iproc) {
        for(unsigned j = 0; j < nDofx2; j++) {
          unsigned xDof  = msh->GetSolutionDof(j, jel, solType_coords);  // global to global mapping between coordinates node and coordinate dof
          for(unsigned k = 0; k < dim; k++) {
            x2[k][j] = (*msh->_topology->_Sol[k])(xDof);  // global extraction and local storage for the element coordinates
          }
        }
        for(unsigned j = 0; j < nDof2; j++) {
          unsigned jDof  = msh->GetSolutionDof(j, jel, solType);  // global to global mapping between coordinates node and coordinate dof
          solu2[j] = (*sol->_Sol[soluIndex])(jDof);  // global extraction and local storage for the element coordinates
        }
      }
      for(unsigned k = 0; k < dim; k++) {
        MPI_Bcast(& x2[k][0], nDofx2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }
      MPI_Bcast(& solu2[0], nDof2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);

      // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      if(kproc == iproc) {
        geom_element_jel.set_coords_at_dofs_and_geom_type(jel, solType_coords);
      }
      for(unsigned k = 0; k < dim; k++) {
        MPI_Bcast(& geom_element_jel.get_coords_at_dofs()[k][0], nDofx2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }
      for(unsigned k = 0; k < space_dim; k++) {
        MPI_Bcast(& geom_element_jel.get_coords_at_dofs_3d()[k][0], nDofx2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }
      // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// --- geometry and solution


// --- geom_el type
      if(kproc == iproc) {   ielGeom2 = msh->GetElementType(jel);  }
      MPI_Bcast(&ielGeom2, 1, MPI_UNSIGNED_SHORT, kproc, MPI_COMM_WORLD);
// --- geom_el type

//       const unsigned jgNumber = msh->_finiteElement[ielGeom2][solType]->GetGaussPointNumber();
      const unsigned jgNumber = ml_prob.GetQuadratureRule(ielGeom2).GetGaussPointsNumber();

      vector < vector < double > > xg2(jgNumber);
      vector <double> weight2(jgNumber);
      vector < vector <double> > phi2(jgNumber);  // local test function
      std::vector< double > solY(jgNumber, 0.);

      
      for(unsigned jg = 0; jg < jgNumber; jg++) {

//         msh->_finiteElement[ielGeom2][solType]->Jacobian(x2, jg, weight2[jg], phi2[jg], phi_x);

        elem_all[ielGeom2][solType_coords]->JacJacInv(/*x2*/geom_element_jel.get_coords_at_dofs_3d(), jg, Jac_qp, JacI_qp, detJac_qp, space_dim);
        weight2[jg] = detJac_qp * ml_prob.GetQuadratureRule(ielGeom2).GetGaussWeightsPointer()[jg];
        elem_all[ielGeom2][solType]->shape_funcs_current_elem(jg, JacI_qp, phi2[jg], phi_x /*boost::none*/, boost::none /*phi_u_xx*/, space_dim);



//--- geometry and solution, at qp
        xg2[jg].assign(dim, 0.);
        solY[jg] = 0.;

        for(unsigned j = 0; j < nDof2; j++) {
          solY[jg] += solu2[j] * phi2[jg][j];
          for(unsigned k = 0; k < dim; k++) {
            xg2[jg][k] += x2[k][j] * phi2[jg][j];
          }
        }
      }       
//--- geometry and solution, at qp
        
     
       
       unsigned counter_verify = 0;
       
       for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
                
           geom_element.set_coords_at_dofs_and_geom_type(iel, solType_coords);
              
           geom_element.set_elem_center(iel, solType_coords);


   //************ set control flag *********************
  int control_el_flag = 0;
        control_el_flag = ControlDomainFlag_bdry(geom_element.get_elem_center());
  std::vector<int> control_node_flag(Sol_n_el_dofs_Mat[pos_mat_ctrl], 0);
 //*************************************************** 
      
	// Perform face loop over elements that contain some control face
	if (control_el_flag == 1) {
        
        
        	  double tau=0.;
	  std::vector<double> normal(space_dim, 0.);
	       
	  // loop on faces of the current element

	  for(unsigned jface=0; jface < msh->GetElementFaceNumber(iel); jface++) {
          
       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, jface);    
       
       std::vector<unsigned int> Sol_el_n_dofs_current_face(n_quantities); ///@todo the active flag is not an unknown!

       for (unsigned  k = 0; k < Sol_el_n_dofs_current_face.size(); k++) {
                 if (SolFEType_quantities[k] < 3) Sol_el_n_dofs_current_face[k] = msh->GetElementFaceDofNumber(iel, jface, SolFEType_quantities[k]);  ///@todo fix this absence
       }
       
       const unsigned nDof_max_bdry = ElementJacRes<double>::compute_max_n_dofs(Sol_el_n_dofs_current_face);
       
       geom_element.set_coords_at_dofs_bdry_3d(iel, jface, solType_coords);
 
       geom_element.set_elem_center_bdry_3d();

	    // look for boundary faces
            const int bdry_index = msh->el->GetFaceElementIndex(iel, jface);
            
	    if( bdry_index < 0) {
	      const unsigned int face_in_rectangle_domain = -( msh->el->GetFaceElementIndex(iel,jface)+1);
		
// 	      if( !ml_sol->_SetBoundaryConditionFunction(xx,"U",tau,face,0.) && tau!=0.){
	      if(  face_in_rectangle_domain == FACE_FOR_CONTROL) { //control face

        counter_verify++;
              
              //Quadrature loop
                      const unsigned n_gauss_bdry = ml_prob.GetQuadratureRule(ielGeom_bdry).GetGaussPointsNumber();
      
     //**** Evaluating coarse FE functions on Quadrature Points of the "sub-elements"
                      
                      
     //**** Evaluating coarse FE functions on Quadrature Points of the "sub-elements"
                      
    
		for(unsigned ig_bdry = 0; ig_bdry < n_gauss_bdry; ig_bdry++) {
            
            
      //============ Adaptive quadrature for iel == jel ==================
          if(iel == jel) {
                                          
            if(Nsplit != 0) {
                
                
                
                
                
                
            }  //end if Nsplit != 0
                                          
        } //iel == jel                
    //============ Adaptive quadrature for iel == jel ==================
             
                
                
            
            
            
            
        }
              
              
              
          }
        }
      }
        
	  
    } //end control elem flag
           
           
        short unsigned ielGeom1 = msh->GetElementType(iel);
        unsigned nDof1  = msh->GetElementDofNumber(iel, solType);
        unsigned nDofx1 = msh->GetElementDofNumber(iel, solType_coords);

// --- l2GMap
        l2GMap1.resize(nDof1);
        for(unsigned i = 0; i < nDof1; i++) {
          l2GMap1[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);
        }
// --- l2GMap


// --- geometry        
        for(int k = 0; k < dim; k++) {
          x1[k].resize(nDofx1);
        }
        
        // local storage of coordinates
        for(unsigned i = 0; i < nDofx1; i++) {
          unsigned xDof  = msh->GetSolutionDof(i, iel, xType);
          for(unsigned k = 0; k < dim; k++) {
            x1[k][i] = (*msh->_topology->_Sol[k])(xDof);
          }
        }
// --- geometry        
        
        
// --- solution        
       solu1.resize(nDof1);
        
        for(unsigned i = 0; i < nDof1; i++) {
          unsigned iDof  = msh->GetSolutionDof(i, iel, solType);  // global to global mapping between coordinates node and coordinate dof
          solu1[i] = (*sol->_Sol[soluIndex])(iDof);  // global extraction and local storage for the element coordinates
        }
// --- solution        
                  
  //****** matrix resizing ******
//           CC_local.assign(nDof1 * nDof2, 0.);   //resize
        CC_nonlocal_II.assign(nDof1 * nDof2, 0.);   //resize
        CC_nonlocal_IJ.assign(nDof1 * nDof2, 0.);   //resize
        CC_nonlocal_JI.assign(nDof1 * nDof2, 0.);   //resize
        CC_nonlocal_JJ.assign(nDof1 * nDof2, 0.);   //resize
//         Res_nonlocal.assign(nDof1, 0);    //resize
        Res_nonlocalI.assign(nDof1, 0);    //resize
        Res_nonlocalJ.assign(nDof1, 0);    //resize

        Res_mixed.assign(nDof1, 0);    //resize
        KK_mixed.assign(nDof1 * nDof1, 0.);

        if(iel == jel) {
          Res_local.assign(nDof1, 0);    //resize
          KK_local.assign(nDof1 * nDof1, 0.);
          if(Nsplit != 0) {
//             Vectors and matrices for adaptive quadrature
            Res_local_refined.assign(nDof1, 0);    //resize
            CClocal_refined.assign(nDof1 * nDof1, 0.);
          }
        }
  //****** matrix resizing ******

  
  const unsigned igNumber = msh->_finiteElement[ielGeom1][solType]->GetGaussPointNumber();

         double weight1;
        vector < double > phi1;  // local test function

        double weight3;
        vector < double > phi3;  // local test function

        double solX = 0.;
        std::vector<double> sol_u_x(space_dim);
        std::fill(sol_u_x.begin(), sol_u_x.end(), 0.);


        std::vector < std::vector < std::vector <double > > > aP(3);
        if(Nsplit > 0) {
          for(unsigned j_fe_type = 0; j_fe_type < solType + 1; j_fe_type++) {
            ProjectNodalToPolynomialCoefficients(aP[j_fe_type], x1, ielGeom1, j_fe_type) ;
          }
        }                   
                    
                    
        for(unsigned ig = 0; ig < igNumber; ig++) {
  
            
          msh->_finiteElement[ielGeom1][solType]->Jacobian(x1, ig, weight1, phi1, phi_x);

          // evaluate the solution, the solution derivatives and the coordinates in the gauss point
          vector < double > xg1(dim, 0.);
          solX = 0.;

          for(unsigned i = 0; i < nDof1; i++) {
            solX += solu1[i] * phi1[i];
            for(unsigned d = 0; d < sol_u_x.size(); d++)   sol_u_x[d] += solu1[i] * phi_x[i * dim + d];
            for(unsigned k = 0; k < dim; k++) {
              xg1[k] += x1[k][i] * phi1[i];
            }
          }
          
//============ Adaptive quadrature for iel == jel ==================
          if(iel == jel) {
                                          
            if(Nsplit != 0) {
                                         
            std::cout.precision(14);
              std::vector< std::vector<std::vector<double>>> x3;

              for(unsigned split = 0; split <= Nsplit; split++) {

//                 unsigned size_part;
//                 if(dim == 1) size_part = 2;
//                 else size_part = (split != Nsplit) ? 12 : 4;

                if(dim == 1) GetElementPartition1D(xg1, x1, split, Nsplit, x3);
                else if(dim == 2) {
                  //GetElementPartition2D(xg1, x1, split, Nsplit, x3);
                  GetElementPartitionQuad(xg1, x1, split, Nsplit, x3);
                }

                //for(unsigned r = 0; r < size_part; r++) {
                for(unsigned r = 0; r < x3.size(); r++) {


                  for(unsigned jg = 0; jg < igNumber; jg++) {


                    msh->_finiteElement[ielGeom1][solType]->Jacobian(x3[r], jg, weight3, phi3, phi_x);

                    vector < double > xg3(dim, 0.);

                    for(unsigned i = 0; i < nDof1; i++) {
                      for(unsigned k = 0; k < dim; k++) {
                        xg3[k] += x3[r][k][i] * phi3[i];
                      }
                    }

                    std::vector<double> xi3(dim, 0.);

                    GetClosestPointInReferenceElement(x1, xg3, ielGeom1, xi3);
                    GetInverseMapping(solType, ielGeom1, aP, xg3, xi3, 1000);

                    msh->_finiteElement[ielGeom1][solType]->GetPhi(phi3, xi3);

                    double solY3 = 0.;
                    for(unsigned i = 0; i < nDof1; i++) {
                      solY3 += solu1[i] * phi3[i];
                    }

                    double dist_xyz3 = 0;
                    for(unsigned k = 0; k < dim; k++) {
                      dist_xyz3 += (xg1[k] - xg3[k]) * (xg1[k] - xg3[k]);
                    }

                    const double denom3 = pow(dist_xyz3, (double)((dim / 2.) + s_frac));

                    for(unsigned i = 0; i < nDof1; i++) {

                      Res_local_refined[ i ]    +=      - (C_ns / 2.) * OP_Hhalf * check_limits *
                                                        ((solX - solY3) * (phi1[i] - phi3[i]) * weight3 / denom3
                                                        ) * weight1 ;

                      for(unsigned j = 0; j < nDof2; j++) {
                        CClocal_refined[ i * nDof2 + j ] += (C_ns / 2.) * OP_Hhalf * check_limits *
                                                            ((phi1[j] - phi3[j]) * (phi1[i] - phi3[i]) * weight3 / denom3
                                                            ) * weight1 ;

                      }
                    }
                                          
                  } //end jg
                } //end r
              }  //end split
              
            }  //end if Nsplit != 0

                                          
        } //iel == jel
  //============ Adaptive quadrature for iel == jel ==================

                            
                            
                            
                            
  } //end ig
                
                
//============ add to global ==================
         if(iel == jel) {
          KK->add_matrix_blocked(KK_local, l2GMap1, l2GMap1);
          RES->add_vector_blocked(Res_local, l2GMap1);

          if(Nsplit != 0) {
            KK->add_matrix_blocked(CClocal_refined, l2GMap1, l2GMap1);
            RES->add_vector_blocked(Res_local_refined, l2GMap1);
          }
        }
//        KK->add_matrix_blocked(CClocal, l2GMap1, l2GMap2);
        KK->add_matrix_blocked(KK_mixed, l2GMap1, l2GMap1);
        RES->add_vector_blocked(Res_mixed, l2GMap1);

        KK->add_matrix_blocked(CC_nonlocal_II, l2GMap1, l2GMap1);
//         KK->add_matrix_blocked(CC_nonlocal_JI, l2GMap2, l2GMap1);  ///@todo
//         KK->add_matrix_blocked(CC_nonlocal_IJ, l2GMap1, l2GMap2);  ///@todo
        KK->add_matrix_blocked(CC_nonlocal_JJ, l2GMap2, l2GMap2);

//        RES->add_vector_blocked(Res_nonlocal, l2GMap1);
        RES->add_vector_blocked(Res_nonlocalI, l2GMap1);
        RES->add_vector_blocked(Res_nonlocalJ, l2GMap2);
//============ add to global - end ==================
        
        
    } //end iel
    
    } //end jel
    } //end kproc
    
    
// AAA do not close this because later they will be filled with the rest of the system!!!      
//    KK->close();
//   RES->close();

   //print JAC and RES to files
KK->close();KK->zero();
const unsigned nonlin_iter = 0/*mlPdeSys->GetNonlinearIt()*/;
    assemble_jacobian< double, double >::print_global_jacobian(assembleMatrix, ml_prob, KK, nonlin_iter);
//     assemble_jacobian< double, double >::print_global_residual(ml_prob, RES, nonlin_iter);
abort();
//   std::cout << "****************************" << std::endl;
  
  }
    
  //********** END FRAC CONTROL *****************************************
  
  
    
  void control_eqn_bdry(const unsigned iproc,
                        MultiLevelProblem &    ml_prob,
                        MultiLevelSolution*    ml_sol,
                        const Solution*        sol,
                        const Mesh * msh,
                        const  LinearEquationSolver* pdeSys,
                        //-----------
                        CurrentElem < double > & geom_element,
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
                        std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > elem_all,
                        std::vector < std::vector < double > >  Jac_qp_bdry,
                        std::vector < std::vector < double > >  JacI_qp_bdry,
                        double detJac_qp_bdry,
                        double weight_bdry,
                        vector <double> phi_ctrl_bdry,
                        vector <double> phi_ctrl_x_bdry, 
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
                        const double beta     
                       ) {

 
    for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
        
            geom_element.set_coords_at_dofs_and_geom_type(iel, solType_coords);

            const short unsigned ielGeom = geom_element.geom_type();

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

   geom_element.set_elem_center(iel, solType_coords);


   //************ set control flag *********************
  int control_el_flag = 0;
        control_el_flag = ControlDomainFlag_bdry(geom_element.get_elem_center());
  std::vector<int> control_node_flag(Sol_n_el_dofs_Mat[pos_mat_ctrl],0);
 //*************************************************** 
      
	// Perform face loop over elements that contain some control face
	if (control_el_flag == 1) {
	  
	  double tau=0.;
	  std::vector<double> normal(space_dim, 0.);
	       
	  // loop on faces of the current element

	  for(unsigned jface=0; jface < msh->GetElementFaceNumber(iel); jface++) {
          
       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, jface);    
       
       std::vector<unsigned int> Sol_el_n_dofs_current_face(n_quantities); ///@todo the active flag is not an unknown!

       for (unsigned  k = 0; k < Sol_el_n_dofs_current_face.size(); k++) {
                 if (SolFEType_quantities[k] < 3) Sol_el_n_dofs_current_face[k] = msh->GetElementFaceDofNumber(iel, jface, SolFEType_quantities[k]);  ///@todo fix this absence
       }
       
       const unsigned nDof_max_bdry = ElementJacRes<double>::compute_max_n_dofs(Sol_el_n_dofs_current_face);
       
       geom_element.set_coords_at_dofs_bdry_3d(iel, jface, solType_coords);
 
       geom_element.set_elem_center_bdry_3d();

	    // look for boundary faces
            const int bdry_index = msh->el->GetFaceElementIndex(iel, jface);
            
	    if( bdry_index < 0) {
	      const unsigned int face_in_rectangle_domain = -( msh->el->GetFaceElementIndex(iel,jface)+1);
		
// 	      if( !ml_sol->_SetBoundaryConditionFunction(xx,"U",tau,face,0.) && tau!=0.){
	      if(  face_in_rectangle_domain == FACE_FOR_CONTROL) { //control face
              
 //=================================================== 
		//we use the dirichlet flag to say: if dirichlet = true, we set 1 on the diagonal. if dirichlet = false, we put the boundary equation
	      bool  dir_bool = ml_sol->GetBdcFunction()(geom_element.get_elem_center_bdry(), Solname_Mat[pos_mat_ctrl].c_str(), tau, face_in_rectangle_domain, 0.);

 //=================================================== 
        
 
//========= initialize gauss quantities on the boundary ============================================
                double sol_ctrl_bdry_gss = 0.;
                std::vector<double> sol_ctrl_x_bdry_gss(space_dim);   std::fill(sol_ctrl_x_bdry_gss.begin(), sol_ctrl_x_bdry_gss.end(), 0.);

//========= initialize gauss quantities on the boundary ============================================
		
        const unsigned n_gauss_bdry = ml_prob.GetQuadratureRule(ielGeom_bdry).GetGaussPointsNumber();
        
    
		for(unsigned ig_bdry = 0; ig_bdry < n_gauss_bdry; ig_bdry++) {
    
    elem_all[ielGeom_bdry][solType_coords]->JacJacInv(geom_element.get_coords_at_dofs_bdry_3d(), ig_bdry, Jac_qp_bdry, JacI_qp_bdry, detJac_qp_bdry, space_dim);
	elem_all[ielGeom_bdry][solType_coords]->compute_normal(Jac_qp_bdry, normal);
    
    weight_bdry = detJac_qp_bdry * ml_prob.GetQuadratureRule(ielGeom_bdry).GetGaussWeightsPointer()[ig_bdry];

    elem_all[ielGeom_bdry][SolFEType_quantities[pos_sol_ctrl]] ->shape_funcs_current_elem(ig_bdry, JacI_qp_bdry, phi_ctrl_bdry, phi_ctrl_x_bdry, boost::none, space_dim);
  
//========== compute gauss quantities on the boundary ===============================================
		  sol_ctrl_bdry_gss = 0.;
                  std::fill(sol_ctrl_x_bdry_gss.begin(), sol_ctrl_x_bdry_gss.end(), 0.);
		      for (int i_bdry = 0; i_bdry < Sol_n_el_dofs_quantities[pos_sol_ctrl]; i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
			
			sol_ctrl_bdry_gss +=  sol_eldofs_Mat[pos_mat_ctrl][i_vol] * phi_ctrl_bdry[i_bdry];
                            for (int d = 0; d < space_dim; d++) {
			      sol_ctrl_x_bdry_gss[d] += sol_eldofs_Mat[pos_mat_ctrl][i_vol] * phi_ctrl_x_bdry[i_bdry * space_dim + d];
			    }
		      }
		      
	      
//========== compute gauss quantities on the boundary ================================================

		  // *** phi_i loop ***
		  for(unsigned i_bdry=0; i_bdry < nDof_max_bdry; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);

                 double lap_rhs_dctrl_ctrl_bdry_gss_i = 0.;
                 for (unsigned d = 0; d < space_dim; d++) {
                       if ( i_vol < Sol_n_el_dofs_Mat[pos_mat_ctrl] )  lap_rhs_dctrl_ctrl_bdry_gss_i +=  phi_ctrl_x_bdry[i_bdry * space_dim + d] * sol_ctrl_x_bdry_gss[d];
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
                                                                                (     ( 1 - is_block_dctrl_ctrl_inside_bdry ) * alpha * phi_ctrl_bdry[i_bdry] * sol_ctrl_bdry_gss
							                           +  ( 1 - is_block_dctrl_ctrl_inside_bdry ) * beta * lap_rhs_dctrl_ctrl_bdry_gss_i 
							                         );  //boundary optimality condition
//============ Bdry Residuals ==================    
		    
		    for(unsigned j_bdry=0; j_bdry < nDof_max_bdry; j_bdry ++) {
		         unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, jface, j_bdry);

//============ Bdry Jacobians ==================	
//============ Bdry Jacobians ==================	



// SECOND BLOCK ROW
//=========== boundary control eqn =============	    

//========block delta_control / control ========
              double  lap_mat_dctrl_ctrl_bdry_gss = 0.;
		      for (unsigned d = 0; d < space_dim; d++) {  lap_mat_dctrl_ctrl_bdry_gss += phi_ctrl_x_bdry[i_bdry * space_dim + d] * phi_ctrl_x_bdry[j_bdry * space_dim + d];    }

          
              Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat, sum_Sol_n_el_dofs, pos_mat_ctrl, pos_mat_ctrl, i_vol, j_vol) ] 
			+=  control_node_flag[i_vol] *  weight_bdry * ( ( 1 - is_block_dctrl_ctrl_inside_bdry ) * alpha * phi_ctrl_bdry[i_bdry] * phi_ctrl_bdry[j_bdry] 
			                                              + ( 1 - is_block_dctrl_ctrl_inside_bdry ) * beta *  lap_mat_dctrl_ctrl_bdry_gss);   
    
		   
//============ End Bdry Jacobians ==================	
//============ End Bdry Jacobians ==================	
				
	      }  //end j loop
	      
		  }  //end i loop
		}  //end ig_bdry loop
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
