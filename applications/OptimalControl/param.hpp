#ifndef ELLIPTIC_PARAMETERS
#define ELLIPTIC_PARAMETERS



#include "Mesh.hpp"
#include "CurrentElem.hpp"
#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"
#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"
#include "MED_IO.hpp"

#include "fractional_functions.hpp"


#include <functional>

#include <boost/mpi.hpp>



//*********************** Sets Number of refinements *****************************************
#define N_UNIFORM_LEVELS 4
#define N_ERASED_LEVELS   N_UNIFORM_LEVELS - 1


//*********************** Sets Number of subdivisions in X and Y direction *****************************************

#define NSUB_X  2
#define NSUB_Y  2
#define NSUB_Z  2


//*********************** Sets the regularization parameters *******************************************************
// for pure boundary approaches
#define ALPHA_CTRL_BDRY 0.001 
#define BETA_CTRL_BDRY   ALPHA_CTRL_BDRY

// for lifting approaches (both internal and external)
#define ALPHA_CTRL_VOL 0.005 
#define BETA_CTRL_VOL ALPHA_CTRL_VOL


//*********************** Control boundary extremes *******************************************************

#define GAMMA_CONTROL_LOWER 0.25
#define GAMMA_CONTROL_UPPER 0.75


//*********************** Lifting internal extension *******************************************************
#define LIFTING_INTERNAL_DEPTH  0.25
#define LIFTING_INTERNAL_WIDTH_LOWER  GAMMA_CONTROL_LOWER
#define LIFTING_INTERNAL_WIDTH_UPPER  GAMMA_CONTROL_UPPER




//*********************** Control box constraints *******************************************************
#define  INEQ_FLAG 1.
#define  C_COMPL 1.


namespace femus {


 double InequalityConstraint(const std::vector<double> & dof_obj_coord, const bool upper) {

     double constr_value = 0.;
     double constr_value_upper =  1000.;// dof_obj_coord[1]*(1. - dof_obj_coord[1]);
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

    const double target_line_position_along_coordinate = 0.5;
 //***** set target domain flag ******
  int target_flag = 0; //set 0 to 1 to get the entire domain
  
  const double offset_to_include_line = 1.e-5;
   
  const int  target_line_sign = target_line_sign_func(/*FACE_FOR_CONTROL*/FACE_FOR_TARGET);
  
  const unsigned int axis_dir = axis_direction_target_reg(/*FACE_FOR_CONTROL*/FACE_FOR_TARGET);
   
   const double target_line = target_line_position_along_coordinate + target_line_sign * offset_to_include_line; 
   
   
   
      if ((  target_line_sign * elem_center[axis_dir] < target_line_sign * target_line ) && 
          (  target_line_sign * elem_center[axis_dir] > - target_line_position_along_coordinate + target_line_sign * (target_line_position_along_coordinate - target_line_sign * offset_to_include_line)))
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
  
  const double control_domain_depth = LIFTING_INTERNAL_DEPTH;
  
  const double control_domain_width_lower = LIFTING_INTERNAL_WIDTH_LOWER;
  const double control_domain_width_upper = LIFTING_INTERNAL_WIDTH_UPPER;
   
   const int  target_line_sign = target_line_sign_func(FACE_FOR_CONTROL);

   const double extreme_pos = extreme_position(FACE_FOR_CONTROL);

   const unsigned int axis_dir = axis_direction_Gamma_control(FACE_FOR_CONTROL);

   
   if ( ( target_line_sign * elem_center[1 - axis_dir] <   target_line_sign * ( extreme_pos + target_line_sign * control_domain_depth ) )
       && ( elem_center[axis_dir] > control_domain_width_lower - offset_to_include_line ) 
       && ( elem_center[axis_dir] < control_domain_width_upper + offset_to_include_line ) )
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

// // //         std::cout << " mu dofs " << std::endl;
// // //                 for (unsigned i = 0; i < sol_actflag.size(); i++) {
// // //                 std::cout << sol_eldofs[pos_mu][i] << " ";
// // //         }
// // //         
// // //         std::cout << std::endl;
        
// // //         std::cout << " ctrl dofs " << std::endl;
// // //         for (unsigned i = 0; i < sol_actflag.size(); i++) {
// // //                 std::cout << sol_eldofs[pos_ctrl][i] << " ";
// // //         }
// // //         std::cout << std::endl;
        
        for (unsigned i = 0; i < sol_actflag.size(); i++) {
            std::vector<double> node_coords_i(dim, 0.);
            for (unsigned d = 0; d < dim; d++) node_coords_i[d] = coords_at_dofs[d][i];
            
            ctrl_lower[i] = InequalityConstraint(node_coords_i, false);
            ctrl_upper[i] = InequalityConstraint(node_coords_i, true);
            
            const double lower_test_value = sol_eldofs[pos_mu][i] + c_compl * ( sol_eldofs[pos_ctrl][i] - ctrl_lower[i] );
            const double upper_test_value = sol_eldofs[pos_mu][i] + c_compl * ( sol_eldofs[pos_ctrl][i] - ctrl_upper[i] );

            if      ( lower_test_value < 0. )  {
                std::cout << "Found active node below" << std::endl;
                std::cout << "The current value of mu is " <<  sol_eldofs[pos_mu][i] << std::endl;
                   sol_actflag[i] = 1;
            }
            else if ( upper_test_value > 0. )  {
                std::cout << "Found active node above" << std::endl;
                std::cout << "The current value of mu is " <<  sol_eldofs[pos_mu][i] << std::endl;
                sol_actflag[i] = 2;
            }
        }

//************** local to global act flag ***************************

        for (unsigned i = 0; i < sol_actflag.size(); i++) {
            unsigned solDof_mu = msh->GetSolutionDof(i, iel, solFEType_act_flag);
            (sol->_Sol[solIndex_act_flag])->set(solDof_mu, sol_actflag[i]);
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
     
		const unsigned nve_bdry = msh->GetElementFaceDofNumber(iel, iface, solFEType_act_flag);
        
        const unsigned dim = coords_at_dofs.size();
        
         // 0: inactive; 1: active_a; 2: active_b
        assert(Sol_n_el_dofs[pos_mu] == Sol_n_el_dofs[pos_ctrl]);///@todo More appropriately, 
            sol_actflag.resize(nve_bdry/*nDof_mu*/);
            ctrl_lower.resize(nve_bdry/*nDof_mu*/);
            ctrl_upper.resize(nve_bdry/*nDof_mu*/);
           std::fill(sol_actflag.begin(), sol_actflag.end(), 0);
           std::fill(ctrl_lower.begin(), ctrl_lower.end(), 0.);
           std::fill(ctrl_upper.begin(), ctrl_upper.end(), 0.);

      for (unsigned int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
        std::vector<double> node_coords_i(dim, 0.);
        for (unsigned d = 0; d < dim; d++) node_coords_i[d] = coords_at_dofs[d][i_bdry];
        
        ctrl_lower[i_bdry] = InequalityConstraint(node_coords_i, false);
        ctrl_upper[i_bdry] = InequalityConstraint(node_coords_i, true);

        const double lower_test_value = sol_eldofs[pos_mu][i_vol] + c_compl * ( sol_eldofs[pos_ctrl][i_vol] - ctrl_lower[i_bdry] );
        const double upper_test_value = sol_eldofs[pos_mu][i_vol] + c_compl * ( sol_eldofs[pos_ctrl][i_vol] - ctrl_upper[i_bdry] );
        
        if      ( lower_test_value < 0. )  sol_actflag[i_bdry] = 1;
        else if ( upper_test_value > 0. )  sol_actflag[i_bdry] = 2;
            }
            
        //************** act flag **************************** 
      for (unsigned int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
	    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
      unsigned solDof_actflag = msh->GetSolutionDof(i_vol, iel, solFEType_act_flag); 
      (sol->_Sol[solIndex_act_flag])->set(solDof_actflag, sol_actflag[i_bdry]);     
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
//       std::vector<double> Res_mu (Sol_n_el_dofs[pos_mu]);       std::fill(Res_mu.begin(),Res_mu.end(), 0.);
      
      for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
	    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
        
      if (sol_actflag[i_bdry] == 0) {  //inactive
//          Res_mu [i_vol]      = - ineq_flag * ( 1. * sol_eldofs[pos_mu][i_vol] - 0. ); 
         Res_mu_bdry[i_bdry] = - ineq_flag * ( 1. * sol_eldofs[pos_mu][i_vol] - 0. ); 
      }
      else if (sol_actflag[i_bdry] == 1) {  //active_a 
// 	 Res_mu [i_vol]      = - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl][i_vol] - c_compl * ctrl_lower[i_bdry]);
     Res_mu_bdry[i_bdry] = - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl][i_vol] - c_compl * ctrl_lower[i_bdry]);
          
    }
      else if (sol_actflag[i_bdry] == 2) {  //active_b 
// 	Res_mu [i_vol]      =  - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl][i_vol] - c_compl * ctrl_upper[i_bdry]);
    Res_mu_bdry[i_bdry] =  - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl][i_vol] - c_compl * ctrl_upper[i_bdry]);
      }
    }

    
//     RES->insert(Res_mu,  L2G_dofmap[pos_mu]);    
    RES->insert(Res_mu_bdry,  L2G_dofmap_mu_bdry);    
 //============= delta_mu row - end ===============================
    
 //============= delta_mu-delta_ctrl row ===============================
 //auxiliary volume vector for act flag
//  unsigned nDof_actflag_vol  = msh->GetElementDofNumber(iel, solFEType_act_flag);
//  std::vector<double> sol_actflag_vol(nDof_actflag_vol); 


 for (unsigned i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++) if (sol_actflag[i_bdry] != 0 ) sol_actflag[i_bdry] = ineq_flag * c_compl;    
 
//  std::fill(sol_actflag_vol.begin(), sol_actflag_vol.end(), 0.);
//     for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
//        unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
//        sol_actflag_vol[i_vol] = sol_actflag[i_bdry];
//     }
 
//  KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap[pos_mu], L2G_dofmap[pos_ctrl], sol_actflag_vol);
 if (assembleMatrix) { KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap_mu_bdry, L2G_dofmap_ctrl_bdry, sol_actflag); }
 //============= delta_mu-delta_ctrl row - end ===============================

 //============= delta_mu-delta_mu row ===============================
 // Attention: this equation goes in contrast with \mu = 0 on \Omega \setminus \Gamma_c
 // In fact, here we shouldn't insert all VOLUME values, but only the BOUNDARY ones
 // The best way is to then do a L2G_map ON THE BOUNDARY only
 
  for (unsigned i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++) sol_actflag[i_bdry] =  ineq_flag * (1 - sol_actflag[i_bdry]/c_compl)  + (1-ineq_flag) * 1.;  //can do better to avoid division, maybe use modulo operator 

//  std::fill(sol_actflag_vol.begin(), sol_actflag_vol.end(), 0.);
//     for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
//        unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
//        sol_actflag_vol[i_vol] = sol_actflag[i_bdry];
//     }
  
//   KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap[pos_mu], L2G_dofmap[pos_mu], sol_actflag_vol );
  if (assembleMatrix) { KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap_mu_bdry, L2G_dofmap_mu_bdry, sol_actflag);  }
 //============= delta_mu-delta_mu row - end ===============================
  

}




 void node_insertion(const unsigned int iel,
                         const   Mesh* msh,
                         const     vector < vector < int > > & L2G_dofmap_Mat,
                         const unsigned int pos_mat_mu,
                         const unsigned int pos_mat_ctrl,
                         const std::vector < std::vector < double > > & sol_eldofs_Mat,
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
     
     
     
      //*************************************************** 
    std::vector < int > l2GMap_mu(sol_actflag.size());
    std::vector < int > l2GMap_ctrl(sol_actflag.size());
    for (unsigned i = 0; i < sol_actflag.size(); i++) {
      l2GMap_mu[i]   = L2G_dofmap_Mat[pos_mat_mu][i];   //pdeSys->GetSystemDof(solIndex_mu, solPdeIndex_mu, i, iel);
      l2GMap_ctrl[i] = L2G_dofmap_Mat[pos_mat_ctrl][i]; //pdeSys->GetSystemDof(solIndex_ctrl, solPdeIndex_ctrl, i, iel);
    } 
 //*************************************************** 

 //============= delta_mu row ===============================
      std::vector<double> Res_mu (sol_actflag.size()); std::fill(Res_mu.begin(),Res_mu.end(), 0.);
      
    for (unsigned i = 0; i < sol_actflag.size(); i++) {
      if (sol_actflag[i] == 0){  //inactive
         Res_mu [i] = - ineq_flag * ( 1. * sol_eldofs_Mat[pos_mat_mu][i] - 0. ); 
// 	 Res_mu [i] = Res[nDof_u + nDof_ctrl + nDof_adj + i]; 
      }
      else if (sol_actflag[i] == 1){  //active_a 
	 Res_mu [i] = - ineq_flag * ( c_compl *  sol_eldofs_Mat[pos_mat_ctrl][i] - c_compl * ctrl_lower[i]);
      }
      else if (sol_actflag[i] == 2){  //active_b 
	Res_mu [i]  =  - ineq_flag * ( c_compl *  sol_eldofs_Mat[pos_mat_ctrl][i] - c_compl * ctrl_upper[i]);
      }
    }
//          Res[nDof_u + nDof_ctrl + nDof_adj + i]  = c_compl * (  (2 - sol_actflag[i]) * (ctrl_lower[i] - sol_ctrl[i]) + ( sol_actflag[i] - 1 ) * (ctrl_upper[i] - sol_ctrl[i])  ) ;
//          Res_mu [i] = Res[nDof_u + nDof_ctrl + nDof_adj + i] ;

    
    RES->insert(Res_mu, l2GMap_mu);
//     RES->insert(Res_ctrl, l2GMap_ctrl);
//     RES->insert(Res_u, l2GMap_u);
//     RES->insert(Res_adj, l2GMap_adj);
    
//  //============= delta_state-delta_state row ===============================
//  KK->matrix_set_off_diagonal_values_blocked(l2GMap_u, l2GMap_u, 1.);

//  //============= delta_ctrl-delta_ctrl row ===============================
//  KK->matrix_set_off_diagonal_values_blocked(l2GMap_ctrl, l2GMap_ctrl, 1.);
 
//  //============= delta_adj-delta_adj row ===============================
//  KK->matrix_set_off_diagonal_values_blocked(l2GMap_adj, l2GMap_adj, 1.);
  
 //============= delta_mu-delta_ctrl row ===============================
 for (unsigned i = 0; i < sol_actflag.size(); i++) if (sol_actflag[i] != 0 ) sol_actflag[i] = ineq_flag * c_compl;    
  
  if (assembleMatrix) { KK->matrix_set_off_diagonal_values_blocked(l2GMap_mu, l2GMap_ctrl, sol_actflag); }

 //============= delta_mu-delta_mu row ===============================
  for (unsigned i = 0; i < sol_actflag.size(); i++) sol_actflag[i] =   ineq_flag * (1 - sol_actflag[i]/c_compl)  + (1-ineq_flag) * 1.;  //can do better to avoid division, maybe use modulo operator 

  if (assembleMatrix) {    KK->matrix_set_off_diagonal_values_blocked(l2GMap_mu, l2GMap_mu, sol_actflag );  }
  

     
     



 }



///@todo This is being added to a weak form?
///this is the same for volume or boundary
 void add_one_times_mu_res_ctrl(const unsigned iproc,
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
 
 
void el_dofs_quantities_vol(const Solution*                sol,
                        const Mesh * msh,
                        const unsigned int iel,
                        const    vector < unsigned > & SolFEType,
                        vector < unsigned > & Sol_n_el_dofs 
     ) {
    
        for (unsigned  k = 0; k < Sol_n_el_dofs.size(); k++) {
            unsigned  ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
            Sol_n_el_dofs[k] = ndofs_unk;
        }   
    
}


 
void el_dofs_unknowns_vol(const Solution*                sol,
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
  

//=============== construct control node flag field  =========================================    
  std::vector< std::vector< int > > is_dof_associated_to_boundary_control_equation(
      const Mesh * msh,
      /*const*/ MultiLevelSolution * ml_sol,
         const MultiLevelProblem *    ml_prob,
      const unsigned iel,
      CurrentElem < double > & geom_element_iel,
      const unsigned solType_coords,
      std::vector < std::string > Solname_Mat,      
      std::vector < unsigned > SolFEType_Mat,    
      std::vector < unsigned > Sol_n_el_dofs_Mat,    
      const unsigned pos_mat_ctrl,
      const unsigned n_components_ctrl

) {
	      /* For every component:
           * (control_node_flag[c][i])       picks nodes on \Gamma_c
           * (1 - control_node_flag[c][i])   picks nodes on \Omega \setminus \Gamma_c
	       */
          
       std::vector< std::vector< int > > control_node_flag(n_components_ctrl);
       
	  for (unsigned c = 0; c < n_components_ctrl; c++) {
              control_node_flag[c].resize(Sol_n_el_dofs_Mat[pos_mat_ctrl + c]);
              std::fill(control_node_flag[c].begin(), control_node_flag[c].end(), 0);   
         }
       
          
	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
          
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();

          
	      const unsigned int face_in_rectangle_domain = - ( msh->el->GetFaceElementIndex(iel,iface) + 1);
	     double tau = 0.;
         
         
         
	  for (unsigned c = 0; c < n_components_ctrl; c++) {
          
	      const bool  dir_bool = ml_sol->GetBdcFunctionMLProb()(ml_prob, geom_element_iel.get_elem_center_bdry_3d(), Solname_Mat[pos_mat_ctrl + c].c_str(), tau, face_in_rectangle_domain, 0.);

	      if (dir_bool == false) { 

          const unsigned ndofs_ctrl_bdry = msh->GetElementFaceDofNumber(iel, iface, SolFEType_Mat[pos_mat_ctrl + c]);
		  for(unsigned i_bdry = 0; i_bdry < ndofs_ctrl_bdry; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
		//we use the dirichlet flag to say: if dirichlet = true, we set 1 on the diagonal. if dirichlet = false, we put the boundary equation
		
// 		std::cout << " found boundary control nodes ==== " << std::endl;
			for(unsigned k = 0; k < control_node_flag[c].size(); k++) {
				  control_node_flag[c][i_vol] = 1;
			    }
              }
        }
      } 
        
        
          
    }
    
    return   control_node_flag;
    
  }
  
  //=============== construct control node flag field  =========================================    
  
  
  
  

  
  bool volume_elem_contains_a_boundary_control_face( const std::vector<double> & elem_center ) {

      int control_flag_jel = 0;
        control_flag_jel = ControlDomainFlag_bdry(elem_center);
        
        return (control_flag_jel == 1);

      
  }
  
  
  
  bool face_is_a_boundary_control_face(/*const*/ elem * el, const unsigned jel, const unsigned jface) {
      
  	    // look for boundary faces
            const int bdry_index_j = el->GetFaceElementIndex(jel, jface);
	    // look for face number equal to control face
	      const unsigned int face_in_rectangle_domain_j = - ( el->GetFaceElementIndex(jel,jface) + 1);

	    // look for boundary faces && look for control faces
		
   return ( bdry_index_j < 0 && face_in_rectangle_domain_j == FACE_FOR_CONTROL );
      
  }
  
  
  bool check_if_same_elem(const unsigned iel, const unsigned jel) {
      
   return (iel == jel);
      
  }
  
  bool check_if_same_elem_bdry(const unsigned iel, const unsigned jel, const unsigned iface, const unsigned jface) {

   return (iel == jel && iface == jface);
      
  }


  
  void  set_dense_pattern_for_unknowns(NonLinearImplicitSystemWithPrimalDualActiveSetMethod  & system, const std::vector < Unknown > unknowns)  {
  
///     system.init();   ///@todo Understand why I cannot put this here but it has to be in the main(), there must be some objects that get destroyed, passing the reference is not enough

  const MultiLevelProblem &  ml_prob = system.GetMLProb();
  const MultiLevelMesh *  ml_mesh = ml_prob.GetMLMesh();
  
  unsigned n_levels = ml_mesh->GetNumberOfLevels();
  std::ostringstream sp_out_base; sp_out_base << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "sp_";
  system._LinSolver[n_levels - 1]->sparsity_pattern_print_nonzeros(sp_out_base.str(), "on");
  system._LinSolver[n_levels - 1]->sparsity_pattern_print_nonzeros(sp_out_base.str(), "off");

  
  const Mesh* msh = ml_mesh->GetLevel(n_levels - 1);
  unsigned nprocs = msh->n_processors();
  unsigned iproc = msh->processor_id();

  
    for(int ivar = 0; ivar < unknowns.size(); ivar++) {
  
        if (unknowns[ivar]._is_sparse == false ) { 
        
  unsigned variable_index = system.GetSolPdeIndex(unknowns[ivar]._name.c_str());
  
  
  
  unsigned n_dofs_var_all_procs = 0;
  for (int ip = 0; ip < nprocs; ip++) {
     n_dofs_var_all_procs += system._LinSolver[n_levels - 1]->KKoffset[variable_index + 1][ip] - system._LinSolver[n_levels - 1]->KKoffset[variable_index][ip];
  // how does this depend on the number of levels and the number of processors? 
  // For the processors I summed over them and it seems to work fine
  // For the levels... should I pick the coarsest level instead of the finest one, or is it the same?
   } 

  unsigned n_vars = system.GetSolPdeIndex().size();   //assume all variables are dense: we should sum 3 sparse + 1 dense... or better n_components_ctrl * dense and the rest sparse
  //check that the dofs are picked correctly, it doesn't seem so 
  
  system.SetSparsityPatternMinimumSize (n_dofs_var_all_procs * n_vars, unknowns[ivar]._name);  ///@todo this is like AddSolution: it increases a vector
        }
    }
    
    
    
  return; 
  
  }
  
  
  
  
  
  
  
  void mixed_integral(const unsigned UNBOUNDED,
                      const unsigned dim,
                      const unsigned dim_bdry,
                      const std::vector < std::vector < double > > & ex_control,
                      const unsigned div, 
                      const double weight_iqp_bdry,
                      const std::vector < double > & x_iqp_bdry,
                      const std::vector < double > & phi_ctrl_iel_bdry_iqp_bdry,
                      const double sol_ctrl_iqp_bdry,
                      const double s_frac,
                      const double check_limits,
                      const double C_ns,
                      const unsigned int OP_Hhalf,
                      const double beta,
                      const unsigned int nDof_vol_iel,
                      std::vector < double > & KK_local_iel,
                      std::vector < double > & Res_local_iel,
                      std::vector < double > & KK_local_iel_mixed_num,
                      std::vector < double > & Res_local_iel_mixed_num,
                      const Mesh * msh,
                      const Solution *    sol,
                      const MultiLevelSolution *    ml_sol,
                      const int iel,
                      const int jel,
                      const unsigned iface,
                      std::vector <int> bdry_bdry,
                      CurrentElem < double > & geom_element_jel,
                      const unsigned jelGeom_bdry,
                      unsigned solType_coords
                     ) {
      
     
  const unsigned  sol_node_flag_index =  ml_sol->GetIndex("node_based_bdry_flag");
  const unsigned  group_salome = 2;   ///@todo fix here, maybe pass it in the args
      
      
  if(UNBOUNDED == 1) {
      
      //============ Mixed Integral 1D - Analytical ==================      
      if (dim_bdry == 1) {
      
//               double ex_1 = EX_1;
//               double ex_2 = EX_2;
//               std::vector < double > ex_1_vec(dim);
//               std::vector < double > ex_2_vec(dim);
//               ex_1_vec[0] = EX_1;
//               ex_1_vec[1] = 0.;
//               ex_2_vec[0] = EX_2;
//               ex_2_vec[1] = 0.;
// //               ex_1_vec[0] = 1.;
// //               ex_1_vec[1] = EX_1;
// //               ex_2_vec[0] = 1.;
// //               ex_2_vec[1] = EX_2;
              
              double dist2_1 = 0.;
              double dist2_2 = 0.;
              double dist_1 = 0.;  //distance from node to extreme 1
              double dist_2 = 0.;  //distance from node to extreme 2
              
//               for(int d = 0; d < dim; d++) {
//                 dist2_1 += (x_iqp_bdry[d] - ex_1_vec[d]) * (x_iqp_bdry[d] - ex_1_vec[d]);
//                 dist2_2 += (x_iqp_bdry[d] - ex_2_vec[d]) * (x_iqp_bdry[d] - ex_2_vec[d]);
//               }
//               // ex_control Built as  ( x1 y1 ) ( x2 y2 ) 
//               
              for(int d = 0; d < dim; d++) {
                dist2_1 += (x_iqp_bdry[d] - ex_control[d][0]) * (x_iqp_bdry[d] - ex_control[d][0]);
                dist2_2 += (x_iqp_bdry[d] - ex_control[d][1]) * (x_iqp_bdry[d] - ex_control[d][1]);
//                 std::cout<< ex_control[d][0] << "  " << ex_control[d][1] << "\n";
              }
              
              
              dist_1 = sqrt( dist2_1 );
              dist_2 = sqrt( dist2_2 );
              
              double mixed_term = pow(dist_1, -2. * s_frac) + pow(dist_2, - 2. * s_frac);

              for(unsigned i_bdry = 0; i_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); i_bdry++) {
                unsigned int i_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
                for(unsigned j_bdry = 0; j_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); j_bdry++) {
                  unsigned int j_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, j_bdry);
                  KK_local_iel[ i_vol_iel * nDof_vol_iel + j_vol_iel ] += 0.5 * C_ns * check_limits * (1. / s_frac) * OP_Hhalf * beta * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * phi_ctrl_iel_bdry_iqp_bdry[j_bdry] * weight_iqp_bdry * mixed_term;
                }
                Res_local_iel[ i_vol_iel ] += - 0.5 * C_ns * check_limits * (1. / s_frac) * OP_Hhalf  * beta * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * sol_ctrl_iqp_bdry * weight_iqp_bdry * mixed_term;
              }
   
      }
      
    //============ Mixed Integral 2D - Numerical ==================      
      else if (dim_bdry == 2) {           
          
            unsigned iel_geom_type = msh->GetElementType(iel);
            unsigned iel_geom_type_face = msh->GetElementFaceType(iel, iface);

            unsigned f_n_faces_faces  =  msh->el->GetNFC(iel_geom_type, iel_geom_type_face); /* ElementFaceFaceNumber */

         double mixed_term1 = 0.;
         
            // *** Face Gauss point loop (Integral over 2d object) ***
            for(unsigned e_bdry_bdry = 0; e_bdry_bdry < f_n_faces_faces/*bdry_bdry.size()*/; e_bdry_bdry++) {  ///@todo I think this has to be fixed

              // look for boundary faces
            

//               unsigned n_dofs_bdry_bdry = msh->el->GetNFC(LINE, solType_coords);
              
              unsigned n_dofs_bdry_bdry =  msh->el->GetNFACENODES(iel_geom_type_face/*jelGeom_bdry*/, e_bdry_bdry, solType_coords); //TODO ///@todo this is only taking linear nodes, do we want to do that for the coordinates too?

              vector  < vector  <  double> > delta_coordinates_bdry_bdry(dim);    // A matrix holding the face coordinates rowwise.
              for(int k = 0; k < dim; k++) {
                delta_coordinates_bdry_bdry[k].resize(n_dofs_bdry_bdry);
              }
              
                  std::vector < int > nodes_face_face_flags(n_dofs_bdry_bdry, 0); 
              
              for(unsigned i_bdry_bdry = 0; i_bdry_bdry < n_dofs_bdry_bdry; i_bdry_bdry++) {
                  
                unsigned inode_bdry = msh->el->GetIG(jelGeom_bdry, e_bdry_bdry/*bdry_bdry[e_bdry_bdry]*/, i_bdry_bdry); // face-to-element local node mapping. TODO: verify jelGeom_bdry
                unsigned inode_vol    = msh->el->GetIG(iel_geom_type, iface, inode_bdry); //from n-1 to n

                unsigned node_global = msh->el->GetElementDofIndex(jel, inode_vol);
                
                nodes_face_face_flags[i_bdry_bdry] = (*sol->_Sol[sol_node_flag_index])(node_global);
                
              // delta coords  -----
                for(unsigned k = 0; k < dim; k++) {
                  delta_coordinates_bdry_bdry[k][i_bdry_bdry] = geom_element_jel.get_coords_at_dofs_3d()[k][inode_vol] - x_iqp_bdry[k];  ///@todo// TODO We extract the local coordinates on the face from local coordinates on the element.
                }
              }
              
              
                bool is_face_bdry_bdry  =  MED_IO::boundary_of_boundary_3d_check_face_of_face_via_nodes( nodes_face_face_flags, group_salome);

              if(is_face_bdry_bdry) {
                
              // delta coords - refinement -----
              
// // //               std::vector  <  double > delta_coordinates_bdry_bdry_refined( (div + 1) * dim);
              
              vector  < vector  <  double > > delta_coordinates_bdry_bdry_refined(dim);
              for(int k = 0; k < dim; k++) {
                delta_coordinates_bdry_bdry_refined[k].resize(div + 1); // set "4" as a parameter
              }
              
              for(unsigned n = 0; n <= div; n++) {
                for(int k = 0; k < dim; k++) {
// // //                   delta_coordinates_bdry_bdry_refined[n + k * div] = delta_coordinates_bdry_bdry[k][0] + n * (delta_coordinates_bdry_bdry[k][1] - delta_coordinates_bdry_bdry[k][0]) /  div ;
                  delta_coordinates_bdry_bdry_refined[k][n] = delta_coordinates_bdry_bdry[k][0] + n * (delta_coordinates_bdry_bdry[k][1] - delta_coordinates_bdry_bdry[k][0]) /  div ;
                }
              }
              for(unsigned n = 0; n < div; n++) {
                  
                const unsigned dir_x_for_atan = ( ( (FACE_FOR_CONTROL - 1) / 2 ) + 1 ) % 3;  ///@todo I think needs to be changed
                const unsigned dir_y_for_atan = ( dir_x_for_atan + 1 ) % 3 ;  ///@todo I think needs to be changed
// // //                 double teta2 = atan2(delta_coordinates_bdry_bdry_refined[(n+1) + dir_y_for_atan * div], delta_coordinates_bdry_bdry_refined[(n+1) + dir_x_for_atan * div]);
// // //                 double teta1 = atan2(delta_coordinates_bdry_bdry_refined[n + dir_y_for_atan * div], delta_coordinates_bdry_bdry_refined[n + dir_x_for_atan * div]);
                double teta2 = atan2(delta_coordinates_bdry_bdry_refined[dir_y_for_atan][n + 1], delta_coordinates_bdry_bdry_refined[dir_x_for_atan][n + 1]);
                double teta1 = atan2(delta_coordinates_bdry_bdry_refined[dir_y_for_atan][n], delta_coordinates_bdry_bdry_refined[dir_x_for_atan][n]);

// // //                 double delta_teta = 0.;
// // //                 if(teta2 < teta1) delta_teta = std::min(teta1 - teta2, 2. * M_PI + teta2 - teta1);
// // //                 else delta_teta = std::min(teta2 - teta1, 2. * M_PI + teta1 - teta2);
                if(teta2 < teta1) {
                    teta2 += 2. * M_PI;
                }
                
                double delta_teta = teta2 - teta1;

                vector <double> mid_point;
                mid_point.resize(dim);
                for(unsigned k = 0; k < dim; k++) {
// // //                   mid_point[k] = (delta_coordinates_bdry_bdry_refined[(n+1) + k * div] + delta_coordinates_bdry_bdry_refined[n + k * div]) * 0.5;
                  mid_point[k] = (delta_coordinates_bdry_bdry_refined[k][n + 1] + delta_coordinates_bdry_bdry_refined[k][n]) * 0.5;
                }
                double dist2 = 0;
                for(int k = 0; k < dim; k++) {
                  dist2 += mid_point[k] * mid_point[k];
                }
                double dist = sqrt(dist2);
                
                mixed_term1 += 2. * pow(dist, -  2. * s_frac) * (1. / (2. * s_frac))  * delta_teta;
              }
//               delta coords - refinement -----

          
              }
              
            }
            
//             for(unsigned i = 0; i < phi_ctrl_iel_bdry_iqp_bdry.size(); i++) {
//               for(unsigned j = 0; j < phi_ctrl_iel_bdry_iqp_bdry.size(); j++) {
//                 KK_local_iel[ i * nDof_vol_iel + j ] += 0.5 * C_ns * check_limits * OP_Hhalf * beta * weight_iqp_bdry * phi_ctrl_iel_bdry_iqp_bdry[i] * phi_ctrl_iel_bdry_iqp_bdry[j] * mixed_term1;
//               }
//               Res_local_iel[ i ] += - 0.5 * C_ns * check_limits * OP_Hhalf * beta * weight_iqp_bdry * phi_ctrl_iel_bdry_iqp_bdry[i] * sol_ctrl_iqp_bdry * mixed_term1;
//             }
             for(unsigned i_bdry = 0; i_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); i_bdry++) {
                unsigned int i_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
                for(unsigned j_bdry = 0; j_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); j_bdry++) {
                  unsigned int j_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, j_bdry);
                  KK_local_iel_mixed_num[ i_vol_iel * nDof_vol_iel + j_vol_iel ] += 0.5 * C_ns * check_limits * OP_Hhalf * beta * weight_iqp_bdry * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * phi_ctrl_iel_bdry_iqp_bdry[j_bdry] * mixed_term1;
                }
                Res_local_iel_mixed_num[ i_vol_iel ] += - 0.5 * C_ns * check_limits * OP_Hhalf * beta * weight_iqp_bdry * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * sol_ctrl_iqp_bdry * mixed_term1;
              }
            
            
            
          
      }          
      
  }
  
  
 }
  
  
  
  
    //********** FRAC CONTROL - BEGIN *****************************************

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
                        const unsigned int max_size,
                        //-----------
                         const unsigned int n_quantities,
                        vector < unsigned > SolFEType_quantities,
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
                        //---- Control unknown -------
                        const unsigned int n_components_ctrl,
                        const unsigned int pos_mat_ctrl,
                        const unsigned int pos_sol_ctrl,
                        const unsigned int is_block_dctrl_ctrl_inside_main_big_assembly,
                        //-----------
                        SparseMatrix*  KK,
                        NumericVector* RES,
                        const bool assembleMatrix,
                        //-----------
                        const double alpha,
                        const double beta,
                        //-----------
                        const double s_frac,
                        const double check_limits,
                        const double C_ns,
                        const unsigned int OP_Hhalf,
                        const unsigned int OP_L2,
                        const double RHS_ONE,
                        const unsigned int UNBOUNDED,
                        //--- Control domain -------
                        const std::vector < std::vector < double > > & ex_control,
                        //--- Quadrature --------
                        const unsigned qrule_i,
                        const unsigned qrule_j,
                        const unsigned qrule_k,
                        const unsigned int Nsplit,
                        const unsigned int Quadrature_split_index,
                        const unsigned int N_div_unbounded
                       ) {
      
      
      
  unsigned solType =   SolFEType_Mat[pos_mat_ctrl];
  unsigned soluIndex = SolIndex_Mat[pos_mat_ctrl];
  unsigned soluPdeIndex = SolPdeIndex[pos_mat_ctrl];
  
  std::vector < double > sol_ctrl_iel;
  std::vector < double > sol_ctrl_jel;
    
 
 
  //-------- local to global mappings --------------
  vector< int > l2gMap_iel;  l2gMap_iel.reserve(max_size);
  vector< int > l2gMap_jel;  l2gMap_jel.reserve(max_size);


  //-------- Local matrices and rhs --------------
  vector < double > Res_local_iel; Res_local_iel.reserve(max_size);
  vector < double > KK_local_iel;  KK_local_iel.reserve(max_size * max_size);

//   Local matrices and rhs for adaptive quadrature (iel == jel)
  vector < double > Res_local_iel_refined; Res_local_iel_refined.reserve(max_size);
  vector < double > KK_local_iel_refined;   KK_local_iel_refined.reserve(max_size * max_size);

//   Local matrices and rhs for the mixed internal-external integral term (both adaptive and non-adaptive)
  vector < double > Res_local_iel_mixed_num;  Res_local_iel_mixed_num.reserve(max_size);
  vector < double > KK_local_iel_mixed_num;   KK_local_iel_mixed_num.reserve(max_size * max_size);

//   Non local matrices and vectors for H^s laplacian operator
  vector< double >         Res_nonlocal_iel;  Res_nonlocal_iel.reserve(max_size);
  vector< double >         Res_nonlocal_jel;  Res_nonlocal_jel.reserve(max_size);

  vector < double > KK_nonlocal_iel_iel;  KK_nonlocal_iel_iel.reserve(max_size * max_size);
  vector < double > KK_nonlocal_iel_jel;  KK_nonlocal_iel_jel.reserve(max_size * max_size);
  vector < double > KK_nonlocal_jel_iel;  KK_nonlocal_jel_iel.reserve(max_size * max_size);
  vector < double > KK_nonlocal_jel_jel;  KK_nonlocal_jel_jel.reserve(max_size * max_size); 
 
 //----------------------
  KK->zero();
  RES->zero(); 
 //----------------------

//   phi_x as unused input of certain functions
  vector < double > phi_x;
  phi_x.reserve(max_size * dim);
  
 
///   boost::mpi::communicator world(MPI_COMM_WORLD, boost::mpi::comm_attach);  /// @todo future solution: broadcast whole class instances


  unsigned count_visits_of_boundary_faces = 0;
  
   for(int kproc = 0; kproc < nprocs; kproc++) {
       
  const int proc_to_bcast_from = kproc;
  
       
    for(int jel = msh->_elementOffset[kproc]; jel < msh->_elementOffset[kproc + 1]; jel++) {

// --- geometry
   // all these little vectors are filled in one proc and broadcast to all       

        // --- 
        unsigned nDof_jel_coords;
         
      if (iproc == kproc) {
        nDof_jel_coords = msh->GetElementDofNumber(jel, solType_coords);
      }
     MPI_Bcast(& nDof_jel_coords, 1, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
    // ---       
        
    // --- 
         unsigned short  jel_geommm;
      if (iproc == kproc) {
       jel_geommm = msh->el->GetElementType(jel);
//         std::cout  << " current_proc " << iproc << " from external_proc " << kproc << " elem " << jel << " (before bcast) "  << jel_geommm << std::endl;
//         geom_element_jel.set_geom_type(jel);
//         jel_geom = geom_element_jel.geom_type();
     }
      MPI_Bcast(& jel_geommm, 1, MPI_UNSIGNED_SHORT, proc_to_bcast_from, MPI_COMM_WORLD);
//         std::cout << " current_proc " << iproc << " from external_proc " << kproc << " elem " << jel << " (after  bcast) "  << jel_geommm << std::endl;
        // --- 

        // --- coords - other way
        geom_element_jel.allocate_coords_at_dofs_3d(jel, nDof_jel_coords, solType_coords);
        
      if(kproc == iproc) {
        geom_element_jel.fill_coords_at_dofs_3d(jel, solType_coords);
      }
      for(unsigned k = 0; k < dim; k++) {
        MPI_Bcast(& geom_element_jel.get_coords_at_dofs()[k][0], nDof_jel_coords, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
      }
      for(unsigned k = 0; k < space_dim; k++) {
        MPI_Bcast(& geom_element_jel.get_coords_at_dofs_3d()[k][0], nDof_jel_coords, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
      }
// --- coords - other way

      if(kproc == iproc) {
        geom_element_jel.set_elem_center_3d(jel, solType_coords);
      }
        MPI_Bcast(& geom_element_jel.get_elem_center_3d()[0], space_dim, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);

// --- geometry        

        
// // // all of this is not used right now in this routine        
// // //       if(kproc == iproc) {
// // //  //***************************************************
// // //    el_dofs_unknowns_vol(sol, msh, pdeSys, jel,
// // //                         SolFEType_Mat,
// // //                         SolIndex_Mat,
// // //                         SolPdeIndex,
// // //                         Sol_n_el_dofs_Mat, 
// // //                         sol_eldofs_Mat,  
// // //                         L2G_dofmap_Mat);  //all unknowns here, perhaps we could restrict it to the ctrl components only
// // //       }
// // //         MPI_Bcast(& SolFEType_Mat[0], , MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
// // //         MPI_Bcast(& SolIndex_Mat[0], , MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
// // //         MPI_Bcast(& SolPdeIndex[0], , MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
// // //         MPI_Bcast(& Sol_n_el_dofs_Mat[0], , MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
// // //   //***************************************************

   
	// Perform face loop over elements that contain some control face
        
        
	if ( volume_elem_contains_a_boundary_control_face(geom_element_jel.get_elem_center_3d()) ) {

      
// ***************************************
// ******* jel-related stuff - BEGIN *************
// ***************************************

// --- 2 - solution -----------------
      unsigned nDof_jel;

      if(kproc == iproc) {
        nDof_jel  = msh->GetElementDofNumber(jel, solType);    // number of solution element dofs
      }

      MPI_Bcast(&nDof_jel, 1, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
      
      
      sol_ctrl_jel.resize(nDof_jel);
      
      if(kproc == iproc) {
        for(unsigned j = 0; j < nDof_jel; j++) {
          unsigned jDof  = msh->GetSolutionDof(j, jel, solType);
          sol_ctrl_jel[j] = (*sol->_Sol[soluIndex])(jDof);
        }
      }
      
      MPI_Bcast(& sol_ctrl_jel[0], nDof_jel, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
// --- 2 - solution -----------------

      
// --- 3 - l2GMap -----------------
        l2gMap_jel.resize(nDof_jel);

      // local storage of global mapping and solution ********************
      if(kproc == iproc) {
        for(unsigned j = 0; j < nDof_jel; j++) {
          l2gMap_jel[j] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, j, jel);  // global to global mapping between solution node and pdeSys dof
        }
      }
      MPI_Bcast(&l2gMap_jel[0], nDof_jel, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
      // ******************************************************************
// --- 3 - l2GMap -----------------

// ***************************************
// ******* jel-related stuff - END *************
// ***************************************


      
//------------------------------------        
//------------ jface opening ---------        
//------------------------------------        
      
// --- 
      unsigned n_faces_jel;
      if(kproc == iproc) {
          n_faces_jel = msh->GetElementFaceNumber(jel); 
    }
      MPI_Bcast(& n_faces_jel, 1, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
// ---
      
	  // loop on faces of the current element
	  for(unsigned jface = 0; jface < n_faces_jel; jface++) {
          
          
// --- geometry        
       unsigned jelGeom_bdry;
       if(kproc == iproc) {
           jelGeom_bdry = msh->GetElementFaceType(jel, jface);
        }  
      MPI_Bcast(& jelGeom_bdry, 1, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);

        unsigned nve_bdry;
       if(kproc == iproc) {
          nve_bdry  = msh->GetElementFaceDofNumber(jel, jface, solType_coords);
        }  
      MPI_Bcast(& nve_bdry, 1, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);

      
       geom_element_jel.allocate_coords_at_dofs_bdry_3d(jel, jface, nve_bdry);
      unsigned coords_at_dofs_bdry_3d_size = 0;
      
       if(kproc == iproc) {
       geom_element_jel.fill_coords_at_dofs_bdry_3d(jel, jface, solType_coords);
       coords_at_dofs_bdry_3d_size = geom_element_jel.get_coords_at_dofs_bdry_3d()[0].size();  ///@todo all coords have the same ndofs
        }  
      MPI_Bcast(& coords_at_dofs_bdry_3d_size, 1, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
      for(unsigned k = 0; k < space_dim; k++) {
         MPI_Bcast(& geom_element_jel.get_coords_at_dofs_bdry_3d()[k][0], coords_at_dofs_bdry_3d_size, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }
      
      
      
           std::fill(geom_element_jel.get_elem_center_bdry_3d().begin(), geom_element_jel.get_elem_center_bdry_3d().end(), 0.);

       if(kproc == iproc) {
       geom_element_jel.set_elem_center_bdry_3d();
        }  
      for(unsigned k = 0; k < space_dim; k++) {
        MPI_Bcast(& geom_element_jel.get_elem_center_bdry_3d()[0], space_dim, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
      }
// --- geometry        

     /*bool*/int jface_is_a_boundary_control;
       if(kproc == iproc) {
           jface_is_a_boundary_control = face_is_a_boundary_control_face(msh->el, jel, jface);
       }
      MPI_Bcast(& jface_is_a_boundary_control, 1, MPI_INTEGER, proc_to_bcast_from, MPI_COMM_WORLD);

      
	    if( jface_is_a_boundary_control ) {
//------------------------------------        
//------------ jface opening ---------    
//------------------------------------        

            

// // // //---- Quadrature in jqp_bdry, preparation right before iel - BEGIN ------- 

// Wait a second... If I want to prepare the jqp_bdry loop, I must be inside jface as well...
// So now it seems to me that I have to do jel jface iel iface instead...
// Previously it was jel iel iqp jqp
// Now, it has to be jel jface  - iel iface - iqp_bdry jqp_bdry
// The two quadrature loops must be the innermost. In this way you exclude all non-needed volume elements and all non-needed faces, so that you minimize the number of inner ifs. You keep them outside as much as possible
// There will be a storage of jqp_bdry

            
      const unsigned n_jqp_bdry = ml_prob.GetQuadratureRuleMultiple(qrule_j, jelGeom_bdry).GetGaussPointsNumber();

      std::vector < std::vector < double > > x_jqp_bdry(n_jqp_bdry);
      std::vector < double > weight_jqp_bdry(n_jqp_bdry);
      std::vector < std::vector < double > > phi_ctrl_jel_bdry_jqp_bdry(n_jqp_bdry);
      std::vector < std::vector < double > > phi_ctrl_x_jel_bdry_jqp_bdry(n_jqp_bdry);
      std::vector < std::vector < double > > phi_coords_jel_bdry_jqp_bdry(n_jqp_bdry);
      std::vector < std::vector < double > > phi_coords_x_jel_bdry_jqp_bdry(n_jqp_bdry);
      std::vector< double > sol_ctrl_jqp_bdry(n_jqp_bdry, 0.);


      
         for(unsigned jqp_bdry = 0; jqp_bdry < n_jqp_bdry; jqp_bdry++) {

    elem_all[qrule_j][jelGeom_bdry][solType_coords]->JacJacInv(geom_element_jel.get_coords_at_dofs_bdry_3d(), jqp_bdry, Jac_jel_bdry_jqp_bdry, JacI_jel_bdry_jqp_bdry, detJac_jel_bdry_jqp_bdry, space_dim);
    
    weight_jqp_bdry[jqp_bdry] = detJac_jel_bdry_jqp_bdry * ml_prob.GetQuadratureRuleMultiple(qrule_j, jelGeom_bdry).GetGaussWeightsPointer()[jqp_bdry];

    elem_all[qrule_j][jelGeom_bdry][SolFEType_quantities[pos_sol_ctrl]] ->shape_funcs_current_elem(jqp_bdry, JacI_jel_bdry_jqp_bdry, phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry], phi_ctrl_x_jel_bdry_jqp_bdry[jqp_bdry], boost::none, space_dim);
            
    elem_all[qrule_j][jelGeom_bdry][solType_coords] ->shape_funcs_current_elem(jqp_bdry, JacI_jel_bdry_jqp_bdry, phi_coords_jel_bdry_jqp_bdry[jqp_bdry], phi_coords_x_jel_bdry_jqp_bdry[jqp_bdry], boost::none, space_dim);

//========== compute gauss quantities on the boundary ===============================================
//--- geom
           x_jqp_bdry[jqp_bdry].assign(dim, 0.);
          
//        if(kproc == iproc) {
            for(unsigned d = 0; d < dim; d++) {
	      for (int j_bdry = 0; j_bdry < geom_element_jel.get_coords_at_dofs_bdry_3d()[d].size(); j_bdry++)  {
			
              x_jqp_bdry[jqp_bdry][d] += geom_element_jel.get_coords_at_dofs_bdry_3d()[d][j_bdry] * phi_coords_jel_bdry_jqp_bdry[jqp_bdry][j_bdry];

		      }
            }
//        }
//       MPI_Bcast(& x_jqp_bdry[jqp_bdry][0], dim, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
//--- geom
    
//--- solution
    sol_ctrl_jqp_bdry[jqp_bdry] = 0.;
//        if(kproc == iproc) {
	      for (int j_bdry = 0; j_bdry < phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry].size()/*Sol_n_el_dofs_quantities[pos_sol_ctrl]*/; j_bdry++)  {
		    unsigned int j_vol = msh->GetLocalFaceVertexIndex_PassElemType(jel_geommm, jface, j_bdry);
			
			sol_ctrl_jqp_bdry[jqp_bdry] +=  /*sol_eldofs_Mat[pos_mat_ctrl]*/sol_ctrl_jel[j_vol] * phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][j_bdry];

		      }
//        }
//       MPI_Bcast(& sol_ctrl_jqp_bdry[dim], 1, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
//--- solution
//========== compute gauss quantities on the boundary ================================================


        }  //jqp_bdry

// // //         //we can do the broadcast after the loop, faster
// // //         for(unsigned jqp_bdry = 0; jqp_bdry < n_jqp_bdry; jqp_bdry++) {
// // //             MPI_Bcast(& x_jqp_bdry[jqp_bdry][0], dim, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
// // //         }
// // //             MPI_Bcast(& sol_ctrl_jqp_bdry[0], n_jqp_bdry, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
  
        
// // // //---- Quadrature in jqp_bdry, preparation right before iel - END ------- 

    
// ---- boundary faces in jface: compute and broadcast - BEGIN ----
// This is only needed for when the boundary is a 2D face. We'll look at it later on
// look for what face of jface are on the boundary of the domain
// I believe we have to see deeply how we can extend this to the boundary case
// TODO: instead of faceIndex we want to have the new group condition that we impose 
// on the boundary of the boundary.


       std::vector< int > bdry_bdry(0);
///       unsigned n_faces;
/// 
///       if(iproc == kproc) {
///         for(unsigned j_bd_face = 0; j_bd_face < msh->GetElementFaceNumber_PassElemType(jelGeom_bdry); j_bd_face++) {
///           int faceIndex = msh->el->GetBoundaryIndex(jface, j_bd_face); /// TODO find a new condition and correct msh->GetElementFaceNumber ///@todo this is wrong
/// 
///           // look for boundary faces of the boundary
///           if(faceIndex >= 1) {
///             unsigned i = bdry_bdry.size();
///             bdry_bdry.resize(i + 1);
///             bdry_bdry[i] = jface;
///           }
///         }
///         n_faces = bdry_bdry.size();
///       }
/// 
///       MPI_Bcast(& n_faces, 1, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
/// 
///       bdry_bdry.resize(n_faces);
///       MPI_Bcast(& bdry_bdry[0], n_faces, MPI_INT, proc_to_bcast_from, MPI_COMM_WORLD);  

// ---- boundary faces in jface: compute and broadcast - END ----    


              
       for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
           
           
// --- geometry        
        geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);

        short unsigned ielGeom = geom_element_iel.geom_type();
              
        geom_element_iel.set_elem_center_3d(iel, solType_coords);
// --- geometry        

      
	// Perform face loop over elements that contain some control face
	if ( volume_elem_contains_a_boundary_control_face( geom_element_iel.get_elem_center_3d() ) ) {
        
// ***************************************
// ******* iel-related stuff - BEGIN *************
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
// ******* iel-related stuff - END *************
// ***************************************
        


        
//------------ iface opening ---------        
	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
          
// --- geom          
       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, iface);    
       
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();
// --- geom          

   
	    if( face_is_a_boundary_control_face(msh->el, iel, iface) ) {
//------------ iface opening ---------        
		
//                 count_visits_of_boundary_faces++;


             
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
          for(unsigned fe_type = 0; fe_type < /*solType*/solType_coords + 1; fe_type++) { //loop up to the FE type + 1 of the unknown
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
          std::vector < double > x_iqp_bdry(dim, 0.);  ///@todo is this dim or dim_bdry?

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
                KK_local_iel[ l_vol * nDof_iel + m_vol ] += OP_L2 * alpha * phi_ctrl_iel_bdry_iqp_bdry[l_bdry] * phi_ctrl_iel_bdry_iqp_bdry[m_bdry] * weight_iqp_bdry;
              }
              double mass_res_i = phi_ctrl_iel_bdry_iqp_bdry[l_bdry] * sol_ctrl_iqp_bdry ;
              Res_local_iel[ l_vol ] += OP_L2 * alpha * weight_iqp_bdry * mass_res_i ;
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
                
          /*const*/ short unsigned kelGeom_bdry = ielGeom_bdry;
           
          const unsigned n_kqp_bdry = ml_prob.GetQuadratureRuleMultiple(qrule_k, kelGeom_bdry).GetGaussPointsNumber();
                
                
          std::cout.precision(14);
          std::vector< std::vector< std::vector<double> > > x3;

          for(unsigned split = 0; split <= Nsplit; split++) {

            
            if (dim_bdry/*dim*/ == 1) GetElementPartition1D(x_iqp_bdry, geom_element_iel.get_coords_at_dofs_bdry_3d(), split, Nsplit, x3, space_dim); //TODO space_dim or dim?
            else if (dim_bdry/*dim*/ == 2) {
              //TODO need to be checked !!!
              //GetElementPartition2D(x_iqp_bdry, geom_element_iel.get_coords_at_dofs_bdry_3d(), split, Nsplit, x3);
              GetElementPartitionQuad(x_iqp_bdry, geom_element_iel.get_coords_at_dofs_bdry_3d(), split, Nsplit, x3);
            }

            //for(unsigned r = 0; r < size_part; r++) {
            for(unsigned r = 0; r < x3.size(); r++) {


              for(unsigned k_qp_bdry = 0; k_qp_bdry < n_kqp_bdry; k_qp_bdry++) {

// ********* PREPARATION PART - BEGIN ***************
                elem_all[qrule_k][kelGeom_bdry][solType_coords]->JacJacInv(x3[r]/*geom_element_iel.get_coords_at_dofs_bdry_3d()*/, k_qp_bdry, Jac_kel_bdry_kqp_bdry, JacI_kel_bdry_kqp_bdry, detJac_kel_bdry_kqp_bdry, space_dim);
    
//                 weight_kqp_bdry = detJac_kel_bdry_kqp_bdry * ml_prob.GetQuadratureRuleMultiple(qrule_k, kelGeom_bdry).GetGaussWeightsPointer()[k_qp_bdry];
                
                msh->_finiteElement[kelGeom_bdry][solType]->Jacobian(x3[r], k_qp_bdry, weight_kqp_bdry, phi_coords_kel_bdry_kqp_bdry, phi_x);
                
//                 elem_all[qrule_k][kelGeom_bdry][SolFEType_quantities[pos_sol_ctrl]] ->shape_funcs_current_elem(k_qp_bdry, JacI_kel_bdry_kqp_bdry, phi_ctrl_kel_bdry_kqp_bdry, phi_ctrl_x_kel_bdry_kqp_bdry, boost::none, space_dim);
            
//                 elem_all[qrule_k][kelGeom_bdry][solType_coords] ->shape_funcs_current_elem(k_qp_bdry, JacI_kel_bdry_kqp_bdry, phi_coords_kel_bdry_kqp_bdry, phi_coords_x_kel_bdry_kqp_bdry, boost::none, space_dim);

//--- geom
                vector < double > x_kqp_bdry(dim, 0.);  ///@todo is this dim or dim_bdry?

                for(unsigned d = 0; d < x_kqp_bdry.size(); d++) {
                  for (int k_bdry = 0; k_bdry < /*geom_element_iel.get_coords_at_dofs_bdry_3d()[d].size()*/phi_coords_kel_bdry_kqp_bdry.size(); k_bdry++)  {
			
                    x_kqp_bdry[d] += x3[r][d][k_bdry]/*geom_element_iel.get_coords_at_dofs_bdry_3d()[d][i_vol]*/ * phi_coords_kel_bdry_kqp_bdry[k_bdry];

                  }
                }
//--- geom
    
                std::vector<double> xi3(dim, 0.);

                GetClosestPointInReferenceElement(geom_element_iel.get_coords_at_dofs_bdry_3d(), x_kqp_bdry, kelGeom_bdry, xi3);
//                 GetInverseMapping(solType_coords, kelGeom_bdry, aP, x_kqp_bdry, xi3, 1000);  ///@todo generalize to rectangular Jacobian TODO is needed?

                msh->_finiteElement[kelGeom_bdry][/*solType*/ solType_coords]->GetPhi(phi_ctrl_kel_bdry_kqp_bdry, xi3); //TODO solType or solType_coords?

                double solY3 = 0.;
                for(unsigned i_bdry = 0; i_bdry < phi_ctrl_kel_bdry_kqp_bdry.size()/*nDof_iel*/; i_bdry++) {
                  unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
                  solY3 += sol_ctrl_iel[i_vol] * phi_ctrl_kel_bdry_kqp_bdry[i_bdry];
                }
    
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
                double dist_xyz3 = 0;
                for(unsigned k = 0; k < x_iqp_bdry.size(); k++) {
                  dist_xyz3 += (x_iqp_bdry[k] - x_kqp_bdry[k]) * (x_iqp_bdry[k] - x_kqp_bdry[k]);
                }

                const double denom_ik = pow(dist_xyz3, (double)( 0.5 * dim_bdry + s_frac));
                
                const double common_weight =  0.5 * C_ns * OP_Hhalf * beta * check_limits * weight_iqp_bdry * weight_kqp_bdry  / denom_ik;

//                 for(unsigned i = 0; i < nDof_iel; i++) {
                for(unsigned l_bdry = 0; l_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); l_bdry++) { //dofs of test function
                  unsigned int l_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, l_bdry);

                  Res_local_iel_refined[ l_vol_iel ]    +=      - common_weight * (sol_ctrl_iqp_bdry - solY3) * (phi_ctrl_iel_bdry_iqp_bdry[l_bdry] - phi_ctrl_kel_bdry_kqp_bdry[l_bdry]);

                for(unsigned m_bdry = 0; m_bdry < phi_ctrl_kel_bdry_kqp_bdry.size(); m_bdry++) { //dofs of unknown function
                    unsigned int m_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, m_bdry);
                    
                    KK_local_iel_refined[ l_vol_iel * nDof_jel + m_vol_iel ] += common_weight * (phi_ctrl_iel_bdry_iqp_bdry[m_bdry] - phi_ctrl_kel_bdry_kqp_bdry[m_bdry]) * 
                                                        (phi_ctrl_iel_bdry_iqp_bdry[l_bdry] - phi_ctrl_kel_bdry_kqp_bdry[l_bdry]);


                  }
                }
// ********* BOUNDED PART - END ***************
// ********* UNBOUNDED PART - BEGIN ***************
                if( iqp_bdry == Quadrature_split_index ) { ///@todo is there a way to put this outside of the quadrature loop?
                    
              mixed_integral(UNBOUNDED,
                              dim,
                              dim_bdry,
                              ex_control,
                              N_div_unbounded,
                              weight_kqp_bdry,
                              x_kqp_bdry,
                              phi_ctrl_kel_bdry_kqp_bdry,
                              sol_ctrl_kqp_bdry,
                              s_frac,
                              check_limits,
                              C_ns,
                              OP_Hhalf,
                              beta,
                              nDof_jel, ///@todo
//                               KK_local_iel_refined,
//                               Res_local_iel_refined,
                              KK_local_iel,
                              Res_local_iel,
                              KK_local_iel_mixed_num,
                              Res_local_iel_mixed_num,
                              msh,
                              sol,
                              ml_sol,
                              iel,
                             jel,
                              iface,
                              bdry_bdry,
                              geom_element_jel,
                              jelGeom_bdry,
                              solType_coords
                             ); 
                }


// ********* UNBOUNDED PART - END ***************
                                         
                  } //end k_qp_bdry
                } //end r
              }  //end split
              
              
                
                
                
            }  //end iel == jel && Nsplit != 0
      //============ Same elem, && Adaptive quadrature - END ==================
              
             
      //============ Either different elements, or lack of adaptivity (so all elements) - BEGIN ==================
        else {  //  if(iel != jel || Nsplit == 0) 
            
// ********* UNBOUNDED PART - BEGIN ***************
//           if(check_if_same_elem_bdry(iel, jel, iface, jface)) { //TODO I removed this since we don't want iel==jel here
              
               mixed_integral(UNBOUNDED,
                              dim,
                              dim_bdry,
                              ex_control,
                              N_div_unbounded,
                              weight_iqp_bdry,
                              x_iqp_bdry,
                              phi_ctrl_iel_bdry_iqp_bdry,
                              sol_ctrl_iqp_bdry,
                              s_frac,
                              check_limits,
                              C_ns,
                              OP_Hhalf,
                              beta,
                              nDof_iel,
                              KK_local_iel,
                              Res_local_iel,
                              KK_local_iel_mixed_num,
                              Res_local_iel_mixed_num,
                              msh,
                              sol,
                              ml_sol,
                              iel,
                              jel,
                              iface,
                              bdry_bdry,
                              geom_element_jel,
                              jelGeom_bdry,
                              solType_coords
                             );
               
//           }
              
// ********* UNBOUNDED PART - END ***************
            
// ********* BOUNDED PART - BEGIN ***************
            for(unsigned jqp_bdry = 0; jqp_bdry < n_jqp_bdry; jqp_bdry++) {

              double dist_xyz = 0.;
              for(unsigned d = 0; d < x_iqp_bdry.size(); d++) {
                dist_xyz += (x_iqp_bdry[d] - x_jqp_bdry[jqp_bdry][d]) * (x_iqp_bdry[d] - x_jqp_bdry[jqp_bdry][d]);
              }

              const double denom = pow(dist_xyz, (double)(  0.5 * /*dim*/dim_bdry + s_frac));
              
              const double common_weight = (0.5 * C_ns) * OP_Hhalf * beta * check_limits * weight_iqp_bdry * weight_jqp_bdry[jqp_bdry]  / denom;

//               for(unsigned i = 0; i < nDof_iel; i++) {
           for(unsigned l_bdry = 0; l_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); l_bdry++) { //dofs of test function
               		    unsigned int l_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, l_bdry);
               		    unsigned int l_vol_jel = msh->el->GetIG(jel_geommm, jface, l_bdry)/*msh->GetLocalFaceVertexIndex(jel, jface, l_bdry)*/;

                Res_nonlocal_iel[ l_vol_iel ]      +=      - common_weight * (sol_ctrl_iqp_bdry - sol_ctrl_jqp_bdry[jqp_bdry]) * (phi_ctrl_iel_bdry_iqp_bdry[l_bdry]);

                Res_nonlocal_jel[ l_vol_jel ]      +=      - common_weight * (sol_ctrl_iqp_bdry - sol_ctrl_jqp_bdry[jqp_bdry]) * (- phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][l_bdry]);

               
//                 for(unsigned j = 0; j < nDof_jel; j++) {
           for(unsigned m_bdry = 0; m_bdry < phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry].size(); m_bdry++) { //dofs of unknown function
               		    unsigned int m_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, m_bdry);
               		    unsigned int m_vol_jel = msh->el->GetIG(jel_geommm, jface, m_bdry)/*msh->GetLocalFaceVertexIndex(jel, jface, m_bdry)*/;

             /*  u(x) v(x)*/     KK_nonlocal_iel_iel[ l_vol_iel * nDof_jel + m_vol_iel ] += common_weight *          phi_ctrl_iel_bdry_iqp_bdry[m_bdry]            *    phi_ctrl_iel_bdry_iqp_bdry[l_bdry];

             /*- u(y) v(x)*/     KK_nonlocal_iel_jel[ l_vol_iel * nDof_jel + m_vol_jel ] += common_weight * (- 1.) * phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][m_bdry]  *    phi_ctrl_iel_bdry_iqp_bdry[l_bdry];

             /*- u(x) v(y)*/     KK_nonlocal_jel_iel[ l_vol_jel * nDof_jel + m_vol_iel ] += common_weight * (- 1.) * phi_ctrl_iel_bdry_iqp_bdry[m_bdry]            *   phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][l_bdry];

             /*  u(y) v(y)*/     KK_nonlocal_jel_jel[ l_vol_jel * nDof_jel + m_vol_jel ] += common_weight *          phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][m_bdry]  *    phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][l_bdry];


                  }
                }
                
              } //endl jqp_bdry loop
// ********* BOUNDED PART - END ***************
            
            
         } //end if(iel != jel || Nsplit == 0)
      //============ Either different elements, or lack of adaptivity (so all elements) - END ==================

        
         } //OP_Hhalf != 0              
      //============  Fractional assembly - END ==================
            
      }   //end iqp_bdry
      
               std::vector<unsigned> Sol_n_el_dofs_Mat_vol2(1, nDof_jel);
//          assemble_jacobian<double,double>::print_element_residual(iel, Res, Sol_n_el_dofs_Mat_vol, 10, 5);
//          assemble_jacobian<double,double>::print_element_jacobian(iel, KK_local_iel_mixed_num, Sol_n_el_dofs_Mat_vol2, 10, 5);
      
      
              
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


 if( check_if_same_elem(iel, jel)/*check_if_same_elem_bdry(iel, jel, iface, jface) if you put it inside the face loops */ ) {
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
        KK->add_matrix_blocked(KK_nonlocal_iel_jel, l2gMap_iel, l2gMap_jel);
        KK->add_matrix_blocked(KK_nonlocal_jel_iel, l2gMap_jel, l2gMap_iel);
        KK->add_matrix_blocked(KK_nonlocal_jel_jel, l2gMap_jel, l2gMap_jel);
// Since 1 is dense and 3 are sparse, and the dense dofs are 30, we should have at most 3x9 + 30 = 57, but in the sparsity print it shows 30. That's the problem


        RES->add_vector_blocked(Res_nonlocal_iel, l2gMap_iel);
        RES->add_vector_blocked(Res_nonlocal_jel, l2gMap_jel);
//============ add to global - END ==================

        
          
// //              if (print_algebra_local) {
//          std::vector<unsigned> Sol_n_el_dofs_Mat_vol2(1, nDof_jel);
// //          assemble_jacobian<double,double>::print_element_residual(iel, Res, Sol_n_el_dofs_Mat_vol, 10, 5);
//          assemble_jacobian<double,double>::print_element_jacobian(iel, KK_local_iel_mixed_num, Sol_n_el_dofs_Mat_vol2, 10, 5);
// //      }
         

        
//----- iel ---        
    } //end control elem flag i (control_flag_iel == 1)
  } //end iel
//----- iel ---        


//----- jface ---        
     
} //end if(bdry_index_j < 0)//end if(face_in_rectangle_domain_j == FACE_FOR_CONTROL)


      } //end jface

 //----- jface ---   
 
     
    }  //end control elem flag jel
    



   } //end jel
//----- jel ---        
   
 } //end kproc
    
    
    
  
  }
    
  //**********FRAC CONTROL - END *****************************************
  
  
    
  void control_eqn_bdry(const unsigned iproc,
                        MultiLevelProblem &    ml_prob,
                        MultiLevelSolution*    ml_sol,
                        const Solution*        sol,
                        const Mesh * msh,
                        const  LinearEquationSolver* pdeSys,
                        //----- Geom Element ------
                        CurrentElem < double > & geom_element_iel,
                        const unsigned int solType_coords,
                        const unsigned int space_dim,
                        //----- Mat ------
                        const unsigned int n_unknowns,
                        const    vector < std::string > & Solname_Mat,
                        const    vector < unsigned > & SolFEType_Mat,
                        const vector < unsigned > & SolIndex_Mat,
                        const vector < unsigned > & SolPdeIndex,
                        vector < unsigned > & Sol_n_el_dofs_Mat, 
                        vector < vector < double > > & sol_eldofs_Mat,  
                        vector < vector < int > > & L2G_dofmap_Mat,
                        //--- Equation, local --------
                        const unsigned max_size,
                        //----- Sol ------
                        const unsigned int n_quantities,
                        vector < unsigned > SolFEType_quantities,
                        //---- Quadrature - FE Evaluations -------
                        std::vector < std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > > elem_all,
                        //---- Quadrature ------
                        std::vector < std::vector < double > >  Jac_iel_bdry_iqp_bdry,
                        std::vector < std::vector < double > >  JacI_iel_bdry_iqp_bdry,
                        double detJac_iel_bdry_iqp_bdry,
                        double weight_iqp_bdry,
                        vector <double> phi_ctrl_iel_bdry_iqp_bdry,
                        vector <double> phi_ctrl_x_iel_bdry_iqp_bdry, 
                        //---- Control -------
                        const unsigned int n_components_ctrl,
                        const unsigned int pos_mat_ctrl,
                        const unsigned int pos_sol_ctrl,
                        //-----------
                        const unsigned int is_block_dctrl_ctrl_inside_main_big_assembly,
                        //-----------
                        SparseMatrix *  KK,
                        NumericVector * RES,
                        const bool assembleMatrix,
                        //-----------
                        const double alpha,
                        const double beta,
                        const double RHS_ONE,
                        //-----------
                        const unsigned qrule_i
                       ) {
      
   
   const unsigned  elem_dof_size_max = n_components_ctrl * max_size;
   
   std::vector < double >  Res;      Res.reserve( elem_dof_size_max );                         //should have Mat order
   std::vector < double >  Jac;      Jac.reserve( elem_dof_size_max * elem_dof_size_max);   //should have Mat order

 
    for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
        
// --- geometry        
            geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);

            const short unsigned ielGeom = geom_element_iel.geom_type();

           geom_element_iel.set_elem_center_3d(iel, solType_coords);
// --- geometry        

           
 //***************************************************
   el_dofs_unknowns_vol(sol, msh, pdeSys, iel,
                        SolFEType_Mat,
                        SolIndex_Mat,
                        SolPdeIndex,
                        Sol_n_el_dofs_Mat, 
                        sol_eldofs_Mat,  
                        L2G_dofmap_Mat);  //all unknowns here, perhaps we could restrict it to the ctrl components only
  //***************************************************
      
 //***************************************************
 //***************************************************
   //extract a subvector containing only the control components, starting from zero      
      
   std::vector< unsigned >::const_iterator first  = Sol_n_el_dofs_Mat.begin() + pos_mat_ctrl;
   std::vector< unsigned >::const_iterator last   = Sol_n_el_dofs_Mat.begin() + pos_mat_ctrl + n_components_ctrl;
   std::vector< unsigned > Sol_n_el_dofs_Mat_ctrl_only(first, last);
   
   unsigned int sum_Sol_n_el_dofs_ctrl_only = ElementJacRes< double >::compute_sum_n_dofs(Sol_n_el_dofs_Mat_ctrl_only);
 //***************************************************
 //***************************************************
   //create a subvector containing only the control components, starting from zero  
   std::vector< unsigned > L2G_dofmap_Mat_ctrl_only; 
   L2G_dofmap_Mat_ctrl_only.resize(0);
      for (unsigned  k = 0; k < n_components_ctrl; k++)     L2G_dofmap_Mat_ctrl_only.insert(L2G_dofmap_Mat_ctrl_only.end(), L2G_dofmap_Mat[pos_mat_ctrl + k].begin(), L2G_dofmap_Mat[pos_mat_ctrl + k].end());
 //***************************************************
 //***************************************************

   
 //***************************************************
    const unsigned int res_length =  n_components_ctrl * Sol_n_el_dofs_Mat[pos_mat_ctrl];

    Res.resize(res_length);                 std::fill(Res.begin(), Res.end(), 0.);

    Jac.resize(res_length * res_length);    std::fill(Jac.begin(), Jac.end(), 0.);
 //***************************************************

  
 //************ set control flag *********************
  std::vector< std::vector< int > > control_node_flag = 
       is_dof_associated_to_boundary_control_equation(msh, ml_sol, & ml_prob, iel, geom_element_iel, solType_coords, Solname_Mat, SolFEType_Mat, Sol_n_el_dofs_Mat, pos_mat_ctrl, n_components_ctrl);
  //*************************************************** 
      

    if ( volume_elem_contains_a_boundary_control_face(geom_element_iel.get_elem_center_3d()) ) {
        
	  
	  // loop on faces of the current element
	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
          
// ---
       std::vector<unsigned int> Sol_el_n_dofs_current_face(n_quantities); ///@todo the active flag is not an unknown!

       for (unsigned  k = 0; k < Sol_el_n_dofs_current_face.size(); k++) {
                 if (SolFEType_quantities[k] < 3) Sol_el_n_dofs_current_face[k] = msh->GetElementFaceDofNumber(iel, iface, SolFEType_quantities[k]);  ///@todo fix this absence
       }
       
// ---     
       
// --- geometry        
       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, iface);    
       
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();

// --- geometry        
         
	    if( face_is_a_boundary_control_face(msh->el, iel, iface) ) {
              

//========= initialize gauss quantities on the boundary ============================================
                std::vector< double > sol_ctrl_iqp_bdry(n_components_ctrl);
                std::vector< std::vector< double > > sol_ctrl_x_iqp_bdry(n_components_ctrl);
                
           for (unsigned c = 0; c < n_components_ctrl; c++) {
               sol_ctrl_iqp_bdry[c] = 0.;                
               sol_ctrl_x_iqp_bdry[c].assign(space_dim, 0.); 
            }
//========= initialize gauss quantities on the boundary ============================================
		
        const unsigned n_qp_bdry = ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussPointsNumber();
        
    
		for(unsigned iqp_bdry = 0; iqp_bdry < n_qp_bdry; iqp_bdry++) {
    
    elem_all[qrule_i][ielGeom_bdry][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_bdry_3d(), iqp_bdry, Jac_iel_bdry_iqp_bdry, JacI_iel_bdry_iqp_bdry, detJac_iel_bdry_iqp_bdry, space_dim);
// 	elem_all[qrule_i][ielGeom_bdry][solType_coords]->compute_normal(Jac_qp_bdry, normal);
    
    weight_iqp_bdry = detJac_iel_bdry_iqp_bdry * ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussWeightsPointer()[iqp_bdry];

    
   for (unsigned c = 0; c < n_components_ctrl; c++) {
    elem_all[qrule_i][ielGeom_bdry][SolFEType_quantities[pos_sol_ctrl + c]] ->shape_funcs_current_elem(iqp_bdry, JacI_iel_bdry_iqp_bdry, phi_ctrl_iel_bdry_iqp_bdry, phi_ctrl_x_iel_bdry_iqp_bdry, boost::none, space_dim);
   }  
   
//========== compute gauss quantities on the boundary ===============================================
   for (unsigned c = 0; c < n_components_ctrl; c++) {
          sol_ctrl_iqp_bdry[c] = 0.;
                  std::fill(sol_ctrl_x_iqp_bdry[c].begin(), sol_ctrl_x_iqp_bdry[c].end(), 0.);
		      for (int i_bdry = 0; i_bdry < phi_ctrl_iel_bdry_iqp_bdry.size()/*Sol_n_el_dofs_quantities[pos_sol_ctrl + c]*/; i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
			
			sol_ctrl_iqp_bdry[c] +=  sol_eldofs_Mat[pos_mat_ctrl + c][i_vol] * phi_ctrl_iel_bdry_iqp_bdry[i_bdry];
                            for (int d = 0; d < space_dim; d++) {
			      sol_ctrl_x_iqp_bdry[c][d] += sol_eldofs_Mat[pos_mat_ctrl + c][i_vol] * phi_ctrl_x_iel_bdry_iqp_bdry[i_bdry * space_dim + d];
			    }
		      }
   }
//========== compute gauss quantities on the boundary ================================================

//equation ************
   for (unsigned c = 0; c < n_components_ctrl; c++) {

		  // *** phi_i loop ***
		  for(unsigned i_bdry = 0; i_bdry < Sol_el_n_dofs_current_face[pos_sol_ctrl + c] ; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);

                 double lap_rhs_dctrl_ctrl_bdry_gss_i_c = 0.;
                 for (unsigned d = 0; d < space_dim; d++) {
                       /*if ( i_vol < Sol_n_el_dofs_Mat[pos_mat_ctrl] )*/  lap_rhs_dctrl_ctrl_bdry_gss_i_c +=  phi_ctrl_x_iel_bdry_iqp_bdry[i_bdry * space_dim + d] * sol_ctrl_x_iqp_bdry[c][d];
                 }
                 
		 
//============ Bdry Residuals - BEGIN ==================	
           const unsigned res_pos = assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_ctrl_only, /*pos_mat_ctrl +*/ c, i_vol);
           
                Res[ res_pos ]  +=  - control_node_flag[c][i_vol] *  weight_iqp_bdry *
                                         (     ( 1 - is_block_dctrl_ctrl_inside_main_big_assembly ) * alpha * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * sol_ctrl_iqp_bdry[c]
							                +  ( 1 - is_block_dctrl_ctrl_inside_main_big_assembly ) * beta  * lap_rhs_dctrl_ctrl_bdry_gss_i_c
							                         );  //boundary optimality condition
                Res[ res_pos ] += - RHS_ONE * weight_iqp_bdry * (phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * (-1.) /** ( sin(2 * acos(0.0) * x1[0][l_bdry])) * ( sin(2 * acos(0.0) * x1[1][l_bdry]))*/);
                          
//============ Bdry Residuals - END ==================    
		    
//============ Bdry Jacobians - BEGIN ==================	
   for (unsigned e = 0; e < n_components_ctrl; e++) {
        if (e == c) {
		    for(unsigned j_bdry = 0; j_bdry < Sol_el_n_dofs_current_face[pos_sol_ctrl + e]; j_bdry++) {
		         unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, iface, j_bdry);

//========block delta_control / control ========
              double  lap_mat_dctrl_ctrl_bdry_gss = 0.;
		      for (unsigned d = 0; d < space_dim; d++) {  
                  lap_mat_dctrl_ctrl_bdry_gss += phi_ctrl_x_iel_bdry_iqp_bdry[i_bdry * space_dim + d] * phi_ctrl_x_iel_bdry_iqp_bdry[j_bdry * space_dim + d];    
               }

          
           const unsigned jac_pos = assemble_jacobian< double, double >::jac_row_col_index(Sol_n_el_dofs_Mat_ctrl_only, sum_Sol_n_el_dofs_ctrl_only, /*pos_mat_ctrl +*/ c, /*pos_mat_ctrl +*/ e, i_vol, j_vol);
              Jac[ jac_pos ]   +=  control_node_flag[c][i_vol] *  weight_iqp_bdry * (
                                    ( 1 - is_block_dctrl_ctrl_inside_main_big_assembly ) * alpha * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * phi_ctrl_iel_bdry_iqp_bdry[j_bdry] 
			                      + ( 1 - is_block_dctrl_ctrl_inside_main_big_assembly ) * beta  * lap_mat_dctrl_ctrl_bdry_gss);   
				
	        }   //end j loop
	      }
   }
//============ Bdry Jacobians - END ==================	
	      
		  }  //end i loop
		  
   }  //end components ctrl	  
		  
		}  //end iqp_bdry loop
	  }    //end if control face
	      
	  }    //end loop over faces
	  
	} //end if control element flag

	
  /*is_res_control_only*/
                    RES->add_vector_blocked(Res, L2G_dofmap_Mat_ctrl_only);
              
                  if (assembleMatrix) {
                     KK->add_matrix_blocked(Jac, L2G_dofmap_Mat_ctrl_only, L2G_dofmap_Mat_ctrl_only);
                  }   
               
      
    }  //iel

    
}


  
  
  
  
  
} //end namespace
 
#endif
