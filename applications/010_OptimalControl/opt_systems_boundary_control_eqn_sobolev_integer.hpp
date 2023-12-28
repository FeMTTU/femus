#ifndef __opt_systems_boundary_control_eqn_hpp__
#define __opt_systems_boundary_control_eqn_hpp__



#include "Assemble_useful_functions.hpp"
#include "Assemble_unknown_jacres.hpp"



namespace femus {




   
namespace ctrl {
    


  
template < class LIST_OF_CTRL_FACES, class DOMAIN_CONTAINING_CTRL_FACES >
class Gamma_control_equation_integer_sobolev_differentiability_index {
  
public: 
 
 static void control_eqn_bdry(const unsigned iproc,
                        MultiLevelProblem &    ml_prob,
                        MultiLevelSolution*    ml_sol,
                        const Solution*        sol,
                        const Mesh * msh,
                        //----- Eqn ------
                        const  LinearEquationSolver* pdeSys,
                        //----- Geom Element ------
                        CurrentElem < double > & geom_element_iel,
                        const unsigned int solType_coords,
                        const unsigned int space_dim,
                        //----- Mat ------
                        const unsigned int n_unknowns,
                        const    std::vector < std::string > & Solname_Mat,
                        const    std::vector < unsigned > & SolFEType_Mat,
                        const std::vector < unsigned > & SolIndex_Mat,
                        const std::vector < unsigned > & SolPdeIndex,
                        std::vector < unsigned > & Sol_n_el_dofs_Mat, 
                        std::vector < std::vector < double > > & sol_eldofs_Mat,  
                        std::vector < std::vector < int > > & L2G_dofmap_Mat,
                        //--- Equation, local --------
                        const unsigned max_size,
                        //----- Sol ------
                        const unsigned int n_quantities,
                        std::vector < unsigned > SolFEType_quantities,
                        //---- Quadrature - FE Evaluations -------
                        std::vector < std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > > elem_all,
                        //---- Quadrature, Geometry ------
                        std::vector < std::vector < double > >  Jac_iel_bdry_iqp_bdry,
                        std::vector < std::vector < double > >  JacI_iel_bdry_iqp_bdry,
                        double detJac_iel_bdry_iqp_bdry,
                        double weight_iqp_bdry,
                        //---- Quadrature, Control ------
                        std::vector <double> phi_ctrl_iel_bdry_iqp_bdry,
                        std::vector <double> phi_ctrl_x_iel_bdry_iqp_bdry, 
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
                        const double rhs_one,
                        //-----------
                        const unsigned int operator_L2,
                        const unsigned int operator_H1semi,
                        //-----------
                        const unsigned qrule_i,
                        //-----------
                        const bool print_algebra_local
                       ) {
      
   
   const unsigned  elem_dof_size_max = n_components_ctrl * max_size;
   
   std::vector < double >  Res_ctrl_only;      Res_ctrl_only.reserve( elem_dof_size_max );                         //should have Mat order
   std::vector < double >  Jac_ctrl_only;      Jac_ctrl_only.reserve( elem_dof_size_max * elem_dof_size_max);   //should have Mat order

   
// integral - BEGIN ************
  double integral_alpha  = 0.;
  double integral_beta   = 0.;
// integral - END ************

 
    for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {
        
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

    Res_ctrl_only.resize(res_length);                 std::fill(Res_ctrl_only.begin(), Res_ctrl_only.end(), 0.);

    Jac_ctrl_only.resize(res_length * res_length);    std::fill(Jac_ctrl_only.begin(), Jac_ctrl_only.end(), 0.);
 //***************************************************


    if ( DOMAIN_CONTAINING_CTRL_FACES ::volume_elem_contains_a_Gamma_control_face(sol, msh, iel) ) {
        
  
 //************ set control flag - BEGIN *********************
  std::vector< std::vector< int > > control_node_flag_iel_all_faces = 
       femus::is_dof_associated_to_Gamma_control_equation(msh /*ok*/, ml_sol /*ok*/, & ml_prob /*ok*/, iel /*ok*/, ielGeom, geom_element_iel /*ok*/, solType_coords /*ok*/, Solname_Mat /*ok*/, SolFEType_Mat /*ok*/, Sol_n_el_dofs_Mat, pos_mat_ctrl /*ok*/, n_components_ctrl /*ok*/);
       
       ///@todo here I have to do it "on the go", for each boundary dof!!!
 //************ set control flag - END *********************
      
	  
	  // loop on faces of the current element
	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
          
// ---
       std::vector<unsigned int> Sol_el_n_dofs_current_face(n_quantities); ///@todo the active flag is not an unknown!

       for (unsigned  k = 0; k < Sol_el_n_dofs_current_face.size(); k++) {
                 if (SolFEType_quantities[k] < 3) Sol_el_n_dofs_current_face[k] = msh->GetElementFaceDofNumber(iel, iface, SolFEType_quantities[k]);  ///@todo fix this absence
       }
       
// ---     
       
// --- geometry - BEGIN        
       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, iface);    
       
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();

// --- geometry - END        
         
        std::pair< int, unsigned int > pair_control_iface = femus::face_is_a_Gamma_control_face_of_some_index< LIST_OF_CTRL_FACES >(msh->GetMeshElements(), iel, iface);

         int   iface_is_a_boundary_control  = pair_control_iface.first;

	    if( iface_is_a_boundary_control ) {
              

//========= initialize gauss quantities on the boundary - BEGIN  ============================================
                std::vector< double > sol_ctrl_iqp_bdry(n_components_ctrl);
                std::vector< std::vector< double > > sol_ctrl_x_iqp_bdry(n_components_ctrl);
                
           for (unsigned c = 0; c < n_components_ctrl; c++) {
               sol_ctrl_iqp_bdry[c] = 0.;                
               sol_ctrl_x_iqp_bdry[c].assign(space_dim, 0.); 
            }
//========= initialize gauss quantities on the boundary - END ============================================
		
        const unsigned n_qp_bdry = ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussPointsNumber();
        
    
		for(unsigned iqp_bdry = 0; iqp_bdry < n_qp_bdry; iqp_bdry++) {
    
    elem_all[qrule_i][ielGeom_bdry][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_bdry_3d(), iqp_bdry, Jac_iel_bdry_iqp_bdry, JacI_iel_bdry_iqp_bdry, detJac_iel_bdry_iqp_bdry, space_dim);
// 	elem_all[qrule_i][ielGeom_bdry][solType_coords]->compute_normal(Jac_qp_bdry, normal);
    
    weight_iqp_bdry = detJac_iel_bdry_iqp_bdry * ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussWeightsPointer()[iqp_bdry];

    
   for (unsigned c = 0; c < n_components_ctrl; c++) {
    elem_all[qrule_i][ielGeom_bdry][SolFEType_quantities[pos_sol_ctrl + c]] ->shape_funcs_current_elem(iqp_bdry, JacI_iel_bdry_iqp_bdry, phi_ctrl_iel_bdry_iqp_bdry, phi_ctrl_x_iel_bdry_iqp_bdry, boost::none, space_dim);
   }  
   
//========== compute gauss quantities on the boundary - BEGIN ===============================================
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
//========== compute gauss quantities on the boundary - END ================================================

// integral - BEGIN ************
   for (unsigned c = 0; c < n_components_ctrl; c++) {
                  integral_alpha += operator_L2 * weight_iqp_bdry * sol_ctrl_iqp_bdry[c] * sol_ctrl_iqp_bdry[c];
                 for (unsigned d = 0; d < space_dim; d++) {
                  integral_beta  +=  operator_H1semi * weight_iqp_bdry * sol_ctrl_x_iqp_bdry[c][d] * sol_ctrl_x_iqp_bdry[c][d];
      }
   }
  // integral - END ************


//equation - BEGIN ************
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
           
                Res_ctrl_only[ res_pos ]  +=  - control_node_flag_iel_all_faces[c][i_vol] *  weight_iqp_bdry *
                                         (     ( 1 - is_block_dctrl_ctrl_inside_main_big_assembly ) * operator_L2 * alpha * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * sol_ctrl_iqp_bdry[c]
							                +  ( 1 - is_block_dctrl_ctrl_inside_main_big_assembly ) * operator_H1semi * beta  * lap_rhs_dctrl_ctrl_bdry_gss_i_c
							                         );  //boundary optimality condition
                Res_ctrl_only[ res_pos ] += - rhs_one * weight_iqp_bdry * (phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * (-1.) /** ( sin(2 * acos(0.0) * x1[0][l_bdry])) * ( sin(2 * acos(0.0) * x1[1][l_bdry]))*/);
                          
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
              Jac_ctrl_only[ jac_pos ]   +=  control_node_flag_iel_all_faces[c][i_vol] *  weight_iqp_bdry * (
                                    ( 1 - is_block_dctrl_ctrl_inside_main_big_assembly ) * operator_L2 * alpha * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * phi_ctrl_iel_bdry_iqp_bdry[j_bdry]
			                      + ( 1 - is_block_dctrl_ctrl_inside_main_big_assembly ) * operator_H1semi * beta  * lap_mat_dctrl_ctrl_bdry_gss);
				
	        }   //end j loop
	      }
   }
//============ Bdry Jacobians - END ==================	
	      
		  }  //end i loop
		  
   }  //end components ctrl	  
//equation - END ************
		  
		}  //end iqp_bdry loop
	  }    //end if control face
	      
	  }    //end loop over faces
	  
	} //end if control element flag

	
  /*is_res_control_only*/
                     RES->add_vector_blocked(Res_ctrl_only, L2G_dofmap_Mat_ctrl_only);
              
                  if (assembleMatrix) {
                     KK->add_matrix_blocked(Jac_ctrl_only, L2G_dofmap_Mat_ctrl_only, L2G_dofmap_Mat_ctrl_only);
                  }   
         
         
         
        if (print_algebra_local) {
  
         assemble_jacobian<double,double>::print_element_residual(iel, Res_ctrl_only, Sol_n_el_dofs_Mat_ctrl_only, 10, 5);
         assemble_jacobian<double,double>::print_element_jacobian(iel, Jac_ctrl_only, Sol_n_el_dofs_Mat_ctrl_only, 10, 5);
        }
     
    }  //iel

    
// integral - BEGIN ************
  std::cout << std::endl;
  std::cout <<  "&&&&&& Integrals from previous iteration &&&&&&&&&&" << std::endl;
  double integral_alpha_parallel = 0.; MPI_Allreduce( &integral_alpha, &integral_alpha_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  std::cout << "The value of the integral_alpha  is " << std::setw(11) << std::setprecision(10) << 0.5 * integral_alpha_parallel << std::endl;
  double integral_beta_parallel = 0.;  MPI_Allreduce( &integral_beta, &integral_beta_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  std::cout << "The value of the integral_beta   is " << std::setw(11) << std::setprecision(10) << 0.5 * integral_beta_parallel << std::endl;
  std::cout <<  "&&&&&& Integrals from previous iteration &&&&&&&&&&" << std::endl;
  std::cout << std::endl;
// integral - END ************
   
    
    
}




};
 

 
   
 
} //end namespace ctrl




  
  
  
} //end namespace femus
 






#endif
