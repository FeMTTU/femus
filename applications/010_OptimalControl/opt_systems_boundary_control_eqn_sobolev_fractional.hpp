#ifndef __opt_systems_boundary_control_eqn_sobolev_fractional_hpp__
#define __opt_systems_boundary_control_eqn_sobolev_fractional_hpp__


#include "Assemble_useful_functions.hpp"
#include "opt_common.hpp"

//for reading additional fields from MED file (based on MED ordering)
#include "Files.hpp"
#include "MED_IO.hpp"

#include "PolynomialBases.hpp"


#include "fractional_functions.hpp"


#define TEST_JEL_SINGLE_FOR_LOOP  /*0*/iel

//*******************************************************************************************
//*********************** Mesh independent - BEGIN *****************************************
//*******************************************************************************************

namespace femus {

namespace ctrl {




template < class LIST_OF_CTRL_FACES, class DOMAIN_CONTAINING_CTRL_FACES >
 class Gamma_control_equation_fractional_sobolev_differentiability_index {
//*********************** Mesh independent, ALMOST - BEGIN *****************************************
 
 public:
  
  static void unbounded_integral_over_exterior_of_boundary_control_face(
//
//                       int & count,
//
                      const unsigned unbounded,
                      const unsigned dim,
                      const unsigned dim_bdry,
//////////
                      const unsigned int operator_Hhalf,
                      const double s_frac,
                      const double check_limits,
                      const double C_ns,
                      const double beta,
//////////
                      const Mesh * msh,
                      const Solution *    sol,
                      const MultiLevelSolution *    ml_sol,
//////////
//////////   Both 2D and 3D - BEGIN
                      const int iel,
                      const unsigned iface,
                      unsigned int i_element_face_index,
                      std::vector< std::vector< int > > control_node_flag_iel_all_faces,
                      std::vector< std::vector< int > > control_node_flag_jel_all_faces,
                      const std::vector<unsigned int> nDof_vol_iel, 
//////////
                      const double weight_qp_of_iface,
                      const std::vector < double > & x_qp_of_iface,
                      const std::vector < double > & phi_ctrl_iel_bdry_qp_of_iface,
                      const std::vector<double> sol_ctrl_qp_of_iface,
//////////   Both 2D and 3D - END
  //////////   2D only - BEGIN
                      std::vector < double > & KK_local_iel_unbounded_integral_analytical_both_ref_and_non_ref,
                      std::vector < double > & Res_local_iel_unbounded_integral_analytical_both_ref_and_non_ref,
  //////////   2D only - END
  //////////   3D only - BEGIN
                      const int jel,
                      const unsigned jface,
                      const unsigned int j_element_face_index,
                      CurrentElem < double > & geom_element_jel,
                      std::string node_based_bdry_bdry_in,
                      const unsigned n_divisions_face_of_face, 
                      std::vector < double > & KK_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref,
                      std::vector < double > & Res_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref,
                      unsigned solType_coords,
  //////////   3D only - END
                      double & integral
                     ) {
//       count++;
      
     
  const unsigned int n_components_ctrl = nDof_vol_iel.size();
  
  unsigned nDof_iel_vec = 0;
        for (unsigned c = 0; c < n_components_ctrl; c++) {nDof_iel_vec  +=  nDof_vol_iel[c];    }

  
      
  if(unbounded == 1) {
      
      //============ Mixed Integral 1D - Analytical ==================      
      if (dim_bdry == 1 && check_if_same_elem( TEST_JEL_SINGLE_FOR_LOOP , jel) ) {

// ---------- extreme coordinates, initialize - BEGIN
  std::vector < double > extreme_coords_Gamma_c_along_abscissa(2); //number of extremes of the 1D ctrl segment

  std::vector < std::vector < double > > extreme_coords_Gamma_c_xy_components(dim); // control zone extremes. Built as: 2D) ( x1 y1 ) ( x2 y2 )    3D) ( x1 y1 z1 ) ( x2 y2 z2 ) ( x3 y3 z3 )

  for(unsigned d = 0; d < extreme_coords_Gamma_c_xy_components.size(); d++) {
      extreme_coords_Gamma_c_xy_components[d].resize( extreme_coords_Gamma_c_along_abscissa.size() );
  }
// ---------- extreme coordinates, initialize - END
  
  
  //loop over element of the class 'list of face for control' to find the t-index of the current Gamma_c - BEGIN
  
  
  const unsigned int number_of_face_for_controls = LIST_OF_CTRL_FACES :: _face_with_extremes_index_size; // get the number of Gamma_c
  
  for(int t = 0; t < number_of_face_for_controls; t++) { 
     if( i_element_face_index == LIST_OF_CTRL_FACES :: _face_with_extremes_index[t] ) {

         for(int k = 0; k < extreme_coords_Gamma_c_along_abscissa.size(); k++){
             extreme_coords_Gamma_c_along_abscissa[k] = LIST_OF_CTRL_FACES :: _face_with_extremes_extremes_on_tang_surface[t][LIST_OF_CTRL_FACES ::  _num_of_tang_components_per_face -1][k];
         }

         int control_xyz                    = (LIST_OF_CTRL_FACES :: _face_with_extremes_index[t] - 1) / 2;
         bool is_face_on_maximal_coordinate = (LIST_OF_CTRL_FACES :: _face_with_extremes_index[t] - 1) % 2;
         for(unsigned d = 0; d < extreme_coords_Gamma_c_xy_components.size(); d++) {
                        for(unsigned n_e = 0; n_e < extreme_coords_Gamma_c_xy_components[d].size(); n_e++) {
                          if(control_xyz == d) extreme_coords_Gamma_c_xy_components[d][n_e] =  LIST_OF_CTRL_FACES :: normal_coordinate(is_face_on_maximal_coordinate);
                          else                 extreme_coords_Gamma_c_xy_components[d][n_e] = extreme_coords_Gamma_c_along_abscissa[n_e];
                        }
                      }
      }
    }
    
  //loop over element of the class 'list of face for control' to find the t-index of the current Gamma_c - END ---------




//--- Denominator - BEGIN -------

              
              double dist2_1 = 0.;
              double dist2_2 = 0.;
              double dist_1 = 0.;  //distance from node to extreme 1
              double dist_2 = 0.;  //distance from node to extreme 2
              
//               for(int d = 0; d < dim; d++) {
//                 dist2_1 += (x_qp_of_iface[d] - ex_1_vec[d]) * (x_qp_of_iface[d] - ex_1_vec[d]);
//                 dist2_2 += (x_qp_of_iface[d] - ex_2_vec[d]) * (x_qp_of_iface[d] - ex_2_vec[d]);
//               }
//               // extreme_coords_Gamma_c_xy_components Built as  ( x1 y1 ) ( x2 y2 ) 
//               
              for(int d = 0; d < dim; d++) {
                dist2_1 += (x_qp_of_iface[d] - extreme_coords_Gamma_c_xy_components[d][0]) * (x_qp_of_iface[d] - extreme_coords_Gamma_c_xy_components[d][0]);
                dist2_2 += (x_qp_of_iface[d] - extreme_coords_Gamma_c_xy_components[d][1]) * (x_qp_of_iface[d] - extreme_coords_Gamma_c_xy_components[d][1]);
//                 std::cout<< extreme_coords_Gamma_c_xy_components[d][0] << "  " << extreme_coords_Gamma_c_xy_components[d][1] << "\n";
              }
              
              
              dist_1 = sqrt( dist2_1 );
              dist_2 = sqrt( dist2_2 );
              
              double mixed_denominator = pow(dist_1, -2. * s_frac) + pow(dist_2, - 2. * s_frac);
//--- Denominator - END -------


              
//--- Integral - BEGIN -------
              for (unsigned c = 0; c < n_components_ctrl; c++) {
                  integral +=  0.5 * C_ns * check_limits * operator_Hhalf  * beta * sol_ctrl_qp_of_iface[c] * sol_ctrl_qp_of_iface[c] * weight_qp_of_iface  * mixed_denominator * (1. / s_frac);
              }   
//--- Integral - END -------
              
              
              
//--- Equation - BEGIN -------
              for (unsigned c = 0; c < n_components_ctrl; c++) {

              for(unsigned i_bdry = 0; i_bdry < phi_ctrl_iel_bdry_qp_of_iface.size(); i_bdry++) {
                unsigned int i_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
                
                const unsigned res_pos = assemble_jacobian<double,double>::res_row_index(nDof_vol_iel, c, i_vol_iel);

                Res_local_iel_unbounded_integral_analytical_both_ref_and_non_ref[ res_pos/*i_vol_iel*/ ] += - 0.5 *
                 control_node_flag_iel_all_faces[c][i_vol_iel] *
                C_ns * check_limits * operator_Hhalf  * beta  * phi_ctrl_iel_bdry_qp_of_iface[i_bdry] * sol_ctrl_qp_of_iface[c] * weight_qp_of_iface * mixed_denominator * (1. / s_frac);
                
                for (unsigned e = 0; e < n_components_ctrl; e++) {
                  if (e == c) { 
                     for(unsigned j_bdry = 0; j_bdry < phi_ctrl_iel_bdry_qp_of_iface.size(); j_bdry++) {
                  unsigned int j_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, j_bdry);
                  const unsigned jac_pos = assemble_jacobian< double, double >::jac_row_col_index(nDof_vol_iel, nDof_iel_vec, c, e, i_vol_iel, j_vol_iel);         

                  KK_local_iel_unbounded_integral_analytical_both_ref_and_non_ref[ jac_pos /*i_vol_iel * nDof_vol_iel + j_vol_iel*/ ] += 0.5 *
                   control_node_flag_iel_all_faces[c][i_vol_iel] * control_node_flag_jel_all_faces[c][j_vol_iel] *
                  C_ns * check_limits * operator_Hhalf * beta * phi_ctrl_iel_bdry_qp_of_iface[i_bdry] * phi_ctrl_iel_bdry_qp_of_iface[j_bdry] * weight_qp_of_iface * mixed_denominator * (1. / s_frac);
                     }
                
                  }
                }
                
              }
           }
//--- Equation - END -------
           
           
      }
      
    //============ Mixed Integral 2D - Numerical ==================      
      else if (dim_bdry == 2) {           
          
          
  // ctrl face to node-node association
 std::map<unsigned int, unsigned int >  ctrl_faces_VS_their_nodes = LIST_OF_CTRL_FACES :: from_ctrl_faces_to_their_boundaries();
 
 
          const unsigned t_face_that_I_want = j_element_face_index;

  //t_med_flag_of_node_bdry_bdry_for_control_face BEGIN
          
            unsigned  t_med_flag_of_node_bdry_bdry_for_control_face;
          
    for (unsigned m = 0; m < ml_sol->GetMLMesh()->_group_info_all_meshes.size(); m++)   {
    
                //find its corresponding node_node group
                 unsigned gr_node_node;
          for ( unsigned gr = 0; gr < ml_sol->GetMLMesh()->_group_info_all_meshes[m].size(); gr++) {
                      if ( ml_sol->GetMLMesh()->_group_info_all_meshes[m][gr]._user_defined_flag  == ctrl_faces_VS_their_nodes[ t_face_that_I_want ]) {
                          gr_node_node = gr;
                      }
          }
                
                t_med_flag_of_node_bdry_bdry_for_control_face = ml_sol->GetMLMesh()->_group_info_all_meshes[m][gr_node_node]._med_flag;
                
     }    
  //t_med_flag_of_node_bdry_bdry_for_control_face END

  
  // bdry bdry flag
  const unsigned  sol_node_flag_index =  ml_sol->GetIndex( node_based_bdry_bdry_in.c_str() );
 
//--- Denominator, unbounded integral, numerical - BEGIN -------
            unsigned jel_geom_type = msh->GetElementType(jel);
            unsigned jel_geom_type_face = msh->GetElementFaceType(jel, jface);

            unsigned jel_n_faces_faces  =  msh->GetMeshElements()->GetNFC(jel_geom_type, jel_geom_type_face); /* ElementFaceFaceNumber */

         double mixed_denominator_numerical = 0.;
         

         for(unsigned e_bdry_bdry = 0; e_bdry_bdry < jel_n_faces_faces; e_bdry_bdry++) {

              // look for boundary faces - BEGIN
              
              unsigned jel_n_dofs_bdry_bdry =  msh->GetMeshElements()->GetNFACENODES(jel_geom_type_face, e_bdry_bdry, solType_coords);

              // look for boundary of boundary faces - BEGIN
                  std::vector < int > nodes_face_face_flags(jel_n_dofs_bdry_bdry, 0); 
              
              for(unsigned jdof_bdry_bdry = 0; jdof_bdry_bdry < jel_n_dofs_bdry_bdry; jdof_bdry_bdry++) { 
                  
                unsigned jnode_bdry_bdry     = msh->GetMeshElements()->GetIG(jel_geom_type_face, e_bdry_bdry, jdof_bdry_bdry); // face-to-element local node mapping.
                unsigned jnode_bdry_bdry_vol = msh->GetMeshElements()->GetIG(jel_geom_type, jface, jnode_bdry_bdry);
                
              // nodes_face_face_flags  -----
                unsigned node_global = msh->GetMeshElements()->GetElementDofIndex(jel, jnode_bdry_bdry_vol);
                nodes_face_face_flags[jdof_bdry_bdry] = (*sol->_Sol[sol_node_flag_index])(node_global);
              }
              
              
                bool is_face_bdry_bdry  =  MED_IO::boundary_of_boundary_3d_check_face_of_face_via_nodes( nodes_face_face_flags, t_med_flag_of_node_bdry_bdry_for_control_face);
              // look for boundary of boundary faces - END
                

              if(is_face_bdry_bdry) {
                  
              // delta coords - BEGIN
              std::vector < std::vector <  double> > radius_centered_at_x_qp_of_iface_bdry_bdry(dim);    // A matrix holding the face coordinates rowwise.
              for(int k = 0; k < dim; k++) {
                radius_centered_at_x_qp_of_iface_bdry_bdry[k].resize(jel_n_dofs_bdry_bdry);
              }
              
              for(unsigned jdof_bdry_bdry = 0; jdof_bdry_bdry < jel_n_dofs_bdry_bdry; jdof_bdry_bdry++) { 
                  
                unsigned jnode_bdry_bdry     = msh->GetMeshElements()->GetIG(jel_geom_type_face, e_bdry_bdry, jdof_bdry_bdry); // face-to-element local node mapping.
                unsigned jnode_bdry_bdry_vol = msh->GetMeshElements()->GetIG(jel_geom_type, jface, jnode_bdry_bdry);
                
                for(unsigned k = 0; k < dim; k++) {
                  radius_centered_at_x_qp_of_iface_bdry_bdry[k][jdof_bdry_bdry] = geom_element_jel.get_coords_at_dofs_3d()[k][jnode_bdry_bdry_vol] - x_qp_of_iface[k];
                }
                
              }
              
                  
              // delta coords - END
                
              // delta coords - refinement - BEGIN -----
              std::vector < std::vector <  double > > radius_centered_at_x_qp_of_iface_bdry_bdry_refined(dim);
              for(int k = 0; k < dim; k++) {
                radius_centered_at_x_qp_of_iface_bdry_bdry_refined[k].resize(n_divisions_face_of_face + 1);
              }
              
                constexpr unsigned point_along_face_of_face_first = 0;
                constexpr unsigned point_along_face_of_face_last  = 1;

              
              for(unsigned n = 0; n <= n_divisions_face_of_face; n++) {
                for(int k = 0; k < dim; k++) {
                  const double increment_in_current_dim = 
                  (radius_centered_at_x_qp_of_iface_bdry_bdry[k][ point_along_face_of_face_last ] - 
                   radius_centered_at_x_qp_of_iface_bdry_bdry[k][ point_along_face_of_face_first]) /  n_divisions_face_of_face;
                  radius_centered_at_x_qp_of_iface_bdry_bdry_refined[k][n] = radius_centered_at_x_qp_of_iface_bdry_bdry[k][0] + n * increment_in_current_dim ;
                }
              }
              // delta coords - refinement - END -----
              
              
              // compute unbounded integral - BEGIN -----
              for(unsigned n = 0; n < n_divisions_face_of_face; n++) {
                
///@todooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
              // what are the global axes that are consistent with outgoing normal - BEGIN -----
                const unsigned current_jctrl_face = j_element_face_index /*LIST_OF_CTRL_FACES :: _face_with_extremes_index[0]*/; /*FACE_FOR_CONTROL*/

                std::vector< unsigned > global_dirs_for_atan(dim_bdry);
                constexpr unsigned global_dir_first = 0;
                constexpr unsigned global_dir_second = 1;
                global_dirs_for_atan[global_dir_first ] = ( ( ( current_jctrl_face  - 1) / 2 ) + 1 ) % 3;
                global_dirs_for_atan[global_dir_second] = ( global_dirs_for_atan[0] + 1 ) % 3 ;
                
                if ( (current_jctrl_face % 2) == 1 ) {  std::reverse(global_dirs_for_atan.begin(), global_dirs_for_atan.end()); } 
              //  what are the global axes that are consistent with outgoing normal - END -----
                
              // theta's - BEGIN -----
                std::vector< double > theta_first_and_last_radius(2);
                constexpr unsigned theta_of_radius_first = 0;
                constexpr unsigned theta_of_radius_second = 1;

                for(unsigned p = 0; p < theta_first_and_last_radius.size(); p++) {
                   theta_first_and_last_radius[p] = atan2(radius_centered_at_x_qp_of_iface_bdry_bdry_refined[ global_dirs_for_atan[global_dir_second] ][n + p],
                                                          radius_centered_at_x_qp_of_iface_bdry_bdry_refined[ global_dirs_for_atan[global_dir_first] ][n + p]);
                              }

// // //                 double delta_theta = 0.;
// // //                 if(theta_first_and_last_radius[1] < theta_first_and_last_radius[0]) delta_theta = std::min(theta_first_and_last_radius[0] - theta_first_and_last_radius[1], 2. * M_PI + theta_first_and_last_radius[1] - theta_first_and_last_radius[0]);
// // //                 else delta_theta = std::min(theta_first_and_last_radius[1] - theta_first_and_last_radius[0], 2. * M_PI + theta_first_and_last_radius[0] - theta_first_and_last_radius[1]);
                if(theta_first_and_last_radius[ theta_of_radius_second ] < theta_first_and_last_radius[ theta_of_radius_first ]) {
                    theta_first_and_last_radius[theta_of_radius_second ] += 2. * M_PI;
                }
                
                double delta_theta = theta_first_and_last_radius[ theta_of_radius_second ] -
                                     theta_first_and_last_radius[ theta_of_radius_first ];
              // theta's - END -----

               // distance |x - y| at midpoint - BEGIN -----
               std::vector <double> mid_point(dim, 0.);
                for(unsigned k = 0; k < dim; k++) {
                  mid_point[k] = (radius_centered_at_x_qp_of_iface_bdry_bdry_refined[k][n + 1] + radius_centered_at_x_qp_of_iface_bdry_bdry_refined[k][n]) * 0.5;
                }
                double dist2 = 0;
                for(int k = 0; k < dim; k++) {
                  dist2 += mid_point[k] * mid_point[k];
                }
                double dist = sqrt(dist2);
               // distance |x - y| at midpoint - END -----
                
               // integral - BEGIN -----
                mixed_denominator_numerical += pow(dist, -  2. * s_frac)  * delta_theta;
               // integral - END -----
                
              }
              // compute unbounded integral - END -----

          
              }
              
            }
            
            
            /// @todo ONLY DIFFERENCES: the mixed_denominator is numerical, and so also the corresponding Res and Jac. It could be done with a single function
//--- Denominator, unbounded integral, numerical - END -------
            
//--- Integral - BEGIN -------
              for (unsigned c = 0; c < n_components_ctrl; c++) {
                 integral +=  0.5 * C_ns * check_limits * operator_Hhalf  * beta  * sol_ctrl_qp_of_iface[c] * sol_ctrl_qp_of_iface[c] * weight_qp_of_iface * mixed_denominator_numerical * (1. / s_frac);
              }
//--- Integral - END -------
              
              
              
//--- Equation - BEGIN -------
              for (unsigned c = 0; c < n_components_ctrl; c++) {

              for(unsigned row_bdry = 0; row_bdry < phi_ctrl_iel_bdry_qp_of_iface.size(); row_bdry++) {
                unsigned int row_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, row_bdry);
                
                const unsigned res_pos = assemble_jacobian<double,double>::res_row_index(nDof_vol_iel, c, row_vol_iel);

                Res_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref[ res_pos/*row_vol_iel*/ ] += - 0.5 *
                 control_node_flag_iel_all_faces[c][row_vol_iel] *
                C_ns * check_limits * operator_Hhalf  * beta * phi_ctrl_iel_bdry_qp_of_iface[row_bdry] * sol_ctrl_qp_of_iface[c] * weight_qp_of_iface * mixed_denominator_numerical * (1. / s_frac);
                
                for (unsigned e = 0; e < n_components_ctrl; e++) {
                  if (e == c) { 
                     for(unsigned col_bdry = 0; col_bdry < phi_ctrl_iel_bdry_qp_of_iface.size(); col_bdry++) {
                  unsigned int col_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, col_bdry);
                  const unsigned jac_pos = assemble_jacobian< double, double >::jac_row_col_index(nDof_vol_iel, nDof_iel_vec, c, e, row_vol_iel, col_vol_iel);         

                  KK_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref[ jac_pos /*row_vol_iel * nDof_vol_iel + col_vol_iel*/ ] += 0.5 *
                   control_node_flag_iel_all_faces[c][col_vol_iel] *
                  C_ns * check_limits * operator_Hhalf * beta * phi_ctrl_iel_bdry_qp_of_iface[row_bdry] * phi_ctrl_iel_bdry_qp_of_iface[col_bdry] * weight_qp_of_iface * mixed_denominator_numerical * (1. / s_frac);
                     }
                
                  }
                }
                
              }
           }
//--- Equation - END -------
            
          
      }          
      
  }
  
  
 }
  
//********** BDRY_BDRY_FLAG_COPY_AND_DELETE- BEGIN *****************************************

  static   void   bdry_bdry_flag_copy_and_delete(   MultiLevelProblem & ml_prob,
                                             MultiLevelSolution & ml_sol,
                                           const MultiLevelMesh & ml_mesh, 
                                              const unsigned erased_levels,
                                             MultiLevelSolution *  ml_sol_bdry_bdry_flag,
                                             const std::string node_based_bdry_bdry_flag_name,
                                             const unsigned  steady_flag,
                                             const bool      is_an_unknown_of_a_pde,
                                             const FEFamily node_bdry_bdry_flag_fe_fam,
                                             const FEOrder node_bdry_bdry_flag_fe_ord) {
        
 ml_sol.AddSolution(node_based_bdry_bdry_flag_name.c_str(), node_bdry_bdry_flag_fe_fam, node_bdry_bdry_flag_fe_ord, steady_flag, is_an_unknown_of_a_pde);
  ml_sol.Initialize(node_based_bdry_bdry_flag_name.c_str());
  // copy ml_sol_bdry_bdry_flag at the non-removed levels into ml_sol
  for(unsigned l = 0; l < ml_mesh.GetNumberOfLevels(); l++) {
      *(ml_sol.GetSolutionLevel(l)->_Sol[ ml_sol.GetIndex(node_based_bdry_bdry_flag_name.c_str()) ]) =
      *(ml_sol_bdry_bdry_flag->GetSolutionLevel(l + erased_levels)->_Sol[ ml_sol_bdry_bdry_flag->GetIndex(node_based_bdry_bdry_flag_name.c_str()) ]);
  }
  
  delete ml_sol_bdry_bdry_flag;
     
     }  
//********** BDRY_BDRY_FLAG_COPY_AND_DELETE- END*****************************************

//********** BDRY_BDRY_FLAG- BEGIN *****************************************

  static    MultiLevelSolution *  bdry_bdry_flag(  Files & files,
                                            MultiLevelMesh & ml_mesh, 
                                               const std::string infile,
                                               std::vector < unsigned > & node_mapping_from_mesh_file_to_new,
                                             const std::string node_based_bdry_bdry_flag_name,
                                             const unsigned  steady_flag,
                                             const bool      is_an_unknown_of_a_pde,
                                             const FEFamily node_bdry_bdry_flag_fe_fam,
                                             const FEOrder node_bdry_bdry_flag_fe_ord) {
     
     
  MultiLevelSolution * ml_sol_bdry_bdry_flag = new MultiLevelSolution(&ml_mesh);
  ml_sol_bdry_bdry_flag->SetWriter(VTK);
  ml_sol_bdry_bdry_flag->GetWriter()->SetDebugOutput(true);
  

  ml_sol_bdry_bdry_flag->AddSolution(node_based_bdry_bdry_flag_name.c_str(), node_bdry_bdry_flag_fe_fam, node_bdry_bdry_flag_fe_ord, steady_flag, is_an_unknown_of_a_pde);
  ml_sol_bdry_bdry_flag->Initialize(node_based_bdry_bdry_flag_name.c_str());

  // ======= COARSE READING and REFINEMENT ========================
  ml_sol_bdry_bdry_flag->GetSolutionLevel(0)->GetSolutionByName( node_based_bdry_bdry_flag_name ) = MED_IO(*ml_mesh.GetLevel(0)).node_based_flag_read_from_file(infile, node_mapping_from_mesh_file_to_new);

  ml_mesh.GetLevelZero(0)->deallocate_node_mapping(node_mapping_from_mesh_file_to_new);

  for(unsigned l = 1; l < ml_mesh.GetNumberOfLevels(); l++) {
     ml_sol_bdry_bdry_flag->RefineSolution(l);
  }
  
  std::vector < std::string > variablesToBePrinted_aux;
  variablesToBePrinted_aux.push_back(node_based_bdry_bdry_flag_name/*"all"*/);
  for(unsigned l = 0; l < ml_mesh.GetNumberOfLevels(); l++) {
  ml_sol_bdry_bdry_flag->GetWriter()->Write("aux", files.GetOutputPath(), "", "biquadratic", variablesToBePrinted_aux, l+1);
   }
 
 
 return ml_sol_bdry_bdry_flag;
 
}
//********** BDRY_BDRY_FLAG- END*****************************************


  
  

  
//********** FRAC CONTROL - BEGIN *****************************************

 static  void control_eqn_bdry(const unsigned iproc,
                                   const unsigned nprocs,
                        MultiLevelProblem &    ml_prob,
                        MultiLevelSolution*    ml_sol,
                        const Solution*        sol,
                        /*const*/ Mesh * msh,
                        const  LinearEquationSolver* pdeSys,
                        //-----------
                        CurrentElem < double > & geom_element_iel,
                        const unsigned int solType_coords,
                        const unsigned int dim,
                        const unsigned int space_dim,
                        const unsigned int dim_bdry,
                        //-----------
                        const unsigned int n_unknowns,
                        const    std::vector < std::string > & Solname_Mat,
                        const    std::vector < unsigned > & SolFEType_Mat,
                        const std::vector < unsigned > & SolIndex_Mat,
                        const std::vector < unsigned > & SolPdeIndex,
                        std::vector < unsigned > & Sol_n_el_dofs_Mat, 
                        std::vector < std::vector < double > > & sol_eldofs_Mat,  
                        std::vector < std::vector < int > > & L2G_dofmap_Mat,
                        const unsigned int max_size,
                        //-----------
                         const unsigned int n_quantities,
                        std::vector < unsigned > SolFEType_quantities,
                        //-----------
                        std::vector < std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > > elem_all,
                        //-----------
                        std::vector < std::vector < double > >  Jac_iel_bdry_qp_of_iface,
                        std::vector < std::vector < double > >  JacI_iel_bdry_qp_of_iface,
                        double detJac_iel_bdry_qp_of_iface,
                        double weight_qp_of_iface,
                        std::vector <double> phi_ctrl_iel_bdry_qp_of_iface,
                        std::vector <double> phi_ctrl_x_iel_bdry_qp_of_iface, 
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
                        const unsigned use_Cns,
                        const unsigned int operator_Hhalf,
                        const unsigned int operator_L2,
                        const double rhs_one,
                        const unsigned int unbounded,
                        const std::string node_based_bdry_bdry_in,
                        //--- Quadrature --------
                        const unsigned qrule_i,
                        const unsigned qrule_j,
                        const unsigned qrule_k,
                        const unsigned int integration_num_split,
                        const unsigned int integration_split_index,
                        const unsigned int N_div_face_of_face,
                        //-----------
                        const bool print_algebra_local
                       ) {


// --- Fractional - BEGIN
const double C_ns =    compute_C_ns(dim_bdry, s_frac, use_Cns);  
// --- Fractional - END
      

// --- Geometry - BEGIN
  CurrentElem < double > geom_element_jel(dim, msh);            // must be adept if the domain is moving, otherwise double
// --- Geometry - END
      
// --- Quadrature - BEGIN
     std::vector < std::vector < double > >  JacI_jel_bdry_qp_of_jface(space_dim);
     std::vector < std::vector < double > >  Jac_jel_bdry_qp_of_jface(dim-1);
    for (unsigned d = 0; d < Jac_jel_bdry_qp_of_jface.size(); d++) {   Jac_jel_bdry_qp_of_jface[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_jel_bdry_qp_of_jface.size(); d++) { JacI_jel_bdry_qp_of_jface[d].resize(dim-1); }
    
    double detJac_jel_bdry_qp_of_jface;
// --- Quadrature - END


 //********************* bdry cont *******************
 //*************************************************** 
  std::vector <double> phi_coords_iel_bdry_qp_of_iface;  
  std::vector <double> phi_coords_x_iel_bdry_qp_of_iface; 

  phi_coords_iel_bdry_qp_of_iface.reserve(max_size);
  phi_coords_x_iel_bdry_qp_of_iface.reserve(max_size * space_dim);

 //*************************************************** 

  
  std::vector <unsigned> solType_ctrl(n_components_ctrl);
  std::vector <unsigned> solIndex_ctrl(n_components_ctrl);
  std::vector <unsigned> solPdeIndex_ctrl(n_components_ctrl);
  
  for (unsigned c = 0; c < n_components_ctrl; c++) {
          solType_ctrl[c]     = SolFEType_Mat[pos_mat_ctrl + c];
          solIndex_ctrl[c]    = SolIndex_Mat[pos_mat_ctrl + c];
          solPdeIndex_ctrl[c] = SolPdeIndex[pos_mat_ctrl + c];
  }
  
  std::vector < std::vector < double > > sol_ctrl_iel(n_components_ctrl);  ///@todo 1
  std::vector < std::vector < double > > sol_ctrl_jel(n_components_ctrl);
    
    for (unsigned c = 0; c < n_components_ctrl; c++) {
        sol_ctrl_iel[c].reserve(max_size);
        sol_ctrl_jel[c].reserve(max_size);
    }
    
    
   const unsigned  elem_dof_size_max = n_components_ctrl * max_size;

 
  //-------- local to global mappings - BEGIN --------------
  std::vector < std::vector < int > > l2gMap_iel(n_components_ctrl); 
  std::vector < std::vector < int > > l2gMap_jel(n_components_ctrl);

    for (unsigned c = 0; c < n_components_ctrl; c++) {
        l2gMap_iel[c].reserve(max_size);
        l2gMap_jel[c].reserve(max_size);
    }
  
  std::vector < int > l2gMap_iel_vec;  l2gMap_iel_vec.reserve(elem_dof_size_max);
  std::vector < int > l2gMap_jel_vec;  l2gMap_jel_vec.reserve(elem_dof_size_max);
  //-------- local to global mappings - END --------------


  //-------- Local matrices and rhs - BEGIN --------------
  std::vector < double > Res_local_iel_integer_operators; Res_local_iel_integer_operators.reserve(elem_dof_size_max);
  std::vector < double > KK_local_iel_integer_operators;  KK_local_iel_integer_operators.reserve(elem_dof_size_max * elem_dof_size_max);

//   Local matrices and rhs for adaptive quadrature (iel == jel)
  std::vector < double > Res_local_iel_only_refined; Res_local_iel_only_refined.reserve(elem_dof_size_max);
  std::vector < double > KK_local_iel_only_refined;   KK_local_iel_only_refined.reserve(elem_dof_size_max * elem_dof_size_max);

//   (both adaptive and non-adaptive)
  std::vector < double > Res_local_iel_unbounded_integral_analytical_both_ref_and_non_ref; Res_local_iel_unbounded_integral_analytical_both_ref_and_non_ref.reserve(elem_dof_size_max);
  std::vector < double > KK_local_iel_unbounded_integral_analytical_both_ref_and_non_ref;  KK_local_iel_unbounded_integral_analytical_both_ref_and_non_ref.reserve(elem_dof_size_max * elem_dof_size_max);

//   Local matrices and rhs for the mixed internal-external integral term (both adaptive and non-adaptive)
  std::vector < double > Res_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref;  Res_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref.reserve(elem_dof_size_max);
  std::vector < double > KK_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref;   KK_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref.reserve(elem_dof_size_max * elem_dof_size_max);

//   Non local matrices and vectors for H^s laplacian operator
  std::vector < double >         Res_nonlocal_iel_bounded_integral;  Res_nonlocal_iel_bounded_integral.reserve(elem_dof_size_max);
  std::vector < double >         Res_nonlocal_jel_bounded_integral;  Res_nonlocal_jel_bounded_integral.reserve(elem_dof_size_max);

  std::vector < double > KK_nonlocal_iel_iel_bounded_integral;  KK_nonlocal_iel_iel_bounded_integral.reserve(elem_dof_size_max * elem_dof_size_max);
  std::vector < double > KK_nonlocal_iel_jel_bounded_integral;  KK_nonlocal_iel_jel_bounded_integral.reserve(elem_dof_size_max * elem_dof_size_max);
  std::vector < double > KK_nonlocal_jel_iel_bounded_integral;  KK_nonlocal_jel_iel_bounded_integral.reserve(elem_dof_size_max * elem_dof_size_max);
  std::vector < double > KK_nonlocal_jel_jel_bounded_integral;  KK_nonlocal_jel_jel_bounded_integral.reserve(elem_dof_size_max * elem_dof_size_max); 
  //-------- Local matrices and rhs - END --------------

  
  //-------- Global matrices and rhs - BEGIN --------------
  KK->zero();
  RES->zero(); 
  //-------- Global matrices and rhs - END --------------

  
//   phi_x_placeholder_unused as unused input of certain functions
  std::vector < double > phi_x_placeholder_unused;
  phi_x_placeholder_unused.reserve(max_size * dim);
  
 
///   boost::mpi::communicator world(MPI_COMM_WORLD, boost::mpi::comm_attach);  /// @todo future solution: broadcast whole class instances


  unsigned count_visits_of_boundary_faces = 0;


// integral - BEGIN ************
  double integral  = 0.;
// integral - END ************
      

  
   for(int jproc = 0; jproc < nprocs; jproc++) {
       
  const int proc_to_bcast_from = jproc;
  
// --- jel opening - BEGIN
    for(int jel = msh->GetElementOffset(jproc); jel < msh->GetElementOffset(jproc + 1); jel++) {


// ------- Sol_n_el_dofs_Mat - BEGIN
    if (iproc == jproc) {
        for (unsigned  k = 0; k < Sol_n_el_dofs_Mat.size(); k++) {
                    unsigned  ndofs_unk = msh->GetElementDofNumber(jel, SolFEType_Mat[k]);
                    Sol_n_el_dofs_Mat[k] = ndofs_unk;
        }
    }
    MPI_Bcast(& Sol_n_el_dofs_Mat[0], n_unknowns , MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
// ------ Sol_n_el_dofs_Mat - END


// --- geometry - BEGIN 
   // all these little vectors are filled in one proc and broadcast to all       

// --- - BEGIN  
        unsigned nDof_jel_coords;
         
      if (iproc == jproc) {
        nDof_jel_coords = msh->GetElementDofNumber(jel, solType_coords);
      }
     MPI_Bcast(& nDof_jel_coords, 1, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
// ---  - END       
        
// ---  - BEGIN
         unsigned short  jel_geommm;
      if (iproc == jproc) {
       jel_geommm = msh->GetElementType(jel);
//         std::cout  << " current_proc " << iproc << " from external_proc " << jproc << " elem " << jel << " (before bcast) "  << jel_geommm << std::endl;
//         geom_element_jel.set_geom_type(jel);
//         jel_geom = geom_element_jel.geom_type();
     }
      MPI_Bcast(& jel_geommm, 1, MPI_UNSIGNED_SHORT, proc_to_bcast_from, MPI_COMM_WORLD);
//         std::cout << " current_proc " << iproc << " from external_proc " << jproc << " elem " << jel << " (after  bcast) "  << jel_geommm << std::endl;
// ---  - END

// --- coords - other way - BEGIN
        geom_element_jel.allocate_coords_at_dofs_3d(jel, nDof_jel_coords, solType_coords);
        
      if(jproc == iproc) {
        geom_element_jel.fill_coords_at_dofs_3d(jel, solType_coords);
      }
      for(unsigned k = 0; k < dim; k++) {
        MPI_Bcast(& geom_element_jel.get_coords_at_dofs()[k][0], nDof_jel_coords, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
      }
      for(unsigned k = 0; k < space_dim; k++) {
        MPI_Bcast(& geom_element_jel.get_coords_at_dofs_3d()[k][0], nDof_jel_coords, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
      }
// --- coords - other way - END

// --- coords - elem center - BEGIN
      if(jproc == iproc) {
        geom_element_jel.set_elem_center_3d(jel, solType_coords);
      }
        MPI_Bcast(& geom_element_jel.get_elem_center_3d()[0], space_dim, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
// --- coords - elem center - END

// --- geometry - END

// comment: el_dofs_unknowns_vol - BEGIN
// // // all of this is not used right now in this routine        
// // //       if(jproc == iproc) {
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
// comment: el_dofs_unknowns_vol - END
   
	// Perform face loop over elements that contain some control face
        bool vol_elem_contains_control;
        if(jproc == iproc) {
           vol_elem_contains_control =  DOMAIN_CONTAINING_CTRL_FACES :: volume_elem_contains_a_Gamma_control_face(sol, msh, jel);
        }
      MPI_Bcast(& vol_elem_contains_control, 1, MPI_C_BOOL, proc_to_bcast_from, MPI_COMM_WORLD);
        
        
        
	if ( vol_elem_contains_control ) {
// --- jel opening - END

//************ set control flag - BEGIN *********************
       std::vector< std::vector< int > > control_node_flag_jel_all_faces =
       femus::is_dof_associated_to_Gamma_control_equation(msh /*ok*/, ml_sol /*ok*/, & ml_prob /*ok*/, jel /*ok*/, jel_geommm, geom_element_jel /*ok*/, solType_coords /*ok*/, Solname_Mat /*ok*/, SolFEType_Mat /*ok*/, Sol_n_el_dofs_Mat, pos_mat_ctrl /*ok*/, n_components_ctrl /*ok*/);
//************ set control flag - END *********************

// ***************************************
// ******* jel-related stuff - BEGIN *************
// ***************************************

// --- 2 - solution - BEGIN -----------------
      std::vector < unsigned > nDof_jel(n_components_ctrl);

      if(jproc == iproc) {
        for (unsigned c = 0; c < n_components_ctrl; c++) {
           nDof_jel[c]  = msh->GetElementDofNumber(jel, solType_ctrl[c]);    // number of solution element dofs
           }
    }

      MPI_Bcast(& nDof_jel[0], n_components_ctrl, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
      
      
      unsigned nDof_jel_vec = 0;
        for (unsigned c = 0; c < n_components_ctrl; c++) { nDof_jel_vec  +=  nDof_jel[c];    }

      
      
        for (unsigned c = 0; c < n_components_ctrl; c++) {
           sol_ctrl_jel[c].resize(nDof_jel[c]);
        }
        
      if(jproc == iproc) {
          
      for (unsigned c = 0; c < n_components_ctrl; c++) {
       for(unsigned j = 0; j < nDof_jel[c]; j++) {
          unsigned jDof  = msh->GetSolutionDof(j, jel, solType_ctrl[c]);
          sol_ctrl_jel[c][j] = (*sol->_Sol[solIndex_ctrl[c]])(jDof);
           }
         }
      }
      
      for (unsigned c = 0; c < n_components_ctrl; c++) {
         MPI_Bcast(& sol_ctrl_jel[c][0], nDof_jel[c], MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
      }
// --- 2 - solution - END -----------------

      
// --- 3 - l2GMap - BEGIN -----------------
      for (unsigned c = 0; c < n_components_ctrl; c++) {
          l2gMap_jel[c].resize(nDof_jel[c]);
      }
           
           
      // local storage of global mapping and solution ********************
      if(jproc == iproc) {
       for (unsigned c = 0; c < n_components_ctrl; c++) {
         for(unsigned j = 0; j < nDof_jel[c]; j++) {
          l2gMap_jel[c][j] = pdeSys->GetSystemDof(solIndex_ctrl[c], solPdeIndex_ctrl[c], j, jel);  // global to global mapping between solution node and pdeSys dof
         }
       }
    }
    
      for (unsigned c = 0; c < n_components_ctrl; c++) {
         MPI_Bcast(& l2gMap_jel[c][0], nDof_jel[c], MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
      }
      // ******************************************************************
// --- 3 - l2GMap - END -----------------



// --- l2GMapVec - BEGIN -----------------
    l2gMap_jel_vec.resize(0);
      for (unsigned c = 0; c < n_components_ctrl; c++) {
          l2gMap_jel_vec.insert(l2gMap_jel_vec.end(), l2gMap_jel[c].begin(), l2gMap_jel[c].end());
      }
// --- l2GMapVec - END -----------------



// ***************************************
// ******* jel-related stuff - END *************
// ***************************************


      
//------------------------------------        
//------------ jface opening - BEGIN ---------        
//------------------------------------        
      
// --- - BEGIN  
      unsigned n_faces_jel;
      if(jproc == iproc) {
          n_faces_jel = msh->GetElementFaceNumber(jel); 
    }
      MPI_Bcast(& n_faces_jel, 1, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
// --- - END
      
	  // loop on faces of the current element
	  for(unsigned jface = 0; jface < n_faces_jel; jface++) {

// // // //            const unsigned int number_of_face_for_controls = LIST_OF_CTRL_FACES :: _face_with_extremes_index_size; // get the number of Gamma_c
// // // //            for(int t = 0; t < number_of_face_for_controls; t++) {

          
// --- geometry about jface - BEGIN         
       unsigned jelGeom_bdry;
       if(jproc == iproc) {
           jelGeom_bdry = msh->GetElementFaceType(jel, jface);
        }  
      MPI_Bcast(& jelGeom_bdry, 1, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);

        unsigned nve_bdry;
       if(jproc == iproc) {
          nve_bdry  = msh->GetElementFaceDofNumber(jel, jface, solType_coords);
        }
      MPI_Bcast(& nve_bdry, 1, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);


       geom_element_jel.allocate_coords_at_dofs_bdry_3d(jel, jface, nve_bdry);
      unsigned coords_at_dofs_bdry_3d_size = 0;

       if(jproc == iproc) {
         geom_element_jel.fill_coords_at_dofs_bdry_3d(jel, jface, solType_coords);
         coords_at_dofs_bdry_3d_size = geom_element_jel.get_coords_at_dofs_bdry_3d()[0].size();  ///@todo all coords have the same ndofs
        }
      MPI_Bcast(& coords_at_dofs_bdry_3d_size, 1, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);

      for(unsigned k = 0; k < space_dim; k++) {
         MPI_Bcast(& geom_element_jel.get_coords_at_dofs_bdry_3d()[k][0], coords_at_dofs_bdry_3d_size, MPI_DOUBLE, jproc, MPI_COMM_WORLD);
      }


           std::fill(geom_element_jel.get_elem_center_bdry_3d().begin(), geom_element_jel.get_elem_center_bdry_3d().end(), 0.);

       if(jproc == iproc) {
       geom_element_jel.set_elem_center_bdry_3d();
        }
      for(unsigned k = 0; k < space_dim; k++) {
        MPI_Bcast(& geom_element_jel.get_elem_center_bdry_3d()[0], space_dim, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
      }
// --- geometry about jface - END

// --- jface is control - BEGIN
//        std::make_pair jface_is_a_boundary_control

     /*bool*/int    jface_is_a_boundary_control;
       unsigned int jface_boundary_control_index;

       if(jproc == iproc) {
           
           std::pair< int, unsigned int > pair_control_jface = femus::face_is_a_Gamma_control_face_of_some_index< LIST_OF_CTRL_FACES >(msh->GetMeshElements(), jel, jface);
           
           jface_is_a_boundary_control  = pair_control_jface.first;
           jface_boundary_control_index = pair_control_jface.second;

       }
      MPI_Bcast(& jface_is_a_boundary_control, 1, MPI_INTEGER, proc_to_bcast_from, MPI_COMM_WORLD);
      MPI_Bcast(& jface_boundary_control_index, 1, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
// --- jface is control - END

      
	    if( jface_is_a_boundary_control /*jface_boundary_control_index== LIST_OF_CTRL_FACES :: _face_with_extremes_index[t] */  ) {

            int count_unbounded = 0;
            int count_bounded = 0;
//------------------------------------        
//------------ jface opening - END ---------    
//------------------------------------        

            

// // // //----  qp_of_jface, preparation right before iel - BEGIN ------- 

// Wait a second... If I want to prepare the qp_of_jface loop, I must be inside jface as well...
// So now it seems to me that I have to do jel jface iel iface instead...
// Previously it was jel          iel         iqp           jqp
// Now, it has to be jel jface  - iel iface - qp_of_iface - qp_of_jface
// The two quadrature loops must be the innermost. In this way you exclude all non-needed volume elements and all non-needed faces, so that you minimize the number of inner ifs. You keep them outside as much as possible
// There will be a storage of qp_of_jface

            
      const unsigned n_qp_of_jface = ml_prob.GetQuadratureRuleMultiple(qrule_j, jelGeom_bdry).GetGaussPointsNumber();

      std::vector < std::vector < double > > x_qp_of_jface(n_qp_of_jface);
      std::vector < double > weight_qp_of_jface(n_qp_of_jface);
      std::vector < std::vector < double > > phi_coords_jel_bdry_qp_of_jface(n_qp_of_jface);
      std::vector < std::vector < double > > phi_coords_x_jel_bdry_qp_of_jface(n_qp_of_jface);
      
      std::vector < std::vector < double > > phi_ctrl_jel_bdry_qp_of_jface(n_qp_of_jface);   /// @todo assume all ctrl components have the same FE family
      std::vector < std::vector < double > > phi_ctrl_x_jel_bdry_qp_of_jface(n_qp_of_jface);   /// @todo assume all ctrl components have the same FE family

      std::vector< std::vector< double > > sol_ctrl_qp_of_jface(n_components_ctrl);
      
        for (unsigned c = 0; c < n_components_ctrl; c++) {
             sol_ctrl_qp_of_jface[c].resize(n_qp_of_jface);
                 std::fill(sol_ctrl_qp_of_jface[c].begin(), sol_ctrl_qp_of_jface[c].end(), 0.);
            }

      
         for(unsigned qp_of_jface = 0; qp_of_jface < n_qp_of_jface; qp_of_jface++) {

    elem_all[qrule_j][jelGeom_bdry][solType_coords]->JacJacInv(geom_element_jel.get_coords_at_dofs_bdry_3d(), qp_of_jface, Jac_jel_bdry_qp_of_jface, JacI_jel_bdry_qp_of_jface, detJac_jel_bdry_qp_of_jface, space_dim);
    
    weight_qp_of_jface[qp_of_jface] = detJac_jel_bdry_qp_of_jface * ml_prob.GetQuadratureRuleMultiple(qrule_j, jelGeom_bdry).GetGaussWeightsPointer()[qp_of_jface];

    elem_all[qrule_j][jelGeom_bdry][SolFEType_quantities[pos_sol_ctrl]] ->shape_funcs_current_elem(qp_of_jface, JacI_jel_bdry_qp_of_jface, phi_ctrl_jel_bdry_qp_of_jface[qp_of_jface], phi_ctrl_x_jel_bdry_qp_of_jface[qp_of_jface], boost::none, space_dim);
            
    elem_all[qrule_j][jelGeom_bdry][solType_coords] ->shape_funcs_current_elem(qp_of_jface, JacI_jel_bdry_qp_of_jface, phi_coords_jel_bdry_qp_of_jface[qp_of_jface], phi_coords_x_jel_bdry_qp_of_jface[qp_of_jface], boost::none, space_dim);

//========== compute gauss quantities on the boundary - BEGIN ===============================================
//--- geom - BEGIN
           x_qp_of_jface[qp_of_jface].assign(dim, 0.);
          
//        if(jproc == iproc) {
         for(unsigned d = 0; d < dim; d++) {
	      for (int j_bdry = 0; j_bdry < geom_element_jel.get_coords_at_dofs_bdry_3d()[d].size(); j_bdry++)  {
			
              x_qp_of_jface[qp_of_jface][d] += geom_element_jel.get_coords_at_dofs_bdry_3d()[d][j_bdry] * phi_coords_jel_bdry_qp_of_jface[qp_of_jface][j_bdry];

		      }
            }
//        }
//       MPI_Bcast(& x_qp_of_jface[qp_of_jface][0], dim, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
//--- geom - END
    
//--- solution - BEGIN
      for (unsigned c = 0; c < n_components_ctrl; c++) {
         sol_ctrl_qp_of_jface[c][qp_of_jface] = 0.;
//        if(jproc == iproc) {
	      for (int j_bdry = 0; j_bdry < phi_ctrl_jel_bdry_qp_of_jface[qp_of_jface].size()/*Sol_n_el_dofs_quantities[pos_sol_ctrl]*/; j_bdry++)  {
		    unsigned int j_vol = msh->GetLocalFaceVertexIndex_PassElemType(jel_geommm, jface, j_bdry);
			
			sol_ctrl_qp_of_jface[c][qp_of_jface] +=  /*sol_eldofs_Mat[pos_mat_ctrl]*/sol_ctrl_jel[c][j_vol] * phi_ctrl_jel_bdry_qp_of_jface[qp_of_jface][j_bdry];

	      }
       }
//        }
//       MPI_Bcast(& sol_ctrl_qp_of_jface[dim], 1, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
//--- solution - END
//========== compute gauss quantities on the boundary - END ================================================


        }  //qp_of_jface

// // //         //we can do the broadcast after the loop, faster
// // //         for(unsigned qp_of_jface = 0; qp_of_jface < n_qp_of_jface; qp_of_jface++) {
// // //             MPI_Bcast(& x_qp_of_jface[qp_of_jface][0], dim, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
// // //         }
// // //             MPI_Bcast(& sol_ctrl_qp_of_jface[0], n_qp_of_jface, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
  
        
// // // //----  qp_of_jface, preparation right before iel - END ------- 

    
   


              
// --- iel opening - BEGIN
       
       for(int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1) ; iel++) {
           
           
// --- geometry - BEGIN        
        geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);

        short unsigned ielGeom = geom_element_iel.geom_type();
              
        geom_element_iel.set_elem_center_3d(iel, solType_coords);
// --- geometry - END        

      
	// Perform face loop over elements that contain some control face
	if ( DOMAIN_CONTAINING_CTRL_FACES ::volume_elem_contains_a_Gamma_control_face(sol, msh, iel ) ) {

// --- iel opening - END

//************ set control flag - BEGIN *********************
        std::vector< std::vector< int > > control_node_flag_iel_all_faces =
        femus::is_dof_associated_to_Gamma_control_equation(msh /*ok*/, ml_sol /*ok*/, & ml_prob /*ok*/, iel /*ok*/, ielGeom, geom_element_iel /*ok*/, solType_coords /*ok*/, Solname_Mat /*ok*/, SolFEType_Mat, Sol_n_el_dofs_Mat, pos_mat_ctrl, n_components_ctrl);
//************ set control flag - END *********************

        

        
// ***************************************
// ******* iel-related stuff - BEGIN *************
// ***************************************
        
      
// --- l2GMap - BEGIN
        std::vector < unsigned > nDof_iel(n_components_ctrl);
        
       for (unsigned c = 0; c < n_components_ctrl; c++) {
         nDof_iel[c]  = msh->GetElementDofNumber(iel, solType_ctrl[c]);
       }
              
      for (unsigned c = 0; c < n_components_ctrl; c++) {
        l2gMap_iel[c].resize(nDof_iel[c]);
        for(unsigned i = 0; i < l2gMap_iel[c].size(); i++) {
          l2gMap_iel[c][i] = pdeSys->GetSystemDof(solIndex_ctrl[c], solPdeIndex_ctrl[c], i, iel);
           }
         }
// --- l2GMap - END 

// --- l2GMapVec - BEGIN -----------------
    l2gMap_iel_vec.resize(0);
      for (unsigned c = 0; c < n_components_ctrl; c++) {
          l2gMap_iel_vec.insert(l2gMap_iel_vec.end(), l2gMap_iel[c].begin(), l2gMap_iel[c].end());
      }
// --- l2GMapVec - END -----------------


// --- element matrix and vector resizing - BEGIN
unsigned nDof_iel_vec = 0;
        for (unsigned c = 0; c < n_components_ctrl; c++) {  nDof_iel_vec  +=  nDof_iel[c];    }


        KK_nonlocal_iel_iel_bounded_integral.assign(nDof_iel_vec * nDof_iel_vec, 0.);   //resize
        KK_nonlocal_iel_jel_bounded_integral.assign(nDof_iel_vec * nDof_jel_vec, 0.);   //resize
        KK_nonlocal_jel_iel_bounded_integral.assign(nDof_jel_vec * nDof_iel_vec, 0.);   //resize
        KK_nonlocal_jel_jel_bounded_integral.assign(nDof_jel_vec * nDof_jel_vec, 0.);   //resize
        Res_nonlocal_iel_bounded_integral.assign(nDof_iel_vec, 0.);    //resize
        Res_nonlocal_jel_bounded_integral.assign(nDof_jel_vec, 0.);    //resize

        Res_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref.assign(nDof_iel_vec, 0.);    //resize
        KK_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref.assign(nDof_iel_vec * nDof_iel_vec, 0.);

        if( check_if_same_elem(TEST_JEL_SINGLE_FOR_LOOP /*iel*/, jel) ) {
            
          Res_local_iel_integer_operators.assign(nDof_iel_vec, 0.);    //resize
          KK_local_iel_integer_operators.assign(nDof_iel_vec * nDof_iel_vec, 0.);
          
          Res_local_iel_unbounded_integral_analytical_both_ref_and_non_ref.assign(nDof_iel_vec, 0.);    //resize
          KK_local_iel_unbounded_integral_analytical_both_ref_and_non_ref.assign(nDof_iel_vec * nDof_iel_vec, 0.);
          
          if(integration_num_split != 0) {
//             Vectors and matrices for adaptive quadrature
            Res_local_iel_only_refined.assign(nDof_iel_vec, 0.);    //resize
            KK_local_iel_only_refined.assign(nDof_iel_vec * nDof_iel_vec, 0.);
          }
          
        }
// --- element matrix and vector resizing - END 


// --- solution - BEGIN        

        for (unsigned c = 0; c < n_components_ctrl; c++) {
           sol_ctrl_iel[c].resize(nDof_iel[c]);

           
        for(unsigned i = 0; i < sol_ctrl_iel[c].size(); i++) {
          unsigned iDof  = msh->GetSolutionDof(i, iel, solType_ctrl[c]);  // global to global mapping between coordinates node and coordinate dof
          sol_ctrl_iel[c][i] = (*sol->_Sol[solIndex_ctrl[c]])(iDof);  // global extraction and local storage for the element coordinates
          }
        }
// --- solution - END        

// ***************************************
// ******* iel-related stuff - END *************
// ***************************************
        


        
//------------ iface opening - BEGIN  ---------        
	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {

// //             const unsigned int number_of_face_for_controls = LIST_OF_CTRL_FACES :: _face_with_extremes_index_size; // get the number of Gamma_c
// //            for(int w = 0; w < number_of_face_for_controls; w++) {


// --- geom - BEGIN
       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, iface);    
       
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();
// --- geom - END           

       int iface_is_a_boundary_control;
       unsigned int iface_boundary_control_index;

           std::pair< int, unsigned int > pair_control_iface = femus::face_is_a_Gamma_control_face_of_some_index< LIST_OF_CTRL_FACES >(msh->GetMeshElements(), iel, iface);
           
           iface_is_a_boundary_control  = pair_control_iface.first;
           iface_boundary_control_index = pair_control_iface.second;

	    if( iface_is_a_boundary_control && (iface_boundary_control_index == jface_boundary_control_index) ) {  //the integration is along the same Boundary face, not cross integration along different boundary faces
//------------ iface opening - END ---------        
		
//                 count_visits_of_boundary_faces++;


             
     //**** Adaptive preparation - BEGIN ******** 
     std::vector < std::vector < double > >  Jac_kel_bdry_kqp_bdry;
     std::vector < std::vector < double > >  JacI_kel_bdry_kqp_bdry;
     double detJac_kel_bdry_kqp_bdry;
      std::vector < double >  x_kqp_bdry;
      double  weight_kqp_bdry;
      std::vector < double >  phi_coords_kel_bdry_kqp_bdry;
      std::vector < double >  phi_coords_x_kel_bdry_kqp_bdry;
      
      std::vector < double >  phi_ctrl_kel_bdry_kqp_bdry; /// @todo assume all ctrl components have the same FE family
      std::vector < double >  phi_ctrl_x_kel_bdry_kqp_bdry;
      
       std::vector< double > sol_ctrl_kqp_bdry(n_components_ctrl);
            
        for (unsigned c = 0; c < n_components_ctrl; c++) {   sol_ctrl_kqp_bdry[c] = 0.;   }


      
//         double weight3;
//         std::vector < double > phi3;

//        Evaluating coarse FE functions on Quadrature Points of the "sub-elements"
       std::vector < std::vector < std::vector <double > > > aP(3);  //[NFE_FAMS][DIM==3][N_DOFS]
        if(integration_num_split != 0) {
          for(unsigned fe_type = 0; fe_type < /*solType*/solType_coords + 1; fe_type++) { //loop up to the FE type + 1 of the unknown
            ProjectNodalToPolynomialCoefficients(aP[fe_type], geom_element_iel.get_coords_at_dofs_bdry_3d(), ielGeom_bdry, fe_type) ;         ///@todo check this!!!!!!!!!!!!
          }
        }                      
     //**** Adaptive preparation - END ********  
                      
              //qp_of_iface initialization - BEGIN
        const unsigned n_qp_of_iface = ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussPointsNumber();
        
        
       std::vector< double > sol_ctrl_qp_of_iface(n_components_ctrl);
        for (unsigned c = 0; c < n_components_ctrl; c++) {   sol_ctrl_qp_of_iface[c] = 0.;   }
              //qp_of_iface initialization - END
    
         
//------------ qp_of_iface opening - BEGIN  ---------        
		for(unsigned qp_of_iface = 0; qp_of_iface < n_qp_of_iface; qp_of_iface++) {
//------------ qp_of_iface opening - END  ---------        

            
//========== qp_of_iface FE shape - BEGIN ===============================================
    elem_all[qrule_i][ielGeom_bdry][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_bdry_3d(), qp_of_iface, Jac_iel_bdry_qp_of_iface, JacI_iel_bdry_qp_of_iface, detJac_iel_bdry_qp_of_iface, space_dim);
    
    weight_qp_of_iface = detJac_iel_bdry_qp_of_iface * ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussWeightsPointer()[qp_of_iface];

    elem_all[qrule_i][ielGeom_bdry][solType_coords] ->shape_funcs_current_elem(qp_of_iface, JacI_iel_bdry_qp_of_iface, phi_coords_iel_bdry_qp_of_iface, phi_coords_x_iel_bdry_qp_of_iface, boost::none, space_dim);

    elem_all[qrule_i][ielGeom_bdry][SolFEType_quantities[pos_sol_ctrl]] ->shape_funcs_current_elem(qp_of_iface, JacI_iel_bdry_qp_of_iface, phi_ctrl_iel_bdry_qp_of_iface, phi_ctrl_x_iel_bdry_qp_of_iface, boost::none, space_dim);
//========== qp_of_iface FE shape - END ===============================================
            
//========== qp_of_iface quantities - BEGIN ===============================================
//--- geom - BEGIN 
          std::vector < double > x_qp_of_iface(dim, 0.);  ///@todo is this dim or dim_bdry?

            for(unsigned d = 0; d < x_qp_of_iface.size(); d++) {
	      for (int i_bdry = 0; i_bdry < geom_element_iel.get_coords_at_dofs_bdry_3d()[d].size(); i_bdry++)  {
			
              x_qp_of_iface[d] += geom_element_iel.get_coords_at_dofs_bdry_3d()[d][i_bdry] * phi_coords_iel_bdry_qp_of_iface[i_bdry];

		      }
            }
//--- geom - END
    
//--- solution - BEGIN 
        for (unsigned c = 0; c < n_components_ctrl; c++) {
          sol_ctrl_qp_of_iface[c] = 0.;
	        for (int i_bdry = 0; i_bdry < phi_ctrl_iel_bdry_qp_of_iface.size(); i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
			
			sol_ctrl_qp_of_iface[c] +=  /*sol_eldofs_Mat[pos_mat_ctrl]*/sol_ctrl_iel[c][i_vol] * phi_ctrl_iel_bdry_qp_of_iface[i_bdry];

		      }
       }
//--- solution - END
//========== qp_of_iface quantities - END ================================================


  
      //============  Non-fractional assembly - BEGIN ==================
          
            if( check_if_same_elem_bdry(iel, jel, iface, jface) /*&& (t == w)*/ ) {  //in principle, "jel" and "jface" could be any values in their range,
                                                                                    //but "jel == iel" and "jface == iface" is the only value that always guarantees that we reach this point,
                                                                                    //because we are under additional constraints "jface_is_a_boundary_control" from above
              
                
       //============  Mass assembly (local operator) - BEGIN ==================
    for (unsigned c = 0; c < n_components_ctrl; c++) {
        
        integral += operator_L2 * alpha * weight_qp_of_iface * sol_ctrl_qp_of_iface[c] * sol_ctrl_qp_of_iface[c];
    }
  
    
    
    for (unsigned c = 0; c < n_components_ctrl; c++) {
          for(unsigned l_bdry = 0; l_bdry < phi_ctrl_iel_bdry_qp_of_iface.size(); l_bdry++) {
               		    unsigned int l_vol = msh->GetLocalFaceVertexIndex(iel, iface, l_bdry);

              double mass_res_i = phi_ctrl_iel_bdry_qp_of_iface[l_bdry] * sol_ctrl_qp_of_iface[c];
              const unsigned res_pos = assemble_jacobian<double,double>::res_row_index(nDof_iel, c, l_vol);
              Res_local_iel_integer_operators[ res_pos/*l_vol*/ ] += (-1.) *
              control_node_flag_iel_all_faces[c][l_vol]  * operator_L2 * alpha * weight_qp_of_iface * mass_res_i ;
              Res_local_iel_integer_operators[ res_pos/*l_vol*/ ] += (-1.) *
              control_node_flag_iel_all_faces[c][l_vol]  * rhs_one * weight_qp_of_iface * (phi_ctrl_iel_bdry_qp_of_iface[l_bdry] * (-1.) /** ( sin(2 * acos(0.0) * x1[0][l_bdry])) * ( sin(2 * acos(0.0) * x1[1][l_bdry]))*/);
              
              for (unsigned e = 0; e < n_components_ctrl; e++) {
                  if (e == c) { 
            for(unsigned m_bdry = 0; m_bdry < phi_ctrl_iel_bdry_qp_of_iface.size(); m_bdry++) {
               		    unsigned int m_vol = msh->GetLocalFaceVertexIndex(iel, iface, m_bdry);
                const unsigned jac_pos = assemble_jacobian< double, double >::jac_row_col_index(nDof_iel, nDof_iel_vec, c, e, l_vol, m_vol);         
                KK_local_iel_integer_operators[ jac_pos/*l_vol * nDof_iel + m_vol*/ ] +=
                control_node_flag_iel_all_faces[c][l_vol] *
                operator_L2 * alpha * phi_ctrl_iel_bdry_qp_of_iface[l_bdry] * phi_ctrl_iel_bdry_qp_of_iface[m_bdry] * weight_qp_of_iface;
                 }
               }
             }
             
            } 

       }
       //============  Mass assembly (local operator) - END ==================
         
      //============  Laplacian assembly - BEGIN ==================
      //============  Laplacian assembly - END ==================
            }
    
      //============  Non-fractional assembly - END ==================
             
             
      //============  Fractional assembly - BEGIN ==================
        if(operator_Hhalf != 0) {

      //============ Or different elements, or lack of adaptivity (so all elements) - BEGIN ==================
        if( !( check_if_same_elem_bdry(iel, jel, iface, jface) && integration_num_split != 0 ) ) {  //  if(iel != jel || integration_num_split == 0)

// // //            if ( iface_boundary_control_index == jface_boundary_control_index )  {
// ********* BOUNDED PART - BEGIN ***************

//------------ qp_of_jface opening - BEGIN  ---------
            for(unsigned qp_of_jface = 0; qp_of_jface < n_qp_of_jface; qp_of_jface++) {
//------------ qp_of_jface opening - END  ---------

              double dist_xyz = 0.;
              for(unsigned d = 0; d < x_qp_of_iface.size(); d++) {
                dist_xyz += (x_qp_of_iface[d] - x_qp_of_jface[qp_of_jface][d]) * (x_qp_of_iface[d] - x_qp_of_jface[qp_of_jface][d]);
              }

              const double denom = pow(dist_xyz, (double)(  0.5 * /*dim*/dim_bdry + s_frac));

              const double common_weight = (0.5 * C_ns) * operator_Hhalf * beta * check_limits * weight_qp_of_iface * weight_qp_of_jface[qp_of_jface]  / denom;

              for (unsigned c = 0; c < n_components_ctrl; c++) {

                  integral +=  common_weight * (sol_ctrl_qp_of_iface[c] - sol_ctrl_qp_of_jface[c][qp_of_jface]) * sol_ctrl_qp_of_iface[c];
                  integral +=  common_weight * (sol_ctrl_qp_of_iface[c] - sol_ctrl_qp_of_jface[c][qp_of_jface]) * (- sol_ctrl_qp_of_jface[c][qp_of_jface]);

              }


              for (unsigned c = 0; c < n_components_ctrl; c++) {
                 for(unsigned l_bdry = 0; l_bdry < phi_ctrl_iel_bdry_qp_of_iface.size(); l_bdry++) { //dofs of test function
               		    unsigned int l_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, l_bdry);
               		    unsigned int l_vol_jel = msh->GetMeshElements()->GetIG(jel_geommm, jface, l_bdry)/*msh->GetLocalFaceVertexIndex(jel, jface, l_bdry)*/;

              const unsigned res_pos_iel = assemble_jacobian<double,double>::res_row_index(nDof_iel, c, l_vol_iel);
              const unsigned res_pos_jel = assemble_jacobian<double,double>::res_row_index(nDof_jel, c, l_vol_jel);
                Res_nonlocal_iel_bounded_integral[ res_pos_iel /*l_vol_iel*/ ]      +=      -
                 control_node_flag_iel_all_faces[c][l_vol_iel] * control_node_flag_jel_all_faces[c][l_vol_jel] *
                common_weight * (sol_ctrl_qp_of_iface[c] - sol_ctrl_qp_of_jface[c][qp_of_jface]) * (phi_ctrl_iel_bdry_qp_of_iface[l_bdry]);

                Res_nonlocal_jel_bounded_integral[ res_pos_jel /*l_vol_jel*/ ]      +=      -
                 control_node_flag_iel_all_faces[c][l_vol_iel] * control_node_flag_jel_all_faces[c][l_vol_jel]*
                common_weight * (sol_ctrl_qp_of_iface[c] - sol_ctrl_qp_of_jface[c][qp_of_jface]) * (- phi_ctrl_jel_bdry_qp_of_jface[qp_of_jface][l_bdry]);


//                 for(unsigned j = 0; j < nDof_jel; j++) {
                for (unsigned e = 0; e < n_components_ctrl; e++) {
                  if (e == c) {
                      for(unsigned m_bdry = 0; m_bdry < phi_ctrl_jel_bdry_qp_of_jface[qp_of_jface].size(); m_bdry++) { //dofs of unknown function
               		    unsigned int m_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, m_bdry);
               		    unsigned int m_vol_jel = msh->GetMeshElements()->GetIG(jel_geommm, jface, m_bdry)/*msh->GetLocalFaceVertexIndex(jel, jface, m_bdry)*/;

                    const unsigned jac_pos_iel_iel = assemble_jacobian< double, double >::jac_row_col_index(nDof_iel, nDof_iel_vec, c, e, l_vol_iel, m_vol_iel);
                    const unsigned jac_pos_iel_jel = assemble_jacobian< double, double >::jac_row_col_index(nDof_iel, nDof_iel_vec, c, e, l_vol_iel, m_vol_jel);
                    const unsigned jac_pos_jel_iel = assemble_jacobian< double, double >::jac_row_col_index(nDof_iel, nDof_iel_vec, c, e, l_vol_jel, m_vol_iel);
                    const unsigned jac_pos_jel_jel = assemble_jacobian< double, double >::jac_row_col_index(nDof_iel, nDof_iel_vec, c, e, l_vol_jel, m_vol_jel);

             /*  u(x) v(x)*/     KK_nonlocal_iel_iel_bounded_integral[ jac_pos_iel_iel /*l_vol_iel * nDof_jel + m_vol_iel*/ ] +=
              control_node_flag_iel_all_faces[c][l_vol_iel] * control_node_flag_jel_all_faces[c][l_vol_jel] *
             common_weight *          phi_ctrl_iel_bdry_qp_of_iface[m_bdry]            *    phi_ctrl_iel_bdry_qp_of_iface[l_bdry];

             /*- u(y) v(x)*/     KK_nonlocal_iel_jel_bounded_integral[ jac_pos_iel_jel /*l_vol_iel * nDof_jel + m_vol_jel*/ ] +=
              control_node_flag_iel_all_faces[c][l_vol_iel] * control_node_flag_jel_all_faces[c][l_vol_jel] *
             common_weight * (- 1.) * phi_ctrl_jel_bdry_qp_of_jface[qp_of_jface][m_bdry]  *    phi_ctrl_iel_bdry_qp_of_iface[l_bdry];

             /*- u(x) v(y)*/     KK_nonlocal_jel_iel_bounded_integral[ jac_pos_jel_iel /*l_vol_jel * nDof_jel + m_vol_iel*/ ] +=
              control_node_flag_iel_all_faces[c][l_vol_iel] * control_node_flag_jel_all_faces[c][l_vol_jel] *
             common_weight * (- 1.) * phi_ctrl_iel_bdry_qp_of_iface[m_bdry]            *   phi_ctrl_jel_bdry_qp_of_jface[qp_of_jface][l_bdry];

             /*  u(y) v(y)*/     KK_nonlocal_jel_jel_bounded_integral[ jac_pos_jel_jel /*l_vol_jel * nDof_jel + m_vol_jel*/ ] +=
              control_node_flag_iel_all_faces[c][l_vol_iel] * control_node_flag_jel_all_faces[c][l_vol_jel] *
             common_weight *          phi_ctrl_jel_bdry_qp_of_jface[qp_of_jface][m_bdry]  *    phi_ctrl_jel_bdry_qp_of_jface[qp_of_jface][l_bdry];


                     }
                   }
                 }
                }

               }

//------------ qp_of_jface closing - BEGIN  ---------
              } //endl qp_of_jface loop
//------------ qp_of_jface closing - END ---------

//  count_bounded++;
// ********* BOUNDED PART - END ***************

// // //             }
//                 if ( iface_boundary_control_index == jface_boundary_control_index  /*t == w*/)  {
// ********* UNBOUNDED PART - BEGIN ***************

               unbounded_integral_over_exterior_of_boundary_control_face(
//
//                               count_unbounded,
//
                              unbounded,
                              dim,
                              dim_bdry,
//////////
                              operator_Hhalf,
                              s_frac,
                              check_limits,
                              C_ns,
                              beta,
//////////
                              msh,
                              sol,
                              ml_sol,
//////////
                              iel,
                              iface,
                              iface_boundary_control_index,
                              control_node_flag_iel_all_faces,
                              control_node_flag_jel_all_faces,
                              nDof_iel,
                              weight_qp_of_iface,
                              x_qp_of_iface,
                              phi_ctrl_iel_bdry_qp_of_iface,
                              sol_ctrl_qp_of_iface,
                              KK_local_iel_unbounded_integral_analytical_both_ref_and_non_ref,
                              Res_local_iel_unbounded_integral_analytical_both_ref_and_non_ref,
//////////
  //////////   3D only - BEGIN
                              jel,
                              jface,
                              jface_boundary_control_index,
                              geom_element_jel,
                              node_based_bdry_bdry_in,
                              N_div_face_of_face,
                              KK_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref,
                              Res_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref,
                              solType_coords,
  //////////   3D only - END
                              integral
                             );


// ********* UNBOUNDED PART - END ***************
//                    }

         } //end if(iel != jel || integration_num_split == 0)
      //============ Or different elements, or lack of adaptivity (so all elements) - END ==================


      //============ Either Same element && Adaptive quadrature - BEGIN ==================
       else  {
                
          /*const*/ short unsigned kelGeom_bdry = ielGeom_bdry;
           
          const unsigned n_kqp_bdry = ml_prob.GetQuadratureRuleMultiple(qrule_k, kelGeom_bdry).GetGaussPointsNumber();
                
                
          std::vector< std::vector< std::vector<double> > > x3;  // Number of element subdivisions x Dimension x Number of dofs

          for(unsigned split = 0; split <= integration_num_split; split++) {

            
            if (dim_bdry/*dim*/ == 1) GetElementPartition1D(x_qp_of_iface, geom_element_iel.get_coords_at_dofs_bdry_3d(), split, integration_num_split, x3, space_dim); //TODO space_dim or dim?
            else if (dim_bdry/*dim*/ == 2) {
              //TODO need to be checked !!!
              //GetElementPartition2D(x_qp_of_iface, geom_element_iel.get_coords_at_dofs_bdry_3d(), split, integration_num_split, x3);
              GetElementPartitionQuad(x_qp_of_iface, geom_element_iel.get_coords_at_dofs_bdry_3d(), split, integration_num_split, x3);
            }

            //for(unsigned r = 0; r < size_part; r++) {
            for(unsigned r = 0; r < x3.size(); r++) {


              for(unsigned k_qp_bdry = 0; k_qp_bdry < n_kqp_bdry; k_qp_bdry++) {

// ********* PREPARATION PART - BEGIN ***************
                elem_all[qrule_k][kelGeom_bdry][solType_coords]->JacJacInv(x3[r]/*geom_element_iel.get_coords_at_dofs_bdry_3d()*/, k_qp_bdry, Jac_kel_bdry_kqp_bdry, JacI_kel_bdry_kqp_bdry, detJac_kel_bdry_kqp_bdry, space_dim);
    
//                 weight_kqp_bdry = detJac_kel_bdry_kqp_bdry * ml_prob.GetQuadratureRuleMultiple(qrule_k, kelGeom_bdry).GetGaussWeightsPointer()[k_qp_bdry];
                
                msh->_finiteElement[kelGeom_bdry][solType_coords]->Jacobian(x3[r], k_qp_bdry, weight_kqp_bdry, phi_coords_kel_bdry_kqp_bdry, phi_x_placeholder_unused);
                
//                 elem_all[qrule_k][kelGeom_bdry][SolFEType_quantities[pos_sol_ctrl]] ->shape_funcs_current_elem(k_qp_bdry, JacI_kel_bdry_kqp_bdry, phi_ctrl_kel_bdry_kqp_bdry, phi_ctrl_x_kel_bdry_kqp_bdry, boost::none, space_dim);
            
//                 elem_all[qrule_k][kelGeom_bdry][solType_coords] ->shape_funcs_current_elem(k_qp_bdry, JacI_kel_bdry_kqp_bdry, phi_coords_kel_bdry_kqp_bdry, phi_coords_x_kel_bdry_kqp_bdry, boost::none, space_dim);

//--- geom - BEGIN
                std::vector < double > x_kqp_bdry(dim, 0.);  ///@todo is this dim or dim_bdry?

                for(unsigned d = 0; d < x_kqp_bdry.size(); d++) {
                  for (int k_bdry = 0; k_bdry < /*geom_element_iel.get_coords_at_dofs_bdry_3d()[d].size()*/phi_coords_kel_bdry_kqp_bdry.size(); k_bdry++)  {
			
                    x_kqp_bdry[d] += x3[r][d][k_bdry]/*geom_element_iel.get_coords_at_dofs_bdry_3d()[d][i_vol]*/ * phi_coords_kel_bdry_kqp_bdry[k_bdry];

                  }
                }
//--- geom - END
    
                std::vector<double> xi3(dim, 0.);

                GetClosestPointInReferenceElement(geom_element_iel.get_coords_at_dofs_bdry_3d(), x_kqp_bdry, kelGeom_bdry, xi3);
//                 GetInverseMapping(solType_coords, kelGeom_bdry, aP, x_kqp_bdry, xi3, 1000);  ///@todo generalize to rectangular Jacobian TODO is needed?

                msh->_finiteElement[kelGeom_bdry][/*solType*/ solType_coords]->GetPhi(phi_ctrl_kel_bdry_kqp_bdry, xi3); //TODO solType or solType_coords?

                 std::vector<double> solY3(n_components_ctrl, 0.);

                for (unsigned c = 0; c < n_components_ctrl; c++) {
     for(unsigned i_bdry = 0; i_bdry < phi_ctrl_kel_bdry_kqp_bdry.size(); i_bdry++) {
                  unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
                  solY3[c] += sol_ctrl_iel[c][i_vol] * phi_ctrl_kel_bdry_kqp_bdry[i_bdry];
                }
               }
// // // 
// // // 
// // //                     msh->_finiteElement[ielGeom][solType]->Jacobian(x3[r], k_qp_bdry, weight_kqp_bdry, phi3, phi_x_placeholder_unused);
// // // 
// // //                     std::vector < double > xg3(dim, 0.);
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
                for(unsigned k = 0; k < x_qp_of_iface.size(); k++) {
                  dist_xyz3 += (x_qp_of_iface[k] - x_kqp_bdry[k]) * (x_qp_of_iface[k] - x_kqp_bdry[k]);
                }

                const double denom_ik = pow(dist_xyz3, (double)( 0.5 * dim_bdry + s_frac));
                
                const double common_weight =  0.5 * C_ns * operator_Hhalf * beta * check_limits * weight_qp_of_iface * weight_kqp_bdry  / denom_ik;


     for (unsigned c = 0; c < n_components_ctrl; c++) {
         integral    +=       common_weight * (sol_ctrl_qp_of_iface[c] - solY3[c]) * (sol_ctrl_qp_of_iface[c] - solY3[c]);
     }
     
                for (unsigned c = 0; c < n_components_ctrl; c++) {
               for(unsigned l_bdry = 0; l_bdry < phi_ctrl_iel_bdry_qp_of_iface.size(); l_bdry++) { //dofs of test function
                  unsigned int l_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, l_bdry);

              const unsigned res_pos = assemble_jacobian<double,double>::res_row_index(nDof_iel, c, l_vol_iel);
                  Res_local_iel_only_refined[res_pos/* l_vol_iel*/ ]    += -
                   control_node_flag_iel_all_faces[c][l_vol_iel] *
                  common_weight * (sol_ctrl_qp_of_iface[c] - solY3[c]) * (phi_ctrl_iel_bdry_qp_of_iface[l_bdry] - phi_ctrl_kel_bdry_kqp_bdry[l_bdry]);

              for (unsigned e = 0; e < n_components_ctrl; e++) {
                  if (e == c) {
                for(unsigned m_bdry = 0; m_bdry < phi_ctrl_kel_bdry_kqp_bdry.size(); m_bdry++) { //dofs of unknown function
                    unsigned int m_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, m_bdry);
                    
                    const unsigned jac_pos = assemble_jacobian< double, double >::jac_row_col_index(nDof_iel, nDof_iel_vec, c, e, l_vol_iel, m_vol_iel);         
                KK_local_iel_only_refined[ jac_pos/*l_vol_iel * nDof_jel + m_vol_iel*/ ] +=
                 control_node_flag_iel_all_faces[c][l_vol_iel] *
                common_weight * (phi_ctrl_iel_bdry_qp_of_iface[m_bdry] - phi_ctrl_kel_bdry_kqp_bdry[m_bdry]) *
                                                        (phi_ctrl_iel_bdry_qp_of_iface[l_bdry] - phi_ctrl_kel_bdry_kqp_bdry[l_bdry]);

                  }
                }
                  
              }
            }
     }
//      count_bounded++;
// ********* BOUNDED PART - END ***************
// // //  if ( /*iface_boundary_control_index == jface_boundary_control_index */ t == w)  {
// ********* UNBOUNDED PART - BEGIN ***************

                if( qp_of_iface == integration_split_index ) { ///@todo is there a way to put this outside of the quadrature loop?

              unbounded_integral_over_exterior_of_boundary_control_face(
                  //
//                   count_unbounded,
                  //
                  unbounded,
                              dim,
                              dim_bdry,
//////////                             
                              operator_Hhalf,
                              s_frac,
                              check_limits,
                              C_ns,
                              beta,
//////////                             
                              msh,
                              sol,
                              ml_sol,
//////////                             
                              iel,
                              iface,
                              iface_boundary_control_index,
                              control_node_flag_iel_all_faces,
                              control_node_flag_jel_all_faces,
                              nDof_iel,
                              weight_kqp_bdry,
                              x_kqp_bdry,
                              phi_ctrl_kel_bdry_kqp_bdry,
                              sol_ctrl_kqp_bdry,
                              KK_local_iel_unbounded_integral_analytical_both_ref_and_non_ref,
                              Res_local_iel_unbounded_integral_analytical_both_ref_and_non_ref,
  //////////   3D only - BEGIN
                             jel,
                             jface,
                             jface_boundary_control_index,
                              geom_element_jel,
                             node_based_bdry_bdry_in,
                              N_div_face_of_face,
                              KK_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref,
                              Res_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref,
                              solType_coords,
  //////////   3D only - END
                              integral
                             ); 
                }


// ********* UNBOUNDED PART - END ***************
// // //  }
                  } //end k_qp_bdry
                } //end r
              }  //end split
              
              
                
                
                
            }  //end iel == jel && integration_num_split != 0
      //============ Either Same element && Adaptive quadrature - END ==================


         } //operator_Hhalf != 0              
      //============  Fractional assembly - END ==================
            
//------------ qp_of_iface closing - BEGIN  ---------        
      }   //end qp_of_iface
       count_bounded++;
//        count_unbounded++;

//------------ qp_of_iface closing - END ---------        
      
      

      
              
//----- iface closing - BEGIN ---        
        } //end face for control

// //            } // end w face
        
      } //end iface
//----- iface closing - END ---        


                
//============ add to global - BEGIN ==================
// // multiply everything by -1.? Don't think so - BEGIN 
// std::transform(KK_local_iel_unbounded_integral_analytical_both_ref_and_non_ref.begin(), KK_local_iel_unbounded_integral_analytical_both_ref_and_non_ref.end(), KK_local_iel_unbounded_integral_analytical_both_ref_and_non_ref.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(Res_local_iel_unbounded_integral_analytical_both_ref_and_non_ref.begin(), Res_local_iel_unbounded_integral_analytical_both_ref_and_non_ref.end(), Res_local_iel_unbounded_integral_analytical_both_ref_and_non_ref.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// 
// std::transform(KK_local_iel_only_refined.begin(), KK_local_iel_only_refined.end(), KK_local_iel_only_refined.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(Res_local_iel_only_refined.begin(), Res_local_iel_only_refined.end(), Res_local_iel_only_refined.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// 
// std::transform(KK_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref.begin(), KK_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref.end(), KK_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(Res_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref.begin(), Res_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref.end(), Res_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref_both_ref_and_non_ref.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// 
// std::transform(KK_nonlocal_iel_iel_bounded_integral.begin(), KK_nonlocal_iel_iel_bounded_integral.end(), KK_nonlocal_iel_iel_bounded_integral.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(KK_nonlocal_iel_jel_bounded_integral.begin(), KK_nonlocal_iel_jel_bounded_integral.end(), KK_nonlocal_iel_jel_bounded_integral.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(KK_nonlocal_jel_iel_bounded_integral.begin(), KK_nonlocal_jel_iel_bounded_integral.end(), KK_nonlocal_jel_iel_bounded_integral.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(KK_nonlocal_jel_jel_bounded_integral.begin(), KK_nonlocal_jel_jel_bounded_integral.end(), KK_nonlocal_jel_jel_bounded_integral.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// 
// std::transform(Res_nonlocal_iel_bounded_integral.begin(), Res_nonlocal_iel_bounded_integral.end(), Res_nonlocal_iel_bounded_integral.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(Res_nonlocal_jel_bounded_integral.begin(), Res_nonlocal_jel_bounded_integral.end(), Res_nonlocal_jel_bounded_integral.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// // multiply everything by -1. - END


 if( check_if_same_elem( TEST_JEL_SINGLE_FOR_LOOP /*iel*/, jel) /*check_if_same_elem_bdry(iel, jel, iface, jface)*/ /* if you put it inside the face loops */ ) {

        // RES & KK of half norm second integral in a 2D domain (double integral over Gamma_c and over R\Gamma_c) - BEGIN

         RES->add_vector_blocked(Res_local_iel_integer_operators, l2gMap_iel_vec);
          KK->add_matrix_blocked(KK_local_iel_integer_operators, l2gMap_iel_vec, l2gMap_iel_vec);
        

         RES->add_vector_blocked(Res_local_iel_unbounded_integral_analytical_both_ref_and_non_ref, l2gMap_iel_vec);
          KK->add_matrix_blocked(KK_local_iel_unbounded_integral_analytical_both_ref_and_non_ref, l2gMap_iel_vec, l2gMap_iel_vec);
        // RES & KK of half norm second integral in a 2D domain (double integral over Gamma_c and over R\Gamma_c) - END

        
          if(integration_num_split != 0) {
            RES->add_vector_blocked(Res_local_iel_only_refined, l2gMap_iel_vec);
            KK->add_matrix_blocked(KK_local_iel_only_refined, l2gMap_iel_vec, l2gMap_iel_vec);
          }

          count_unbounded++;
          
       }

//         count_unbounded--;

        // RES & KK of half norm second integral in a 3D domain (double integral over Gamma_c and over R\Gamma_c) - BEGIN
        RES->add_vector_blocked(Res_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref, l2gMap_iel_vec);
        KK->add_matrix_blocked(KK_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref, l2gMap_iel_vec, l2gMap_iel_vec);
        // RES & KK of half norm second integral in a 3D domain (double integral over Gamma_c and over R\Gamma_c)- END

        // RES & KK of half norm first integral (double integral over Gamma_c) - BEGIN
        RES->add_vector_blocked(Res_nonlocal_iel_bounded_integral, l2gMap_iel_vec);
        RES->add_vector_blocked(Res_nonlocal_jel_bounded_integral, l2gMap_jel_vec);
        
        KK->add_matrix_blocked(KK_nonlocal_iel_iel_bounded_integral, l2gMap_iel_vec, l2gMap_iel_vec);
        KK->add_matrix_blocked(KK_nonlocal_iel_jel_bounded_integral, l2gMap_iel_vec, l2gMap_jel_vec);
        KK->add_matrix_blocked(KK_nonlocal_jel_iel_bounded_integral, l2gMap_jel_vec, l2gMap_iel_vec);
        KK->add_matrix_blocked(KK_nonlocal_jel_jel_bounded_integral, l2gMap_jel_vec, l2gMap_jel_vec);
        // RES & KK of half norm first integral (double integral over Gamma_c) - END

// Since 1 is dense and 3 are sparse, and the dense dofs are 30, we should have at most 3x9 + 30 = 57, but in the sparsity print it shows 30. That's the problem


        
          
// //              if (print_algebra_local) {
//          std::vector<unsigned> Sol_n_el_dofs_Mat_vol2(1, nDof_jel);
// //          assemble_jacobian<double,double>::print_element_residual(iel, Res, Sol_n_el_dofs_Mat_vol, 10, 5);
//          assemble_jacobian<double,double>::print_element_jacobian(iel, KK_nonlocal_iel_unbounded_integral_numerical_both_ref_and_non_ref, Sol_n_el_dofs_Mat_vol2, 10, 5);
// //      }
         
//============ add to global - END ==================


        
//----- iel closing - BEGIN ---        
          } //end control elem flag i (control_flag_iel == 1)
      } //end iel
//----- iel closing - END ---        


//----- jface closing - BEGIN ---        
     
        } //end face for control
// // // // //            } //end t face
      } //end jface

//----- jface closing - END ---   
 
     
//----- jel closing - BEGIN ---        
    }  //end control elem flag jel
   } //end jel
//----- jel closing - END ---        
   
 } //end jproc
    
    
// integral - BEGIN ************
  std::cout << std::endl;
  std::cout <<  "&&&&&& Integrals from previous iteration &&&&&&&&&&" << std::endl;
  double integral_parallel = 0.; MPI_Allreduce( &integral, &integral_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  std::cout << "The value of the integral  is " << std::setw(11) << std::setprecision(10) << 0.5 * integral_parallel << std::endl;
  std::cout <<  "&&&&&& Integrals from previous iteration &&&&&&&&&&" << std::endl;
  std::cout << std::endl;
                  
       std::cout << "SQUARE ROOOOOOOOTTTTTTTTTTTTTTTTTTTT AFTER" << std::endl;
       
// integral - END ************

    
    
    
    
    
  
  }
    
//**********FRAC CONTROL - END *****************************************
  

//*********************** Mesh independent, ALMOST - END *****************************************
    
};

 
 
} //end namespace ctrl




  
  
  
} //end namespace femus
 

 
//*******************************************************************************************
//*********************** Mesh independent - END *****************************************
//*******************************************************************************************




#endif
