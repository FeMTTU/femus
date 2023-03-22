#ifndef __opt_systems_boundary_control_eqn_sobolev_fractional_hpp__
#define __opt_systems_boundary_control_eqn_sobolev_fractional_hpp__



#include "Assemble_useful_functions.hpp"
#include "opt_common.hpp"

//for reading additional fields from MED file (based on MED ordering)
#include "MED_IO.hpp"

#include "fractional_functions.hpp"




//*******************************************************************************************
//*********************** Mesh independent - BEGIN *****************************************
//*******************************************************************************************

namespace femus {

namespace ctrl {




template < class LIST_OF_CTRL_FACES, class DOMAIN_CONTAINING_CTRL_FACES >
 class Gamma_control_equation_fractional_sobolev_differentiability_index {
//*********************** Mesh independent, ALMOST - BEGIN *****************************************
 
 public:
  
  static void mixed_integral(const unsigned unbounded,
                      const unsigned dim,
                      const unsigned dim_bdry,
//                       const std::vector < std::vector < double > > & ex_control,
                      const unsigned div, 
                      const double weight_iqp_bdry,
                      const std::vector < double > & x_iqp_bdry,
                      const std::vector < double > & phi_ctrl_iel_bdry_iqp_bdry,
                      const std::vector<double> sol_ctrl_iqp_bdry,   //////////
                      const double s_frac,
                      const double check_limits,
                      const double C_ns,
                      const unsigned int operator_Hhalf,
                      const double beta,
                      const std::vector<unsigned int> nDof_vol_iel,  //////////
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
                      std::string node_based_bdry_bdry_in,
                      std::vector <int> bdry_bdry,
                      CurrentElem < double > & geom_element_jel,
                      const unsigned jelGeom_bdry,
                      unsigned solType_coords,
                      double & integral
                     ) {
      
     
  const unsigned  sol_node_flag_index =  ml_sol->GetIndex( node_based_bdry_bdry_in.c_str() );
  const unsigned  group_salome = 2;   ///@todo fix here, maybe pass it in the args
  
  const unsigned int n_components_ctrl = nDof_vol_iel.size();
  
  unsigned nDof_iel_vec = 0;
        for (unsigned c = 0; c < n_components_ctrl; c++) {nDof_iel_vec  +=  nDof_vol_iel[c];    }

  
      
  if(unbounded == 1) {
      
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



 //------------ my change BEGIN ---------
  const unsigned int i_element_face_index = - ( msh->el->GetFaceElementIndex(iel, iface) + 1);    //get the face index
  const unsigned int number_of_face_for_controls = LIST_OF_CTRL_FACES :: _face_with_extremes_index_size; // get the number of Gamma_c
  unsigned int n_max = 2 /*pow(2,dim_bdry)*/;

  std::vector < double > extremes(n_max);

  std::vector < std::vector < double > > ex_control(dim); // control zone extremes. Built as: 2D) ( x1 y1 ) ( x2 y2 )    3D) ( x1 y1 z1 ) ( x2 y2 z2 ) ( x3 y3 z3 )
  for(unsigned d = 0; d < dim; d++) {
      ex_control[d].reserve(n_max);
  }

  for(int t = 0; t < number_of_face_for_controls; t++){   //loop over element of the class 'list of face for control' for found the t-index of the current Gamma_c
     if( LIST_OF_CTRL_FACES :: _face_with_extremes_index[t] == i_element_face_index){

         for(int k = 0; k < n_max; k++){
             extremes[k] = LIST_OF_CTRL_FACES :: _face_with_extremes_extremes[t][k] /*EX_1*/;
         }

         int control_xyz = (LIST_OF_CTRL_FACES :: _face_with_extremes_index[t] - 1) / 2;
         bool ctrl_min_max = (LIST_OF_CTRL_FACES :: _face_with_extremes_index[t] - 1) % 2;
         for(unsigned d = 0; d < dim; d++) {
                        for(unsigned n_e = 0; n_e < n_max; n_e++){
                          if(control_xyz == d) ex_control[d][n_e] =  (ctrl_min_max)? DOMAIN_EX_2:DOMAIN_EX_1;
                          else ex_control[d][n_e] = extremes[n_e];
                        }
                      }
                }
    }
  //------------ my change - END---------


//--- Control domain - BEGIN -------
//***************************************************

//   unsigned n_max = pow(2,dim_bdry);
//   std::vector < double > extremes(n_max);
//   extremes[0] = LIST_OF_CTRL_FACES :: _face_with_extremes_extremes[0][0] /*EX_1*/;
//   extremes[1] = LIST_OF_CTRL_FACES :: _face_with_extremes_extremes[0][1] /*EX_2*/;
// //   if(dim_bdry == 2){
// //     extremes[2] = EY_1;
// //     extremes[3] = EY_2;
// //   }

//   std::vector < std::vector < double > > ex_control(dim); // control zone extremes. Built as: 2D) ( x1 y1 ) ( x2 y2 )    3D) ( x1 y1 z1 ) ( x2 y2 z2 ) ( x3 y3 z3 )
//   for(unsigned d = 0; d < dim; d++) {
//       ex_control[d].reserve(n_max);
//   }
//   int control_xyz   = ( LIST_OF_CTRL_FACES :: _face_with_extremes_index[0] /*FACE_FOR_CONTROL*/ - 1) / 2;
//   bool ctrl_min_max = ( LIST_OF_CTRL_FACES :: _face_with_extremes_index[0] /*FACE_FOR_CONTROL*/ - 1) % 2;
//   for(unsigned d = 0; d < dim; d++) {
//     for(unsigned n_e = 0; n_e < n_max; n_e++){
//       if(control_xyz == d) ex_control[d][n_e] =  (ctrl_min_max)? DOMAIN_EX_2:DOMAIN_EX_1;
//       else ex_control[d][n_e] = extremes[n_e];
//     }
//   }
//***************************************************


//--- Control domain - END -------

              
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
              
              double mixed_denominator = pow(dist_1, -2. * s_frac) + pow(dist_2, - 2. * s_frac);

              
              for (unsigned c = 0; c < n_components_ctrl; c++) {
                  integral +=  0.5 * C_ns * check_limits * operator_Hhalf  * beta * sol_ctrl_iqp_bdry[c] * sol_ctrl_iqp_bdry[c] * weight_iqp_bdry * mixed_denominator * (1. / s_frac);
              }   
              
              
              
              for (unsigned c = 0; c < n_components_ctrl; c++) {

              for(unsigned i_bdry = 0; i_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); i_bdry++) {
                unsigned int i_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
                
                const unsigned res_pos = assemble_jacobian<double,double>::res_row_index(nDof_vol_iel, c, i_vol_iel);

                Res_local_iel[ res_pos/*i_vol_iel*/ ] += - 0.5 * C_ns * check_limits * operator_Hhalf  * beta * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * sol_ctrl_iqp_bdry[c] * weight_iqp_bdry * mixed_denominator * (1. / s_frac);
                
                for (unsigned e = 0; e < n_components_ctrl; e++) {
                  if (e == c) { 
                     for(unsigned j_bdry = 0; j_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); j_bdry++) {
                  unsigned int j_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, j_bdry);
                  const unsigned jac_pos = assemble_jacobian< double, double >::jac_row_col_index(nDof_vol_iel, nDof_iel_vec, c, e, i_vol_iel, j_vol_iel);         

                  KK_local_iel[ jac_pos /*i_vol_iel * nDof_vol_iel + j_vol_iel*/ ] += 0.5 * C_ns * check_limits * operator_Hhalf * beta * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * phi_ctrl_iel_bdry_iqp_bdry[j_bdry] * weight_iqp_bdry * mixed_denominator * (1. / s_frac);
                     }
                
                  }
                }
                
              }
           }
      }
      
    //============ Mixed Integral 2D - Numerical ==================      
      else if (dim_bdry == 2) {           
          
            unsigned iel_geom_type = msh->GetElementType(iel);
            unsigned iel_geom_type_face = msh->GetElementFaceType(iel, iface);

            unsigned f_n_faces_faces  =  msh->el->GetNFC(iel_geom_type, iel_geom_type_face); /* ElementFaceFaceNumber */

         double mixed_denominator_numerical = 0.;
         
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
                  
                const unsigned dir_x_for_atan = ( ( (LIST_OF_CTRL_FACES :: _face_with_extremes_index[0] /*FACE_FOR_CONTROL*/ - 1) / 2 ) + 1 ) % 3;  ///@todo I think needs to be changed
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
                
                mixed_denominator_numerical += pow(dist, -  2. * s_frac) * delta_teta;
              }
//               delta coords - refinement -----

          
              }
              
            }
            

            
            
            
//              for(unsigned i_bdry = 0; i_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); i_bdry++) {
//                 unsigned int i_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
//                 
//                 Res_local_iel_mixed_num[ i_vol_iel ] += - 0.5 * C_ns * check_limits * operator_Hhalf * beta * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * sol_ctrl_iqp_bdry * weight_iqp_bdry * mixed_denominator_numerical  * (1. / s_frac);
//                 
//                 for(unsigned j_bdry = 0; j_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); j_bdry++) {
//                   unsigned int j_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, j_bdry);
//                   KK_local_iel_mixed_num[ i_vol_iel * nDof_vol_iel + j_vol_iel ] += 0.5 * C_ns * check_limits * operator_Hhalf * beta * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * phi_ctrl_iel_bdry_iqp_bdry[j_bdry]  * weight_iqp_bdry * mixed_denominator_numerical * (1. / s_frac);
//                 }
//                 
//               }
            
            /// @todo ONLY DIFFERENCES: the mixed_denominator is numerical, and so also the corresponding Res and Jac. It could be done with a single function
            
              for (unsigned c = 0; c < n_components_ctrl; c++) {
                 integral +=  0.5 * C_ns * check_limits * operator_Hhalf  * beta  * sol_ctrl_iqp_bdry[c] * sol_ctrl_iqp_bdry[c] * weight_iqp_bdry * mixed_denominator_numerical * (1. / s_frac);
              }
              
              
              
              for (unsigned c = 0; c < n_components_ctrl; c++) {

              for(unsigned i_bdry = 0; i_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); i_bdry++) {
                unsigned int i_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
                
                const unsigned res_pos = assemble_jacobian<double,double>::res_row_index(nDof_vol_iel, c, i_vol_iel);

                Res_local_iel_mixed_num[ res_pos/*i_vol_iel*/ ] += - 0.5 * C_ns * check_limits * operator_Hhalf  * beta * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * sol_ctrl_iqp_bdry[c] * weight_iqp_bdry * mixed_denominator_numerical * (1. / s_frac);
                
                for (unsigned e = 0; e < n_components_ctrl; e++) {
                  if (e == c) { 
                     for(unsigned j_bdry = 0; j_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); j_bdry++) {
                  unsigned int j_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, j_bdry);
                  const unsigned jac_pos = assemble_jacobian< double, double >::jac_row_col_index(nDof_vol_iel, nDof_iel_vec, c, e, i_vol_iel, j_vol_iel);         

                  KK_local_iel_mixed_num[ jac_pos /*i_vol_iel * nDof_vol_iel + j_vol_iel*/ ] += 0.5 * C_ns * check_limits * operator_Hhalf * beta * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * phi_ctrl_iel_bdry_iqp_bdry[j_bdry] * weight_iqp_bdry * mixed_denominator_numerical * (1. / s_frac);
                     }
                
                  }
                }
                
              }
           }
            
          
      }          
      
  }
  
  
 }
  
  
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
  ml_sol_bdry_bdry_flag->GetSolutionLevel(0)->GetSolutionName(node_based_bdry_bdry_flag_name.c_str()) = MED_IO(*ml_mesh.GetLevel(0)).node_based_flag_read_from_file(infile, node_mapping_from_mesh_file_to_new);

  ml_mesh.GetLevelZero(0)->deallocate_node_mapping(node_mapping_from_mesh_file_to_new);

  for(unsigned l = 1; l < ml_mesh.GetNumberOfLevels(); l++) {
     ml_sol_bdry_bdry_flag->RefineSolution(l);
  }
  
  std::vector < std::string > variablesToBePrinted_aux;
  variablesToBePrinted_aux.push_back("all");
  for(unsigned l = 0; l < ml_mesh.GetNumberOfLevels(); l++) {
  ml_sol_bdry_bdry_flag->GetWriter()->Write(l+1, "aux", files.GetOutputPath(), "", "biquadratic", variablesToBePrinted_aux);
   }
 
 
 return ml_sol_bdry_bdry_flag;
 
}


  
  

  
    //********** FRAC CONTROL - BEGIN *****************************************

 static  void control_eqn_bdry(const unsigned iproc,
                                   const unsigned nprocs,
                        MultiLevelProblem &    ml_prob,
                        MultiLevelSolution*    ml_sol,
                        const Solution*        sol,
                        const Mesh * msh,
                        const  LinearEquationSolver* pdeSys,
                        //-----------
                        CurrentElem < double > & geom_element_iel,
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
                        const unsigned int N_div_unbounded,
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
     std::vector < std::vector < double > >  JacI_jel_bdry_jqp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_jel_bdry_jqp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_jel_bdry_jqp_bdry.size(); d++) {   Jac_jel_bdry_jqp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_jel_bdry_jqp_bdry.size(); d++) { JacI_jel_bdry_jqp_bdry[d].resize(dim-1); }
    
    double detJac_jel_bdry_jqp_bdry;
// --- Quadrature - END


// //--- Control domain - BEGIN -------
// //***************************************************
//   unsigned n_max = pow(2,dim_bdry);
//   std::vector < double > extremes(n_max);
//   extremes[0] = LIST_OF_CTRL_FACES :: _face_with_extremes_extremes[0][0] /*EX_1*/;
//   extremes[1] = LIST_OF_CTRL_FACES :: _face_with_extremes_extremes[0][1] /*EX_2*/;
// //   if(dim_bdry == 2){
// //     extremes[2] = EY_1;
// //     extremes[3] = EY_2;
// //   }
//
//   std::vector < std::vector < double > > ex_control(dim); // control zone extremes. Built as: 2D) ( x1 y1 ) ( x2 y2 )    3D) ( x1 y1 z1 ) ( x2 y2 z2 ) ( x3 y3 z3 )
//   for(unsigned d = 0; d < dim; d++) {
//       ex_control[d].reserve(n_max);
//   }
//   int control_xyz   = ( LIST_OF_CTRL_FACES :: _face_with_extremes_index[0] /*FACE_FOR_CONTROL*/ - 1) / 2;
//   bool ctrl_min_max = ( LIST_OF_CTRL_FACES :: _face_with_extremes_index[0] /*FACE_FOR_CONTROL*/ - 1) % 2;
//   for(unsigned d = 0; d < dim; d++) {
//     for(unsigned n_e = 0; n_e < n_max; n_e++){
//       if(control_xyz == d) ex_control[d][n_e] =  (ctrl_min_max)? DOMAIN_EX_2:DOMAIN_EX_1;
//       else ex_control[d][n_e] = extremes[n_e];
//     }
//   }
// //***************************************************
//
//
//--- Control domain - END -------



 //********************* bdry cont *******************
 //*************************************************** 
  std::vector <double> phi_coords_iel_bdry_iqp_bdry;  
  std::vector <double> phi_coords_x_iel_bdry_iqp_bdry; 

  phi_coords_iel_bdry_iqp_bdry.reserve(max_size);
  phi_coords_x_iel_bdry_iqp_bdry.reserve(max_size * space_dim);

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
  vector< vector< int > > l2gMap_iel(n_components_ctrl); 
  vector< vector< int > > l2gMap_jel(n_components_ctrl);

    for (unsigned c = 0; c < n_components_ctrl; c++) {
        l2gMap_iel[c].reserve(max_size);
        l2gMap_jel[c].reserve(max_size);
    }
  
  vector< int > l2gMap_iel_vec;  l2gMap_iel_vec.reserve(elem_dof_size_max);
  vector< int > l2gMap_jel_vec;  l2gMap_jel_vec.reserve(elem_dof_size_max);
  //-------- local to global mappings - END --------------


  //-------- Local matrices and rhs - BEGIN --------------
  vector < double > Res_local_iel; Res_local_iel.reserve(elem_dof_size_max);
  vector < double > KK_local_iel;  KK_local_iel.reserve(elem_dof_size_max * elem_dof_size_max);

//   Local matrices and rhs for adaptive quadrature (iel == jel)
  vector < double > Res_local_iel_refined; Res_local_iel_refined.reserve(elem_dof_size_max);
  vector < double > KK_local_iel_refined;   KK_local_iel_refined.reserve(elem_dof_size_max * elem_dof_size_max);

//   Local matrices and rhs for the mixed internal-external integral term (both adaptive and non-adaptive)
  vector < double > Res_local_iel_mixed_num;  Res_local_iel_mixed_num.reserve(elem_dof_size_max);
  vector < double > KK_local_iel_mixed_num;   KK_local_iel_mixed_num.reserve(elem_dof_size_max * elem_dof_size_max);

//   Non local matrices and vectors for H^s laplacian operator
  vector< double >         Res_nonlocal_iel;  Res_nonlocal_iel.reserve(elem_dof_size_max);
  vector< double >         Res_nonlocal_jel;  Res_nonlocal_jel.reserve(elem_dof_size_max);

  vector < double > KK_nonlocal_iel_iel;  KK_nonlocal_iel_iel.reserve(elem_dof_size_max * elem_dof_size_max);
  vector < double > KK_nonlocal_iel_jel;  KK_nonlocal_iel_jel.reserve(elem_dof_size_max * elem_dof_size_max);
  vector < double > KK_nonlocal_jel_iel;  KK_nonlocal_jel_iel.reserve(elem_dof_size_max * elem_dof_size_max);
  vector < double > KK_nonlocal_jel_jel;  KK_nonlocal_jel_jel.reserve(elem_dof_size_max * elem_dof_size_max); 
  //-------- Local matrices and rhs - END --------------

  
  //-------- Global matrices and rhs - BEGIN --------------
  KK->zero();
  RES->zero(); 
  //-------- Global matrices and rhs - END --------------

  
//   phi_x_placeholder_unused as unused input of certain functions
  vector < double > phi_x_placeholder_unused;
  phi_x_placeholder_unused.reserve(max_size * dim);
  
 
///   boost::mpi::communicator world(MPI_COMM_WORLD, boost::mpi::comm_attach);  /// @todo future solution: broadcast whole class instances


  unsigned count_visits_of_boundary_faces = 0;


 // integral - BEGIN ************
  double integral  = 0.;
// integral - END ************
      

  
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
        
        
	if ( DOMAIN_CONTAINING_CTRL_FACES :: volume_elem_contains_a_Gamma_control_face(geom_element_jel.get_elem_center_3d()) ) {

      
// ***************************************
// ******* jel-related stuff - BEGIN *************
// ***************************************

// --- 2 - solution -----------------
      vector< unsigned > nDof_jel(n_components_ctrl);

      if(kproc == iproc) {
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
        
      if(kproc == iproc) {
          
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
    // --- 2 - solution -----------------

      
// --- 3 - l2GMap -----------------
      for (unsigned c = 0; c < n_components_ctrl; c++) {
          l2gMap_jel[c].resize(nDof_jel[c]);
      }
           
           
      // local storage of global mapping and solution ********************
      if(kproc == iproc) {
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
// --- 3 - l2GMap -----------------



// --- l2GMapVec -----------------
    l2gMap_jel_vec.resize(0);
      for (unsigned c = 0; c < n_components_ctrl; c++) {
          l2gMap_jel_vec.insert(l2gMap_jel_vec.end(), l2gMap_jel[c].begin(), l2gMap_jel[c].end());
      }
// --- l2GMapVec -----------------



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
           jface_is_a_boundary_control = femus::face_is_a_Gamma_control_face< LIST_OF_CTRL_FACES >(msh->el, jel, jface);
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
      std::vector < std::vector < double > > phi_coords_jel_bdry_jqp_bdry(n_jqp_bdry);
      std::vector < std::vector < double > > phi_coords_x_jel_bdry_jqp_bdry(n_jqp_bdry);
      
      std::vector < std::vector < double > > phi_ctrl_jel_bdry_jqp_bdry(n_jqp_bdry);   /// @todo assume all ctrl components have the same FE family
      std::vector < std::vector < double > > phi_ctrl_x_jel_bdry_jqp_bdry(n_jqp_bdry);   /// @todo assume all ctrl components have the same FE family

      std::vector< std::vector< double > > sol_ctrl_jqp_bdry(n_components_ctrl);
      
        for (unsigned c = 0; c < n_components_ctrl; c++) {
             sol_ctrl_jqp_bdry[c].resize(n_jqp_bdry);
                 std::fill(sol_ctrl_jqp_bdry[c].begin(), sol_ctrl_jqp_bdry[c].end(), 0.);
            }

      
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
      for (unsigned c = 0; c < n_components_ctrl; c++) {
         sol_ctrl_jqp_bdry[c][jqp_bdry] = 0.;
//        if(kproc == iproc) {
	      for (int j_bdry = 0; j_bdry < phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry].size()/*Sol_n_el_dofs_quantities[pos_sol_ctrl]*/; j_bdry++)  {
		    unsigned int j_vol = msh->GetLocalFaceVertexIndex_PassElemType(jel_geommm, jface, j_bdry);
			
			sol_ctrl_jqp_bdry[c][jqp_bdry] +=  /*sol_eldofs_Mat[pos_mat_ctrl]*/sol_ctrl_jel[c][j_vol] * phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][j_bdry];

	      }
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
	if ( DOMAIN_CONTAINING_CTRL_FACES ::volume_elem_contains_a_Gamma_control_face( geom_element_iel.get_elem_center_3d() ) ) {
        
// ***************************************
// ******* iel-related stuff - BEGIN *************
// ***************************************
        
      
// --- l2GMap
        vector< unsigned > nDof_iel(n_components_ctrl);
        
       for (unsigned c = 0; c < n_components_ctrl; c++) {
         nDof_iel[c]  = msh->GetElementDofNumber(iel, solType_ctrl[c]);
       }
              
      for (unsigned c = 0; c < n_components_ctrl; c++) {
        l2gMap_iel[c].resize(nDof_iel[c]);
        for(unsigned i = 0; i < l2gMap_iel[c].size(); i++) {
          l2gMap_iel[c][i] = pdeSys->GetSystemDof(solIndex_ctrl[c], solPdeIndex_ctrl[c], i, iel);
           }
         }
// --- l2GMap

// --- l2GMapVec -----------------
    l2gMap_iel_vec.resize(0);
      for (unsigned c = 0; c < n_components_ctrl; c++) {
          l2gMap_iel_vec.insert(l2gMap_iel_vec.end(), l2gMap_iel[c].begin(), l2gMap_iel[c].end());
      }
// --- l2GMapVec -----------------


// --- element matrix and vector resizing
unsigned nDof_iel_vec = 0;
        for (unsigned c = 0; c < n_components_ctrl; c++) {  nDof_iel_vec  +=  nDof_iel[c];    }


        KK_nonlocal_iel_iel.assign(nDof_iel_vec * nDof_iel_vec, 0.);   //resize
        KK_nonlocal_iel_jel.assign(nDof_iel_vec * nDof_jel_vec, 0.);   //resize
        KK_nonlocal_jel_iel.assign(nDof_jel_vec * nDof_iel_vec, 0.);   //resize
        KK_nonlocal_jel_jel.assign(nDof_jel_vec * nDof_jel_vec, 0.);   //resize
        Res_nonlocal_iel.assign(nDof_iel_vec, 0.);    //resize
        Res_nonlocal_jel.assign(nDof_jel_vec, 0.);    //resize

        Res_local_iel_mixed_num.assign(nDof_iel_vec, 0.);    //resize
        KK_local_iel_mixed_num.assign(nDof_iel_vec * nDof_iel_vec, 0.);

        if( check_if_same_elem(iel, jel) ) {
          Res_local_iel.assign(nDof_iel_vec, 0.);    //resize
          KK_local_iel.assign(nDof_iel_vec * nDof_iel_vec, 0.);
          if(integration_num_split != 0) {
//             Vectors and matrices for adaptive quadrature
            Res_local_iel_refined.assign(nDof_iel_vec, 0.);    //resize
            KK_local_iel_refined.assign(nDof_iel_vec * nDof_iel_vec, 0.);
          }
        }
// --- element matrix and vector resizing 


// --- solution        

        for (unsigned c = 0; c < n_components_ctrl; c++) {
           sol_ctrl_iel[c].resize(nDof_iel[c]);

           
        for(unsigned i = 0; i < sol_ctrl_iel[c].size(); i++) {
          unsigned iDof  = msh->GetSolutionDof(i, iel, solType_ctrl[c]);  // global to global mapping between coordinates node and coordinate dof
          sol_ctrl_iel[c][i] = (*sol->_Sol[solIndex_ctrl[c]])(iDof);  // global extraction and local storage for the element coordinates
          }
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

   
	    if( femus::face_is_a_Gamma_control_face< LIST_OF_CTRL_FACES >(msh->el, iel, iface) ) {
//------------ iface opening ---------        
		
//                 count_visits_of_boundary_faces++;


             
     //**** Adaptive preparation - BEGIN ******** 
     std::vector < std::vector < double > >  Jac_kel_bdry_kqp_bdry;
     std::vector < std::vector < double > >  JacI_kel_bdry_kqp_bdry;
     double detJac_kel_bdry_kqp_bdry;
      vector < double >  x_kqp_bdry;
      double  weight_kqp_bdry;
      vector < double >  phi_coords_kel_bdry_kqp_bdry;
      vector < double >  phi_coords_x_kel_bdry_kqp_bdry;
      
      vector < double >  phi_ctrl_kel_bdry_kqp_bdry; /// @todo assume all ctrl components have the same FE family
      vector < double >  phi_ctrl_x_kel_bdry_kqp_bdry;
      
       std::vector< double > sol_ctrl_kqp_bdry(n_components_ctrl);
            
        for (unsigned c = 0; c < n_components_ctrl; c++) {   sol_ctrl_kqp_bdry[c] = 0.;   }


      
//         double weight3;
//         vector < double > phi3;

//        Evaluating coarse FE functions on Quadrature Points of the "sub-elements"
       std::vector < std::vector < std::vector <double > > > aP(3);  //[NFE_FAMS][DIM==3][N_DOFS]
        if(integration_num_split != 0) {
          for(unsigned fe_type = 0; fe_type < /*solType*/solType_coords + 1; fe_type++) { //loop up to the FE type + 1 of the unknown
            ProjectNodalToPolynomialCoefficients(aP[fe_type], geom_element_iel.get_coords_at_dofs_bdry_3d(), ielGeom_bdry, fe_type) ;         ///@todo check this!!!!!!!!!!!!
          }
        }                      
     //**** Adaptive preparation - END ********  
                      
              //Quadrature loop initialization - BEGIN
        const unsigned n_iqp_bdry = ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussPointsNumber();
        
        
       std::vector< double > sol_ctrl_iqp_bdry(n_components_ctrl);
        for (unsigned c = 0; c < n_components_ctrl; c++) {   sol_ctrl_iqp_bdry[c] = 0.;   }
              //Quadrature loop initialization - END
    
         
		for(unsigned iqp_bdry = 0; iqp_bdry < n_iqp_bdry; iqp_bdry++) {
            
    elem_all[qrule_i][ielGeom_bdry][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_bdry_3d(), iqp_bdry, Jac_iel_bdry_iqp_bdry, JacI_iel_bdry_iqp_bdry, detJac_iel_bdry_iqp_bdry, space_dim);
    
    weight_iqp_bdry = detJac_iel_bdry_iqp_bdry * ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussWeightsPointer()[iqp_bdry];

    elem_all[qrule_i][ielGeom_bdry][solType_coords] ->shape_funcs_current_elem(iqp_bdry, JacI_iel_bdry_iqp_bdry, phi_coords_iel_bdry_iqp_bdry, phi_coords_x_iel_bdry_iqp_bdry, boost::none, space_dim);

    elem_all[qrule_i][ielGeom_bdry][SolFEType_quantities[pos_sol_ctrl]] ->shape_funcs_current_elem(iqp_bdry, JacI_iel_bdry_iqp_bdry, phi_ctrl_iel_bdry_iqp_bdry, phi_ctrl_x_iel_bdry_iqp_bdry, boost::none, space_dim);
            
//========== compute gauss quantities on the boundary - BEGIN ===============================================
//--- geom
          std::vector < double > x_iqp_bdry(dim, 0.);  ///@todo is this dim or dim_bdry?

            for(unsigned d = 0; d < x_iqp_bdry.size(); d++) {
	      for (int i_bdry = 0; i_bdry < geom_element_iel.get_coords_at_dofs_bdry_3d()[d].size(); i_bdry++)  {
			
              x_iqp_bdry[d] += geom_element_iel.get_coords_at_dofs_bdry_3d()[d][i_bdry] * phi_coords_iel_bdry_iqp_bdry[i_bdry];

		      }
            }
//--- geom
    
//--- solution
        for (unsigned c = 0; c < n_components_ctrl; c++) {
          sol_ctrl_iqp_bdry[c] = 0.;
	        for (int i_bdry = 0; i_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
			
			sol_ctrl_iqp_bdry[c] +=  /*sol_eldofs_Mat[pos_mat_ctrl]*/sol_ctrl_iel[c][i_vol] * phi_ctrl_iel_bdry_iqp_bdry[i_bdry];

		      }
       }
//--- solution
//========== compute gauss quantities on the boundary - END ================================================


  
      //============  Non-fractional assembly - BEGIN ==================
          
            if( check_if_same_elem_bdry(iel, jel, iface, jface) ) {
              
                
       //============  Mass assembly - BEGIN ==================
    for (unsigned c = 0; c < n_components_ctrl; c++) {
        
        integral += operator_L2 * alpha * weight_iqp_bdry * sol_ctrl_iqp_bdry[c] * sol_ctrl_iqp_bdry[c];
    }
  
    
    
    for (unsigned c = 0; c < n_components_ctrl; c++) {
          for(unsigned l_bdry = 0; l_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); l_bdry++) {
               		    unsigned int l_vol = msh->GetLocalFaceVertexIndex(iel, iface, l_bdry);

              double mass_res_i = phi_ctrl_iel_bdry_iqp_bdry[l_bdry] * sol_ctrl_iqp_bdry[c];
              const unsigned res_pos = assemble_jacobian<double,double>::res_row_index(nDof_iel, c, l_vol);
              Res_local_iel[ res_pos/*l_vol*/ ] += operator_L2 * alpha * weight_iqp_bdry * mass_res_i ;
              Res_local_iel[ res_pos/*l_vol*/ ] += - rhs_one * weight_iqp_bdry * (phi_ctrl_iel_bdry_iqp_bdry[l_bdry] * (-1.) /** ( sin(2 * acos(0.0) * x1[0][l_bdry])) * ( sin(2 * acos(0.0) * x1[1][l_bdry]))*/);
              
              for (unsigned e = 0; e < n_components_ctrl; e++) {
                  if (e == c) { 
            for(unsigned m_bdry = 0; m_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); m_bdry++) {
               		    unsigned int m_vol = msh->GetLocalFaceVertexIndex(iel, iface, m_bdry);
                const unsigned jac_pos = assemble_jacobian< double, double >::jac_row_col_index(nDof_iel, nDof_iel_vec, c, e, l_vol, m_vol);         
                KK_local_iel[ jac_pos/*l_vol * nDof_iel + m_vol*/ ] += operator_L2 * alpha * phi_ctrl_iel_bdry_iqp_bdry[l_bdry] * phi_ctrl_iel_bdry_iqp_bdry[m_bdry] * weight_iqp_bdry;
                 }
               }
             }
             
            } 
       }
        //============  Mass assembly - END ==================
         
      //============  Laplacian assembly - BEGIN ==================
      //============  Laplacian assembly - END ==================
            }
    
      //============  Non-fractional assembly - END ==================
             
             
      //============  Fractional assembly - BEGIN ==================
        if(operator_Hhalf != 0) {
                 
      //============ Same elem, && Adaptive quadrature - BEGIN ==================
        if( check_if_same_elem_bdry(iel, jel, iface, jface) && integration_num_split != 0) {
                
          /*const*/ short unsigned kelGeom_bdry = ielGeom_bdry;
           
          const unsigned n_kqp_bdry = ml_prob.GetQuadratureRuleMultiple(qrule_k, kelGeom_bdry).GetGaussPointsNumber();
                
                
          std::cout.precision(14);
          std::vector< std::vector< std::vector<double> > > x3;

          for(unsigned split = 0; split <= integration_num_split; split++) {

            
            if (dim_bdry/*dim*/ == 1) GetElementPartition1D(x_iqp_bdry, geom_element_iel.get_coords_at_dofs_bdry_3d(), split, integration_num_split, x3, space_dim); //TODO space_dim or dim?
            else if (dim_bdry/*dim*/ == 2) {
              //TODO need to be checked !!!
              //GetElementPartition2D(x_iqp_bdry, geom_element_iel.get_coords_at_dofs_bdry_3d(), split, integration_num_split, x3);
              GetElementPartitionQuad(x_iqp_bdry, geom_element_iel.get_coords_at_dofs_bdry_3d(), split, integration_num_split, x3);
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
                
                const double common_weight =  0.5 * C_ns * operator_Hhalf * beta * check_limits * weight_iqp_bdry * weight_kqp_bdry  / denom_ik;


     for (unsigned c = 0; c < n_components_ctrl; c++) {
         integral    +=       common_weight * (sol_ctrl_iqp_bdry[c] - solY3[c]) * (sol_ctrl_iqp_bdry[c] - solY3[c]);
     }
     
                for (unsigned c = 0; c < n_components_ctrl; c++) {
               for(unsigned l_bdry = 0; l_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); l_bdry++) { //dofs of test function
                  unsigned int l_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, l_bdry);

              const unsigned res_pos = assemble_jacobian<double,double>::res_row_index(nDof_iel, c, l_vol_iel);
                  Res_local_iel_refined[res_pos/* l_vol_iel*/ ]    +=      - common_weight * (sol_ctrl_iqp_bdry[c] - solY3[c]) * (phi_ctrl_iel_bdry_iqp_bdry[l_bdry] - phi_ctrl_kel_bdry_kqp_bdry[l_bdry]);

              for (unsigned e = 0; e < n_components_ctrl; e++) {
                  if (e == c) {
                for(unsigned m_bdry = 0; m_bdry < phi_ctrl_kel_bdry_kqp_bdry.size(); m_bdry++) { //dofs of unknown function
                    unsigned int m_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, m_bdry);
                    
                    const unsigned jac_pos = assemble_jacobian< double, double >::jac_row_col_index(nDof_iel, nDof_iel_vec, c, e, l_vol_iel, m_vol_iel);         
                KK_local_iel_refined[ jac_pos/*l_vol_iel * nDof_jel + m_vol_iel*/ ] += common_weight * (phi_ctrl_iel_bdry_iqp_bdry[m_bdry] - phi_ctrl_kel_bdry_kqp_bdry[m_bdry]) * 
                                                        (phi_ctrl_iel_bdry_iqp_bdry[l_bdry] - phi_ctrl_kel_bdry_kqp_bdry[l_bdry]);

                  }
                }
                  
              }
            }
     }
// ********* BOUNDED PART - END ***************
// ********* UNBOUNDED PART - BEGIN ***************
                if( iqp_bdry == integration_split_index ) { ///@todo is there a way to put this outside of the quadrature loop?
                    
              mixed_integral(unbounded,
                              dim,
                              dim_bdry,
//                               ex_control,
                              N_div_unbounded,
                              weight_kqp_bdry,
                              x_kqp_bdry,
                              phi_ctrl_kel_bdry_kqp_bdry,
                              sol_ctrl_kqp_bdry,
                              s_frac,
                              check_limits,
                              C_ns,
                              operator_Hhalf,
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
                             node_based_bdry_bdry_in,
                              bdry_bdry,
                              geom_element_jel,
                              jelGeom_bdry,
                              solType_coords,
                              integral
                             ); 
                }


// ********* UNBOUNDED PART - END ***************
                                         
                  } //end k_qp_bdry
                } //end r
              }  //end split
              
              
                
                
                
            }  //end iel == jel && integration_num_split != 0
      //============ Same elem, && Adaptive quadrature - END ==================
              
             
      //============ Either different elements, or lack of adaptivity (so all elements) - BEGIN ==================
        else {  //  if(iel != jel || integration_num_split == 0) 
            
// ********* UNBOUNDED PART - BEGIN ***************
//           if(check_if_same_elem_bdry(iel, jel, iface, jface)) { //TODO I removed this since we don't want iel==jel here
              
               mixed_integral(unbounded,
                              dim,
                              dim_bdry,
//                               ex_control,
                              N_div_unbounded,
                              weight_iqp_bdry,
                              x_iqp_bdry,
                              phi_ctrl_iel_bdry_iqp_bdry,
                              sol_ctrl_iqp_bdry,
                              s_frac,
                              check_limits,
                              C_ns,
                              operator_Hhalf,
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
                              node_based_bdry_bdry_in,
                              bdry_bdry,
                              geom_element_jel,
                              jelGeom_bdry,
                              solType_coords,
                              integral
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
              
              const double common_weight = (0.5 * C_ns) * operator_Hhalf * beta * check_limits * weight_iqp_bdry * weight_jqp_bdry[jqp_bdry]  / denom;

              for (unsigned c = 0; c < n_components_ctrl; c++) {
                  
                  integral +=  common_weight * (sol_ctrl_iqp_bdry[c] - sol_ctrl_jqp_bdry[c][jqp_bdry]) * sol_ctrl_iqp_bdry[c];
                  integral +=  common_weight * (sol_ctrl_iqp_bdry[c] - sol_ctrl_jqp_bdry[c][jqp_bdry]) * (- sol_ctrl_jqp_bdry[c][jqp_bdry]);
                  
              }              
              
              
              for (unsigned c = 0; c < n_components_ctrl; c++) {
                 for(unsigned l_bdry = 0; l_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); l_bdry++) { //dofs of test function
               		    unsigned int l_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, l_bdry);
               		    unsigned int l_vol_jel = msh->el->GetIG(jel_geommm, jface, l_bdry)/*msh->GetLocalFaceVertexIndex(jel, jface, l_bdry)*/;

              const unsigned res_pos_iel = assemble_jacobian<double,double>::res_row_index(nDof_iel, c, l_vol_iel);
              const unsigned res_pos_jel = assemble_jacobian<double,double>::res_row_index(nDof_jel, c, l_vol_jel);
                Res_nonlocal_iel[ res_pos_iel /*l_vol_iel*/ ]      +=      - common_weight * (sol_ctrl_iqp_bdry[c] - sol_ctrl_jqp_bdry[c][jqp_bdry]) * (phi_ctrl_iel_bdry_iqp_bdry[l_bdry]);

                Res_nonlocal_jel[ res_pos_jel /*l_vol_jel*/ ]      +=      - common_weight * (sol_ctrl_iqp_bdry[c] - sol_ctrl_jqp_bdry[c][jqp_bdry]) * (- phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][l_bdry]);

               
//                 for(unsigned j = 0; j < nDof_jel; j++) {
                for (unsigned e = 0; e < n_components_ctrl; e++) {
                  if (e == c) {
                      for(unsigned m_bdry = 0; m_bdry < phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry].size(); m_bdry++) { //dofs of unknown function
               		    unsigned int m_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, m_bdry);
               		    unsigned int m_vol_jel = msh->el->GetIG(jel_geommm, jface, m_bdry)/*msh->GetLocalFaceVertexIndex(jel, jface, m_bdry)*/;

                    const unsigned jac_pos_iel_iel = assemble_jacobian< double, double >::jac_row_col_index(nDof_iel, nDof_iel_vec, c, e, l_vol_iel, m_vol_iel);         
                    const unsigned jac_pos_iel_jel = assemble_jacobian< double, double >::jac_row_col_index(nDof_iel, nDof_iel_vec, c, e, l_vol_iel, m_vol_jel);
                    const unsigned jac_pos_jel_iel = assemble_jacobian< double, double >::jac_row_col_index(nDof_iel, nDof_iel_vec, c, e, l_vol_jel, m_vol_iel);
                    const unsigned jac_pos_jel_jel = assemble_jacobian< double, double >::jac_row_col_index(nDof_iel, nDof_iel_vec, c, e, l_vol_jel, m_vol_jel);
                    
             /*  u(x) v(x)*/     KK_nonlocal_iel_iel[ jac_pos_iel_iel /*l_vol_iel * nDof_jel + m_vol_iel*/ ] += common_weight *          phi_ctrl_iel_bdry_iqp_bdry[m_bdry]            *    phi_ctrl_iel_bdry_iqp_bdry[l_bdry];

             /*- u(y) v(x)*/     KK_nonlocal_iel_jel[ jac_pos_iel_jel /*l_vol_iel * nDof_jel + m_vol_jel*/ ] += common_weight * (- 1.) * phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][m_bdry]  *    phi_ctrl_iel_bdry_iqp_bdry[l_bdry];

             /*- u(x) v(y)*/     KK_nonlocal_jel_iel[ jac_pos_jel_iel /*l_vol_jel * nDof_jel + m_vol_iel*/ ] += common_weight * (- 1.) * phi_ctrl_iel_bdry_iqp_bdry[m_bdry]            *   phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][l_bdry];

             /*  u(y) v(y)*/     KK_nonlocal_jel_jel[ jac_pos_jel_jel /*l_vol_jel * nDof_jel + m_vol_jel*/ ] += common_weight *          phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][m_bdry]  *    phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][l_bdry];


                     }
                   }
                 }
                }
                
               }
               
              } //endl jqp_bdry loop
// ********* BOUNDED PART - END ***************
            
            
         } //end if(iel != jel || integration_num_split == 0)
      //============ Either different elements, or lack of adaptivity (so all elements) - END ==================

        
         } //operator_Hhalf != 0              
      //============  Fractional assembly - END ==================
            
      }   //end iqp_bdry
      
//                std::vector<unsigned> Sol_n_el_dofs_Mat_vol2(1, nDof_jel);
//          assemble_jacobian<double,double>::print_element_residual(iel, Res, Sol_n_el_dofs_Mat_vol, 10, 5);
//          assemble_jacobian<double,double>::print_element_jacobian(iel, KK_local_iel_mixed_num, Sol_n_el_dofs_Mat_vol2, 10, 5);
      
      
              
//----- iface ---        
        } //end face for control
        
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
     

         RES->add_vector_blocked(Res_local_iel, l2gMap_iel_vec);
          KK->add_matrix_blocked(KK_local_iel, l2gMap_iel_vec, l2gMap_iel_vec);
        
        
          if(integration_num_split != 0) {
            RES->add_vector_blocked(Res_local_iel_refined, l2gMap_iel_vec);
            KK->add_matrix_blocked(KK_local_iel_refined, l2gMap_iel_vec, l2gMap_iel_vec);
          }
          
       }

        RES->add_vector_blocked(Res_local_iel_mixed_num, l2gMap_iel_vec);
        KK->add_matrix_blocked(KK_local_iel_mixed_num, l2gMap_iel_vec, l2gMap_iel_vec);
        

        RES->add_vector_blocked(Res_nonlocal_iel, l2gMap_iel_vec);
        RES->add_vector_blocked(Res_nonlocal_jel, l2gMap_jel_vec);
        
        KK->add_matrix_blocked(KK_nonlocal_iel_iel, l2gMap_iel_vec, l2gMap_iel_vec);
        KK->add_matrix_blocked(KK_nonlocal_iel_jel, l2gMap_iel_vec, l2gMap_jel_vec);
        KK->add_matrix_blocked(KK_nonlocal_jel_iel, l2gMap_jel_vec, l2gMap_iel_vec);
        KK->add_matrix_blocked(KK_nonlocal_jel_jel, l2gMap_jel_vec, l2gMap_jel_vec);
// Since 1 is dense and 3 are sparse, and the dense dofs are 30, we should have at most 3x9 + 30 = 57, but in the sparsity print it shows 30. That's the problem

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
     
} //end face for control


      } //end jface

 //----- jface ---   
 
     
    }  //end control elem flag jel
    



   } //end jel
//----- jel ---        
   
 } //end kproc
    
    
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
