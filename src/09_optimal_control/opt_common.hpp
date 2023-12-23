#ifndef __OPT_COMMON_HPP__
#define __OPT_COMMON_HPP__


#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "CurrentElem.hpp"
// #include "Mesh.hpp"
#include "Elem.hpp"


 namespace femus {
     
     
      std::vector< std::vector< int > > is_dof_associated_to_Gamma_control_equation(
      const femus::Mesh * msh,
      /*const*/ femus::MultiLevelSolution * ml_sol,
         const femus::MultiLevelProblem *    ml_prob,
      const unsigned iel,
      const unsigned short iel_geom_elem_type,
      femus::CurrentElem < double > & geom_element_iel,
      const unsigned solType_coords,
      std::vector < std::string > Solname_Mat,
      std::vector < unsigned > SolFEType_Mat,
      std::vector < unsigned > Sol_n_el_dofs_Mat,
      const unsigned pos_mat_ctrl,
      const unsigned n_components_ctrl
            ) {
//=============== construct control node flag field  =========================    
	      /* For every component:
           * (control_node_flag[c][i])       picks nodes on \Gamma_c
           * (1 - control_node_flag[c][i])   picks nodes on \Omega \setminus \Gamma_c
	       */
          
       std::vector< std::vector< int > > control_node_flag_iel_all_faces(n_components_ctrl);
       
	  for (unsigned c = 0; c < n_components_ctrl; c++) {
              control_node_flag_iel_all_faces[c].resize(Sol_n_el_dofs_Mat[pos_mat_ctrl + c]);
              std::fill(control_node_flag_iel_all_faces[c].begin(), control_node_flag_iel_all_faces[c].end(), 0);   
         }
       
          
	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
          
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();

          
	    // look for boundary faces
            const int bdry_index = msh->GetMeshElements()->GetFaceElementIndex(iel, iface);
            
	    if( bdry_index < 0) {
	      const unsigned int face_in_rectangle_domain = - ( bdry_index + 1);
         
         
         
	  for (unsigned c = 0; c < n_components_ctrl; c++) {
          
	     double tau = 0.;
         
	      const bool  dir_bool_c = ml_sol->GetBdcFunctionMLProb()(ml_prob, geom_element_iel.get_elem_center_bdry_3d(), Solname_Mat[pos_mat_ctrl + c].c_str(), tau, face_in_rectangle_domain, 0.);

	      if (dir_bool_c == false) {
              
          const unsigned ndofs_ctrl_bdry = msh->GetElementFaceDofNumber(iel, iface, SolFEType_Mat[pos_mat_ctrl + c]);
		  for(unsigned i_bdry = 0; i_bdry < ndofs_ctrl_bdry; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex_PassElemType(iel_geom_elem_type, iface, i_bdry);
		//we use the dirichlet flag to say: if dirichlet = true, we set 1 on the diagonal. if dirichlet = false, we put the boundary equation
		
				  control_node_flag_iel_all_faces[c][i_vol] = 1;
			    
              }
         }
      } 
        
        
      }   
    }
    
    return   control_node_flag_iel_all_faces;
    
  }
  

  template < class T >
  std::pair< bool, unsigned int > face_is_a_Gamma_control_face_of_some_index(/*const*/ femus::elem * el, const unsigned jel, const unsigned jface) {

  	    // look for boundary faces
            const int bdry_index_j = el->GetFaceElementIndex(jel, jface);
	    // look for face number equal to control face
	      const unsigned int face_index_in_rectangle_domain_j = - ( el->GetFaceElementIndex(jel,jface) + 1);

	    // look for boundary faces && look for control faces

          bool  is_face_for_control = false;

          		  for(unsigned f = 0; f <  T ::_face_with_extremes_index_size; f++) {
                      if (face_index_in_rectangle_domain_j ==  T ::_face_with_extremes_index[f]) { is_face_for_control = true; }
                  }


   return std::pair< bool, unsigned int >/*std::make_pair*/( ( bdry_index_j < 0 && is_face_for_control ), face_index_in_rectangle_domain_j );

  }



  
  
 }

  
#endif
