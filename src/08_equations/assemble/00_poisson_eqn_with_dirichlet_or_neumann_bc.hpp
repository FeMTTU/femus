#ifndef __femus_00_laplacian_hpp__
#define __femus_00_laplacian_hpp__

#include "MultiLevelProblem.hpp"
#include "MultiLevelMesh.hpp"
#include "MultiLevelSolution.hpp"

#include "CurrentElem.hpp"
#include "ElemType_template.hpp"

#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"


namespace femus   {

 class poisson_equation   {
/// This only works with the App Specifics class, because only with that it can be abstract enough
     
 public:    


template < class real_num, class real_num_mov >
static void equation_with_dirichlet_or_neumann_bc(MultiLevelProblem& ml_prob) {

  NonLinearImplicitSystem* mlPdeSys  = &ml_prob.get_system< NonLinearImplicitSystem > (ml_prob.get_app_specs_pointer()->_system_name);  
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);
  elem*                     el = msh->GetMeshElements();

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*             JAC = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

  unsigned    iproc = msh->processor_id();


//***************** Check Number of Unknowns **********************************  
  const int n_vars = mlPdeSys->GetSolPdeIndex().size();
  if (n_vars != 1)  { std::cout << "Function written for only 1 scalar unknown";  abort();  }


  //=============== Geometry ========================================
  unsigned xType = CONTINUOUS_BIQUADRATIC; // the FE for the domain need not be biquadratic
  
  CurrentElem < double > geom_element(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
//***************************************************  


  
 //******************** quadrature *******************************  
  double jacXweight_qp; 

 //********************* quadrature, geometry *********************** 
 //***************************************************  
  std::vector <double> phi_coords;
  std::vector <double> phi_coords_x; 
  std::vector <double> phi_coords_xx;

  phi_coords.reserve(max_size);
  phi_coords_x.reserve(max_size * space_dim);
  phi_coords_xx.reserve(max_size * dim2);
  
  
 //********************* quadrature, unknowns *********************** 
 //***************************************************  
  std::vector <double> phi_u;
  std::vector <double> phi_u_x; 
  std::vector <double> phi_u_xx;

  phi_u.reserve(max_size);
  phi_u_x.reserve(max_size * space_dim);
  phi_u_xx.reserve(max_size * dim2);
  
  
  const std::string solname_u = ml_sol->GetSolName_string_vec()[0];
  
  unsigned solIndex_u = ml_sol->GetIndex(solname_u.c_str()); 
  unsigned solFEType_u = ml_sol->GetSolutionType(solIndex_u); 

  unsigned solPdeIndex_u = mlPdeSys->GetSolPdeIndex(solname_u.c_str());


  std::vector< int > l2GMap_u;    l2GMap_u.reserve(max_size);
  std::vector < double >  sol_u;     sol_u.reserve(max_size);
 //***************************************************  
 //***************************************************  

  
 //***************************************************  
 //********* WHOLE SET OF VARIABLES ****************** 

  std::vector< int > l2GMap_AllVars; l2GMap_AllVars.reserve(n_vars*max_size); // local to global mapping
  std::vector< double >         Res;            Res.reserve(n_vars*max_size);  // local redidual vector
  std::vector < double >        Jac;            Jac.reserve(n_vars*max_size * n_vars*max_size);
 //***************************************************

  RES->zero();
  if (assembleMatrix)  JAC->zero();

  
 //***************************************************  
     std::vector < std::vector < real_num_mov > >  JacI_qp(space_dim);
     std::vector < std::vector < real_num_mov > >  Jac_qp(dim);
    for (unsigned d = 0; d < dim; d++) { Jac_qp[d].resize(space_dim); }
    for (unsigned d = 0; d < space_dim; d++) { JacI_qp[d].resize(dim); }
    
    real_num_mov detJac_qp;
    

  //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<real_num, real_num_mov> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
 //***************************************************  
  
  

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {
      

    geom_element.set_coords_at_dofs_and_geom_type(iel, xType);
        
    const short unsigned ielGeom = geom_element.geom_type();

 //**************** state **************************** 
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solFEType_u);
    sol_u    .resize(nDof_u);
    l2GMap_u.resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
     unsigned solDof_u = msh->GetSolutionDof(i, iel, solFEType_u);
      sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);
      l2GMap_u[i] = pdeSys->GetSystemDof(solIndex_u, solPdeIndex_u, i, iel);
    }
 //***************************************************  
 
 //******************** ALL VARS ********************* 
    unsigned nDof_AllVars = nDof_u; 
    int nDof_max    =  nDof_u;   // TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
    
    Res.resize(nDof_AllVars);                  std::fill(Res.begin(), Res.end(), 0.);
    Jac.resize(nDof_AllVars * nDof_AllVars);   std::fill(Jac.begin(), Jac.end(), 0.);
    l2GMap_AllVars.resize(0);                  l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_u.begin(),l2GMap_u.end());
 //*************************************************** 
    

 //========= BOUNDARY - BEGIN ==================   
    if (dim == 1)   femus::poisson_equation::natural_loop_1d/*<real_num, real_num_mov>*/(& ml_prob, msh, ml_sol,
                      iel, geom_element, xType,
                      solname_u, solFEType_u,
                      Res
                     );

    if (dim == 2 || dim == 3)   femus::poisson_equation::natural_loop_2d3d<real_num, real_num_mov>(& ml_prob, msh, ml_sol,
                      iel, geom_element, xType,
                      solname_u, solFEType_u,
                      Res,
                      elem_all,
                      dim,
                      space_dim,
                      max_size
                     );
 //========= BOUNDARY - END ==================   
 
 //========= VOLUME - BEGIN ==================   
   
 //========= gauss value quantities ==================   
	std::vector<double> sol_u_x_gss(space_dim);     std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
 //===================================================   
    
    //--------------    from giorgiobornia.cpp
 /// @assignment You need to evaluate your manufactured right hand side at the quadrature point qp.
 /// Hence, you need to compute the coordinates of the quadrature point. Let us call them x_qp.
 /// These are obtained just like every quantity at a quadrature point, i.e., by interpolating the values of the quantity at the element nodes.
 /// The interpolation is performed by using the shape functions.
 /// In other words, 
 ///         (x_qp) = summation of (x_nodes) * (shape function of that node, evaluated at qp)
 /// 
 ///   (x_nodes) are obtained from   geom_element.get_coords_at_dofs_3d()  (this is a  std::vector < vector >,  where the outer index is the dimension and the inner index ranges over the nodes) 
 ///   (shape function of that node, evaluated at qp)  is obtained from phi_u  (this is a vector, whose index ranges over the nodes)
 
 
//--------------   
    
    
    
      // *** Quadrature point loop ***
      for (unsigned i_qp = 0; i_qp < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); i_qp++) {
          
        // *** get gauss point weight, test function and test function partial derivatives ***
// 	msh->_finiteElement[ielGeom][solFEType_u]->Jacobian(geom_element.get_coords_at_dofs_3d(),    i_qp, weight,    phi_u,    phi_u_x,    boost::none /*phi_u_xx*/);
          
	elem_all[ielGeom][xType]->JacJacInv(geom_element.get_coords_at_dofs_3d(), i_qp, Jac_qp, JacI_qp, detJac_qp, space_dim);
    jacXweight_qp = detJac_qp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[i_qp];
    elem_all[ielGeom][solFEType_u]->shape_funcs_current_elem(i_qp, JacI_qp, phi_u, phi_u_x, boost::none /*phi_u_xx*/, space_dim);
    elem_all[ielGeom][xType]->shape_funcs_current_elem(i_qp, JacI_qp, phi_coords, phi_coords_x, boost::none /*phi_coords_xx*/, space_dim);


//--------------    
	std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	
	for (unsigned i = 0; i < nDof_u; i++) {
// 	                                                sol_u_gss      += sol_u[i] * phi_u[i];
                   for (unsigned d = 0; d < sol_u_x_gss.size(); d++)   sol_u_x_gss[d] += sol_u[i] * phi_u_x[i * space_dim + d];
          }
//--------------    

//--------------    
 /// @assignment You need to evaluate your manufactured right hand side at the quadrature point qp.
 /// Hence, you need to compute the coordinates of the quadrature point. Let us call them x_qp.
 /// These are obtained just like every quantity at a quadrature point, i.e., by interpolating the values of the quantity at the element nodes.
 /// The interpolation is performed by using the shape functions.
 /// In other words, 
 ///         (x_qp) = summation of (x_nodes) * (shape function of that node, evaluated at qp)
 /// 
 ///   (x_nodes) are obtained from   geom_element.get_coords_at_dofs_3d()  (this is a  std::vector < vector >,  where the outer index is the dimension and the inner index ranges over the nodes) 
 ///   (shape function of that node, evaluated at qp)  is obtained from phi_u  (this is a vector, whose index ranges over the nodes)
 
 std::vector<double> x_qp(dim, 0.);
          
        for (unsigned i = 0; i < phi_coords.size(); i++) {
          	for (unsigned d = 0; d < dim; d++) {
	                                                x_qp[d]    += geom_element.get_coords_at_dofs_3d()[d][i] * phi_coords[i]; // fetch of coordinate points
            }
        }
 

          
//==========FILLING WITH THE EQUATIONS ===========
	// *** phi_i loop ***
        for (unsigned i = 0; i < nDof_max; i++) {
	  
//--------------    
	      double laplace_res_du_u_i = 0.;
              if ( i < nDof_u ) {
                  for (unsigned kdim = 0; kdim < space_dim; kdim++) {
                       laplace_res_du_u_i             +=  phi_u_x   [i * space_dim + kdim] * sol_u_x_gss[kdim];
	             }
              }
//--------------    
              
	      
//======================Residuals=======================
          // FIRST ROW
 /// @assignment for your manufactured right-hand side, implement a function that receives the coordinate of the quadrature point
 /// Put it after the includes, in the top part of this file
 if (i < nDof_u)                      Res[0      + i] += jacXweight_qp * ( phi_u[i] * (  - ml_prob.get_app_specs_pointer()->_assemble_function_for_rhs->laplacian(x_qp)  ) - laplace_res_du_u_i);
//           if (i < nDof_u)                      Res[0      + i] +=  jacXweight_qp * ( phi_u[i] * (  1. ) - laplace_beltrami_res_du_u_i);
//======================Residuals=======================
	      
          if (assembleMatrix) {
	    
            // *** phi_j loop ***
            for (unsigned j = 0; j < nDof_max; j++) {


//--------------    
              double laplace_mat_du_u_i_j = 0.;

              if ( i < nDof_u && j < nDof_u ) {
                  for (unsigned kdim = 0; kdim < space_dim; kdim++) {
                         laplace_mat_du_u_i_j           += phi_u_x   [i * space_dim + kdim] *
                                                           phi_u_x   [j * space_dim + kdim];
	              }
             }
//--------------    



              //============ delta_state row ============================
              //DIAG BLOCK delta_state - state
		  if ( i < nDof_u && j < nDof_u )       Jac[ (0 + i) * nDof_AllVars   + 	(0 + j) ]  += jacXweight_qp * laplace_mat_du_u_i_j;
// 		  if ( i < nDof_u && j < nDof_u )       Jac[ (0 + i) * nDof_AllVars   + 	(0 + j) ]  += jacXweight_qp * laplace_beltrami_mat_du_u_i_j; ///@todo On a flat domain, this must coincide with the standard Laplacian, so we can do a double check with this
            } // end phi_j loop
          } // endif assemble_matrix

        } // end phi_i loop
        
      } // end gauss point loop


 //========= VOLUME - END ==================
 
 
 //========= FROM LOCAL TO GLOBAL - BEGIN ==================   
    RES->add_vector_blocked(Res, l2GMap_AllVars);

    if (assembleMatrix) {
      JAC->add_matrix_blocked(Jac, l2GMap_AllVars, l2GMap_AllVars);
    }
 //========= FROM LOCAL TO GLOBAL - END ==================   
   
   
  } //end element loop for each process

  
  RES->close();

  if (assembleMatrix) JAC->close();

     //print JAC and RES to files
    const unsigned nonlin_iter = 0/*mlPdeSys->GetNonlinearIt()*/;
    assemble_jacobian< double, double >::print_global_jacobian(assembleMatrix, ml_prob, JAC, nonlin_iter);
    assemble_jacobian< double, double >::print_global_residual(ml_prob, RES, nonlin_iter);
  


  return;
}



 private:    
     
static void natural_loop_1d(const MultiLevelProblem *    ml_prob, 
                     const Mesh *                    msh,
                     const MultiLevelSolution *    ml_sol, 
                     const unsigned iel,
                     CurrentElem < double > & geom_element,
                     const unsigned xType,
                     const std::string solname_u,
                     const unsigned solFEType_u,
                     std::vector< double > & Res
                    ) {
    
     double grad_u_dot_n = 0.;
    
    for (unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
        
       geom_element.set_coords_at_dofs_bdry_3d(iel, jface, xType);
       
       geom_element.set_elem_center_bdry_3d();
       
       std::vector <  double > xx_face_elem_center(3, 0.); 
          xx_face_elem_center = geom_element.get_elem_center_bdry_3d();
        
       const int boundary_index = msh->GetMeshElements()->GetFaceElementIndex(iel, jface);
       
       if ( boundary_index < 0) { //I am on the boundary
                  
         unsigned int face = - (boundary_index + 1);
    
         bool is_dirichlet =  ml_sol->GetBdcFunctionMLProb()(ml_prob, xx_face_elem_center, solname_u.c_str(), grad_u_dot_n, face, 0.);                     
         //we have to be careful here, because in GenerateBdc those coordinates are passed as NODE coordinates, 
         //while here we pass the FACE ELEMENT CENTER coordinates. 
         // So, if we use this for enforcing space-dependent Dirichlet or Neumann values, we need to be careful!
         
             if ( !(is_dirichlet)  &&  (grad_u_dot_n != 0.) ) {  //dirichlet == false and nonhomogeneous Neumann
                 
                 
                 
                   unsigned n_dofs_face = msh->GetElementFaceDofNumber(iel, jface, solFEType_u);

                  for (unsigned i = 0; i < n_dofs_face; i++) {
                      
                 unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i);

                 Res[i_vol] +=  grad_u_dot_n /* * phi[node] = 1. */;
                 
                         }
                         
                         
                         
        
                    }
                  
              }

    }
    
}




template < class real_num, class real_num_mov >
static void natural_loop_2d3d(const MultiLevelProblem *    ml_prob, 
                       const Mesh *                    msh,
                       const MultiLevelSolution *    ml_sol, 
                       const unsigned iel,
                       CurrentElem < double > & geom_element,
                       const unsigned solType_coords,
                       const std::string solname_u,
                       const unsigned solFEType_u,
                       std::vector< double > & Res,
                       //-----------
                       std::vector < std::vector < /*const*/ elem_type_templ_base<real_num, real_num_mov> *  > >  elem_all,
                       const unsigned dim,
                       const unsigned space_dim,
                       const unsigned max_size
                    ) {
    
    
    /// @todo - should put these outside the iel loop --
    std::vector < std::vector < double > >  JacI_iqp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_iqp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_iqp_bdry.size(); d++) {   Jac_iqp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp_bdry.size(); d++) { JacI_iqp_bdry[d].resize(dim-1); }
    
    double detJac_iqp_bdry;
  double weight_iqp_bdry = 0.;
// ---
  //boundary state shape functions
  std::vector <double> phi_u_bdry;  
  std::vector <double> phi_u_x_bdry; 

  phi_u_bdry.reserve(max_size);
  phi_u_x_bdry.reserve(max_size * space_dim);
// ---
  

// ---
  std::vector <double> phi_coords_bdry;
  std::vector <double> phi_coords_x_bdry;

  phi_coords_bdry.reserve(max_size);
  phi_coords_x_bdry.reserve(max_size * space_dim);
// ---
   


     double grad_u_dot_n = 0.;
    
    for (unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
        
       geom_element.set_coords_at_dofs_bdry_3d(iel, jface, solType_coords);
       
       geom_element.set_elem_center_bdry_3d();

       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, jface);    
       

       const int boundary_index = msh->GetMeshElements()->GetFaceElementIndex(iel, jface);
       
       if ( boundary_index < 0) { //I am on the boundary
                  
       std::vector <  double > xx_face_elem_center(3, 0.); 
       xx_face_elem_center = geom_element.get_elem_center_bdry_3d();
        
         const unsigned int face = - (boundary_index + 1);
    
         bool is_dirichlet =  ml_sol->GetBdcFunctionMLProb()(ml_prob, xx_face_elem_center, solname_u.c_str(), grad_u_dot_n, face, 0.);                     
         //we have to be careful here, because in GenerateBdc those coordinates are passed as NODE coordinates, 
         //while here we pass the FACE ELEMENT CENTER coordinates. 
         // So, if we use this for enforcing space-dependent Dirichlet or Neumann values, we need to be careful!
         
             if ( !(is_dirichlet) /* &&  (grad_u_dot_n != 0.)*/ ) {  //dirichlet == false and nonhomogeneous Neumann

    unsigned n_dofs_face_u = msh->GetElementFaceDofNumber(iel, jface, solFEType_u);

// dof-based - BEGIN
     std::vector< double > grad_u_dot_n_at_dofs(n_dofs_face_u);


    for (unsigned i_bdry = 0; i_bdry < grad_u_dot_n_at_dofs.size(); i_bdry++) {
        std::vector<double> x_at_node(dim, 0.);
        for (unsigned jdim = 0; jdim < x_at_node.size(); jdim++) x_at_node[jdim] = geom_element.get_coords_at_dofs_bdry_3d()[jdim][i_bdry];

      double grad_u_dot_n_at_dofs_temp = 0.;
      ml_sol->GetBdcFunctionMLProb()(ml_prob, x_at_node, solname_u.c_str(), grad_u_dot_n_at_dofs_temp, face, 0.);
     grad_u_dot_n_at_dofs[i_bdry] = grad_u_dot_n_at_dofs_temp;
      
    }

// dof-based - END
    
               
               

                        const unsigned n_gauss_bdry = ml_prob->GetQuadratureRule(ielGeom_bdry).GetGaussPointsNumber();
        
    
		for(unsigned ig_bdry = 0; ig_bdry < n_gauss_bdry; ig_bdry++) {
  
     elem_all[ielGeom_bdry][solType_coords]->JacJacInv(geom_element.get_coords_at_dofs_bdry_3d(), ig_bdry, Jac_iqp_bdry, JacI_iqp_bdry, detJac_iqp_bdry, space_dim);
//      elem_all[ielGeom_bdry][solType_coords]->compute_normal(Jac_iqp_bdry, normal);
    
    weight_iqp_bdry = detJac_iqp_bdry * ml_prob->GetQuadratureRule(ielGeom_bdry).GetGaussWeightsPointer()[ig_bdry];
                
    elem_all[ielGeom_bdry][solFEType_u ]->shape_funcs_current_elem(ig_bdry, JacI_iqp_bdry, phi_u_bdry, phi_u_x_bdry,  boost::none, space_dim);



//---------------------------------------------------------------------------------------------------------

     elem_all[ielGeom_bdry][solType_coords ]->shape_funcs_current_elem(ig_bdry, JacI_iqp_bdry, phi_coords_bdry, phi_coords_x_bdry,  boost::none, space_dim);

  std::vector<double> x_qp_bdry(dim, 0.);

         for (unsigned i = 0; i < phi_coords_bdry.size(); i++) {
           	for (unsigned d = 0; d < dim; d++) {
 	                                                x_qp_bdry[d]    += geom_element.get_coords_at_dofs_bdry_3d()[d][i] * phi_coords_bdry[i]; // fetch of coordinate points
             }
         }
         
           double grad_u_dot_n_qp = 0.;  ///@todo here we should do a function that provides the gradient at the boundary, and then we do "dot n" with the normal at qp
 
// dof-based
         for (unsigned i_bdry = 0; i_bdry < phi_u_bdry.size(); i_bdry ++) {
           grad_u_dot_n_qp +=  grad_u_dot_n_at_dofs[i_bdry] * phi_u_bdry[i_bdry];
         } 

// quadrature point based         
 // // // ml_sol->GetBdcFunctionMLProb()(ml_prob, x_qp_bdry, solname_u.c_str(), grad_u_dot_n_qp, face, 0.);

//---------------------------------------------------------------------------------------------------------





                  for (unsigned i_bdry = 0; i_bdry < n_dofs_face_u; i_bdry++) {

                 unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);

                 Res[i_vol] +=  weight_iqp_bdry * grad_u_dot_n_qp /*grad_u_dot_n*/  * phi_u_bdry[i_bdry];

                           }
                         
                         
                         
                        }
        
                         
        
                    }
                  
              }
    }
    
}




};


} //end namespace femus



#endif
