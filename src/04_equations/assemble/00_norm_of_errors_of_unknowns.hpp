 
#ifndef __femus_00_norms_of_errors_of_unknowns_hpp__
#define __femus_00_norms_of_errors_of_unknowns_hpp__




namespace femus   {


template < class real_num, class real_num_mov >
void compute_norm_of_errors_of_unknowns(MultiLevelProblem& ml_prob) {

  NonLinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystem> (ml_prob.get_app_specs_pointer()->_system_name);  
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);
  elem*                     el = msh->el;

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));

  unsigned    iproc = msh->processor_id();

  //=============== Geometry ========================================
  unsigned xType = BIQUADR_FE; // the FE for the domain need not be biquadratic
  
  CurrentElem < double > geom_element(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
//***************************************************  


//***************** Check Number of Unknowns **********************************  
  const int n_vars = mlPdeSys->GetSolPdeIndex().size();
  if (n_vars != 1)  { std::cout << "Function written for only 1 scalar unknown";  abort();  }

 //******************** quadrature *******************************  
  double jacXweight_qp; 

 //********************* quadrature, geometry *********************** 
 //***************************************************  
  std::vector <double> phi_coords;
  std::vector <double> phi_coords_x; 
  std::vector <double> phi_coords_xx;

  phi_coords.reserve(maxSize);
  phi_coords_x.reserve(maxSize * space_dim);
  phi_coords_xx.reserve(maxSize * dim2);
  
  
 //********************* quadrature, unknowns *********************** 
 //***************************************************  
  std::vector <double> phi_u;
  std::vector <double> phi_u_x; 
  std::vector <double> phi_u_xx;

  phi_u.reserve(maxSize);
  phi_u_x.reserve(maxSize * space_dim);
  phi_u_xx.reserve(maxSize * dim2);
  
  const std::string solname_u = ml_sol->GetSolName_string_vec()[0];
  unsigned solIndex_u = ml_sol->GetIndex(solname_u.c_str()); 
  unsigned solFEType_u = ml_sol->GetSolutionType(solIndex_u); 

  unsigned solPdeIndex_u;
  solPdeIndex_u = mlPdeSys->GetSolPdeIndex(solname_u.c_str());

  std::vector < double >  sol_u;     sol_u.reserve(maxSize);
  std::vector < double >  sol_u_true;     sol_u_true.reserve(maxSize);
 //***************************************************  
 //***************************************************  


  
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
  
  double  norm_dof_based = 0.;
  double  norm_qp_based = 0.;
  

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
      

    geom_element.set_coords_at_dofs_and_geom_type(iel, xType);
        
    const short unsigned ielGeom = geom_element.geom_type();

 //**************** state **************************** 
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solFEType_u);
    sol_u    .resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
     unsigned solDof_u = msh->GetSolutionDof(i, iel, solFEType_u);
      sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);
    }
    
 //**************** true sol **************************** 
    sol_u_true    .resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
        std::vector< double > xyz_i(dim);
             	for (unsigned d = 0; d < xyz_i.size(); d++) {
                   xyz_i[d] = geom_element.get_coords_at_dofs_3d()[d][i];
                }
               sol_u_true[i]  =   ml_prob.get_app_specs_pointer()->_norm_true_solution( xyz_i);  
    }

    //***************************************************  
 
    
 //*************************************************** 
    

// // //  //========= BOUNDARY ==================   
// // //     if (dim == 1)   ml_prob.get_app_specs_pointer()->_assemble_function_natural_boundary_loop_1d(& ml_prob, msh, ml_sol,
// // //                       iel, geom_element, xType,
// // //                       solname_u, solFEType_u,
// // //                       Res
// // //                      );
// // // 
// // //     if (dim == 2 || dim == 3)   ml_prob.get_app_specs_pointer()->_assemble_function_natural_boundary_loop_2d3d(& ml_prob, msh, ml_sol,
// // //                       iel, geom_element, xType,
// // //                       solname_u, solFEType_u,
// // //                       Res,
// // //                       elem_all,
// // //                       dim,
// // //                       space_dim,
// // //                       maxSize
// // //                      );
 
 //========= VOLUME ==================   
   
 //========= gauss value quantities ==================   
	std::vector<double> sol_u_x_gss(space_dim);     std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
 //===================================================   
    
    
    
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
	
	for (unsigned i = 0; i < sol_u.size(); i++) {
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
 ///   (x_nodes) are obtained from   geom_element.get_coords_at_dofs_3d()  (this is a  vector< vector >,  where the outer index is the dimension and the inner index ranges over the nodes) 
 ///   (shape function of that node, evaluated at qp)  is obtained from phi_u  (this is a vector, whose index ranges over the nodes)
 
 std::vector<double> x_qp(dim, 0.);
          
        for (unsigned i = 0; i < phi_coords.size(); i++) {
          	for (unsigned d = 0; d < dim; d++) {
	                                                x_qp[d]    += geom_element.get_coords_at_dofs_3d()[d][i] * phi_coords[i]; // fetch of coordinate points
            }
        }
 
//--------------    

          
	// *** phi_i loop ***
	double u_qp = 0.;        for (unsigned i = 0; i < phi_u.size(); i++) {            u_qp +=  sol_u[i] * phi_u[i];        }

 	double u_true_qp = 0.;   for (unsigned i = 0; i < phi_u.size(); i++) {            u_true_qp +=  sol_u_true[i] * phi_u[i];        }
       
  
	      
       norm_dof_based +=  jacXweight_qp * (  u_qp - u_true_qp) * (  u_qp - u_true_qp) ;
       norm_qp_based  +=  jacXweight_qp * (  u_qp -  ml_prob.get_app_specs_pointer()->_norm_true_solution( x_qp )) * (  u_qp - ml_prob.get_app_specs_pointer()->_norm_true_solution( x_qp )) ;
//======================Residuals=======================
	      
// // //           if (assembleMatrix) {
// // // 	    
// // //             // *** phi_j loop ***
// // //             for (unsigned j = 0; j < nDof_max; j++) {
// // // 
// // // 
// // // //--------------    
// // //               double laplace_mat_du_u_i_j = 0.;
// // // 
// // //               if ( i < nDof_u && j < nDof_u ) {
// // //                   for (unsigned kdim = 0; kdim < space_dim; kdim++) {
// // //                          laplace_mat_du_u_i_j           += phi_u_x   [i * space_dim + kdim] *
// // //                                                            phi_u_x   [j * space_dim + kdim];
// // // 	              }
// // //              }
// // // //--------------    
// // // 
// // // 
// // // 
// // //               //============ delta_state row ============================
// // //               //DIAG BLOCK delta_state - state
// // // 		  if ( i < nDof_u && j < nDof_u )       Jac[ (0 + i) * nDof_AllVars   + 	(0 + j) ]  += jacXweight_qp * laplace_mat_du_u_i_j;
// // // // 		  if ( i < nDof_u && j < nDof_u )       Jac[ (0 + i) * nDof_AllVars   + 	(0 + j) ]  += jacXweight_qp * laplace_beltrami_mat_du_u_i_j; ///@todo On a flat domain, this must coincide with the standard Laplacian, so we can do a double check with this
// // //             } // end phi_j loop
// // //           } // endif assemble_matrix

        
      } // end gauss point loop


   
  } //end element loop for each process

   double norm_dof_based_parallel = 0.; MPI_Allreduce( &norm_dof_based, &norm_dof_based_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   double norm_qp_based_parallel = 0.; MPI_Allreduce( &norm_qp_based, &norm_qp_based_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
 
  std::cout << std::scientific << std::setw(20) << std::setprecision(15) << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&Norm dof: " << norm_dof_based_parallel << std::endl;
  
  std::cout << std::scientific << std::setw(20) << std::setprecision(15) << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&Norm qp: " << norm_qp_based_parallel << std::endl;



  return;
}




} //end namespace femus




#endif
