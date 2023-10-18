#ifndef __femus_00_norms_of_errors_of_unknowns_hpp__
#define __femus_00_norms_of_errors_of_unknowns_hpp__


#include "Assemble_unknown.hpp"
#include "NonLinearImplicitSystem.hpp"


namespace femus   {


template < class real_num, class real_num_mov >
void compute_L2_norm_of_errors_of_unknowns(MultiLevelProblem& ml_prob) {

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


//***************** Check Number of Unknowns **********************************  
  const int n_vars = mlPdeSys->GetSolPdeIndex().size();
  if (n_vars != 1)  { std::cout << "Function written for only 1 scalar unknown";  abort();  }

//***************** Check that true solution is provided **********************************  
  if (ml_prob.get_app_specs_pointer()->_true_solution == NULL) {  std::cout << "No true solution provided";  abort();  }


  //=============== Geometry ========================================
  unsigned xType = BIQUADR_FE; // the FE for the domain need not be biquadratic
  
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

  unsigned solPdeIndex_u = mlPdeSys->GetSolPdeIndex(solname_u.c_str());

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
    //***************************************************  
    
 //**************** true sol **************************** 
    sol_u_true    .resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
        std::vector< double > xyz_i(dim);
             	for (unsigned d = 0; d < xyz_i.size(); d++) {
                   xyz_i[d] = geom_element.get_coords_at_dofs_3d()[d][i];
                }
               sol_u_true[i]  =   ml_prob.get_app_specs_pointer()->_true_solution( xyz_i);  
    }

    //***************************************************  
 
    
 //*************************************************** 
    

 
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
       norm_qp_based  +=  jacXweight_qp * (  u_qp -  ml_prob.get_app_specs_pointer()->_true_solution( x_qp )) * (  u_qp - ml_prob.get_app_specs_pointer()->_true_solution( x_qp )) ;
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
 
  std::cout << std::scientific << std::setw(20) << std::setprecision(15) << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& L2-Norm of the error, dof-based: " << norm_dof_based_parallel << std::endl;
  
  std::cout << std::scientific << std::setw(20) << std::setprecision(15) << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& L2-Norm of the error, qp-based " << norm_qp_based_parallel << std::endl;



  return;
}




//This function does not seem to require an equation
std::pair < double, double > GetErrorNorm_L2_H1_with_analytical_sol(MultiLevelSolution* ml_sol, 
                                                                    std::string solution_name,
                                                                    double    (* function_value )  (const std::vector<double> & ),
                                                                    void    (* function_gradient)  (const std::vector < double > & , std::vector < double >&  )
                                                                    ) {      
  
  unsigned level = ml_sol->_mlMesh->GetNumberOfLevels() - 1u;
  //  extract pointers to the several objects that we are going to use
  Mesh*     msh = ml_sol->_mlMesh->GetLevel(level);    // pointer to the mesh (level) object
  elem*     el  = msh->el;  // pointer to the elem object in msh (level)
  Solution* sol = ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex = ml_sol->GetIndex( solution_name.c_str() );    // get the position of "u" in the ml_sol object
  unsigned soluType = ml_sol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  vector < double >  solu; // local solution

  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives

  double weight; // gauss point weight

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  solu.reserve(maxSize);

  for (unsigned i = 0; i < dim; i++)
    x[i].reserve(maxSize);

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)

  double seminorm = 0.;
  double l2norm = 0.;

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    // resize local arrays
    solu.resize(nDofu);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofx);
    }

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);  // global extraction and local storage for the element coordinates
      }
    }

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
        
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, boost::none);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      double solu_gss = 0;
      vector < double > gradSolu_gss(dim, 0.);
      vector < double > x_gss(dim, 0.);

      for (unsigned i = 0; i < nDofu; i++) {
        solu_gss += phi[i] * solu[i];

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solu[i];
          x_gss[jdim] += x[jdim][i] * phi[i];
        }
      }

      vector <double> exactGradSol(dim);
      function_gradient(x_gss, exactGradSol);

      for (unsigned j = 0; j < dim ; j++) {
        seminorm   += ((gradSolu_gss[j] - exactGradSol[j]) * (gradSolu_gss[j] - exactGradSol[j])) * weight;
      }

      double exactSol = function_value(x_gss);
      l2norm += (exactSol - solu_gss) * (exactSol - solu_gss) * weight;
    } // end gauss point loop
  } //end element loop for each process

  // add the norms of all processes
  NumericVector* norm_vec;
  norm_vec = NumericVector::build().release();
  norm_vec->init(msh->n_processors(), 1 , false, AUTOMATIC);

  norm_vec->set(iproc, l2norm);
  norm_vec->close();
  l2norm = norm_vec->l1_norm();

  norm_vec->set(iproc, seminorm);
  norm_vec->close();
  seminorm = norm_vec->l1_norm();

  delete norm_vec;

  return std::pair < double, double > (sqrt(l2norm), sqrt(seminorm));

}


std::pair < double, double > GetErrorNorm_L2_H1_with_analytical_sol(MultiLevelSolution* ml_sol, 
                                                                    std::vector< Unknown > & unknowns_vec,
                                                                    double    (* function_value )  (const std::vector<double> & ),
                                                                    void    (* function_gradient)  (const std::vector < double > & , std::vector < double >&  )
                                                                    ){
   if (unknowns_vec.size() != 1) abort();
                                                                     
 return   GetErrorNorm_L2_H1_with_analytical_sol(ml_sol, unknowns_vec[0]._name, function_value, function_gradient) ;                                                                       
                                                                      
                                                                    }



// ||u_i - u_h||/||u_i-u_(h/2)|| = 2^alpha, alpha is order of conv 
std::pair < double, double > GetErrorNorm_L2_H1_multiple_methods(MultiLevelSolution* ml_sol,
                                          Solution* sol_finer,
                                          std::vector< Unknown > & unknowns_vec,
                                                                    double    (* function_value )  (const std::vector<double> & ),
                                                                    void    (* function_gradient)  (const std::vector < double > & , std::vector < double >&  )
                                          ) {
  
  if (unknowns_vec.size() != 1) abort();
  
  unsigned level = ml_sol->_mlMesh->GetNumberOfLevels() - 1u;
  //  extract pointers to the several objects that we are going to use
  Mesh*     msh = ml_sol->_mlMesh->GetLevel(level);    // pointer to the mesh (level) object
  elem*     el  = msh->el;  // pointer to the elem object in msh (level)
  Solution* sol = ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex = ml_sol->GetIndex( unknowns_vec[0]._name.c_str() );    // get the position of "u" in the ml_sol object
  unsigned soluType = ml_sol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  vector < double >  solu; // local solution

  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives

  double weight; // gauss point weight

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  solu.reserve(maxSize);

  vector < double >  solu_finer;   solu_finer.reserve(maxSize);
  
  vector < double >  solu_exact_at_dofs;   solu_exact_at_dofs.reserve(maxSize);

  for (unsigned i = 0; i < dim; i++)
    x[i].reserve(maxSize);

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)

  double seminorm = 0.;
  double l2norm = 0.;
  double seminorm_exact_dofs = 0.;
  double l2norm_exact_dofs = 0.;
  double seminorm_inexact = 0.;
  double l2norm_inexact = 0.;

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    // resize local arrays
    solu.resize(nDofu);
    solu_finer.resize(nDofu);
    solu_exact_at_dofs.resize(nDofu);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofx);
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);  // global extraction and local storage for the element coordinates
      }
    }
    
    const double weird_multigrid_factor = 0.25;  //don't know!
    
         vector<double> x_at_node(dim,0.);
      for (unsigned i = 0; i < nDofu; i++) {
         for (unsigned jdim = 0; jdim < dim; jdim++) {
             x_at_node[jdim] = x[jdim][i];
         }
      }

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
//         std::vector<double> x_at_node(dim,0.);
//         for (unsigned jdim = 0; jdim < dim; jdim++) x_at_node[jdim] = x[jdim][i];
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i]       = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      solu_finer[i] = /*weird_multigrid_factor * */(*sol_finer->_Sol[soluIndex])(solDof);
      solu_exact_at_dofs[i] = function_value(x_at_node);
    }


    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
        
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, boost::none);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      double solu_gss = 0.;
      double solu_finer_gss = 0.;
      double exactSol_from_dofs_gss = 0.;
      vector < double > gradSolu_gss(dim, 0.);
      vector < double > gradSolu_exact_at_dofs_gss(dim, 0.);
      vector < double > gradSolu_finer_gss(dim, 0.);
      vector < double > x_gss(dim, 0.);

      for (unsigned i = 0; i < nDofu; i++) {
        solu_gss       += phi[i] * solu[i];
        solu_finer_gss += phi[i] * solu_finer[i];
        exactSol_from_dofs_gss += phi[i] * solu_exact_at_dofs[i];

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solu[i];
          gradSolu_exact_at_dofs_gss[jdim] += phi_x[i * dim + jdim] * solu_exact_at_dofs[i];
          gradSolu_finer_gss[jdim] += phi_x[i * dim + jdim] * solu_finer[i];
          x_gss[jdim] += x[jdim][i] * phi[i];
        }
      }

      vector <double> exactGradSol(dim);
      function_gradient(x_gss, exactGradSol);

      for (unsigned j = 0; j < dim ; j++) {
        seminorm           += ((gradSolu_gss[j]       - exactGradSol[j]) * (gradSolu_gss[j]       - exactGradSol[j])) * weight;
        seminorm_inexact   += ((gradSolu_gss[j] - gradSolu_finer_gss[j]) * (gradSolu_gss[j] - gradSolu_finer_gss[j])) * weight;
        seminorm_exact_dofs      += ((gradSolu_gss[j]       - gradSolu_exact_at_dofs_gss[j]) * (gradSolu_gss[j]       - gradSolu_exact_at_dofs_gss[j])) * weight;
      }

      double exactSol = function_value(x_gss);
      l2norm              += (solu_gss - exactSol)                * (solu_gss - exactSol)       * weight;
      l2norm_exact_dofs   += (solu_gss - exactSol_from_dofs_gss)  * (solu_gss - exactSol_from_dofs_gss)       * weight;
      l2norm_inexact      += (solu_gss - solu_finer_gss)          * (solu_gss - solu_finer_gss) * weight;
   } // end gauss point loop
  } //end element loop for each process

  // add the norms of all processes
  NumericVector* norm_vec;
  norm_vec = NumericVector::build().release();
  norm_vec->init(msh->n_processors(), 1 , false, AUTOMATIC);

  norm_vec->set(iproc, l2norm);
  norm_vec->close();
  l2norm = norm_vec->l1_norm();

  norm_vec->set(iproc, seminorm);
  norm_vec->close();
  seminorm = norm_vec->l1_norm();

  delete norm_vec;

   // add the norms of all processes
  NumericVector* norm_vec_exact_dofs;
  norm_vec_exact_dofs = NumericVector::build().release();
  norm_vec_exact_dofs->init(msh->n_processors(), 1 , false, AUTOMATIC);

  norm_vec_exact_dofs->set(iproc, l2norm_exact_dofs);
  norm_vec_exact_dofs->close();
  l2norm_exact_dofs = norm_vec_exact_dofs->l1_norm();

  norm_vec_exact_dofs->set(iproc, seminorm_exact_dofs);
  norm_vec_exact_dofs->close();
  seminorm_exact_dofs = norm_vec_exact_dofs->l1_norm();

  delete norm_vec_exact_dofs;

  // add the norms of all processes
  NumericVector* norm_vec_inexact;
  norm_vec_inexact = NumericVector::build().release();
  norm_vec_inexact->init(msh->n_processors(), 1 , false, AUTOMATIC);

  norm_vec_inexact->set(iproc, l2norm_inexact);
  norm_vec_inexact->close();
  l2norm_inexact = norm_vec_inexact->l1_norm();

  norm_vec_inexact->set(iproc, seminorm_inexact);
  norm_vec_inexact->close();
  seminorm_inexact = norm_vec_inexact->l1_norm();

  delete norm_vec_inexact;

  std::pair < double, double > inexact_pair(sqrt(l2norm_inexact), sqrt(seminorm_inexact));
  
  return std::pair < double, double > (sqrt(l2norm), sqrt(seminorm));
  // return std::pair < double, double > (sqrt(l2norm_exact_dofs), sqrt(seminorm_exact_dofs));  ///@todo does not seem to work
  // return std::pair < double, double > (sqrt(l2norm_inexact), sqrt(seminorm_inexact));        ///@todo does not seem to work

}



} //end namespace femus




#endif
