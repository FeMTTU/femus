#ifndef __femus_00_laplacian_with_all_dirichlet_bc_ad_or_not_hpp__
#define __femus_00_laplacian_with_all_dirichlet_bc_ad_or_not_hpp__


#include "MultiLevelProblem.hpp"
#include "Assemble_jacobian.hpp"
#include "LinearEquationSolver.hpp"

#include "adept.h"


namespace femus {
    
    
/**
 * This function assemble the stiffnes matrix Jac and the residual vector Res
 * such that
 *                  Jac w = RES = F - Jac u0,
 * and consequently
 *        u = u0 + w satisfies Jac u = F
 **/

void AssemblePoissonProblem_old_fe_quadrature_nonAD(MultiLevelProblem& ml_prob,
                                                    double    (* right_hand_side )  (const std::vector<double> & ),
                                                    double    (* exact_solution )  (const std::vector<double> & )
                                                    ) {
  //  ml_prob is the global object from/to where get/set all the data

  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use

  
    // Ia) System, OLD WAY
//   LinearImplicitSystem * mlPdeSys  = & ml_prob.get_system< LinearImplicitSystem > ("Poisson");   // pointer to the linear implicit system named "Poisson"
    
    // Ib) System, NEW WAY: that we are currently solving (System that is calling this function)
    const unsigned current_system_number = ml_prob.get_current_system_number();
    LinearImplicitSystem * mlPdeSys  = & ml_prob.get_system< LinearImplicitSystem >(current_system_number);   // pointer to the linear implicit system named "Poisson"
  
    // II) Unknowns of the System
  std::vector< Unknown >   unknowns = ml_prob.get_system< LinearImplicitSystem >(current_system_number).get_unknown_list_for_assembly();

  
  if (unknowns.size() != 1) abort();
  
  
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->GetMeshElements();  // pointer to the elem object in msh (level)

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));  // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex = ml_sol->GetIndex( unknowns[0]._name.c_str() );    // get the position of "u" in the ml_sol object
  unsigned soluType = ml_sol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned soluPdeIndex = mlPdeSys->GetSolPdeIndex(  unknowns[0]._name.c_str() );    // get the position of "u" in the pdeSys object

  std::vector < double >  solu; // local solution
  solu.reserve(maxSize);

  std::vector < double >  solu_exact_at_dofs;  solu_exact_at_dofs.reserve(maxSize);

  
  std::vector < std::vector < double > > x (dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE BI/TRIQUADRATIC)

  for (unsigned i = 0; i < dim; i++) {
    x[i].reserve(maxSize);
  }

  std::vector <double> phi;  // local test function
  std::vector <double> phi_x; // local test function first order partial derivatives
  
  double weight = 0.; // gauss point weight

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  

  std::vector < double > Res; // local redidual vector
  Res.reserve(maxSize);

  std::vector < double > Jac; //local Jacobian matrix
  Jac.reserve(maxSize * maxSize);
  
  std::vector < int > l2GMap; // local to global mapping
  l2GMap.reserve(maxSize);
  

  KK->zero(); // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual Vector

  
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
    

    std::vector<unsigned> Sol_n_el_dofs_Mat_vol(1, nDofu);
  
    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofx);
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->GetTopology()->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }


    // resize local arrays
    l2GMap.resize(nDofu);
    
    solu.resize(nDofu);
    solu_exact_at_dofs.resize(nDofu);

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
        std::vector<double> x_at_node(dim, 0.);
        for (unsigned jdim = 0; jdim < dim; jdim++) x_at_node[jdim] = x[jdim][i];
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // local to global solution mapping
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // local storage of solution
      solu_exact_at_dofs[i] = exact_solution(x_at_node);
      l2GMap[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);   // local to global system solution mapping
    }
    
    


    Res.assign(nDofu, 0.);    //resize and set to zero
    Jac.assign(nDofu * nDofu, 0.);    //resize and set to zero
    

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
        
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, boost::none);
      
      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
    
      std::vector < double > gradSolu_gss(dim, 0.);
      std::vector < double > gradSolu_exact_gss(dim, 0.);
      std::vector < double > x_gss(dim, 0.);

      for (unsigned i = 0; i < nDofu; i++) {
       
        for (unsigned jdim = 0; jdim < dim; jdim++) {
          gradSolu_gss[jdim]       += phi_x[i * dim + jdim] * solu[i];
          gradSolu_exact_gss[jdim] += phi_x[i * dim + jdim] * solu_exact_at_dofs[i];
          x_gss[jdim] += x[jdim][i] * phi[i];
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofu; i++) {

        double weakLaplace = 0.;
        double laplace_weak_exact = 0.;

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          weakLaplace        +=  phi_x[i * dim + jdim] * gradSolu_gss[jdim];
          laplace_weak_exact +=  phi_x[i * dim + jdim] * gradSolu_exact_gss[jdim];
        }
        
        Res[i] += ( -  right_hand_side(x_gss) * phi[i] - weakLaplace) * weight;
        // Res[i] += (laplace_weak_exact                  - weakLaplace) * weight;

        // *** phi_j loop ***
        for (unsigned j = 0; j < nDofu; j++) {
          double weakLaplacej = 0.;

          for (unsigned kdim = 0; kdim < dim; kdim++) {
            weakLaplacej += phi_x[i * dim + kdim] * phi_x[j * dim + kdim];
          }

          Jac[i * nDofu + j] +=  weakLaplacej * weight;
        } // end phi_j loop

      } // end phi_i loop
      
    } // end gauss point loop


    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    RES->add_vector_blocked(Res, l2GMap);

    //store K in the global matrix KK
    KK->add_matrix_blocked(Jac, l2GMap, l2GMap);

    
     constexpr bool print_algebra_local = false;
     if (print_algebra_local) {
         assemble_jacobian<double,double>::print_element_residual(iel, Res, Sol_n_el_dofs_Mat_vol, 10, 5);
         assemble_jacobian<double,double>::print_element_jacobian(iel, Jac, Sol_n_el_dofs_Mat_vol, 10, 5);
     }  
    
  } //end element loop for each process

  RES->close();

  KK->close();

  // ***************** END ASSEMBLY *******************
}


/**
 * This function assemble the stiffnes matrix KK and the residual vector Res
 * Using automatic differentiation for Newton iterative scheme
 *                  J(u0) w =  - F(u0)  ,
 *                  with u = u0 + w
 *                  - F = f(x) - J u = Res
 *                  J = \grad_u F
 *
 * thus
 *                  J w = f(x) - J u0
 **/
void AssemblePoissonProblem_old_fe_quadrature_AD(MultiLevelProblem& ml_prob,
                                                    double    (* right_hand_side )  (const std::vector<double> & ),
                                                    double    (* exact_solution )  (const std::vector<double> & )
                                                    ) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object


  adept::Stack& s = FemusInit::_adeptStack;


      // Ia) System, OLD WAY
//   LinearImplicitSystem * mlPdeSys  = & ml_prob.get_system< LinearImplicitSystem > ("Poisson");   // pointer to the linear implicit system named "Poisson"
    
    // Ib) System, NEW WAY: that we are currently solving (System that is calling this function)
    const unsigned current_system_number = ml_prob.get_current_system_number();
    LinearImplicitSystem * mlPdeSys  = & ml_prob.get_system< LinearImplicitSystem >(current_system_number);   // pointer to the linear implicit system named "Poisson"
  
    // II) Unknowns of the System
  std::vector< Unknown >   unknowns = ml_prob.get_system< LinearImplicitSystem >(current_system_number).get_unknown_list_for_assembly();


  if (unknowns.size() != 1) abort();
  
  
  //  extract pointers to the several objects that we are going to use

 
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->GetMeshElements();  // pointer to the elem object in msh (level)

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex = ml_sol->GetIndex( unknowns[0]._name.c_str() );    // get the position of "u" in the ml_sol object
  unsigned soluType = ml_sol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned soluPdeIndex = mlPdeSys->GetSolPdeIndex(  unknowns[0]._name.c_str() );    // get the position of "u" in the pdeSys object

  std::vector < adept::adouble >  solu; // local solution
  solu.reserve(maxSize);
  std::vector < double >  solu_exact_at_dofs;  solu_exact_at_dofs.reserve(maxSize);

  std::vector < std::vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned i = 0; i < dim; i++) {
    x[i].reserve(maxSize);
  }

  std::vector <double> phi;  // local test function
  std::vector <double> phi_x; // local test function first order partial derivatives

  double weight = 0.; // gauss point weight

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);

  std::vector < adept::adouble > aRes; // local redidual vector
  aRes.reserve(maxSize);

  std::vector < int > l2GMap; // local to global mapping
  l2GMap.reserve(maxSize);
  std::vector < double > Res; // local redidual vector
  Res.reserve(maxSize);
  std::vector < double > Jac;
  Jac.reserve(maxSize * maxSize);

  RES->zero(); // Set to zero all the entries of the Global Residual Vector
  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {
     
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    std::vector<unsigned> Sol_n_el_dofs_Mat_vol(1, nDofu);

    

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofx);
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->GetTopology()->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }


    // resize local arrays
    l2GMap.resize(nDofu);
    solu.resize(nDofu);
    solu_exact_at_dofs.resize(nDofu);

     // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
        std::vector<double> x_at_node(dim, 0.);
        for (unsigned jdim = 0; jdim < dim; jdim++) x_at_node[jdim] = x[jdim][i];
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      solu_exact_at_dofs[i] = exact_solution(x_at_node);
      l2GMap[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }


   aRes.resize(nDofu);
    std::fill(aRes.begin(), aRes.end(), 0.);    //set aRes to zero

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
        
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, boost::none);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      std::vector < adept::adouble > gradSolu_gss(dim, 0.);
      std::vector < adept::adouble > gradSolu_exact_gss(dim, 0.);
      
      std::vector < double > x_gss(dim, 0.);

      for (unsigned i = 0; i < nDofu; i++) {

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          gradSolu_gss[jdim]       += phi_x[i * dim + jdim] * solu[i];
          gradSolu_exact_gss[jdim] += phi_x[i * dim + jdim] * solu_exact_at_dofs[i];
          x_gss[jdim] += x[jdim][i] * phi[i];
        }
      }
      
      

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofu; i++) {

        adept::adouble laplace = 0.;
        adept::adouble laplace_weak_exact = 0.;

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          laplace            +=  phi_x[i * dim + jdim] * gradSolu_gss[jdim];
          laplace_weak_exact +=  phi_x[i * dim + jdim] * gradSolu_exact_gss[jdim];
        }

        aRes[i] +=  ( right_hand_side(x_gss) * phi[i] + laplace ) * weight;    ///@todo check the sign here
        // aRes[i] += (- laplace_weak_exact             + laplace) * weight;   ///@todo check the sign here

      } // end phi_i loop
      
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    Res.resize(nDofu);    //resize

    for (int i = 0; i < nDofu; i++) {
      Res[i] = - aRes[i].value();
    }

    RES->add_vector_blocked(Res, l2GMap);



    // define the dependent variables
    s.dependent(&aRes[0], nDofu);

    // define the independent variables
    s.independent(&solu[0], nDofu);

    // get the jacobian matrix (ordered by row major )
    Jac.resize(nDofu * nDofu);    //resize
    s.jacobian(&Jac[0], true);

    //store K in the global matrix KK
    KK->add_matrix_blocked(Jac, l2GMap, l2GMap);

    s.clear_independents();
    s.clear_dependents();

    
  constexpr bool print_algebra_local = false;
     if (print_algebra_local) {
         assemble_jacobian<double,double>::print_element_residual(iel, Res, Sol_n_el_dofs_Mat_vol, 10, 5);
         assemble_jacobian<double,double>::print_element_jacobian(iel, Jac, Sol_n_el_dofs_Mat_vol, 10, 5);
     }
     
     
  } //end element loop for each process

  RES->close();

  KK->close();

  // ***************** END ASSEMBLY *******************
}





}

#endif
