/** tutorial/Ex2
 * This example shows how to set and solve the weak form of the Poisson problem
 *                    $$ \Delta u = \Delta u_exact \text{ on }\Omega, $$
 *                    $$ u=0 \text{ on } \Gamma, $$
 * on a square domain $\Omega$ with boundary $\Gamma$;
 * all the coarse-level meshes are removed;
 * a multilevel problem and an equation system are initialized;
 * a direct solver is used to solve the problem.
 **/

#include "FemusInit.hpp"
#include "Files.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "LinearImplicitSystem.hpp"

#include "Assemble_jacobian.hpp"


#include "00_poisson_eqn_with_all_dirichlet_bc.hpp"



#include "../tutorial_common.hpp"

#include "../all_mesh_generation_methods.hpp"


#include "adept.h"


using namespace femus;








double GetExactSolutionValue(const std::vector < double >& x) {
  double pi = acos(-1.);
  return cos(pi * x[0]) * cos(pi * x[1]);
}


void GetExactSolutionGradient(const std::vector < double >& x, vector < double >& solGrad) {
  double pi = acos(-1.);
  solGrad[0]  = -pi * sin(pi * x[0]) * cos(pi * x[1]);
  solGrad[1] = -pi * cos(pi * x[0]) * sin(pi * x[1]);
}


double GetExactSolutionLaplace(const std::vector < double >& x) {
  double pi = acos(-1.);
  return -pi * pi * cos(pi * x[0]) * cos(pi * x[1]) - pi * pi * cos(pi * x[0]) * cos(pi * x[1]);
}




bool SetBoundaryCondition(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;

  return dirichlet;
}

void AssemblePoissonProblem_old_fe_quadrature(MultiLevelProblem & ml_prob);
void AssemblePoissonProblem_old_fe_quadrature_AD(MultiLevelProblem & ml_prob);

std::pair < double, double > GetErrorNorm(MultiLevelSolution* ml_sol, std::vector< Unknown > & unknowns_vec);






int main(int argc, char** args) {

  // ======= Init ========================
  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // ======= Problem  ==================
  MultiLevelProblem ml_prob;

  
  // ======= Files - BEGIN  ========================
  const bool use_output_time_folder = false;
  const bool redirect_cout_to_file = true;
  Files files; 
        files.CheckIODirectories(use_output_time_folder);
        files.RedirectCout(redirect_cout_to_file);

  // ======= Problem, Files ========================
  ml_prob.SetFilesHandler(&files);
  // ======= Files - END  ========================

  
  
  // ======= Quad Rule - BEGIN ========================
   std::string fe_quad_rule("seventh");
    // ======= Quad Rule - END ========================
 
   
  // ======= Mesh file types (function, salome, gambit) ========================
    for (unsigned  mesh_file_type = 0; mesh_file_type < 3; mesh_file_type++) {  
  
  // define multilevel mesh
  MultiLevelMesh ml_mesh;

  // read coarse level mesh
   std::string mesh_name = Domain_square_m05p05::quad9_all_mesh_generation_methods_Structured(mesh_file_type, ml_mesh, fe_quad_rule);


   
  unsigned dim = ml_mesh.GetDimension();
  unsigned maxNumberOfMeshes;

  if (dim == 1 || dim == 2) {
    maxNumberOfMeshes = 6;
  } else {
    maxNumberOfMeshes = 4;
  }

  vector < vector < double > > l2Norm;
  l2Norm.resize(maxNumberOfMeshes);

  vector < vector < double > > semiNorm;
  semiNorm.resize(maxNumberOfMeshes);

      
  // ======= Assemble methods (AD or NON-AD) ========================
    std::vector < std::pair <femus::System::AssembleFunctionType, std::string> >  assemble_pointer_vec(2);
    
    assemble_pointer_vec[0].first = AssemblePoissonProblem_old_fe_quadrature;
    assemble_pointer_vec[0].second = "non-automatic_diff";
    
    assemble_pointer_vec[1].first = AssemblePoissonProblem_old_fe_quadrature_AD;
    assemble_pointer_vec[1].second = "automatic_diff";

    
    
   for (unsigned func = 0; func < assemble_pointer_vec.size(); func++) {   // loop on the mesh level

  std::cout << std::endl;
  std::cout << std::endl;
     std::cout << "===================================" << std::endl;
     std::cout << "===================================" << std::endl;
     std::cout << assemble_pointer_vec[func].second << std::endl;
     std::cout << "===================================" << std::endl;
     std::cout << "===================================" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;


  
  // ======= Mesh refinements ========================
  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {   // loop on the mesh level

    unsigned numberOfUniformLevels = i + 1;
    unsigned numberOfSelectiveLevels = 0;
    ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

    // erase all the coarse mesh levels
    ml_mesh.EraseCoarseLevels(numberOfUniformLevels - 1);

    // print mesh info
    ml_mesh.PrintInfo();

    std::vector< Unknown > unknowns = systems__generate_list_of_scalar_unknowns_for_each_FE_family_lagrangian();
    
    l2Norm[i].resize(unknowns.size());
    semiNorm[i].resize(unknowns.size());

      
      
  // ======= FE SPACES ========================
    for (unsigned int u = 0; u < unknowns.size(); u++) {

        // define the multilevel solution and attach the ml_mesh object to it
      MultiLevelSolution ml_sol(&ml_mesh);

      ml_sol.SetWriter(VTK);


      // add variables to ml_sol
      ml_sol.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order, unknowns[u]._time_order, unknowns[u]._is_pde_unknown);

      ml_sol.Initialize("All");

      // attach the boundary condition function and generate boundary data
      ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
      ml_sol.GenerateBdc( unknowns[u]._name.c_str() );

      // ======= Problem, Mesh and Solution  ==================
      // attach the ml_sol object
       ml_prob.SetMultiLevelMeshAndSolution(& ml_sol);
       

      ml_prob.get_systems_map().clear();
       ml_prob.set_current_system_number(0/*u*/);               //way to communicate to the assemble function, which doesn't belong to any class
      
    // ======= System - BEGIN ========================
       std::string sys_name = "Poisson" + unknowns[u]._name;
      // add system Poisson in ml_prob as a Linear Implicit System
      LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > (sys_name);


        // add solution "u" to system
      system.AddSolutionToSystemPDE( unknowns[u]._name.c_str() );

        std::vector< Unknown > unknowns_vec(1);
        unknowns_vec[0] = unknowns[u]; //need to turn this into a vector

        system.set_unknown_list_for_assembly(unknowns_vec); //way to communicate to the assemble function, which doesn't belong to any class


        // attach the assembling function to system
      system.SetAssembleFunction( assemble_pointer_vec[func].first );

      // initialize and solve the system
      system.init();
      
      system.SetOuterSolver(PREONLY);
      system.MGsolve();
    // ======= System - END ========================

      std::pair< double , double > norm = GetErrorNorm(& ml_sol, unknowns_vec);
      l2Norm[i][u]  = norm.first;
      semiNorm[i][u] = norm.second;

  // ======= Print - BEGIN  ========================
      std::vector < std::string > variablesToBePrinted;
      variablesToBePrinted.push_back("All");
            
      std::string  output_name = mesh_name + "_" + assemble_pointer_vec[func].second + "_" + unknowns[u]._name;
//       ml_sol.GetWriter()->SetGraphVariable ("u");
      ml_sol.GetWriter()->Write(output_name.c_str(), DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, i);
  // ======= Print - END  ========================

    } //end FE  space

      
      
} //end mesh level


  // ======= H1 - BEGIN  ========================
  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "SEMINORM ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < semiNorm[i].size(); j++) {
      std::cout << semiNorm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < maxNumberOfMeshes - 1) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < semiNorm[i].size(); j++) {
        std::cout << log(semiNorm[i][j] / semiNorm[i + 1][j]) / log(2.) << "\t\t\t";
      }

      std::cout << std::endl;
    }

  }
  // ======= H1 - END  ========================
  
  
  // ======= L2 - BEGIN  ========================
  std::cout << std::endl;
  std::cout << std::endl;

  // print the seminorm of the error and the order of convergence between different levels
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "l2 ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < l2Norm[i].size(); j++) {
      std::cout << l2Norm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < maxNumberOfMeshes - 1) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < l2Norm[i].size(); j++) {
        std::cout << log(l2Norm[i][j] / l2Norm[i + 1][j]) / log(2.) << "\t\t\t";
      }

      std::cout << std::endl;
    }

  }
  // ======= L2 - END  ========================

    } //end assemble func loop

    
    } //end mesh file type
    
    
  return 0;
}





/**
 * This function assemble the stiffnes matrix Jac and the residual vector Res
 * such that
 *                  Jac w = RES = F - Jac u0,
 * and consequently
 *        u = u0 + w satisfies Jac u = F
 **/

void AssemblePoissonProblem_old_fe_quadrature(MultiLevelProblem& ml_prob) {
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

  
  
  
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

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

  vector < double >  solu; // local solution
  solu.reserve(maxSize);

  vector < vector < double > > x (dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE BI/TRIQUADRATIC)

  for (unsigned i = 0; i < dim; i++) {
    x[i].reserve(maxSize);
  }

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  
  double weight = 0.; // gauss point weight

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  

  vector< double > Res; // local redidual vector
  Res.reserve(maxSize);

  vector < double > Jac; //local Jacobian matrix
  Jac.reserve(maxSize * maxSize);
  
  vector< int > l2GMap; // local to global mapping
  l2GMap.reserve(maxSize);
  

  KK->zero(); // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual Vector

  
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
    

    std::vector<unsigned> Sol_n_el_dofs_Mat_vol(1, nDofu);
  
    // resize local arrays
    solu.resize(nDofu);
    l2GMap.resize(nDofu);

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // local to global solution mapping
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // local storage of solution
      l2GMap[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);   // local to global system solution mapping
    }
    
    

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofx);
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }



    Res.assign(nDofu, 0.);    //resize and set to zero
    Jac.assign(nDofu * nDofu, 0.);    //resize and set to zero
    

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
        
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, boost::none);
      
      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
    
      vector < double > gradSolu_gss(dim, 0.);
      vector < double > x_gss(dim, 0.);

      for (unsigned i = 0; i < nDofu; i++) {
       
        for (unsigned jdim = 0; jdim < dim; jdim++) {
          gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solu[i];
          x_gss[jdim] += x[jdim][i] * phi[i];
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofu; i++) {

        double weakLaplace = 0.;

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          weakLaplace   +=  phi_x[i * dim + jdim] * gradSolu_gss[jdim];
        }
        
        Res[i] += - (   GetExactSolutionLaplace(x_gss) * phi[i] + weakLaplace) * weight;

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
void AssemblePoissonProblem_old_fe_quadrature_AD(MultiLevelProblem& ml_prob) {
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

  

  
  //  extract pointers to the several objects that we are going to use

 
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

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

  vector < adept::adouble >  solu; // local solution
  solu.reserve(maxSize);

  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned i = 0; i < dim; i++) {
    x[i].reserve(maxSize);
  }

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives

  double weight = 0.; // gauss point weight

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);

  vector< adept::adouble > aRes; // local redidual vector
  aRes.reserve(maxSize);

  vector< int > l2GMap; // local to global mapping
  l2GMap.reserve(maxSize);
  vector< double > Res; // local redidual vector
  Res.reserve(maxSize);
  vector < double > Jac;
  Jac.reserve(maxSize * maxSize);

  RES->zero(); // Set to zero all the entries of the Global Residual Vector
  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
     
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    std::vector<unsigned> Sol_n_el_dofs_Mat_vol(1, nDofu);


    // resize local arrays
    l2GMap.resize(nDofu);
    solu.resize(nDofu);

     // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      l2GMap[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofx);
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
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
      vector < adept::adouble > gradSolu_gss(dim, 0.);
      
      vector < double > x_gss(dim, 0.);

      for (unsigned i = 0; i < nDofu; i++) {

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solu[i];
          x_gss[jdim] += x[jdim][i] * phi[i];
        }
      }
      
      

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofu; i++) {

        adept::adouble laplace = 0.;

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          laplace   +=  phi_x[i * dim + jdim] * gradSolu_gss[jdim];
        }

        aRes[i] +=  ( GetExactSolutionLaplace(x_gss) * phi[i] + laplace) * weight;

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




std::pair < double, double > GetErrorNorm(MultiLevelSolution* ml_sol, std::vector< Unknown > & unknowns_vec) {
    
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
      GetExactSolutionGradient(x_gss, exactGradSol);

      for (unsigned j = 0; j < dim ; j++) {
        seminorm   += ((gradSolu_gss[j] - exactGradSol[j]) * (gradSolu_gss[j] - exactGradSol[j])) * weight;
      }

      double exactSol = GetExactSolutionValue(x_gss);
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
