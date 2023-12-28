/** tutorial/Ex3
 * This example shows how to set and solve the weak form of the nonlinear problem
 *                     -\Delta u + < u,u,u > \cdot \nabla u = f(x) \text{ on }\Omega,
 *            u=0 \text{ on } \Gamma,
 * on a box domain $\Omega$ with boundary $\Gamma$;
 * all the coarse-level meshes are removed;
 * a multilevel problem and an equation system are initialized;
 * a direct solver is used to solve the problem.
 **/

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "NumericVector.hpp"
#include "SparseMatrix.hpp"
#include "VTKWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "LinearEquationSolver.hpp"

#include "Solution_functions_over_domains_or_mesh_files.hpp"


#include "adept.h"


#include "FE_convergence.hpp"

#include "../tutorial_common.hpp"

#include "../all_mesh_generation_methods.hpp"




using namespace femus;



double Solution_set_initial_conditions_with_analytical_sol(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char * name) {

Math::Function< double > *  exact_sol =  ml_prob->get_ml_solution()->get_analytical_function(name);

double value = exact_sol->value(x);

   return value;   

}


bool SetBoundaryCondition(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;

  if (facename == 2) {
    dirichlet = false;
  }
  
  return dirichlet;
}

void AssemblePoissonPlusNonlinearAdvection(MultiLevelProblem& ml_prob);

void AssemblePoissonPlusNonlinearAdvection_AD(MultiLevelProblem& ml_prob);




double GetExactSolutionValue(const std::vector < double >& x) {
  double pi = acos(-1.);
  return cos(pi * x[0]) * cos(pi * x[1]);
};


void GetExactSolutionGradient(const std::vector < double >& x, std::vector < double >& solGrad) {
  double pi = acos(-1.);
  solGrad[0] = -pi * sin(pi * x[0]) * cos(pi * x[1]);
  solGrad[1] = -pi * cos(pi * x[0]) * sin(pi * x[1]);
};


double GetExactSolutionLaplace(const std::vector < double >& x) {
  double pi = acos(-1.);
  return -2.*pi * pi * cos(pi * x[0]) * cos(pi * x[1]);       // - pi*pi*cos(pi*x[0])*cos(pi*x[1]);
};





int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);


  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  const std::string relative_path_to_build_directory =  "../../../";
  const std::string mesh_file = relative_path_to_build_directory + Files::mesh_folder_path() + "01_gambit/3d/cube/minus0p5-plus0p5_minus0p5-plus0p5_minus0p5-plus0p5/cube_tet_Two_boundary_groups.neu";

  // const std::string mesh_file = relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/2d/square/minus0p5-plus0p5_minus0p5-plus0p5/square_-0p5-0p5x-0p5-0p5_divisions_2x2.med";

  mlMsh.ReadCoarseMesh(mesh_file.c_str(), "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
    probably in the furure it is not going to be an argument of this function   */
  
  unsigned dim = mlMsh.GetDimension();
  unsigned maxNumberOfMeshes = 6;


  std::vector < std::vector < double > > l2Norm;
  l2Norm.resize(maxNumberOfMeshes);

  std::vector < std::vector < double > > semiNorm;
  semiNorm.resize(maxNumberOfMeshes);

  for (unsigned i = 1; i < maxNumberOfMeshes; i++) {   // loop on the mesh level

    unsigned numberOfUniformLevels = i ;
    unsigned numberOfSelectiveLevels = 0;
    mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

    // erase all the coarse mesh levels
    mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

    // print mesh info
    mlMsh.PrintInfo();

    std::vector< Unknown > unknowns = systems__generate_list_of_scalar_unknowns_for_each_FE_family_lagrangian();


    l2Norm[i].resize( unknowns.size() );
    semiNorm[i].resize( unknowns.size() );

    for (unsigned u = 0; u < unknowns.size(); u++) {   // loop on the FE Order

      // define the multilevel solution and attach the mlMsh object to it
      MultiLevelSolution ml_sol(&mlMsh);
      
      // define the multilevel problem attach the ml_sol object to it
      MultiLevelProblem ml_prob(&ml_sol);
      
      Domains::square_m05p05::Function_Zero_on_boundary_4< double >  analytical_function;

      // add variables to ml_sol
      ml_sol.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order, unknowns[u]._time_order, unknowns[u]._is_pde_unknown);
      ml_sol.set_analytical_function(unknowns[u]._name.c_str(), & analytical_function);   
      ml_sol.Initialize(unknowns[u]._name.c_str(), Solution_set_initial_conditions_with_analytical_sol, & ml_prob);

      // attach the boundary condition function and generate boundary data
      ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
      ml_sol.GenerateBdc( unknowns[u]._name.c_str(), "Steady", & ml_prob);


      ml_prob.get_systems_map().clear();
      ml_prob.set_current_system_number(0/*u*/);               //way to communicate to the assemble function, which doesn't belong to any class

    // ======= System - BEGIN ========================
      // add system Poisson in ml_prob as a Linear Implicit System
      NonLinearImplicitSystem& system = ml_prob.add_system < NonLinearImplicitSystem > ("Poisson");

      // add solution "u" to system
      system.AddSolutionToSystemPDE(unknowns[u]._name.c_str());
      
        // set unknown list
        std::vector< Unknown > unknowns_vec(1);
        unknowns_vec[0] = unknowns[u]; //need to turn this into a vector

        system.set_unknown_list_for_assembly(unknowns_vec); //way to communicate to the assemble function, which doesn't belong to any class

        
      // attach the assembling function to system
      system.SetAssembleFunction(AssemblePoissonPlusNonlinearAdvection_AD);

      // initilaize and solve the system
      system.init();
      system.SetOuterSolver(PREONLY);
      system.MGsolve();
    // ======= System - END ========================
      
      

      std::pair< double , double > norm = GetErrorNorm_L2_H1_with_analytical_sol(&ml_sol, unknowns_vec, GetExactSolutionValue, GetExactSolutionGradient);
      l2Norm[i][u]  = norm.first;
      semiNorm[i][u] = norm.second;
      // print solutions
      std::vector < std::string > variablesToBePrinted;
      variablesToBePrinted.push_back("All");

      VTKWriter vtkIO(&ml_sol);
      vtkIO.SetDebugOutput(true);
      vtkIO.Write(Files::_application_output_directory, fe_fams_for_files[ FILES_CONTINUOUS_LINEAR ], variablesToBePrinted, i + u*10);
      vtkIO.Write(Files::_application_output_directory, fe_fams_for_files[ FILES_CONTINUOUS_QUADRATIC ], variablesToBePrinted, i + u*10);
      vtkIO.Write(Files::_application_output_directory, fe_fams_for_files[ FILES_CONTINUOUS_BIQUADRATIC ], variablesToBePrinted, i + u*10);


    }
    
    
  }

  // ======= L2 - BEGIN  ========================
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "l2 ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (unsigned i = 1; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < 3; j++) {
      std::cout << l2Norm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < maxNumberOfMeshes - 1) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < 3; j++) {
        std::cout << log(l2Norm[i][j] / l2Norm[i + 1][j]) / log(2.) << "\t\t\t";
      }

      std::cout << std::endl;
    }

  }
  // ======= L2 - END  ========================

  // ======= H1 - BEGIN  ========================
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "SEMINORM ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (unsigned i = 1; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < 3; j++) {
      std::cout << semiNorm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < maxNumberOfMeshes - 1) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < 3; j++) {
        std::cout << log(semiNorm[i][j] / semiNorm[i + 1][j]) / log(2.) << "\t\t\t";
      }

      std::cout << std::endl;
    }

  }
  // ======= H1 - END ========================





}



/**
 * Given the non linear problem
 *
 *      - \Delta u + < u, u, u > \cdot \nabla u = f(x)
 * in the unit box centered in the origin with
 *                      f(x) = - \Delta u_e + < u_e, u_e, u_e > \cdot \nabla u_e
 *                    u_e = \cos ( \pi * x ) * \cos( \pi * y ),
 * the following function assembles (explicitely) the Jacobian matrix J(u^i) and the residual vector Res(u^i)
 * for the Newton iteration, i.e.
 *
 *                  J(u^i) w = Res(u^i) = f(x) - ( - \Delta u^i + < u^i, u^i, u^i > \cdot \nabla u^i ),
 *        u^{i+1} = u^i + w
 *        where
 *        J(u^i) w = - \Delta w  + < w , w , w >  \cdot \nabla u^i + < u^i , u^i , u^i >  \cdot \nabla w
 *
 **/

void AssemblePoissonPlusNonlinearAdvection(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled



  
    const unsigned current_system_number = ml_prob.get_current_system_number();
    NonLinearImplicitSystem * mlPdeSys  = & ml_prob.get_system< NonLinearImplicitSystem >(current_system_number); 
    
    // II) Unknowns of the System
  std::vector< Unknown >   unknowns = ml_prob.get_system< NonLinearImplicitSystem >(current_system_number).get_unknown_list_for_assembly();
  
  
   const unsigned level = mlPdeSys->GetLevelToAssemble();
 
  
  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->GetMeshElements();  // pointer to the elem object in msh (level)

  MultiLevelSolution*  ml_sol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = ml_sol->GetIndex( unknowns[0]._name.c_str() );    // get the position of "u" in the ml_sol object
  unsigned soluType = ml_sol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex(  unknowns[0]._name.c_str() );    // get the position of "u" in the pdeSys object

  std::vector < double >  solu; // local solution

  std::vector < std::vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector < int > sysDof; // local to global pdeSys dofs
  std::vector <double> phi;  // local test function
  std::vector <double> phi_x; // local test function first order partial derivatives
  std::vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  std::vector < double > Res; // local redidual vector
  std::vector < double > J; // local Jacobian matrix

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  solu.reserve(maxSize);

  for (unsigned i = 0; i < dim; i++)
    x[i].reserve(maxSize);

  sysDof.reserve(maxSize);
  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  phi_xx.reserve(maxSize * dim2);
  Res.reserve(maxSize);
  J.reserve(maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {


    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofs2 = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
    // resize local arrays
    sysDof.resize(nDofs);
    solu.resize(nDofs);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofs2);
    }

    Res.resize(nDofs);    //resize
    std::fill(Res.begin(), Res.end(), 0);    //set Res to zero

    J.resize(nDofs * nDofs);    //resize
    std::fill(J.begin(), J.end(), 0);    //set K to zero

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      sysDof[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofs2; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->GetTopology()->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      double soluGauss = 0;
      std::vector < double > soluGauss_x(dim, 0.);
      std::vector < double > xGauss(dim, 0.);

      for (unsigned i = 0; i < nDofs; i++) {
        soluGauss += phi[i] * solu[i];

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          soluGauss_x[jdim] += phi_x[i * dim + jdim] * solu[i];
          xGauss[jdim] += x[jdim][i] * phi[i];
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {

        double nonLinearTerm = 0.;
        double mLaplace = 0.;

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          mLaplace   +=  phi_x[i * dim + jdim] * soluGauss_x[jdim];
          nonLinearTerm += soluGauss * soluGauss_x[jdim] * phi[i];
        }

        double exactSolValue = GetExactSolutionValue(xGauss);
        std::vector < double > exactSolGrad(dim);
        GetExactSolutionGradient(xGauss , exactSolGrad);
        double exactSolLaplace = GetExactSolutionLaplace(xGauss);


        double f = (- exactSolLaplace + exactSolValue * (exactSolGrad[0] + exactSolGrad[1])) * phi[i] ;
        Res[i] += (f - (mLaplace + nonLinearTerm)) * weight;

        // *** phi_j loop ***
        for (unsigned j = 0; j < nDofs; j++) {
          mLaplace = 0.;
          nonLinearTerm = 0.;

          for (unsigned kdim = 0; kdim < dim; kdim++) {
            mLaplace += (phi_x[i * dim + kdim] * phi_x[j * dim + kdim]);
            nonLinearTerm +=   phi[i] * (phi[j] * soluGauss_x[kdim] + soluGauss * phi_x[j * dim + kdim]);
          }

          J[i * nDofs + j] += (mLaplace + nonLinearTerm) * weight;
        } // end phi_j loop
      } // end phi_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    RES->add_vector_blocked(Res, sysDof);

    KK->add_matrix_blocked(J, sysDof, sysDof);
  } //end element loop for each process

  RES->close();

  KK->close();

  // ***************** END ASSEMBLY *******************
}


/**
 * Given the non linear problem
 *
 *      - \Delta u + < u, u, u > \cdot \nabla u = f(x),
 *
 * in the unit box centered in the origin with
 *
 *                      f(x) = - \Delta u_e + < u_e, u_e, u_e > \cdot \nabla u_e,
 *                    u_e = \cos ( \pi * x ) * \cos( \pi * y ),
 *
 *the following function assembles the residual vector Res(u^i) and using automatic differentiation gets
 *the exact Jacobian matrix J(u^i) for the Newton iteration, i.e.
 *
 *                   J(u^i) w = Res(u^i) = f(x) - ( - \Delta u^i + < u^i, u^i, u^i > \cdot \nabla u^i ),
 *         u^{i+1} = u^i + w,
 *        where
 *        J(u^i) w = - \Delta w  + < w , w , w >  \cdot \nabla u^i + < u^i , u^i , u^i >  \cdot \nabla w.
 *
 **/

void AssemblePoissonPlusNonlinearAdvection_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
    const unsigned current_system_number = ml_prob.get_current_system_number();
    NonLinearImplicitSystem * mlPdeSys  = & ml_prob.get_system< NonLinearImplicitSystem >(current_system_number);
  
    // II) Unknowns of the System
  std::vector< Unknown >   unknowns = ml_prob.get_system< NonLinearImplicitSystem >(current_system_number).get_unknown_list_for_assembly();
  
  
  
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->GetMeshElements();  // pointer to the elem object in msh (level)

  MultiLevelSolution*  ml_sol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = ml_sol->GetIndex(  unknowns[0]._name.c_str() );    // get the position of "u" in the ml_sol object
  unsigned soluType = ml_sol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex(  unknowns[0]._name.c_str() );    // get the position of "u" in the pdeSys object

  std::vector < adept::adouble >  solu; // local solution

  std::vector < std::vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector < int > sysDof; // local to global pdeSys dofs
  std::vector <double> phi;  // local test function
  std::vector <double> phi_x; // local test function first order partial derivatives
  std::vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  std::vector < double > Res; // local redidual vector
  std::vector < adept::adouble > aRes; // local redidual vector


  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  solu.reserve(maxSize);

  for (unsigned i = 0; i < dim; i++)
    x[i].reserve(maxSize);

  sysDof.reserve(maxSize);
  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  phi_xx.reserve(maxSize * dim2);
  Res.reserve(maxSize);
  aRes.reserve(maxSize);

  std::vector < double > Jac; // local Jacobian matrix (ordered by column, adept)
  Jac.reserve(maxSize * maxSize);
  std::vector < double > Jact; // local Jacobian matrix (ordered by raw, PETSC)
  Jact.reserve(maxSize * maxSize);



  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {


    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofs2 = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
    // resize local arrays
    sysDof.resize(nDofs);
    solu.resize(nDofs);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofs2);
    }

    Res.resize(nDofs);    //resize
    std::fill(Res.begin(), Res.end(), 0);    //set Res to zero

    aRes.resize(nDofs);    //resize
    std::fill(aRes.begin(), aRes.end(), 0);    //set aRes to zero


    Jact.resize(nDofs * nDofs);
    Jac.resize(nDofs * nDofs);


    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      sysDof[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofs2; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->GetTopology()->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }


    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      adept::adouble soluGauss = 0;
      std::vector < adept::adouble > soluGauss_x(dim, 0.);
      std::vector < double > xGauss(dim, 0.);

      for (unsigned i = 0; i < nDofs; i++) {
        soluGauss += phi[i] * solu[i];

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          soluGauss_x[jdim] += phi_x[i * dim + jdim] * solu[i];
          xGauss[jdim] += x[jdim][i] * phi[i];
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {

        adept::adouble nonLinearTerm = 0.;
        adept::adouble mLaplace = 0.;

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          mLaplace   +=  phi_x[i * dim + jdim] * soluGauss_x[jdim];
          //nonLinearTerm += soluGauss * soluGauss_x[jdim] * phi[i];
        }

        double exactSolValue = GetExactSolutionValue(xGauss);
        std::vector < double > exactSolGrad(dim);
        GetExactSolutionGradient(xGauss , exactSolGrad);
        double exactSolLaplace = GetExactSolutionLaplace(xGauss);


        double f = (- exactSolLaplace + 0*exactSolValue * (exactSolGrad[0] + exactSolGrad[1])) * phi[i] ;
        aRes[i] += (f - (mLaplace + nonLinearTerm)) * weight;

      } // end phi_i loop
    } // end gauss point loop


    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    for (int i = 0; i < nDofs; i++) {
      Res[i] = aRes[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);



    // define the dependent variables
    s.dependent(&aRes[0], nDofs);

    // define the independent variables
    s.independent(&solu[0], nDofs);

    // get the jacobian matrix (ordered by column)
    s.jacobian(&Jac[0]);

    // get the jacobian matrix (ordered by raw, i.e. Jact=Jac^t)
    for (int inode = 0; inode < nDofs; inode++) {
      for (int jnode = 0; jnode < nDofs; jnode++) {
        Jact[inode * nDofs + jnode] = -Jac[jnode * nDofs + inode];
      }
    }

    //store Jact in the global matrix KK
    KK->add_matrix_blocked(Jact, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();

  KK->close();

  // ***************** END ASSEMBLY *******************
}


