/** tutorial/Ex2
 * This example shows how to set and solve the weak form of the Poisson problem
 *                    $$ \Delta u = 1 \text{ on }\Omega, $$
 *          $$ u=0 \text{ on } \Gamma, $$
 * on a square domain $\Omega$ with boundary $\Gamma$;
 * all the coarse-level meshes are removed;
 * a multilevel problem and an equation system are initialized;
 * a direct solver is used to solve the problem.
 **/

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "LinearImplicitSystem.hpp"
#include "adept.h"
// #include "Files.hpp"

// command to view matrices
// ./tutorial_ex2_b -mat_view ::ascii_info_detail


using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0;

//   if (faceName == 2)
//     dirichlet = false;

  return dirichlet;
}

void AssemblePoissonProblem(MultiLevelProblem& ml_prob);

void AssemblePoissonProblem_AD(MultiLevelProblem& ml_prob);

std::pair < double, double > GetErrorNorm(MultiLevelSolution* mlSol, Solution* sol_from_restriction); 
// ||u_i - u_h||/||u_i-u_(h/2)|| = 2^alpha, alpha is order of conv 

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

//   // ======= Files ========================
//   Files files; 
//         files.CheckIODirectories();
//         files.RedirectCout();
 
  // ======= Quad Rule ========================
  std::string fe_quad_rule("seventh");
 /* "seventh" is the order of accuracy that is used in the gauss integration scheme
    In the future it is not going to be an argument of the mesh function   */

 // define multilevel mesh
  MultiLevelMesh mlMsh;
 // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  const unsigned int nsub_x = 2;
  const unsigned int nsub_y = 2;
  const unsigned int nsub_z = 0;
  const std::vector<double> xyz_min = {-0.5,-0.5,0.};
  const std::vector<double> xyz_max = { 0.5, 0.5,0.};
  const ElemType geom_elem_type = QUAD9;
  mlMsh.GenerateCoarseBoxMesh(nsub_x,nsub_y,nsub_z,xyz_min[0],xyz_max[0],xyz_min[1],xyz_max[1],xyz_min[2],xyz_max[2],geom_elem_type,fe_quad_rule.c_str());
//   mlMsh.ReadCoarseMesh("./input/square_quad.neu", "seventh", scalingFactor);

  MultiLevelMesh mlMsh_finest;
  mlMsh_finest.GenerateCoarseBoxMesh(nsub_x,nsub_y,nsub_z,xyz_min[0],xyz_max[0],xyz_min[1],xyz_max[1],xyz_min[2],xyz_max[2],geom_elem_type,fe_quad_rule.c_str());
//   mlMsh_finest.ReadCoarseMesh("./input/square_quad.neu", "seventh", scalingFactor);

  unsigned dim = mlMsh.GetDimension();
  unsigned maxNumberOfMeshes;

  if (dim == 2) {
    maxNumberOfMeshes = 5;
  } else {
    maxNumberOfMeshes = 4;
  }

  const unsigned gap = 2;

  vector < vector < double > > l2Norm;
  l2Norm.resize(maxNumberOfMeshes + 1 - gap);

  vector < vector < double > > semiNorm;
  semiNorm.resize(maxNumberOfMeshes + 1 - gap);

  
  MultiLevelSolution * mlSol_finest;
  
    std::vector< FEOrder > feOrder = {FIRST, SERENDIPITY, SECOND};

    
    for (int i = maxNumberOfMeshes - 1; i >= 0; i--) {   // loop on the mesh level

   unsigned numberOfUniformLevels = i + 1;
    unsigned numberOfSelectiveLevels = 0;
    mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
    // erase all the coarse mesh levels
    mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

    
    if (i == maxNumberOfMeshes - 1) {
        unsigned numberOfUniformLevels_finest = numberOfUniformLevels;
        mlMsh_finest.RefineMesh(numberOfUniformLevels_finest, numberOfUniformLevels_finest + numberOfSelectiveLevels, NULL);
//         mlMsh_finest.EraseCoarseLevels(numberOfUniformLevels_finest - 1 - 1); //I need to keep the structures at all levels here so I can restrict every time
    }
    
    // print mesh info
    mlMsh.PrintInfo();
    
    if (i < l2Norm.size()) {
    l2Norm[i].resize(feOrder.size());
    semiNorm[i].resize(feOrder.size());
    }
    
    for (unsigned j = 0; j < feOrder.size(); j++) {   // loop on the FE Order
      // define the multilevel solution and attach the mlMsh object to it
      MultiLevelSolution mlSol(&mlMsh);
      
      // add variables to mlSol
      mlSol.AddSolution("u", LAGRANGE, feOrder[j]);
      mlSol.Initialize("All");
      
      // attach the boundary condition function and generate boundary data
      mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
      mlSol.GenerateBdc("u");
      
      // define the multilevel problem attach the mlSol object to it
      MultiLevelProblem mlProb(&mlSol);

//       mlProb.SetFilesHandler(&files);
      
      // add system Poisson in mlProb as a Linear Implicit System
      LinearImplicitSystem& system = mlProb.add_system < LinearImplicitSystem > ("Poisson");

      // add solution "u" to system
      system.AddSolutionToSystemPDE("u");

      // attach the assembling function to system
      system.SetAssembleFunction(AssemblePoissonProblem_AD);

      // initilaize and solve the system
      system.init();
//       system.ClearVariablesToBeSolved();
//       system.AddVariableToBeSolved("All");
// 
//       mlSol.SetWriter(VTK);
//       mlSol.GetWriter()->SetDebugOutput(true);
//   
//       system.SetDebugLinear(true);
//       system.SetMaxNumberOfLinearIterations(6);
//       system.SetAbsoluteLinearConvergenceTolerance(1.e-4);

      system.SetOuterSolver(PREONLY);
      system.MGsolve();
      
     
    if (i < maxNumberOfMeshes - gap/*1*/) {
        //restrict the fine solution at the current level  ==================
//     for(unsigned short j = 0; j < mlSol._mlMesh->GetNumberOfLevels(); j++) { //all levels
      const unsigned coarse = i;
      const unsigned fine   = coarse + gap;
    for (unsigned nf = 0; nf < gap; nf++) {
      mlSol_finest->CoarsenSolutionByOneLevel_wrong( fine - nf);
    }
//        }

        //pass the restriction of the fine solution to the function that computes the error  ==================
      Solution* sol = mlSol_finest->GetSolutionLevel(i);    // pointer to the solution (level) object
      std::pair< double , double > norm = GetErrorNorm(&mlSol,sol);

      l2Norm[i][j]  = norm.first;
      semiNorm[i][j] = norm.second;
      
      
        
    }
    else if (i == maxNumberOfMeshes - 1) {
        //store the fine solution  ==================
              mlSol_finest = new MultiLevelSolution (&mlMsh_finest);  //with the declaration outside and a "new" inside it persists outside the loop scopes
              mlSol_finest->AddSolution("u", LAGRANGE, feOrder[j]);
              mlSol_finest->Initialize("All");
              mlSol_finest->AttachSetBoundaryConditionFunction(SetBoundaryCondition);
              mlSol_finest->GenerateBdc("u");
              
//       assert( mlSol._mlMesh->GetNumberOfLevels() == mlSol_finest->_mlMesh->GetNumberOfLevels() - 1 );

//     for(unsigned short k = 0; k < mlSol._mlMesh->GetNumberOfLevels(); k++) { //all levels (only one)
    //       _solution[k]->CopySolutionToOldSolution();  //started from here
              
       const unsigned level_index_current = 0;
      //@todo there is a duplicate function in MLSol: GetSolutionLevel() and GetLevel()
     for(unsigned short j = 0; j <   mlSol.GetSolutionLevel(level_index_current)->_Sol.size(); j++) {    //all variables
               *(mlSol_finest->GetLevel(i)->_Sol[j]) = *(mlSol.GetSolutionLevel(level_index_current)->_Sol[j]);
          }
//        }
      
      
       
       
    }
        
    
    
      // print solutions
      std::vector < std::string > variablesToBePrinted;
      variablesToBePrinted.push_back("All");

      VTKWriter vtkIO(&mlSol);
      vtkIO.Write(/*files.GetOutputPath()*/ DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, i);
      
      VTKWriter vtkIO_finest(mlSol_finest);
      vtkIO_finest.Write(i+1, /*files.GetOutputPath()*/ DEFAULT_OUTPUTDIR, "biquadratic"/*_finest*/, variablesToBePrinted, i+3 + maxNumberOfMeshes);

    }
  }
  
  delete mlSol_finest; 

//   print the seminorm of the error and the order of convergence between different levels
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "l2 ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\tSERENDIPITY\tSECOND\n";

  for (int i = l2Norm.size() - 1; i >= 0; i--) {   // loop on the mesh level
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < feOrder.size(); j++) {
      std::cout << l2Norm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < l2Norm.size() - 2) {
      std::cout.precision(3);
      std::cout << "\t";

      for (unsigned j = 0; j < feOrder.size(); j++) {
        std::cout << log(l2Norm[i][j] / l2Norm[i + 1][j]) / log(2.) << "  \t\t";
      }

      std::cout << std::endl;
    }

  }

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "SEMINORM ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\tSERENDIPITY\tSECOND\n";

  for (int i = l2Norm.size() - 1; i >= 0; i--) {   // loop on the mesh level
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < feOrder.size(); j++) {
      std::cout << semiNorm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < l2Norm.size() - 2) {
      std::cout.precision(3);
      std::cout << "\t";

      for (unsigned j = 0; j < feOrder.size(); j++) {
        std::cout << log(semiNorm[i][j] / semiNorm[i + 1][j]) / log(2.) << "  \t\t";
      }

      std::cout << std::endl;
    }

  }

  return 0;
}

double GetExactSolutionValue(const std::vector < double >& x) {
  double pi = acos(-1.);
  return cos(pi * x[0]) * cos(pi * x[1]);
};

void GetExactSolutionGradient(const std::vector < double >& x, vector < double >& solGrad) {
  double pi = acos(-1.);
  solGrad[0]  = -pi * sin(pi * x[0]) * cos(pi * x[1]);
  solGrad[1] = -pi * cos(pi * x[0]) * sin(pi * x[1]);
};

double GetExactSolutionLaplace(const std::vector < double >& x) {
  double pi = acos(-1.);
  return -pi * pi * cos(pi * x[0]) * cos(pi * x[1]) - pi * pi * cos(pi * x[0]) * cos(pi * x[1]);
};

/**
 * This function assemble the stiffnes matrix Jac and the residual vector Res
 * such that
 *                  Jac w = RES = F - Jac u0,
 * and consequently
 *        u = u0 + w satisfies Jac u = F
 **/
void AssemblePoissonProblem(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data

  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("Poisson");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex("u");    // get the position of "u" in the pdeSys object

  vector < double >  solu; // local solution
  solu.reserve(maxSize);

  vector < double >  solu_exact_at_dofs; // local solution
  solu_exact_at_dofs.reserve(maxSize);

  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned i = 0; i < dim; i++) {
    x[i].reserve(maxSize);
  }

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  phi_xx.reserve(maxSize * dim2);

  vector< double > Res; // local redidual vector
  Res.reserve(maxSize);

  vector< int > l2GMap; // local to global mapping
  l2GMap.reserve(maxSize);
  vector < double > Jac;
  Jac.reserve(maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    // resize local arrays
    l2GMap.resize(nDofu);
    solu.resize(nDofu);
    solu_exact_at_dofs.resize(nDofu);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofx);
    }

    Res.resize(nDofu);    //resize
    std::fill(Res.begin(), Res.end(), 0);    //set Res to zero

    Jac.resize(nDofu * nDofu);    //resize
    std::fill(Jac.begin(), Jac.end(), 0);    //set Jac to zero

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    } 
    
    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
        std::vector<double> x_at_node(dim,0.);
        for (unsigned jdim = 0; jdim < dim; jdim++) x_at_node[jdim] = x[jdim][i];
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      solu_exact_at_dofs[i] = GetExactSolutionValue(x_at_node);
      l2GMap[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }





    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      double solu_gss = 0;
      vector < double > gradSolu_gss(dim, 0.);
      vector < double > gradSolu_exact_gss(dim, 0.);
      vector < double > x_gss(dim, 0.);

      for (unsigned i = 0; i < nDofu; i++) {
        solu_gss += phi[i] * solu[i];

        for (unsigned jdim = 0; jdim < dim; jdim++) {
                gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solu[i];
          gradSolu_exact_gss[jdim] += phi_x[i * dim + jdim] * solu_exact_at_dofs[i];
          x_gss[jdim] += x[jdim][i] * phi[i];
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofu; i++) {

        double laplace = 0.;
        double laplace_weak_exact = 0.;

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          laplace   +=  phi_x[i * dim + jdim] * gradSolu_gss[jdim];
          laplace_weak_exact +=  phi_x[i * dim + jdim] * gradSolu_exact_gss[jdim];
        }

//         double srcTerm = - GetExactSolutionLaplace(x_gss);
//         Res[i] += (srcTerm * phi[i] - laplace) * weight;
        Res[i] += (laplace_weak_exact - laplace) * weight;

        // *** phi_j loop ***
        for (unsigned j = 0; j < nDofu; j++) {
          laplace = 0.;

          for (unsigned kdim = 0; kdim < dim; kdim++) {
            laplace += (phi_x[i * dim + kdim] * phi_x[j * dim + kdim]) * weight;
          }

          Jac[i * nDofu + j] += laplace;
        } // end phi_j loop

      } // end phi_i loop
    } // end gauss point loop


    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    RES->add_vector_blocked(Res, l2GMap);

    //store K in the global matrix KK
    KK->add_matrix_blocked(Jac, l2GMap, l2GMap);

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
void AssemblePoissonProblem_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object


  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("Poisson");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex("u");    // get the position of "u" in the pdeSys object

  vector < adept::adouble >  solu; // local solution
  solu.reserve(maxSize);
  
  vector < adept::adouble >  solu_exact_at_dofs; // local solution
  solu_exact_at_dofs.reserve(maxSize);

  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned i = 0; i < dim; i++) {
    x[i].reserve(maxSize);
  }

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  phi_xx.reserve(maxSize * dim2);

  vector< adept::adouble > aRes; // local redidual vector
  aRes.reserve(maxSize);

  vector< int > l2GMap; // local to global mapping
  l2GMap.reserve(maxSize);
  vector< double > Res; // local redidual vector
  Res.reserve(maxSize);
  vector < double > Jac;
  Jac.reserve(maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    // resize local arrays
    l2GMap.resize(nDofu);
    solu.resize(nDofu);
    solu_exact_at_dofs.resize(nDofu);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofx);
    }

    aRes.resize(nDofu);    //resize
    std::fill(aRes.begin(), aRes.end(), 0);    //set aRes to zero

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    } 
    
    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
        std::vector<double> x_at_node(dim,0.);
        for (unsigned jdim = 0; jdim < dim; jdim++) x_at_node[jdim] = x[jdim][i];
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
                    solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      solu_exact_at_dofs[i] = GetExactSolutionValue(x_at_node);
      l2GMap[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }




    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      adept::adouble solu_gss = 0;
      vector < adept::adouble > gradSolu_gss(dim, 0.);
      vector < adept::adouble > gradSolu_exact_gss(dim, 0.);
      vector < double > x_gss(dim, 0.);

      for (unsigned i = 0; i < nDofu; i++) {
        solu_gss += phi[i] * solu[i];

        for (unsigned jdim = 0; jdim < dim; jdim++) {
                gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solu[i];
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
        

//         double srcTerm = - GetExactSolutionLaplace(x_gss);
//         aRes[i] += (srcTerm * phi[i] - laplace) * weight;
        aRes[i] += (laplace_weak_exact - laplace) * weight;

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

  } //end element loop for each process

  RES->close();

  KK->close();

  // ***************** END ASSEMBLY *******************
}

std::pair < double, double > GetErrorNorm(MultiLevelSolution* mlSol, Solution* sol_finer) {
  unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1u;
  //  extract pointers to the several objects that we are going to use
  Mesh*     msh = mlSol->_mlMesh->GetLevel(level);    // pointer to the mesh (level) object
  elem*     el  = msh->el;  // pointer to the elem object in msh (level)
  Solution* sol = mlSol->GetSolutionLevel(level);    // pointer to the solution (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  vector < double >  solu; // local solution

  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
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
  phi_xx.reserve(maxSize * dim2);

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
      solu_exact_at_dofs[i] = GetExactSolutionValue(x_at_node);
    }


    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

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

      vector <double> exactGradSol(dim);    GetExactSolutionGradient(x_gss, exactGradSol);

      for (unsigned j = 0; j < dim ; j++) {
        seminorm           += ((gradSolu_gss[j]       - exactGradSol[j]) * (gradSolu_gss[j]       - exactGradSol[j])) * weight;
        seminorm_inexact   += ((gradSolu_gss[j] - gradSolu_finer_gss[j]) * (gradSolu_gss[j] - gradSolu_finer_gss[j])) * weight;
        seminorm_exact_dofs      += ((gradSolu_gss[j]       - gradSolu_exact_at_dofs_gss[j]) * (gradSolu_gss[j]       - gradSolu_exact_at_dofs_gss[j])) * weight;
      }

      double exactSol = GetExactSolutionValue(x_gss);
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
  
//   return std::pair < double, double > (sqrt(l2norm), sqrt(seminorm));
//   return std::pair < double, double > (sqrt(l2norm_exact_dofs), sqrt(seminorm_exact_dofs));
  return std::pair < double, double > (sqrt(l2norm_inexact), sqrt(seminorm_inexact));

}
