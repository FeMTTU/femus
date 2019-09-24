/** \file Ex7.cpp
 *  \brief This example shows how to set and solve the weak form
 *   of the Boussinesq appoximation of the Navier-Stokes Equation
 *
 *  \f{eqnarray*}
 *  && \mathbf{V} \cdot \nabla T - \nabla \cdot\alpha \nabla T = 0 \\
 *  && \mathbf{V} \cdot \nabla \mathbf{V} - \nabla \cdot \nu (\nabla \mathbf{V} +(\nabla \mathbf{V})^T)
 *  +\nabla P = \beta T \mathbf{j} \\
 *  && \nabla \cdot \mathbf{V} = 0
 *  \f}
 *  in a unit box domain (in 2D and 3D) with given temperature 0 and 1 on
 *  the left and right walls, respectively, and insulated walls elsewhere.
 *  \author Eugenio Aulisa
 */

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"


using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;
  return dirichlet;
}


bool SetRefinementFlag(const std::vector < double >& x, const int& elemgroupnumber, const int& level) {

  bool refine = 0;

  if (elemgroupnumber == 6 && level < 4) refine = 1;

  if (elemgroupnumber == 7 && level < 5) refine = 1;

  if (elemgroupnumber == 8 && level < 6) refine = 1;

//   if (elemgroupnumber==6 && level<1) refine=1;
//   if (elemgroupnumber==7 && level<2) refine=1;
//   if (elemgroupnumber==8 && level<3) refine=1;

  return refine;

}



void AssemblePoisson_AD(MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  //mlMsh.ReadCoarseMesh("./input/cube_hex.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/square_quad.neu", "seventh", scalingFactor);
  mlMsh.ReadCoarseMesh("./input/quadAMR.neu", "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

//   unsigned numberOfUniformLevels = 3;
//   unsigned numberOfSelectiveLevels = 0;
//   mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  unsigned numberOfUniformLevels = 4;
  unsigned numberOfSelectiveLevels = 3;
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);

  // erase all the coarse mesh levels
  //mlMsh.EraseCoarseLevels(numberOfUniformLevels - 3);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol

  mlSol.AddSolution("V", LAGRANGE, SECOND);

  mlSol.AddSolution("U", LAGRANGE, SECOND);

  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("Poisson");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("U");

  //system.SetLinearEquationSolverType(FEMuS_DEFAULT);
  system.SetLinearEquationSolverType(FEMuS_ASM); // Additive Swartz Method
  // attach the assembling function to system
  system.SetAssembleFunction(AssemblePoisson_AD);

  system.SetMaxNumberOfNonLinearIterations(10);
  system.SetMaxNumberOfLinearIterations(3);
  system.SetAbsoluteLinearConvergenceTolerance(1.e-12);
  system.SetNonLinearConvergenceTolerance(1.e-8);
  system.SetMgType(F_CYCLE); 

  system.SetNumberPreSmoothingStep(0);  
  system.SetNumberPostSmoothingStep(2);
  // initilaize and solve the system
  system.init();

  system.SetSolverFineGrids(GMRES);
  system.SetPreconditionerFineGrids(ILU_PRECOND);
  system.SetTolerances(1.e-3, 1.e-20, 1.e+50, 5);

  system.SetNumberOfSchurVariables(1);
  system.SetElementBlockNumber(4); 
  //system.SetDirichletBCsHandling(ELIMINATION); 
  //system.solve();
  system.MGsolve();

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  GMVWriter gmvIO(&mlSol);
  variablesToBePrinted.push_back("all");
  gmvIO.SetDebugOutput(true);
  gmvIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);


  return 0;
}


void AssemblePoisson_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("Poisson");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*           msh         = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*           el          = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*   mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*   sol             = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys  = mlPdeSys->_LinSolver[level];  // pointer to the equation (level) object
  SparseMatrix*   KK          = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*  RES         = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //solution variable
  unsigned solUIndex;
  solUIndex = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object = 0
  unsigned solUType = mlSol->GetSolutionType(solUIndex);    // get the finite element type for "T"



  unsigned solUPdeIndex;
  solUPdeIndex = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object = 0

  std::cout << solUIndex << " " << solUPdeIndex << std::endl;


  vector < adept::adouble >  solU; // local solution
  vector< adept::adouble > aResU; // local residual vector

  vector < vector < double > > crdX(dim);    // local coordinates
  unsigned crdXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  solU.reserve(maxSize);
  aResU.reserve(maxSize);

  for (unsigned  k = 0; k < dim; k++) {
    crdX[k].reserve(maxSize);
  }

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  phi_xx.reserve(maxSize * dim2);

  double weight; // gauss point weight

  vector< int > sysDof; // local to global pdeSys dofs
  sysDof.reserve(maxSize);

  vector< double > ResU; // local residual vector
  ResU.reserve(maxSize);

  vector < double > Jac;
  Jac.reserve(maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns 
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);      // element geometry type

    unsigned nDofsU = msh->GetElementDofNumber(iel, solUType);      // number of solution element dofs
    unsigned nDofsX = msh->GetElementDofNumber(iel, crdXType);      // number of solution element dofs

    // resize local arrays
    sysDof.resize(nDofsU);
    solU.resize(nDofsU);

    for (unsigned  k = 0; k < dim; k++) {
      crdX[k].resize(nDofsX);
    }

    aResU.assign(nDofsU, 0);

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofsU; i++) {
      unsigned solUDof = msh->GetSolutionDof(i, iel, solUType);    // local to global mapping of the solution U
      solU[i] = (*sol->_Sol[solUIndex])(solUDof);      // value of the solution U in the dofs
      sysDof[i] = pdeSys->GetSystemDof(solUIndex, solUPdeIndex, i, iel);    // local to global mapping between solution U and system
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, crdXType);   // local to global mapping of the coordinate X[dim]

      for (unsigned k = 0; k < dim; k++) {
        crdX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // value of the solution X[dim]  // 
      }
    }

    s.new_recording();

    // *** Gauss point loop *** // 
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solUType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solUType]->Jacobian(crdX, ig, weight, phi, phi_x, phi_xx);

      adept::adouble solUig = 0; // solution U in the gauss point
      vector < adept::adouble > gradSolUig(dim, 0.); // gradient of solution U in the gauss point

      for (unsigned i = 0; i < nDofsU; i++) {
        solUig += phi[i] * solU[i];

        for (unsigned j = 0; j < dim; j++) {
          gradSolUig[j] += phi_x[i * dim + j] * solU[i]; 
        }
      }

      double nu = 1.;

      // *** phiU_i loop ***
      for (unsigned i = 0; i < nDofsU; i++) {
        adept::adouble LaplaceU = 0.;

        for (unsigned j = 0; j < dim; j++) {
          LaplaceU +=  nu * phi_x[i * dim + j] * gradSolUig[j];
        }

        aResU[i] += (phi[i] - LaplaceU) * weight;
		
      } // end phiU_i loop
    } // end gauss point loop

    // } // endif single element not refined or fine grid loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store them in RES
    ResU.resize(nDofsU);    //resize

    for (int i = 0; i < nDofsU; i++) {
      ResU[i] = -aResU[i].value();
    }

    RES->add_vector_blocked(ResU, sysDof);

    //Extarct and store the Jacobian

    Jac.resize(nDofsU * nDofsU);
    // define the dependent variables
    s.dependent(&aResU[0], nDofsU);

    // define the independent variables
    s.independent(&solU[0], nDofsU);

    // get the and store jacobian matrix (row-major)
    s.jacobian(&Jac[0] , true);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();

  KK->close();

  // ***************** END ASSEMBLY *******************
}


