/** \file Ex11.cpp
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
#include "LinearImplicitSystem.hpp"
#include "adept.h"
#include "FieldSplitTree.hpp"


using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;

  if(!strcmp(SolName, "T0")) {
    if(facename == 1) {
      dirichlet = false;
    }

    if(facename == 2 || facename == 3) {  // Dirichlet
      value = 1.;

      if(x[1] > 1.99999) value = 0.;
    }
  }
//   else if(!strcmp(SolName, "Ts")) {
//     if(facename == 1) {
//       dirichlet = false;
//     }
//
//     if(facename == 2 || facename == 3) {  // Dirichlet
//       value = 1.;
//
//       if(x[1] > 1.99999) value = 0.;
//     }
//   }

  else if(!strcmp(SolName, "U")) {
    if(facename == 3) {
      dirichlet = false; //Neumann
    }
  }
  else if(!strcmp(SolName, "V")) {
    if(facename == 1 || facename == 2) {  // Dirichlet 0
      value = 0.;

      if(x[0] > 0.2499 && x[0] < 0.7501)
        value = 0.375;
    }
    else if(facename == 3) {
      dirichlet = false; //Neumann
    }
  }
  else if(!strcmp(SolName, "P")) {
    dirichlet = false;
  }

  return dirichlet;
}


void AssembleNavierStokes(MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );

void AssembleTemperature(MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );


enum PrecType {
  FS_VTp = 1,
  FS_TVp,
  ASM_VTp,
  ASM_TVp,
  ILU_VTp,
  ILU_TVp,
};

int main(int argc, char** args) {

  unsigned precType = 0;

  if(argc >= 2) {
    if(!strcmp("FS_VT", args[1])) precType = FS_VTp;
    else if(!strcmp("FS_TV", args[1])) precType = FS_TVp;
    else if(!strcmp("ASM_VT", args[1])) precType = ASM_VTp;
    else if(!strcmp("ASM_TV", args[1])) precType = ASM_TVp;
    else if(!strcmp("ILU_VT", args[1])) precType = ILU_VTp;

    if(!strcmp("ILU_TV", args[1])) precType = ILU_TVp;

    if(precType == 0) {
      std::cout << "wrong input arguments!" << std::endl;
      abort();
    }
  }
  else {
    std::cout << "No input argument set default preconditioner = NS+T" << std::endl;
    precType == FS_VTp;
  }

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  //mlMsh.ReadCoarseMesh("./input/cube_hex.neu","seventh",scalingFactor);
  mlMsh.ReadCoarseMesh("./input/temp_control.neu", "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 7;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  // erase all the coarse mesh levels
  //mlMsh.EraseCoarseLevels(0);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol

  mlSol.AddSolution("U", LAGRANGE, SECOND);
  mlSol.AddSolution("V", LAGRANGE, SECOND);

  if(dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND);

  mlSol.AddSolution("P",  DISCONTINUOUS_POLYNOMIAL, FIRST);

  mlSol.AssociatePropertyToSolution("P", "Pressure");

  mlSol.AddSolution("Ts", LAGRANGE, SECOND);
  mlSol.AddSolution("Lambda", LAGRANGE, SECOND);
  mlSol.AddSolution("T0", LAGRANGE, SECOND);

  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  //mlSol.FixSolutionAtOnePoint("P");
  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("NS");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");

  if(dim == 3) system.AddSolutionToSystemPDE("W");

  system.AddSolutionToSystemPDE("P");

  system.SetLinearEquationSolverType(FEMuS_ASM); // Additive Swartz preconditioner

  // attach the assembling function to system
  system.SetAssembleFunction(AssembleNavierStokes);

  system.SetMaxNumberOfNonLinearIterations(10);
  system.SetNonLinearConvergenceTolerance(1.e-8);
  system.SetMaxNumberOfResidualUpdatesForNonlinearIteration(10);
  system.SetResidualUpdateConvergenceTolerance(1.e-15);

  //system.SetMaxNumberOfLinearIterations(10);
  //system.SetAbsoluteLinearConvergenceTolerance(1.e-15);

  system.SetMgType(F_CYCLE);

  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);
  // initilaize and solve the system
  system.init();

  system.SetSolverFineGrids(RICHARDSON);
  system.SetPreconditionerFineGrids(ILU_PRECOND);

  system.SetTolerances(1.e-5, 1.e-20, 1.e+50, 20, 20);

  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");
  system.SetNumberOfSchurVariables(1);
  system.SetElementBlockNumber(4);

  system.MGsolve();

  LinearImplicitSystem& system2 = mlProb.add_system < LinearImplicitSystem > ("T");

  // add solution to system2
  system2.AddSolutionToSystemPDE("Ts");
  system2.AddSolutionToSystemPDE("Lambda");
  system2.AddSolutionToSystemPDE("T0");


  //BEGIN buid fieldSplitTree (only for FieldSplitPreconditioner)
  std::vector < unsigned > fieldTs(1);
  fieldTs[0] = system2.GetSolPdeIndex("Ts");
  std::vector < unsigned > solutionTypeTs(1);
  solutionTypeTs[0] = mlSol.GetSolutionType("Ts");
  FieldSplitTree FS_Ts( PREONLY, ASM_PRECOND, fieldTs, solutionTypeTs, "Ts");
  FS_Ts.SetAsmBlockSize(4);
  FS_Ts.SetAsmNumeberOfSchurVariables(0);


//   std::vector < unsigned > fieldLambda(1);
//   fieldLambda[0] = system2.GetSolPdeIndex("Lambda");
//   //fieldLambda[1] = system2.GetSolPdeIndex("T0");
//   std::vector < unsigned > solutionTypeLambda(1);
//   solutionTypeLambda[0] = mlSol.GetSolutionType("Lambda");
//   //solutionTypeLambda[1] = mlSol.GetSolutionType("T0");
//   FieldSplitTree FS_Lambda( PREONLY, ASM_PRECOND, fieldLambda, solutionTypeLambda, "Lambda");
//   FS_Lambda.SetAsmBlockSize(4);
//   FS_Lambda.SetAsmNumeberOfSchurVariables(0);

  std::vector < unsigned > fieldT0(2);
  fieldT0[0] = system2.GetSolPdeIndex("T0");
  fieldT0[1] = system2.GetSolPdeIndex("Lambda");
  std::vector < unsigned > solutionTypeT0(2);
  solutionTypeT0[0] = mlSol.GetSolutionType("T0");
  solutionTypeT0[1] = mlSol.GetSolutionType("Lambda");
  FieldSplitTree FS_T0( PREONLY, ASM_PRECOND, fieldT0, solutionTypeT0, "T0");
  FS_T0.SetAsmBlockSize(4);
  FS_T0.SetAsmNumeberOfSchurVariables(0);

  std::vector < FieldSplitTree *> FS3;
  FS3.reserve(2);
  FS3.push_back(&FS_Ts);
  //FS3.push_back(&FS_Lambda);
  FS3.push_back(&FS_T0);

  FieldSplitTree FsAll( RICHARDSON, FIELDSPLIT_PRECOND, FS3, "Temperature");

  //system2.SetLinearEquationSolverType(FEMuS_FIELDSPLIT);

  system2.SetLinearEquationSolverType(FEMuS_ASM);
  system2.SetAssembleFunction(AssembleTemperature);

  system2.SetMaxNumberOfLinearIterations(10);
  system2.SetAbsoluteLinearConvergenceTolerance(1.e-15);

  system2.SetMgType(F_CYCLE);

  system2.SetNumberPreSmoothingStep(1);
  system2.SetNumberPostSmoothingStep(1);
  // initilaize and solve the system
  system2.init();

  system2.SetSolverFineGrids(RICHARDSON);
  system2.SetPreconditionerFineGrids(ILU_PRECOND);
  //if ( precType == FS_VTp || precType == FS_TVp )
  system2.SetFieldSplitTree(&FsAll);
  system2.SetTolerances(1.e-5, 1.e-20, 1.e+50, 100, 100);

  system2.ClearVariablesToBeSolved();
  system2.AddVariableToBeSolved("All");
  system2.SetNumberOfSchurVariables(0);
  system2.SetElementBlockNumber(4);
  //system.SetElementBlockNumber("All");

  system2.MGsolve();









  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  mlMsh.PrintInfo();

  return 0;
}


void AssembleNavierStokes(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*           msh         = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*           el          = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*   mlSol         = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*   sol         = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level];  // pointer to the equation (level) object

  bool assembleMatrix = mlPdeSys->GetAssembleMatrix();
  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  if(assembleMatrix) s.continue_recording();
  else s.pause_recording();

  SparseMatrix*   KK          = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*  RES         = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //solution variable
  vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object

  if(dim == 3) solVIndex[2] = mlSol->GetIndex("W");       // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"

  vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object

  if(dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W");

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object

  vector < vector < adept::adouble > >  solV(dim);    // local solution
  vector < adept::adouble >  solP; // local solution

  vector< vector < adept::adouble > > aResV(dim);    // local redidual vector
  vector< adept::adouble > aResP; // local redidual vector

  vector < vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for(unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(maxSize);
    aResV[k].reserve(maxSize);
    coordX[k].reserve(maxSize);
  }

  solP.reserve(maxSize);
  aResP.reserve(maxSize);


  vector <double> phiV;  // local test function
  vector <double> phiV_x; // local test function first order partial derivatives
  vector <double> phiV_xx; // local test function second order partial derivatives

  phiV.reserve(maxSize);
  phiV_x.reserve(maxSize * dim);
  phiV_xx.reserve(maxSize * dim2);

  double* phiP;
  double weight; // gauss point weight

  vector< int > sysDof; // local to global pdeSys dofs
  sysDof.reserve((dim + 1) *maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve((dim + 1) *maxSize);

  vector < double > Jac;
  Jac.reserve((dim + 1) *maxSize * (dim + 1) *maxSize);

  if(assembleMatrix)
    KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // element geometry type
    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs

    unsigned nDofsVP = dim * nDofsV + nDofsP;
    // resize local arrays
    sysDof.resize(nDofsVP);

    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsX);
    }

    solP.resize(nDofsP);

    for(unsigned  k = 0; k < dim; k++) {
      aResV[k].resize(nDofsV);    //resize
      std::fill(aResV[k].begin(), aResV[k].end(), 0);    //set aRes to zero
    }

    aResP.resize(nDofsP);    //resize
    std::fill(aResP.begin(), aResP.end(), 0);    //set aRes to zero

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
        sysDof[i + k * nDofsV] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // global to global mapping between solution node and solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);      // global extraction and local storage for the solution
      sysDof[i + dim * nDofsV] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof

      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }


    // start a new recording of all the operations involving adept::adouble variables
    if(assembleMatrix) s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV, phiV_x, phiV_xx);
      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point

      vector < adept::adouble > solV_gss(dim, 0);
      vector < vector < adept::adouble > > gradSolV_gss(dim);

      for(unsigned  k = 0; k < dim; k++) {
        gradSolV_gss[k].resize(dim);
        std::fill(gradSolV_gss[k].begin(), gradSolV_gss[k].end(), 0);
      }

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solV_gss[k] += phiV[i] * solV[k][i];
        }

        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            gradSolV_gss[k][j] += phiV_x[i * dim + j] * solV[k][i];
          }
        }
      }

      adept::adouble solP_gss = 0;

      for(unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }


      // *** phiV_i loop ***
      for(unsigned i = 0; i < nDofsV; i++) {
        vector < adept::adouble > NSV(dim, 0.);

        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            NSV[k]   +=  phiV_x[i * dim + j] * (gradSolV_gss[k][j] + gradSolV_gss[j][k]);
            NSV[k]   +=  phiV[i] * (solV_gss[j] * gradSolV_gss[k][j]);
          }
        }

        for(unsigned  k = 0; k < dim; k++) {
          NSV[k] += -solP_gss * phiV_x[i * dim + k];
        }

        for(unsigned  k = 0; k < dim; k++) {
          aResV[k][i] += - NSV[k] * weight;
        }
      } // end phiV_i loop

      // *** phiP_i loop ***
      for(unsigned i = 0; i < nDofsP; i++) {
        for(int k = 0; k < dim; k++) {
          aResP[i] += - (gradSolV_gss[k][k]) * phiP[i]  * weight;
        }
      } // end phiP_i loop

    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store them in RES
    Res.resize(nDofsVP);    //resize

    for(int i = 0; i < nDofsV; i++) {
      for(unsigned  k = 0; k < dim; k++) {
        Res[ i +  k * nDofsV ] = -aResV[k][i].value();
      }
    }

    for(int i = 0; i < nDofsP; i++) {
      Res[ i +  dim * nDofsV ] = -aResP[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);

    //Extarct and store the Jacobian
    if(assembleMatrix) {
      Jac.resize(nDofsVP * nDofsVP);
      // define the dependent variables

      for(unsigned  k = 0; k < dim; k++) {
        s.dependent(&aResV[k][0], nDofsV);
      }

      s.dependent(&aResP[0], nDofsP);

      // define the independent variables

      for(unsigned  k = 0; k < dim; k++) {
        s.independent(&solV[k][0], nDofsV);
      }

      s.independent(&solP[0], nDofsP);

      // get the and store jacobian matrix (row-major)
      s.jacobian(&Jac[0] , true);
      KK->add_matrix_blocked(Jac, sysDof, sysDof);

      s.clear_independents();
      s.clear_dependents();
    }
  } //end element loop for each process

  RES->close();

  if(assembleMatrix) {
    KK->close();
  }

  // ***************** END ASSEMBLY *******************
}


void AssembleTemperature(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use
  LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<LinearImplicitSystem> ("T");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*           msh         = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*           el          = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*   mlSol         = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*   sol         = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level];  // pointer to the equation (level) object

  bool assembleMatrix = mlPdeSys->GetAssembleMatrix();
  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  if(assembleMatrix) s.continue_recording();
  else s.pause_recording();

  SparseMatrix*   KK          = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*  RES         = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)


  vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object

  if(dim == 3) solVIndex[2] = mlSol->GetIndex("W");       // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //solution variable
  unsigned solTsIndex;
  unsigned solLambdaIndex;
  unsigned solT0Index;


  solTsIndex = mlSol->GetIndex("Ts");    // get the position of "U" in the ml_sol object
  solLambdaIndex = mlSol->GetIndex("Lambda");    // get the position of "V" in the ml_sol object
  solT0Index = mlSol->GetIndex("T0");      // get the position of "V" in the ml_sol object

  unsigned solType = mlSol->GetSolutionType(solTsIndex);    // get the finite element type for "u"

  vector < adept::adouble >  solTs;    // local solution
  vector < adept::adouble >  solLambda;    // local solution
  vector < adept::adouble >  solT0;    // local solution


  vector < adept::adouble > aResTs;    // local redidual vector
  vector < adept::adouble > aResLambda;    // local redidual vector
  vector < adept::adouble > aResT0;    // local redidual vector

  vector < vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for(unsigned  k = 0; k < dim; k++) {
    coordX[k].reserve(maxSize);
  }

  vector < vector < double > >  solV(dim);    // local solution

  for(unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(maxSize);
  }


  unsigned solTsPdeIndex;
  unsigned solLambdaPdeIndex;
  unsigned solT0PdeIndex;
  solTsPdeIndex = mlPdeSys->GetSolPdeIndex("Ts");
  solLambdaPdeIndex = mlPdeSys->GetSolPdeIndex("Lambda");
  solT0PdeIndex = mlPdeSys->GetSolPdeIndex("T0");

  solTs.reserve(maxSize);
  solLambda.reserve(maxSize);
  solT0.reserve(maxSize);

  aResTs.reserve(maxSize);
  aResLambda.reserve(maxSize);
  aResT0.reserve(maxSize);

  vector <double> phiT;  // local test function
  vector <double> phiT_x; // local test function first order partial derivatives
  vector <double> phiT_xx; // local test function second order partial derivatives

  phiT.reserve(maxSize);
  phiT_x.reserve(maxSize * dim);
  phiT_xx.reserve(maxSize * dim2);

  double* phiV;

  double weight; // gauss point weight

  vector< int > sysDof; // local to global pdeSys dofs
  sysDof.reserve(3 * maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve(3 * maxSize);

  vector < double > Jac;
  Jac.reserve(3 * maxSize * 3 * maxSize);

  if(assembleMatrix)
    KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned group = msh->GetElementGroup(iel);
    bool target = (group == 5) ? 1 : 0;

    // element geometry type
    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsT = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs


    unsigned nDofsAll = 3 * nDofsT;
    // resize local arrays
    sysDof.resize(nDofsAll);

    solTs.resize(nDofsT);
    solLambda.resize(nDofsT);
    solT0.resize(nDofsT);

    for(unsigned  k = 0; k < dim; k++) {
      coordX[k].resize(nDofsX);
      solV[k].resize(nDofsV);
    }

    aResTs.resize(nDofsT);    //resize
    std::fill(aResTs.begin(), aResTs.end(), 0);    //set aRes to zero

    aResLambda.resize(nDofsT);    //resize
    std::fill(aResLambda.begin(), aResLambda.end(), 0);    //set aRes to zero

    aResT0.resize(nDofsT);    //resize
    std::fill(aResT0.begin(), aResT0.end(), 0);    //set aRes to zero


    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsT; i++) {
      unsigned solTDof = msh->GetSolutionDof(i, iel, solType);    // global to global mapping between solution node and solution dof

      solTs[i] = (*sol->_Sol[solTsIndex])(solTDof);      // global extraction and local storage for the solution
      sysDof[i + 0 * nDofsT] = pdeSys->GetSystemDof(solTsIndex, solTsPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof

      solLambda[i] = (*sol->_Sol[solLambdaIndex])(solTDof);      // global extraction and local storage for the solution
      sysDof[i + 1 * nDofsT] = pdeSys->GetSystemDof(solLambdaIndex, solLambdaPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof

      solT0[i] = (*sol->_Sol[solT0Index])(solTDof);      // global extraction and local storage for the solution
      sysDof[i + 2 * nDofsT] = pdeSys->GetSystemDof(solT0Index, solT0PdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof

      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
      }
    }


    // start a new recording of all the operations involving adept::adouble variables
    if(assembleMatrix) s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solType]->Jacobian(coordX, ig, weight, phiT, phiT_x, phiT_xx);
      phiV = msh->_finiteElement[ielGeom][solVType]->GetPhi(ig);
      // evaluate the solution, the solution derivatives and the coordinates in the gauss point

      adept::adouble solTs_gss(0);
      adept::adouble solLambda_gss(0);
      adept::adouble solT0_gss(0);
      vector < adept::adouble > gradSolTs_gss(dim, 0);
      vector < adept::adouble > gradSolLambda_gss(dim, 0);
      vector < adept::adouble > gradSolT0_gss(dim, 0);


      for(unsigned i = 0; i < nDofsT; i++) {
        solTs_gss += phiT[i] * solTs[i];
        solLambda_gss += phiT[i] * solLambda[i];
        solT0_gss += phiT[i] * solT0[i];

        for(unsigned j = 0; j < dim; j++) {
          gradSolTs_gss[j] += phiT_x[i * dim + j] * solTs[i];
          gradSolLambda_gss[j] += phiT_x[i * dim + j] * solLambda[i];
          gradSolT0_gss[j] += phiT_x[i * dim + j] * solT0[i];
        }
      }


      vector < double > solV_gss(dim, 0);

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solV_gss[k] += phiV[i] * solV[k][i];
        }
      }

      double Pe = 100.;
      double alpha = 1000000000. * target;
      double beta = 1.;
      double gamma = 1.;
      double Td = 0.9;

      // *** phiT_i loop ***
      for(unsigned i = 0; i < nDofsT; i++) {
        adept::adouble  Ts = 0.;
        adept::adouble  Lambda = 0.;
        adept::adouble  T0 = 0.;

        for(unsigned j = 0; j < dim; j++) {
          Ts += phiT_x[i * dim + j] * gradSolTs_gss[j] / Pe;
          Ts += phiT_x[i * dim + j] * gradSolT0_gss[j] / Pe;
          Ts += phiT[i] * (solV_gss[j] * gradSolTs_gss[j]);
          Ts += phiT[i] * (solV_gss[j] * gradSolT0_gss[j]);
          Lambda += phiT_x[i * dim + j] * gradSolLambda_gss[j] / Pe;
          Lambda += solLambda_gss * (solV_gss[j] * phiT_x[i * dim + j]);
          T0 += phiT_x[i * dim + j] * gradSolT0_gss[j] * gamma;
          T0 += phiT_x[i * dim + j] * gradSolLambda_gss[j] / Pe;
          T0 += solLambda_gss * (solV_gss[j] * phiT_x[i * dim + j]);
        }
        Lambda += alpha * phiT[i] * (solTs_gss + solT0_gss - Td);
        T0 += beta * phiT[i] * solT0_gss;
        T0 += alpha * phiT[i] * (solTs_gss + solT0_gss - Td);

        aResTs[i] += - Ts * weight;
        aResLambda[i] += - Lambda * weight;
        aResT0[i] += - T0 * weight;

      } // end phiV_i loop

    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store them in RES
    Res.resize(nDofsAll);    //resize

    for(int i = 0; i < nDofsT; i++) {
      Res[ i +  0 * nDofsT ] = -aResTs[i].value();
    }

    for(int i = 0; i < nDofsT; i++) {
      Res[ i +  1 * nDofsT ] = -aResLambda[i].value();
    }

    for(int i = 0; i < nDofsT; i++) {
      Res[ i +  2 * nDofsT ] = -aResT0[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);

    //Extarct and store the Jacobian
    if(assembleMatrix) {
      Jac.resize(nDofsAll * nDofsAll);
      // define the dependent variables
      s.dependent(&aResTs[0], nDofsT);
      s.dependent(&aResLambda[0], nDofsT);
      s.dependent(&aResT0[0], nDofsT);

      // define the independent variables
      s.independent(&solTs[0], nDofsT);
      s.independent(&solLambda[0], nDofsT);
      s.independent(&solT0[0], nDofsT);

      // get the and store jacobian matrix (row-major)
      s.jacobian(&Jac[0] , true);
      KK->add_matrix_blocked(Jac, sysDof, sysDof);

      s.clear_independents();
      s.clear_dependents();
    }
  } //end element loop for each process

  RES->close();

  if(assembleMatrix) {
    KK->close();
  }

  // ***************** END ASSEMBLY *******************
}
