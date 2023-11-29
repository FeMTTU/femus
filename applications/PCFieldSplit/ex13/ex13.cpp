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
#include "SparseMatrix.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "LinearEquationSolver.hpp"
#include "FieldSplitTree.hpp"

#include "adept.h"


using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;
  if (!strcmp(SolName, "V")) {
    if (facename == 2) {
      if(x[1]>-0.5 && x[1]<0.5){
	value = 1.;
      }
    }
  } else if (!strcmp(SolName, "P")) {
    dirichlet = false;
  }
  return dirichlet;
}


void AssembleBoussinesqAppoximation_AD(MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;

  const std::string relative_path_to_build_directory =  "../../../";
  const std::string mesh_file = relative_path_to_build_directory + Files::mesh_folder_path() + "01_gambit/02_2d/square/minus0p5-plus0p5_minus0p5-plus0p5/square_2x2_quad_Three_face_groups.neu";

  mlMsh.ReadCoarseMesh(mesh_file.c_str(), "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 7;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  //mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, SECOND);
  mlSol.AddSolution("V", LAGRANGE, SECOND);

  if (dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND);

  //mlSol.AddSolution("P", LAGRANGE, FIRST);
  mlSol.AddSolution("P",  DISCONTINUOUS_POLYNOMIAL, FIRST);

  mlSol.AssociatePropertyToSolution("P", "Pressure");
  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.FixSolutionAtOnePoint("P");
  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("NS");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  system.AddSolutionToSystemPDE("P");

  if (dim == 3) system.AddSolutionToSystemPDE("W");

//   std::vector < unsigned > solutionTypeU(1);
//   solutionTypeU[0] = mlSol.GetSolutionType("U");
//   std::vector < unsigned > fieldU(1);
//   fieldU[0] = system.GetSolPdeIndex("U");
//   FieldSplitTree FS_U( PREONLY, ASM_PRECOND, fieldU, solutionTypeU,  "VelU");
//   FS_U.SetAsmBlockSize(4);
//   FS_U.SetAsmNumeberOfSchurVariables(0);
//
//   std::vector < unsigned > solutionTypeV(1);
//   solutionTypeV[0] = mlSol.GetSolutionType("V");
//   std::vector < unsigned > fieldV(1);
//   fieldV[0] = system.GetSolPdeIndex("V");
//   FieldSplitTree FS_V( PREONLY, ASM_PRECOND, fieldV, solutionTypeV,  "VelV");
//   FS_V.SetAsmBlockSize(4);
//   FS_V.SetAsmNumeberOfSchurVariables(0);
//
//   std::vector < FieldSplitTree *> FS2;
//   FS2.reserve(2);
//   FS2.push_back(&FS_U);
//   FS2.push_back(&FS_V);
//
//   FieldSplitTree FS_UV( PREONLY, FIELDSPLIT_PRECOND, FS2, "Velocity");

  std::vector < unsigned > solutionTypeUV(2);
  solutionTypeUV[0] = mlSol.GetSolutionType("U");
  solutionTypeUV[1] = mlSol.GetSolutionType("V");

  std::vector < unsigned > fieldUV(2);
  fieldUV[0] = system.GetSolPdeIndex("U");
  fieldUV[1] = system.GetSolPdeIndex("V");
  FieldSplitTree FS_UV( PREONLY, ASM_PRECOND, fieldUV, solutionTypeUV,  "Velocity");
  FS_UV.SetAsmBlockSize(4);
  FS_UV.SetAsmNumeberOfSchurVariables(0);

  std::vector < unsigned > solutionTypeP(1);
  solutionTypeP[0] = mlSol.GetSolutionType("P");
  std::vector < unsigned > fieldP(1);
  fieldP[0] = system.GetSolPdeIndex("P");
  //FieldSplitTree FS_P(GMRES, LSC_PRECOND, fieldP, "Pressure");// It works, but it is slower than ILU_PRECOND
  FieldSplitTree FS_P(PREONLY, ASM_PRECOND, fieldP, solutionTypeP, "Pressure");// It works, but it is slower that Vanka-ASM
  FS_P.SetAsmBlockSize(4);
  FS_P.SetAsmNumeberOfSchurVariables(0);


  std::vector < FieldSplitTree *> FS1;
  FS1.reserve(2);
  FS1.push_back(&FS_UV);
  FS1.push_back(&FS_P);
  FieldSplitTree FS_NS(RICHARDSON, FIELDSPLIT_SCHUR_PRECOND, FS1, "Navier-Stokes");


  //system.SetLinearEquationSolverType(FEMuS_DEFAULT);
  //system.SetLinearEquationSolverType(FEMuS_ASM); // Additive Swartz Method
  system.SetLinearEquationSolverType(FEMuS_FIELDSPLIT); // Additive Swartz Method

  // attach the assembling function to system
  system.SetAssembleFunction(AssembleBoussinesqAppoximation_AD);


  system.SetMaxNumberOfNonLinearIterations(10);
  system.SetNonLinearConvergenceTolerance(1.e-8);
  system.SetMaxNumberOfResidualUpdatesForNonlinearIteration(10);
  system.SetResidualUpdateConvergenceTolerance(1.e-12);

  system.SetMgType(F_CYCLE);

  system.SetNumberPreSmoothingStep(0);
  system.SetNumberPostSmoothingStep(2);
  // initilaize and solve the system
  system.init();

  system.SetSolverFineGrids(RICHARDSON);
  system.SetPreconditionerFineGrids(ILU_PRECOND);
  system.SetFieldSplitTree(&FS_NS);
  system.SetTolerances(1.e-10, 1.e-20, 1.e+50, 20, 5);

  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");
  system.SetNumberOfSchurVariables(1);
  system.SetElementBlockNumber(4);

  system.MGsolve();

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);
  vtkIO.Write(Files::_application_output_directory, "biquadratic", variablesToBePrinted);

  return 0;
}


void AssembleBoussinesqAppoximation_AD(MultiLevelProblem& ml_prob) {
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

  SparseMatrix*   KK          = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*  RES         = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27


  bool assembleMatrix = mlPdeSys->GetAssembleMatrix();
  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;
  if( assembleMatrix ) s.continue_recording();
  else s.pause_recording();



  //solution variable

  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object

  if (dim == 3) solVIndex[2] = mlSol->GetIndex("W");      // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"

  std::vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object

  if (dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W");

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object

  std::vector < std::vector < adept::adouble > >  solV(dim);    // local solution
  std::vector < adept::adouble >  solP; // local solution

  std::vector < std::vector < adept::adouble > > aResV(dim);    // local redidual vector
  std::vector < adept::adouble > aResP; // local redidual vector

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(maxSize);
    aResV[k].reserve(maxSize);
    coordX[k].reserve(maxSize);
  }

  solP.reserve(maxSize);
  aResP.reserve(maxSize);

  std::vector <double> phiV;  // local test function
  std::vector <double> phiV_x; // local test function first order partial derivatives
  std::vector <double> phiV_xx; // local test function second order partial derivatives

  phiV.reserve(maxSize);
  phiV_x.reserve(maxSize * dim);
  phiV_xx.reserve(maxSize * dim2);

  double* phiP;
  double weight; // gauss point weight

  std::vector < int > sysDof; // local to global pdeSys dofs
  sysDof.reserve((dim + 1) *maxSize);

  std::vector < double > Res; // local redidual vector
  Res.reserve((dim + 1) *maxSize);

  std::vector < double > Jac;
  Jac.reserve((dim + 1) *maxSize * (dim + 1) *maxSize);

  if(assembleMatrix)
    KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    // element geometry type
    short unsigned ielGeom = msh->GetElementType(iel);
    
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs

    unsigned nDofsVP =  dim * nDofsV + nDofsP;
    // resize local arrays
    sysDof.resize(nDofsVP);

    for (unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsX);
    }
    solP.resize(nDofsP);

    for (unsigned  k = 0; k < dim; k++) {
      aResV[k].resize(nDofsV);    //resize
      std::fill(aResV[k].begin(), aResV[k].end(), 0);    //set aRes to zero
    }

    aResP.resize(nDofsP);    //resize
    std::fill(aResP.begin(), aResP.end(), 0);    //set aRes to zero

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
        sysDof[i + k * nDofsV] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    for (unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // global to global mapping between solution node and solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);      // global extraction and local storage for the solution
      sysDof[i + dim * nDofsV] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }


      // start a new recording of all the operations involving adept::adouble variables
      if(assembleMatrix) s.new_recording();

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {
        // *** get gauss point weight, test function and test function partial derivatives ***
        msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV, phiV_x, phiV_xx);
        phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

        // evaluate the solution, the solution derivatives and the coordinates in the gauss point

        std::vector < adept::adouble > solV_gss(dim, 0);
        std::vector < std::vector < adept::adouble > > gradSolV_gss(dim);

        for (unsigned  k = 0; k < dim; k++) {
          gradSolV_gss[k].resize(dim);
          std::fill(gradSolV_gss[k].begin(), gradSolV_gss[k].end(), 0);
        }

        for (unsigned i = 0; i < nDofsV; i++) {
          for (unsigned  k = 0; k < dim; k++) {
            solV_gss[k] += phiV[i] * solV[k][i];
          }

          for (unsigned j = 0; j < dim; j++) {
            for (unsigned  k = 0; k < dim; k++) {
              gradSolV_gss[k][j] += phiV_x[i * dim + j] * solV[k][i];
            }
          }
        }

        adept::adouble solP_gss = 0;

        for (unsigned i = 0; i < nDofsP; i++) {
          solP_gss += phiP[i] * solP[i];
        }

        double nu = 1.0;
        double alpha = 1.;
        double beta = 2000.;


        // *** phiV_i loop ***
        for (unsigned i = 0; i < nDofsV; i++) {
          std::vector < adept::adouble > NSV(dim, 0.);

          for (unsigned j = 0; j < dim; j++) {
            for (unsigned  k = 0; k < dim; k++) {
              NSV[k]   +=  nu * phiV_x[i * dim + j] * (gradSolV_gss[k][j]);// + gradSolV_gss[j][k]);
              NSV[k]   +=  phiV[i] * (solV_gss[j] * gradSolV_gss[k][j]);
            }
          }

          for (unsigned  k = 0; k < dim; k++) {
            NSV[k] += -solP_gss * phiV_x[i * dim + k];
          }

          for (unsigned  k = 0; k < dim; k++) {
            aResV[k][i] += - NSV[k] * weight;
          }
        } // end phiV_i loop

        // *** phiP_i loop ***
        for (unsigned i = 0; i < nDofsP; i++) {
          for (int k = 0; k < dim; k++) {
            aResP[i] += - (gradSolV_gss[k][k]) * phiP[i]  * weight;
          }
        } // end phiP_i loop

      } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store them in RES
    Res.resize(nDofsVP);    //resize

    for (int i = 0; i < nDofsV; i++) {
      for (unsigned  k = 0; k < dim; k++) {
        Res[ i + k * nDofsV ] = -aResV[k][i].value();
      }
    }

    for (int i = 0; i < nDofsP; i++) {
      Res[ i + dim * nDofsV ] = -aResP[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);

    if(assembleMatrix){
      //Extarct and store the Jacobian

      Jac.resize(nDofsVP * nDofsVP);
      // define the dependent variables

      for (unsigned  k = 0; k < dim; k++) {
        s.dependent(&aResV[k][0], nDofsV);
      }

      s.dependent(&aResP[0], nDofsP);

      // define the independent variables

      for (unsigned  k = 0; k < dim; k++) {
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

  if(assembleMatrix)
    KK->close();

  // ***************** END ASSEMBLY *******************
}
