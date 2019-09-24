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
#include "adept.h"
#include "FieldSplitTree.hpp"
#include <stdlib.h>

double Re = 0.01;
double Rem = 1.0;
double coeffS = 1.0;
double pi = acos(-1.);
// double Mu = 0.01;
// double Miu = 0.001;  int c0=2; int cn=6; //Re=1000;
// int counter = 0 ;

using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time)
{
  bool dirichlet = true; //dirichlet
  value = 0.;

  if (!strcmp(SolName, "P")) {
    dirichlet = false;
//     if(facename == 1) value = -(2.0 * x[1] - 1.);
//     else if (facename == 2) value = -(2.0 * x[0] - 1.);
//     else if (facename == 3) value = 2.0 * x[1] - 1.;
//     else if (facename == 4) value = 2.0 * x[0] - 1.;
  }
  else if (!strcmp(SolName, "B1")) {
    if(facename == 2) value = sin(pi * x[0]);
    else if (facename == 4) value = -sin(pi * x[0]);
  }
  else if (!strcmp(SolName, "B2")) {
    if(facename == 1) value = - sin(pi * x[1]);
    else if (facename == 3) value = sin(pi * x[1]);
  }
  else if (!strcmp(SolName, "R")) {
    dirichlet = false;
  }
  return dirichlet;
}

void PrintConvergenceInfo(char *stdOutfile, char* infile, const unsigned &numofrefinements);
void AssembleBoussinesqAppoximation(MultiLevelProblem& ml_prob);
// vector <double> GetExactSolutionValue(const std::vector < double >& x);
int main(int argc, char** args)
{

//   unsigned precType = 0;
//
//   if(argc >= 2) {
//     if(!strcmp("FS_VT", args[1])) precType = FS_VTp;
//     else if(!strcmp("FS_TV", args[1])) precType = FS_TVp;
//     else if(!strcmp("ASM_VT", args[1])) precType = ASM_VTp;
//     else if(!strcmp("ASM_TV", args[1])) precType = ASM_TVp;
//     else if(!strcmp("ILU_VT", args[1])) precType = ILU_VTp;
//
//     if(!strcmp("ILU_TV", args[1])) precType = ILU_TVp;
//
//     if(precType == 0) {
//       std::cout << "wrong input arguments!" << std::endl;
//       abort();
//     }
//   }
//   else {
//     std::cout << "No input argument set default preconditioner = NS+T" << std::endl;
//     precType = FS_VTp;
//   }
//
//   if(argc >= 3) {
//     Prandtl = strtod(args[2], NULL);
//     std::cout << Prandtl<<std::endl;
//   }
//
//
//   if(argc >= 4) {
//     Rayleigh = strtod(args[3], NULL);
//     std::cout << Rayleigh <<std::endl;
//   }

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  //mlMsh.ReadCoarseMesh("./input/cube_hex.neu","seventh",scalingFactor);
  mlMsh.ReadCoarseMesh("./input/quad_square.neu", "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 6;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(3);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol

  mlSol.AddSolution("B1", LAGRANGE, SECOND);
  mlSol.AddSolution("B2", LAGRANGE, SECOND);
//  mlSol.AddSolution("R",  LAGRANGE, FIRST);
  mlSol.AddSolution("R",  DISCONTINUOUS_POLYNOMIAL, FIRST);
  mlSol.AssociatePropertyToSolution("R", "Pressure");

  mlSol.AddSolution("U", LAGRANGE, SECOND);
  mlSol.AddSolution("V", LAGRANGE, SECOND);
  //mlSol.AddSolution("P", LAGRANGE, FIRST);
  mlSol.AddSolution("P",  DISCONTINUOUS_POLYNOMIAL, FIRST);
  mlSol.AssociatePropertyToSolution("P", "Pressure");
  mlSol.Initialize("All");

//   mlSol.Initialize("T",InitalValueT);
  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.FixSolutionAtOnePoint("R");
  mlSol.FixSolutionAtOnePoint("P");
  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);
  // add system Poisson in mlProb as a Linear Implicit System
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("NS");

  system.AddSolutionToSystemPDE("B1");
  system.AddSolutionToSystemPDE("B2");
  system.AddSolutionToSystemPDE("R");

  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  system.AddSolutionToSystemPDE("P");

  //BEGIN buid fieldSplitTree (only for FieldSplitPreconditioner)
  std::vector <unsigned> fieldB1B2R(3);
  fieldB1B2R[0] = system.GetSolPdeIndex("B1");
  fieldB1B2R[1] = system.GetSolPdeIndex("B2");
  fieldB1B2R[2] = system.GetSolPdeIndex("R");

  std::vector <unsigned> solutionTypeB1B2R(3);
  solutionTypeB1B2R[0] = mlSol.GetSolutionType("B1");
  solutionTypeB1B2R[1] = mlSol.GetSolutionType("B2");
  solutionTypeB1B2R[2] = mlSol.GetSolutionType("R");

  FieldSplitTree FS_MF(PREONLY, ASM_PRECOND, fieldB1B2R, solutionTypeB1B2R, "Magnetic Field");
  FS_MF.SetAsmBlockSize(4);

  FS_MF.SetAsmBlockPreconditioner(MLU_PRECOND);
  FS_MF.SetAsmNumeberOfSchurVariables(1);

  std::vector < unsigned > fieldUVP(3);
  fieldUVP[0] = system.GetSolPdeIndex("U");
  fieldUVP[1] = system.GetSolPdeIndex("V");
  fieldUVP[2] = system.GetSolPdeIndex("P");

  std::vector < unsigned > solutionTypeUVP(3);
  solutionTypeUVP[0] = mlSol.GetSolutionType("U");
  solutionTypeUVP[1] = mlSol.GetSolutionType("V");
  solutionTypeUVP[2] = mlSol.GetSolutionType("P");

  FieldSplitTree FS_NS(PREONLY, ASM_PRECOND, fieldUVP, solutionTypeUVP, "Navier-Stokes");
  FS_NS.SetAsmBlockSize(4);

//   FS_NS.SetAsmBlockPreconditioner(MLU_PRECOND);
  FS_NS.SetAsmNumeberOfSchurVariables(1);

  std::vector < FieldSplitTree *> FS2;
  FS2.reserve(2);
  FS2.push_back(&FS_MF); // Magnetic Field  First
  FS2.push_back(&FS_NS); // Navier-stokes block Second

  FieldSplitTree FS_MHD(RICHARDSON, FIELDSPLIT_PRECOND, FS2, "MHD");
  FS_MHD.SetRichardsonScaleFactor(.6);

  //END buid fieldSplitTree
  system.SetLinearEquationSolverType(FEMuS_FIELDSPLIT); // Field-Split preconditioner
  // system.SetLinearEquationSolverType(FEMuS_ASM);  // Additive Swartz preconditioner
  // ILU preconditioner

  // attach the assembling function to system
  system.SetAssembleFunction(AssembleBoussinesqAppoximation);

  system.SetMaxNumberOfNonLinearIterations(10);
  system.SetNonLinearConvergenceTolerance(1.e-8);
  //system.SetMaxNumberOfResidualUpdatesForNonlinearIteration(2);
  //system.SetResidualUpdateConvergenceTolerance(1.e-12);

  //system.SetMaxNumberOfLinearIterations(10);
  //system.SetAbsoluteLinearConvergenceTolerance(1.e-15);

  system.SetMaxNumberOfLinearIterations(1);
  system.SetAbsoluteLinearConvergenceTolerance(1.e-15);


  system.SetMgType(F_CYCLE);

  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);
  // initilaize and solve the system
  system.init();

  system.SetSolverFineGrids(RICHARDSON);
  system.SetPreconditionerFineGrids(ILU_PRECOND);
  system.SetRichardsonScaleFactor(.6);

  system.SetFieldSplitTree(&FS_MHD);
  system.SetTolerances(1.e-5, 1.e-8, 1.e+50, 10, 10); //GMRES tolerances
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");
//   system.SetNumberOfSchurVariables(1);
//   system.SetElementBlockNumber(4);
  system.MGsolve();

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);


  vtkIO.SetDebugOutput(true);

  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  mlMsh.PrintInfo();

  //char infile[256]="FSVTtrueResidual.txt";
  //char stdOutfile[256]="output.txt";

  return 0;
}

void GetExactSolutionValue(const std::vector < double >& x, vector < double >& U,
                           double &P, vector < double >& B)
{
  U[0] = x[0] * x[0] * (x[0] - 1.0) * (x[0] - 1.0) * x[1] * (x[1] - 1.) * (2 * x[1] - 1.);
  U[1] = - x[0] * (x[0] - 1.0) * (2 * x[0] - 1) * x[1] * x[1] * (x[1] - 1.) * (x[1] - 1.);
  P = (2.0 * x[0] - 1.) * (2.0 * x[1] - 1.);
  B[0] = sin(pi * x[0]) * cos(pi * x[1]);
  B[1] = -sin(pi * x[1]) * cos(pi * x[0]);
}

void GetGradExactSolutionValue(const std::vector < double >& x, vector < double >& GradU, vector < double >& GradV,
                               vector < double >& GradB1, vector < double >& GradB2)
{
  GradU [0] = 2. * x[0] * (x[0] - 1.) * (2 * x[0] - 1.) * x[1] * (x[1] - 1.) * (2. * x[1] - 1.);
  GradU [1] = x[0] * x[0] * (x[0] - 1.) * (x[0] - 1.) * (6. * x[1] * x[1] - 6. * x[1] + 1.);
  GradV [0] = - (6. * x[0] * x[0] - 6. * x[0] + 1.) * x[1] * x[1] * (x[1] - 1.) * (x[1] - 1.);
  GradV [1] = - GradU [0];

  GradB1 [0] = pi * cos(pi * x[0]) * cos(pi * x[1]);
  GradB1 [1] = -pi * sin(pi * x[0]) * sin(pi * x[1]);
  GradB2 [0] = pi * sin(pi * x[0]) * sin(pi * x[1]);
  GradB2 [1] = - GradB1 [0];



  /*
    GradUx [1] = x[0] * x[0] * (x[0] - 1.0) * (x[0] - 1.0) * (x[1] - 1.0) * (2.0 * x[1] -1. )
  	      + x[0] * x[0] * (x[0] - 1.0) * (x[0] - 1.0) * (x[1] - 1.0) * (2.0 * x[1] -1. )
  	     + 2.0 * x[0] * x[0] * (x[0] - 1.0) * x[1] * (x[1] - 1.0) * (2.0 * x[1] -1. );



      USol[0]= x[0] * x[0] * (x[0]-1.0) * (x[0] - 1.0) * x[1] * (x[1] - 1.) * (2*x[1] -1.);
      Usol[1] = - x[0] * (x[0] - 1.0) * (2*x[0]-1) * x[1] * x[1] * (x[1] - 1.) * (x[1] -1.);
      Psol = (2.0 * x[0] - 1.) * (2.0 * x[1] - 1.);
      Bsol[0] = sin(pi*x[0]) * cos(pi*x[1]);
      Bsol[1] = -sin(pi*x[1]) * cos(pi*x[0]);*/
}

void AssembleBoussinesqAppoximation(MultiLevelProblem& ml_prob)
{
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


//   double Mu;
//   if(counter < c0 ) Mu = 2. * Miu;
//   else if ( counter <= cn ) {
//     Mu = 2*Miu*(cn-counter)/(cn-c0) + Miu*(counter-c0)/(cn-c0);
//   }
//   else{
//     Mu = Miu;
//   }
//   std::cout << counter << " " << Mu <<std::endl;
//   counter++;

//   Re = 1.0/Mu;
//

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

  SparseMatrix*   KK          = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*  RES         = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); 	// get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1)); // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //solution variable
  vector < unsigned > solBIndex(dim);
  solBIndex[0] = mlSol->GetIndex("B1");
  solBIndex[1] = mlSol->GetIndex("B2");
  unsigned solBType = mlSol -> GetSolutionType(solBIndex[0]);

  unsigned solRIndex;
  solRIndex = mlSol->GetIndex("R");
  unsigned solRType = mlSol->GetSolutionType(solRIndex);

  vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object
  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"

  vector < unsigned > solBPdeIndex(dim);
  solBPdeIndex[0] = mlPdeSys->GetSolPdeIndex("B1");    // get the position of "B1" in the pdeSys object
  solBPdeIndex[1] = mlPdeSys->GetSolPdeIndex("B2");    // get the position of "B2" in the pdeSys object
  unsigned solRPdeIndex;
  solRPdeIndex = mlPdeSys->GetSolPdeIndex("R");    // get the position of "R" in the pdeSys object

  vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object
  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object

  vector <vector <double> > solB(dim); // local solution for B1, B2
  vector <double> solR; // local solution for langrange for r
  vector < vector < double > >  solV(dim);    // local solution for U, V
  vector < double >  solP; // local solution for P

  vector < vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for(unsigned  k = 0; k < dim; k++) {
    solB[k].reserve(maxSize);
    coordX[k].reserve(maxSize);
  }
  solR.reserve(maxSize);

  for(unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(maxSize);
  }
  solP.reserve(maxSize);

  vector <double> phiB;  // local test function
  vector <double> phiB_x; // local test function first order partial derivatives
  vector <double> phiB_xx; // local test function second order partial derivatives

  phiB.reserve(maxSize);
  phiB_x.reserve(maxSize * dim);
  phiB_xx.reserve(maxSize * dim2);
  double* phiR;

  vector <double> phiV;  // local test function
  vector <double> phiV_x; // local test function first order partial derivatives
  vector <double> phiV_xx; // local test function second order partial derivatives

  phiV.reserve(maxSize);
  phiV_x.reserve(maxSize * dim);
  phiV_xx.reserve(maxSize * dim2);

  double* phiP;
  double weight; // gauss point weight

  vector< int > sysDof; // local to global pdeSys dofs
  sysDof.reserve((dim * 2 + 2) *maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve((dim * 2 + 2) *maxSize);

  vector < double > Jac; //local Jacobian
  Jac.reserve((dim * 2 + 2) *maxSize * (dim * 2 + 2) *maxSize);

  if(assembleMatrix)
    KK->zero(); // Set to zero all the entries of the Global Matrix

  //BEGIN element loop
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    //BEGIN local dof number extraction
    unsigned nDofsB = msh->GetElementDofNumber(iel, solBType);  //Magnetic
    unsigned nDofsR = msh->GetElementDofNumber(iel, solRType);  //pressure

    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);  //velocity
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);  //pressure

    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType); // coordinates
    unsigned nDofsBRVP = dim * nDofsB + nDofsR + dim * nDofsV + nDofsP; // all solutions
    //END local dof number extraction

    //BEGIN memory allocation
    Res.resize(nDofsBRVP);
    std::fill(Res.begin(), Res.end(), 0.0);
    Jac.resize(nDofsBRVP * nDofsBRVP);
    std::fill(Jac.begin(), Jac.end(), 0.0);
    sysDof.resize(nDofsBRVP);

    for(unsigned  k = 0; k < dim; k++) {
      solB[k].resize(nDofsB);
      coordX[k].resize(nDofsX);
    }
    solR.resize(nDofsR);

    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
    }
    solP.resize(nDofsP);
    //END memory allocation

    //BEGIN global to local extraction
    for(unsigned i = 0; i < nDofsB; i++) { //Magnetic field
      unsigned solBDof = msh->GetSolutionDof(i, iel, solBType);  //local to global solution dof
      for(unsigned  k = 0; k < dim; k++) {
        solB[k][i] = (*sol->_Sol[solBIndex[k]])(solBDof);  //global to local solution value
        sysDof[i + k * nDofsB] = pdeSys->GetSystemDof(solBIndex[k], solBPdeIndex[k], i, iel);  //local to global system dof
      }
    }

    for(unsigned i = 0; i < nDofsR; i++) { //Lagrange
      unsigned solRDof = msh->GetSolutionDof(i, iel, solRType);  //local to global solution dof
      solR[i] = (*sol->_Sol[solRIndex])(solRDof);  //global to local solution value
      sysDof[i + dim * nDofsB] = pdeSys->GetSystemDof(solRIndex, solRPdeIndex, i, iel);  //local to global system dof
    }

    for(unsigned i = 0; i < nDofsV; i++) { //velocity
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);  //local to global solution dof
      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);  //global to local solution value
        sysDof[i + dim * nDofsB + nDofsR + k * nDofsV] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);  //local to global system dof
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) { //pressure
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);  //local to global solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);  //global to local solution value
      sysDof[i + dim * nDofsB + nDofsR + dim * nDofsV] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);  //local to global system dof
    }

    for(unsigned i = 0; i < nDofsX; i++) { //coordinates
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);  //local to global coordinate dof

      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);  //global to local coordinate value
      }
    }

    //END global to local extraction

    //BEGIN Gauss point loop
    short unsigned ielGeom = msh->GetElementType(iel);

    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solBType]->Jacobian(coordX, ig, weight, phiB, phiB_x, phiB_xx);
      phiR = msh->_finiteElement[ielGeom][solRType]->GetPhi(ig);

      msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV, phiV_x, phiV_xx);
      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      vector <double> solB_gss(dim, 0);
      vector <vector <double> > gradSolB_gss(dim);
      for(unsigned j = 0; j < dim; j ++) {
        gradSolB_gss[j].resize(dim);
        std::fill(gradSolB_gss[j].begin(), gradSolB_gss[j].end(), 0.0);
      }

      for(unsigned i = 0; i < nDofsB; i++) {
        for(unsigned j = 0; j < dim; j++) { // define the solutions B1, B2 for solB_gss
          solB_gss[j] += phiB[i] * solB[j][i];
        }

        for(unsigned j = 0; j < dim; j++) { // define the soluitons B1, B2 for gradSolB_gss
          for(unsigned k = 0; k < dim; k++) { // define the direction of derivatives X, Y for gradSolB_gss
            gradSolB_gss[j][k] += phiB_x[i * dim + k] * solB[j][i];
          }
        }
      }

      double solR_gss = 0.0;
      for(unsigned i = 0; i < nDofsR; i++) {
        solR_gss += phiR[i] * solR[i];
      }

      vector < double > solV_gss(dim, 0);
      vector < vector < double > > gradSolV_gss(dim);
      for(unsigned j = 0; j < dim; j++) {
        gradSolV_gss[j].resize(dim);
        std::fill(gradSolV_gss[j].begin(), gradSolV_gss[j].end(), 0.0);
      }

      for(unsigned i = 0; i < nDofsV; i++) {
        for (unsigned j = 0; j < dim; j++) { // define the solutions U,V for solV_gass
          solV_gss[j] += phiV[i] * solV[j][i];
        }

        for(unsigned j = 0; j < dim; j++) { // define the solutions U,V for gradSolV_gss
          for(unsigned k = 0; k < dim; k++) { // define the direction of derivatives X, Y for gradSolV_gss
            gradSolV_gss[j][k] += phiV_x[i * dim + k] * solV[j][i];
          }
        }
      }

      double solP_gss = 0;
      for(unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }

      vector <double> x_gss(dim, 0);
      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned j = 0; j < dim; j ++) {
          x_gss[j] += phiV[i] *  coordX[j][i];
        }
      }


      std::vector < double > Ue(dim);
      double Pe;
      std::vector < double > Be(dim);

      std::vector < std::vector < double > > gradUe(dim);
      std::vector < std::vector < double > > gradBe(dim);
      for(unsigned j = 0; j < dim; j++) {
        gradUe[j].resize(dim);
        gradBe[j].resize(dim);
      }

      GetExactSolutionValue(x_gss, Ue, Pe, Be);
      GetGradExactSolutionValue(x_gss, gradUe[0], gradUe[1], gradBe[0], gradBe[1]);




      // Begin phiB_i loop: Momentum balance
      unsigned kk, ll;
      double coeffSignk, coeffSignl;
      for(unsigned i = 0; i < nDofsB; i ++) {
        for(unsigned k = 0; k < dim; k ++) {
          unsigned irow = k * nDofsB + i;
          for (unsigned l = 0; l < dim; l++) {
            if (k == 0) {
              kk = 1;
              coeffSignk = -1.0;
            }
            else {
              kk = 0;
              coeffSignk = 1.0;
            }
            if (l == 0) {
              ll = 1;
              coeffSignl = -1.0;
            }
            else {
              ll = 0;
              coeffSignl = 1.0;
            }
            Res[irow] -= coeffS * 1.0 / Rem * coeffSignk * phiB_x[i * dim + kk] * coeffSignl * (gradSolB_gss[l][ll] - gradBe[l][ll])  * weight;
            Res[irow] += coeffS * coeffSignk * phiB_x[i * dim + kk] * coeffSignl * (solB_gss[l] * solV_gss[ll] - Be[l] * Ue[ll])  * weight;
          }
          Res[irow] += phiB_x[i * dim + k] * solR_gss * weight;

          if (assembleMatrix) {
            unsigned irowMat = irow * nDofsBRVP;
            for(unsigned l = 0; l < dim; l++) {
              if (k == 0) {
                kk = 1;
                coeffSignk = -1.0;
              }
              else {
                kk = 0;
                coeffSignk = 1.0;
              }
              if (l == 0) {
                ll = 1;
                coeffSignl = -1.0;
              }
              else {
                ll = 0;
                coeffSignl = 1.0;
              }
              for(unsigned j = 0; j < nDofsB; j++) {
                unsigned jcol1 = l * nDofsB + j;
                unsigned jcol2 = dim * nDofsB + nDofsR + ll * nDofsV + j;
                Jac[irowMat + jcol1] += coeffS * 1.0 / Rem * coeffSignk * phiB_x[i * dim + kk] * coeffSignl * phiB_x[j * dim + ll] * weight;
                Jac[irowMat + jcol1] -= coeffS * coeffSignk * phiB_x[i * dim + kk] * coeffSignl * phiB[j] * solV_gss[ll] * weight;
                Jac[irowMat + jcol2] -= coeffS * coeffSignk * phiB_x[i * dim + kk] * coeffSignl * phiV[j] * solB_gss[l] * weight;
              }
            }
            for (unsigned j = 0; j < nDofsR; j++) {
              unsigned jcol = dim * nDofsB + j;
              Jac[irowMat + jcol] -= phiB_x[i * dim + k] * phiR[j] * weight;
            }
          }
        }
      }
      //END phiB_i loop

      // Begin phiR_i loop: Mass balance
      for (unsigned i = 0; i < nDofsR; i++) {
        unsigned irow = dim * nDofsB + i;
        for (unsigned k = 0; k < dim; k++) {
          Res[irow] += phiR[i] * gradSolB_gss[k][k]  * weight;
          if (assembleMatrix) {
            unsigned irowMat = irow * nDofsBRVP;
            for (unsigned j = 0; j < nDofsB; j++) {
              unsigned jcol = k * nDofsB + j;
              Jac[irowMat + jcol] += -phiR[i] * phiB_x[j * dim + k] * weight;
            }
          }
        }
      }
      // END phiR_i loop

      //BEGIN phiV_i loop: Momentum balance
      for (unsigned i = 0; i < nDofsV; i++) {
        for (unsigned k = 0; k < dim; k++) {
          unsigned irow = dim * nDofsB + nDofsR + k * nDofsV + i;
          for (unsigned l = 0; l < dim; l++) {
            if (k == 0) {
              kk = 1;
              coeffSignk = -1.0;
            }
            else {
              kk = 0;
              coeffSignk = 1.0;
            }
            if (l == 0) {
              ll = 1;
              coeffSignl = -1.0;
            }
            else {
              ll = 0;
              coeffSignl = 1.0;
            }
            Res[irow] +=  -1.0 / Re * phiV_x[i * dim + l] * (gradSolV_gss[k][l] + gradSolV_gss[l][k] - gradUe[k][l] - gradUe[l][k]) * weight;
            Res[irow] +=  -phiV[i] * ( solV_gss[l] * gradSolV_gss[k][l] - Ue[l] * gradUe[k][l]) * weight;
// 	    Res[irow] += coeffS * coeffSignk * phiV[i] * coeffSignl * gradSolB_gss[l][ll] * solB_gss[kk]* weight;
          }
          Res[irow] += (solP_gss - Pe) * phiV_x[i * dim + k] * weight;

          if (assembleMatrix) {
            unsigned irowMat = irow * nDofsBRVP;

            for (unsigned l = 0; l < dim; l++) {
              if (k == 0) {
                kk = 1;
                coeffSignk = -1.0;
              }
              else {
                kk = 0;
                coeffSignk = 1.0;
              }
              if (l == 0) {
                ll = 1;
                coeffSignl = -1.0;
              }
              else {
                ll = 0;
                coeffSignl = 1.0;
              }
              for (unsigned j = 0; j < nDofsV; j++) {
                unsigned jcol1 = dim * nDofsB + nDofsR + k * nDofsV + j;
                unsigned jcol2 = dim * nDofsB + nDofsR + l * nDofsV + j;
                unsigned jcol3 = ll * nDofsB + j;
                unsigned jcol4 = l * nDofsB + j;
                Jac[ irowMat + jcol1] += 1.0 / Re * phiV_x[i * dim + l] * phiV_x[j * dim + l] * weight;
                Jac[ irowMat + jcol2] += 1.0 / Re * phiV_x[i * dim + l] * phiV_x[j * dim + k] * weight;
                Jac[ irowMat + jcol1] += phiV[i] * solV_gss[l] * phiV_x[j * dim + l] * weight;
                Jac[ irowMat + jcol2] += phiV[i] * phiV[j] * gradSolV_gss[k][l] * weight;
// 		Jac[irowMat + jcol3] -= coeffS * coeffSignk * phiV[i] * phiB[j] * coeffSignl * gradSolB_gss[l][ll] * weight;
// 		Jac[irowMat + jcol4] -= coeffS * coeffSignk * phiV[i] * solB_gss[ll] * coeffSignl * phiB_x[j*dim+ll] * weight;
              }
            }

            for (unsigned j = 0; j < nDofsP; j++) {
              unsigned jcol =  dim * nDofsB + nDofsR + dim * nDofsV + j;
              Jac[ irowMat + jcol] += - phiV_x[i * dim + k] * phiP[j] * weight;
            }
          }
        }
      }
      //END phiV_i loop

      //BEGIN phiP_i loop: mass balance
      for (unsigned i = 0; i < nDofsP; i++) {
        unsigned irow = dim * nDofsB + nDofsR + dim * nDofsV + i;

        for (int k = 0; k < dim; k++) {
          Res[irow] += gradSolV_gss[k][k] * phiP[i]  * weight;

          if (assembleMatrix) {
            unsigned irowMat = nDofsBRVP * irow;
            for (unsigned j = 0; j < nDofsV; j++) {
              unsigned jcol = dim * nDofsB + nDofsR + k * nDofsV + j;
              Jac[ irowMat + jcol ] += - phiP[i] * phiV_x[j * dim + k] * weight;
            }
          }
        }
      }

      //END phiP_i loop
    }
    //END Gauss point loop

    //BEGIN local to global Matrix/Vector assembly
    RES->add_vector_blocked(Res, sysDof);

    if (assembleMatrix) {
      KK->add_matrix_blocked(Jac, sysDof, sysDof);
    }

    //END local to global Matrix/Vector assembly

  }

  //END element loop

  RES->close();

  if (assembleMatrix) {
    KK->close();
  }
  // ***************** END ASSEMBLY *******************
}

