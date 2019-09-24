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
#include "PetscMatrix.hpp"

double Prandtl = 0.02;
double Rayleigh = 10000.;

unsigned counter = 0;

using namespace femus;

// double InitalValueT0(const std::vector < double >& x){
//   return (x[0]+0.5);
// }

bool SetBoundaryCondition0(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;

  if(!strcmp(SolName, "T0")) {
    if(facename == 2) {
      value = 1.;
    } 
    else if(facename == 3) {
      dirichlet = false; //Neumann
    }
  } 
  else if(!strcmp(SolName, "P0")) {
    dirichlet = false;
  }
  else if(!strcmp(SolName, "T")) {
    if(facename == 2) {
      value = 0.;
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

void AssembleBoussinesqAppoximation0(MultiLevelProblem& ml_prob);
void AssembleBoussinesqAppoximation(MultiLevelProblem& ml_prob);

unsigned preconditioner = 0;
enum PrecType {
  FS_VTp = 1,
  FS_TVp,
  ASM_VTp,
  ASM_TVp,
  ILU_VTp,
  ILU_TVp
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
    precType = FS_VTp;
  }
  
  if(argc >= 3) {
    Prandtl = strtod(args[2], NULL);
    std::cout << Prandtl<<std::endl;
  }
  
  
  if(argc >= 4) {
    Rayleigh = strtod(args[3], NULL);
    std::cout << Rayleigh <<std::endl;
  }
    
  preconditioner = precType; // add by guoyi ke to denote the type of preconditioner in the assembling 
  
  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  //mlMsh.ReadCoarseMesh("./input/cube_hex.neu","seventh",scalingFactor);
  mlMsh.ReadCoarseMesh("./input/square_16quads.neu", "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 3;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(1);
  numberOfUniformLevels -= 1;
  // print mesh info
  mlMsh.PrintInfo();
  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("T0", LAGRANGE, SECOND);
  mlSol.AddSolution("U0", LAGRANGE, SECOND);
  mlSol.AddSolution("V0", LAGRANGE, SECOND);
  if(dim == 3) mlSol.AddSolution("W0", LAGRANGE, SECOND);

  mlSol.AddSolution("P0",  DISCONTINUOUS_POLYNOMIAL, FIRST);
  mlSol.AssociatePropertyToSolution("P0", "Pressure");
  
  mlSol.AddSolution("T", LAGRANGE, SECOND);
  mlSol.AddSolution("U", LAGRANGE, SECOND);
  mlSol.AddSolution("V", LAGRANGE, SECOND);
  if(dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND);

  mlSol.AddSolution("P",  DISCONTINUOUS_POLYNOMIAL, FIRST);
  mlSol.AssociatePropertyToSolution("P", "Pressure");
  mlSol.Initialize("All");
  //mlSol.Initialize("T0",InitalValueT0);

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition0);
  mlSol.FixSolutionAtOnePoint("P0");
  mlSol.FixSolutionAtOnePoint("P");
  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);
  // add system Poisson in mlProb as a Linear Implicit System
  NonLinearImplicitSystem& system0 = mlProb.add_system < NonLinearImplicitSystem > ("NS0");
  if(precType == FS_TVp || precType == ASM_TVp || precType == ILU_TVp)
    system0.AddSolutionToSystemPDE("T0");

  // add solution "u" to system
  system0.AddSolutionToSystemPDE("U0");
  system0.AddSolutionToSystemPDE("V0");
  if(dim == 3) system0.AddSolutionToSystemPDE("W0");
  
  if(precType == ASM_VTp) system0.AddSolutionToSystemPDE("T0");
  system0.AddSolutionToSystemPDE("P0");

  if(precType == FS_VTp || precType == ILU_VTp) system0.AddSolutionToSystemPDE("T0");
  //BEGIN buid fieldSplitTree (only for FieldSplitPreconditioner)
  std::vector < unsigned > fieldUVP0(3);
  fieldUVP0[0] = system0.GetSolPdeIndex("U0");
  fieldUVP0[1] = system0.GetSolPdeIndex("V0");
  fieldUVP0[2] = system0.GetSolPdeIndex("P0");

  std::vector < unsigned > solutionTypeUVP0(3);
  solutionTypeUVP0[0] = mlSol.GetSolutionType("U0");
  solutionTypeUVP0[1] = mlSol.GetSolutionType("V0");
  solutionTypeUVP0[2] = mlSol.GetSolutionType("P0");

  FieldSplitTree FS_NS0(PREONLY, ASM_PRECOND, fieldUVP0, solutionTypeUVP0, "Navier-Stokes");
  FS_NS0.SetAsmBlockSize(4); // size of vanka block is 4^4
  FS_NS0.SetAsmNumeberOfSchurVariables(1);

  std::vector < unsigned > fieldT0(1);
  fieldT0[0] = system0.GetSolPdeIndex("T0");

  std::vector < unsigned > solutionTypeT0(1);
  solutionTypeT0[0] = mlSol.GetSolutionType("T0");

  FieldSplitTree FS_T0( PREONLY, ASM_PRECOND, fieldT0, solutionTypeT0, "Temperature");

  FS_T0.SetAsmBlockSize(4); //size of vanka block is 4^4
  FS_T0.SetAsmNumeberOfSchurVariables(0);

  std::vector < FieldSplitTree *> FS20;
  FS20.reserve(2);
  if(precType == FS_VTp) FS20.push_back(&FS_NS0);   //Navier-Stokes block first
  FS20.push_back(&FS_T0);
  if(precType == FS_TVp) FS20.push_back(&FS_NS0);   //Navier-Stokes block last

  FieldSplitTree FS_NST0(RICHARDSON, FIELDSPLIT_PRECOND, FS20, "Benard");
  FS_NST0.SetRichardsonScaleFactor(.6);

  //END buid fieldSplitTree
  if(precType == FS_VTp || precType == FS_TVp) system0.SetLinearEquationSolverType(FEMuS_FIELDSPLIT);    // Field-Split preconditioned
  else if(precType == ASM_VTp || precType == ASM_TVp) system0.SetLinearEquationSolverType(FEMuS_ASM);  // Additive Swartz preconditioner
  else if(precType == ILU_VTp || precType == ILU_TVp) system0.SetLinearEquationSolverType(FEMuS_DEFAULT);
  // attach the assembling function to system
  system0.SetAssembleFunction(AssembleBoussinesqAppoximation0);

  system0.SetMaxNumberOfNonLinearIterations(10);
  system0.SetNonLinearConvergenceTolerance(1.e-8);
  //system.SetMaxNumberOfResidualUpdatesForNonlinearIteration(10);
  //system.SetResidualUpdateConvergenceTolerance(1.e-15);
  
  system0.SetMaxNumberOfLinearIterations(1);
  system0.SetAbsoluteLinearConvergenceTolerance(1.e-15);
  
  system0.SetMgType(F_CYCLE);

  system0.SetNumberPreSmoothingStep(1);
  system0.SetNumberPostSmoothingStep(1);
  // initilaize and solve the system
  system0.init();

  system0.SetSolverFineGrids(RICHARDSON);
  system0.SetPreconditionerFineGrids(ILU_PRECOND);
  system0.SetRichardsonScaleFactor(.6);

  if(precType == FS_VTp || precType == FS_TVp) system0.SetFieldSplitTree(&FS_NST0);

  system0.SetTolerances(1.e-5, 1.e-8, 1.e+50, 30, 30); //GMRES tolerances
  system0.ClearVariablesToBeSolved();
  system0.AddVariableToBeSolved("All");
  
  if(precType == ASM_VTp || precType == ASM_TVp){
    system0.SetNumberOfSchurVariables(1);
    system0.SetElementBlockNumber(4); // The number of block element 4^4
  }
  else if(precType == ILU_VTp || precType == ILU_TVp){
    system0.SetElementBlockNumber("All");
  }
  system0.MGsolve();
  
  
  std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"<<std::endl;
 // Prandtl = 0.005;
  
  ///////////////////////////////////////////////////////////////
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("NST");
  if(precType == FS_TVp || precType == ASM_TVp || precType == ILU_TVp)
    system.AddSolutionToSystemPDE("T");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  if(dim == 3) system.AddSolutionToSystemPDE("W");
  if(precType == ASM_VTp) system.AddSolutionToSystemPDE("T");
  system.AddSolutionToSystemPDE("P");

  if(precType == FS_VTp || precType == ILU_VTp) system.AddSolutionToSystemPDE("T");

  //BEGIN buid fieldSplitTree (only for FieldSplitPreconditioner)
  std::vector < unsigned > fieldUVP(3);
  fieldUVP[0] = system.GetSolPdeIndex("U");
  fieldUVP[1] = system.GetSolPdeIndex("V");
  fieldUVP[2] = system.GetSolPdeIndex("P");

  std::vector < unsigned > solutionTypeUVP(3);
  solutionTypeUVP[0] = mlSol.GetSolutionType("U");
  solutionTypeUVP[1] = mlSol.GetSolutionType("V");
  solutionTypeUVP[2] = mlSol.GetSolutionType("P");

  FieldSplitTree FS_NS(PREONLY, ASM_PRECOND, fieldUVP, solutionTypeUVP, "Navier-Stokes");
  FS_NS.SetAsmBlockSize(2);
  FS_NS.SetAsmNumeberOfSchurVariables(1);

  std::vector < unsigned > fieldT(1);
  fieldT[0] = system.GetSolPdeIndex("T");

  std::vector < unsigned > solutionTypeT(1);
  solutionTypeT[0] = mlSol.GetSolutionType("T");

  FieldSplitTree FS_T( PREONLY, ASM_PRECOND, fieldT, solutionTypeT, "Temperature");
  FS_T.SetAsmBlockSize(2);
  FS_T.SetAsmNumeberOfSchurVariables(0);

  std::vector < FieldSplitTree *> FS2;
  FS2.reserve(2);
  if(precType == FS_VTp) FS2.push_back(&FS_NS);   //Navier-Stokes block first
  
  FS2.push_back(&FS_T);
  if(precType == FS_TVp) FS2.push_back(&FS_NS);   //Navier-Stokes block last

  FieldSplitTree FS_NST(RICHARDSON, FIELDSPLIT_PRECOND, FS2, "Benard");
  FS_NST.SetRichardsonScaleFactor(.6);

  //END buid fieldSplitTree
  if(precType == FS_VTp || precType == FS_TVp) system.SetLinearEquationSolverType(FEMuS_FIELDSPLIT);    // Field-Split preconditioned
  else if(precType == ASM_VTp || precType == ASM_TVp) system.SetLinearEquationSolverType(FEMuS_ASM);  // Additive Swartz preconditioner
  else if(precType == ILU_VTp || precType == ILU_TVp) system.SetLinearEquationSolverType(FEMuS_DEFAULT);

  // attach the assembling function to system
  system.SetAssembleFunction(AssembleBoussinesqAppoximation);
  system.SetMaxNumberOfNonLinearIterations(1);
  system.SetNonLinearConvergenceTolerance(1.e-8);
  //system.SetMaxNumberOfResidualUpdatesForNonlinearIteration(10);
  //system.SetResidualUpdateConvergenceTolerance(1.e-15);
  
  system.SetMaxNumberOfLinearIterations(1);
  system.SetAbsoluteLinearConvergenceTolerance(1.e-15);	
  
  system.SetMgType(V_CYCLE);

  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);
  // initilaize and solve the system
  system.init();

  system.SetSolverFineGrids(RICHARDSON);
  system.SetPreconditionerFineGrids(ILU_PRECOND);
  
  system.SetRichardsonScaleFactor(.6);

  if(precType == FS_VTp || precType == FS_TVp) system.SetFieldSplitTree(&FS_NST);
  
  system.SetTolerances(1.e-5, 1.e-8, 1.e+50, 1, 1); //GMRES tolerances
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");
  
  if(precType == ASM_VTp || precType == ASM_TVp){
    system.SetNumberOfSchurVariables(1);
    system.SetElementBlockNumber(2);
  }
  else if(precType == ILU_VTp || precType == ILU_TVp){
    system.SetElementBlockNumber("All");
  }
  
 
  //////////////////////////////////////////////////////////////////////
  //solution variable
  unsigned solTIndex;
  solTIndex = mlSol.GetIndex("T");    // get the position of "T" in the ml_sol object
  unsigned solTType = mlSol.GetSolutionType(solTIndex);    // get the finite element type for "T"

  vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol.GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol.GetIndex("V");    // get the position of "V" in the ml_sol object
  if(dim == 3) solVIndex[2] = mlSol.GetIndex("W");       // get the position of "V" in the ml_sol object
  unsigned solVType = mlSol.GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  unsigned solPIndex;
  solPIndex = mlSol.GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol.GetSolutionType(solPIndex);    // get the finite element type for "u"
     
  Mesh* msh = mlMsh.GetLevel(numberOfUniformLevels-1);
  unsigned nprocs = msh->n_processors();
  unsigned sizeT = msh->_dofOffset[solTType][nprocs];
  unsigned sizeU = msh->_dofOffset[solVType][nprocs];
  unsigned sizeV = msh->_dofOffset[solVType][nprocs];
  unsigned sizeP = msh->_dofOffset[solPType][nprocs];
   
  unsigned sizeTUVP = sizeT + sizeU +sizeV +sizeP;
    
  Solution* sol = mlSol.GetLevel(numberOfUniformLevels-1);
    
  for(unsigned i=0; i< sizeTUVP; i++){
    system.SetOuterSolver(PREONLY);  
    system.MGsolve();
    mlSol.GenerateBdc("All");
    
    std::ofstream fout;
    if(i==0){
      fout.open( "preconditioner.txt");
    }
    else{
       fout.open( "preconditioner.txt",  std::ofstream::app);
    }
    if( !fout.is_open() ) {
      std::cout << std::endl << " The output file preconditioner cannot be opened.\n";
      abort();
    }
    
    if(precType == FS_TVp || precType == ASM_TVp || precType == ILU_TVp){
      for(unsigned j = 0; j < sizeT; j++ ){
	 fout << (*sol->_Sol[solTIndex])(j)<< " ";
      } 
    } 
      
    for(unsigned j = 0; j < sizeU; j++ ){
      fout << (*sol->_Sol[solVIndex[0]])(j) << " ";
    }
    
    for(unsigned j = 0; j < sizeV; j++ ){
      fout << (*sol->_Sol[solVIndex[1]])(j) << " ";
    }
    
    if (precType == ASM_VTp){
      for(unsigned j = 0; j < sizeT; j++ ){
	fout << (*sol->_Sol[solTIndex])(j)<< " ";
      }
    }
    
    for(unsigned j = 0; j < sizeP; j++ ){
      fout << (*sol->_Sol[solPIndex])(j)<< " ";
    }
    
    if (precType == FS_VTp || precType == ILU_VTp){
      for(unsigned j = 0; j < sizeT; j++ ){
	fout << (*sol->_Sol[solTIndex])(j)<< " ";
      }
    }
// std::cout << sizeT <<"AAA" << sizeU <<"BBB"<<sizeV<<"CCC" << sizeP<<"DDD"<<std::endl;   
    fout<<std::endl;
    fout.close();
    
    // print solutions
    std::vector < std::string > variablesToBePrinted;
    variablesToBePrinted.push_back("All");

    VTKWriter vtkIO(&mlSol);
    vtkIO.SetDebugOutput( true );
    vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, counter);
  }
  /////////////////////////////////////ultiLevelProb/////////////////////////
  
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput( true );
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);
  mlMsh.PrintInfo();
  
  return 0;
}


void AssembleBoussinesqAppoximation0(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the Mlem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("NS0");   // pointer to the linear implicit system named "Poisson"
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

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //solution variable
  unsigned solTIndex;
  solTIndex = mlSol->GetIndex("T0");    // get the position of "T" in the ml_sol object
  unsigned solTType = mlSol->GetSolutionType(solTIndex);    // get the finite element type for "T"

  vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U0");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V0");    // get the position of "V" in the ml_sol object

  if(dim == 3) solVIndex[2] = mlSol->GetIndex("W0");       // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("P0");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"

  unsigned solTPdeIndex;
  solTPdeIndex = mlPdeSys->GetSolPdeIndex("T0");    // get the position of "T" in the pdeSys object

  // std::cout << solTIndex <<" "<<solTPdeIndex<<std::endl;


  vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U0");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V0");    // get the position of "V" in the pdeSys object

  if(dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W0");

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("P0");    // get the position of "P" in the pdeSys object

  vector < double >  solT; // local solution
  vector < vector < double > >  solV(dim);    // local solution
  vector < double >  solP; // local solution

  vector < vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  solT.reserve(maxSize);

  for(unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(maxSize);
    coordX[k].reserve(maxSize);
  }

  solP.reserve(maxSize);

  vector <double> phiV;  // local test function
  vector <double> phiV_x; // local test function first order partial derivatives
  vector <double> phiV_xx; // local test function second order partial derivatives

  phiV.reserve(maxSize);
  phiV_x.reserve(maxSize * dim);
  phiV_xx.reserve(maxSize * dim2);

  vector <double> phiT;  // local test function
  vector <double> phiT_x; // local test function first order partial derivatives
  vector <double> phiT_xx; // local test function second order partial derivatives

  phiT.reserve(maxSize);
  phiT_x.reserve(maxSize * dim);
  phiT_xx.reserve(maxSize * dim2);

  double* phiP;
  double weight; // gauss point weight

  vector< int > sysDof; // local to global pdeSys dofs
  sysDof.reserve((dim + 2) *maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve((dim + 2) *maxSize);

  vector < double > Jac;
  Jac.reserve((dim + 2) *maxSize * (dim + 2) *maxSize);

  if(assembleMatrix)
    KK->zero(); // Set to zero all the entries of the Global Matrix

  //BEGIN element loop
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    //BEGIN local dof number extraction
    unsigned nDofsT = msh->GetElementDofNumber(iel, solTType);  //temperature
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);  //velocity
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);  //pressure
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType); // coordinates
    unsigned nDofsTVP = nDofsT + dim * nDofsV + nDofsP; // all solutions
    //END local dof number extraction

    //BEGIN memory allocation
    Res.resize(nDofsTVP);
    std::fill(Res.begin(), Res.end(), 0);
    Jac.resize(nDofsTVP * nDofsTVP);
    std::fill(Jac.begin(), Jac.end(), 0);
    sysDof.resize(nDofsTVP);

    solT.resize(nDofsV);

    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsX);
    }

    solP.resize(nDofsP);
    //END memory allocation

    //BEGIN global to local extraction
    for(unsigned i = 0; i < nDofsT; i++) { //temperature
      unsigned solTDof = msh->GetSolutionDof(i, iel, solTType);  //local to global solution dof
      solT[i] = (*sol->_Sol[solTIndex])(solTDof);  //global to local solution value
      sysDof[i] = pdeSys->GetSystemDof(solTIndex, solTPdeIndex, i, iel);  //local to global system dof
    }

    for(unsigned i = 0; i < nDofsV; i++) { //velocity
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);  //local to global solution dof

      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);  //global to local solution value
        sysDof[i + nDofsT + k * nDofsV] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);  //local to global system dof
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) { //pressure
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);  //local to global solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);  //global to local solution value
      sysDof[i + nDofsT + dim * nDofsV] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);  //local to global system dof
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
      msh->_finiteElement[ielGeom][solTType]->Jacobian(coordX, ig, weight, phiT, phiT_x, phiT_xx);
      msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV, phiV_x, phiV_xx);
      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      double solT_gss = 0;
      vector < double > gradSolT_gss(dim, 0.);

      for(unsigned i = 0; i < nDofsT; i++) {
        solT_gss += phiT[i] * solT[i];

        for(unsigned j = 0; j < dim; j++) {
          gradSolT_gss[j] += phiT_x[i * dim + j] * solT[i];
        }
      }

      vector < double > solV_gss(dim, 0);
      vector < vector < double > > gradSolV_gss(dim);

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

      double solP_gss = 0;

      for(unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }


      double alpha = 1.;
      double beta = 1.;//40000.;


      double Pr = Prandtl;
      double Ra = Rayleigh;

      //BEGIN phiT_i loop: Energy balance
      for(unsigned i = 0; i < nDofsT; i++) {
        unsigned irow = i;

        for(unsigned k = 0; k < dim; k++) {
          Res[irow] +=  -alpha / sqrt(Ra * Pr) * phiT_x[i * dim + k] * gradSolT_gss[k] * weight;
          Res[irow] +=  -phiT[i] * solV_gss[k] * gradSolT_gss[k] * weight;

          if(assembleMatrix) {
            unsigned irowMat = irow * nDofsTVP;

            for(unsigned j = 0; j < nDofsT; j++) {
              Jac[ irowMat + j ] +=  alpha / sqrt(Ra * Pr) * phiT_x[i * dim + k] * phiT_x[j * dim + k] * weight;
              Jac[ irowMat + j ] +=  phiT[i] * solV_gss[k] * phiT_x[j * dim + k] * weight;
            }

            for(unsigned j = 0; j < nDofsV; j++) {
              unsigned jcol = nDofsT + k * nDofsV + j;
              Jac[ irowMat + jcol ] += phiT[i] * phiV[j] * gradSolT_gss[k] * weight;
            }
          }

        }
      }

      //END phiT_i loop

      //BEGIN phiV_i loop: Momentum balance
      for(unsigned i = 0; i < nDofsV; i++) {

        for(unsigned k = 0; k < dim; k++) {
          unsigned irow = nDofsT + k * nDofsV + i;

          for(unsigned l = 0; l < dim; l++) {
            Res[irow] +=  -sqrt(Pr / Ra) * phiV_x[i * dim + l] * (gradSolV_gss[k][l] + gradSolV_gss[l][k]) * weight;
            Res[irow] +=  -phiV[i] * solV_gss[l] * gradSolV_gss[k][l] * weight;
          }

          Res[irow] += solP_gss * phiV_x[i * dim + k] * weight;

          if(k == 1) {
            Res[irow] += beta * solT_gss * phiV[i] * weight;
          }

          if(assembleMatrix) {
            unsigned irowMat = nDofsTVP * irow;

            for(unsigned l = 0; l < dim; l++) {
              for(unsigned j = 0; j < nDofsV; j++) {
                unsigned jcol1 = (nDofsT + k * nDofsV + j);
                unsigned jcol2 = (nDofsT + l * nDofsV + j);
                Jac[ irowMat + jcol1] += sqrt(Pr / Ra) * phiV_x[i * dim + l] * phiV_x[j * dim + l] * weight;
                Jac[ irowMat + jcol2] += sqrt(Pr / Ra) * phiV_x[i * dim + l] * phiV_x[j * dim + k] * weight;
                Jac[ irowMat + jcol1] += phiV[i] * solV_gss[l] * phiV_x[j * dim + l] * weight;
                Jac[ irowMat + jcol2] += phiV[i] * phiV[j] * gradSolV_gss[k][l] * weight;
              }
            }

            for(unsigned j = 0; j < nDofsP; j++) {
              unsigned jcol = (nDofsT + dim * nDofsV) + j;
              Jac[ irowMat + jcol] += - phiV_x[i * dim + k] * phiP[j] * weight;
            }

            if(k == 1) {
              for(unsigned j = 0; j < nDofsT; j++) {
                Jac[ irowMat + j ] += -beta * phiV[i] * phiT[j] * weight;
              }
            }
          }

        }
      }

      //END phiV_i loop

      //BEGIN phiP_i loop: mass balance
      for(unsigned i = 0; i < nDofsP; i++) {
        unsigned irow = nDofsT + dim * nDofsV + i;

        for(int k = 0; k < dim; k++) {
          Res[irow] += (gradSolV_gss[k][k]) * phiP[i]  * weight;

          if(assembleMatrix) {
            unsigned irowMat = nDofsTVP * irow;

            for(unsigned j = 0; j < nDofsV; j++) {
              unsigned jcol = (nDofsT + k * nDofsV + j);
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

    if(assembleMatrix) {
      KK->add_matrix_blocked(Jac, sysDof, sysDof);
    }

    //END local to global Matrix/Vector assembly

  }

  //END element loop

  RES->close();

  if(assembleMatrix) {
    KK->close();
  }

  // ***************** END ASSEMBLY *******************
}

void AssembleBoussinesqAppoximation(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("NST");   // pointer to the linear implicit system named "Poisson"
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

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //old solution variable
  unsigned solT0Index;
  solT0Index = mlSol->GetIndex("T0");    // get the position of "T" in the ml_sol object
  vector < unsigned > solV0Index(dim);
  solV0Index[0] = mlSol->GetIndex("U0");    // get the position of "U" in the ml_sol object
  solV0Index[1] = mlSol->GetIndex("V0");    // get the position of "V" in the ml_sol object
  if(dim == 3) solV0Index[2] = mlSol->GetIndex("W0");       // get the position of "V" in the ml_sol object
  
  //solution variable
  unsigned solTIndex;
  solTIndex = mlSol->GetIndex("T");    // get the position of "T" in the ml_sol object
  unsigned solTType = mlSol->GetSolutionType(solTIndex);    // get the finite element type for "T"

  vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object
  if(dim == 3) solVIndex[2] = mlSol->GetIndex("W");       // get the position of "V" in the ml_sol object
  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"

  unsigned solTPdeIndex;
  solTPdeIndex = mlPdeSys->GetSolPdeIndex("T");    // get the position of "T" in the pdeSys object

  vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object
  if(dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W");

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object

  vector < double >  solT0; // local solution
  vector < vector < double > >  solV0(dim);    // local solution
  
  vector < double >  solT; // local solution
  vector < vector < double > >  solV(dim);    // local solution
  
  vector < double >  solP; // local solution

  vector < vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  solT0.reserve(maxSize);
  solT.reserve(maxSize);
  
  for(unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(maxSize);
    solV0[k].reserve(maxSize);
    coordX[k].reserve(maxSize);
  }

  solP.reserve(maxSize);

  //RHS vectors
  vector < double >  fT; // local solution
  fT.reserve(maxSize);
  vector < vector < double > >  fV(dim);    // local solution
  for(unsigned  k = 0; k < dim; k++) {
    fV[k].reserve(maxSize);
  }
  vector < double >  fP; // local solution
  fP.reserve(maxSize);
  
  vector <double> phiV;  // local test function
  vector <double> phiV_x; // local test function first order partial derivatives
  vector <double> phiV_xx; // local test function second order partial derivatives

  phiV.reserve(maxSize);
  phiV_x.reserve(maxSize * dim);
  phiV_xx.reserve(maxSize * dim2);

  vector <double> phiT;  // local test function
  vector <double> phiT_x; // local test function first order partial derivatives
  vector <double> phiT_xx; // local test function second order partial derivatives

  phiT.reserve(maxSize);
  phiT_x.reserve(maxSize * dim);
  phiT_xx.reserve(maxSize * dim2);

  double* phiP;
  double weight; // gauss point weight

  vector< int > sysDof; // local to global pdeSys dofs
  sysDof.reserve((dim + 2) *maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve((dim + 2) *maxSize);

  vector < double > Jac;
  Jac.reserve((dim + 2) *maxSize * (dim + 2) *maxSize);

  if(counter == 10){ 
    KK->print_matlab("matrix.txt", "ascii");
    Mat KKp = (static_cast< PetscMatrix* >(KK))->mat();  
    PetscViewer    viewer;
    PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,300,300,&viewer);
    MatView(KKp,viewer);	
  }   
  
  if(assembleMatrix)
    KK->zero(); // Set to zero all the entries of the Global Matrix  
  
  sol->_Sol[solVIndex[0]]->zero();  
  sol->_Sol[solVIndex[1]]->zero();  
  sol->_Sol[solPIndex]->zero();  
  sol->_Sol[solTIndex]->zero();  
  
  unsigned nprocs = msh->n_processors();  
  unsigned sizeU = msh->_dofOffset[solVType][nprocs];
  unsigned sizeUV = sizeU + msh->_dofOffset[solVType][nprocs];
  unsigned sizeP = msh->_dofOffset[solPType][nprocs];
  unsigned sizeT = msh->_dofOffset[solTType][nprocs];
  
  if(preconditioner == FS_TVp || preconditioner == ASM_TVp || preconditioner == ILU_TVp){
    if (counter < sizeT){
       sol->_Sol[solTIndex]->set(counter, 1.);
       sol->_Sol[solTIndex]->close();
    }
    else if ( counter < sizeT + sizeU ){
       sol->_Sol[solVIndex[0]]->set(counter - sizeT, 1.);
       sol->_Sol[solVIndex[0]]->close();
    }
    else if ( counter < sizeT + sizeUV){
      sol->_Sol[solVIndex[1]]->set(counter - sizeT - sizeU, 1.);
      sol->_Sol[solVIndex[1]]->close();
    }
    else{
      sol->_Sol[solPIndex]->set(counter - sizeT - sizeUV , 1.);
      sol->_Sol[solPIndex]->close();
    }
  }
  
  if(preconditioner == ASM_VTp ){
    if ( counter <  sizeU ){
      sol->_Sol[solVIndex[0]]->set(counter, 1.);
      sol->_Sol[solVIndex[0]]->close();
    }
    else if(counter < sizeUV ){
      sol->_Sol[solVIndex[1]]->set(counter - sizeU, 1.);
      sol->_Sol[solVIndex[1]]->close();
    }
    else if (counter < sizeUV + sizeT){
      sol->_Sol[solTIndex]->set(counter - sizeUV, 1.);
      sol->_Sol[solTIndex]->close();
    }
    else{
      sol->_Sol[solPIndex]->set(counter - sizeUV - sizeT, 1.);
      sol->_Sol[solPIndex]->close();
    }
  }
  
  if (preconditioner == FS_VTp || preconditioner == ILU_VTp){
     if ( counter <  sizeU ){
      sol->_Sol[solVIndex[0]]->set(counter, 1.);
      sol->_Sol[solVIndex[0]]->close();
     }
     else if(counter < sizeUV ){
      sol->_Sol[solVIndex[1]]->set(counter - sizeU, 1.);
      sol->_Sol[solVIndex[1]]->close();
     }
     else if (counter < sizeUV + sizeP ){
      sol->_Sol[solPIndex]->set(counter - sizeUV, 1.);
      sol->_Sol[solPIndex]->close();
    }
    else{
      sol->_Sol[solTIndex]->set(counter- sizeUV - sizeP, 1.);
      sol->_Sol[solTIndex]->close();
    }
  }
   
//   /*if ( counter <  sizeU ){
//     sol->_Sol[solVIndex[0]]->set(counter, 1.);
//     sol->_Sol[solVIndex[0]]->close();
//   }	
//   else if ( counter < sizeUV ){
//     sol->_Sol[solVIndex[1]]->set(counter - sizeU, 1.);
//     sol->_Sol[solVIndex[1]]->close();
//   }
//   else if ( counter < sizeUVP ){
//     sol->_Sol[solPIndex]->set(counter - sizeUV, 1.);
//     sol->_Sol[solPIndex]->close();
//   }
//   else {
//     sol->_Sol[solTIndex]->set(counter- sizeUVP, 1.);
//     sol->_Sol[solTIndex]->close();
//   }*/

  //BEGIN element loop
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    //BEGIN local dof number extraction
    unsigned nDofsT = msh->GetElementDofNumber(iel, solTType);  //temperature
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);  //velocity
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);  //pressure
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType); // coordinates
    unsigned nDofsTVP = nDofsT + dim * nDofsV + nDofsP; // all solutions
    //END local dof number extraction

    //BEGIN memory allocation
    Res.resize(nDofsTVP);
    std::fill(Res.begin(), Res.end(), 0);
    Jac.resize(nDofsTVP * nDofsTVP);
    std::fill(Jac.begin(), Jac.end(), 0);
    sysDof.resize(nDofsTVP);

    solT.resize(nDofsV);
    solT0.resize(nDofsV);

    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      solV0[k].resize(nDofsV);
      coordX[k].resize(nDofsX);
    }

    solP.resize(nDofsP);
    //END memory allocation
    fT.assign(nDofsT,0.);
    fP.assign(nDofsP,0.);
    for(unsigned  k = 0; k < dim; k++) {
      fV[k].assign(nDofsV,0.);
    }
    
    //BEGIN global to local extraction
    for(unsigned i = 0; i < nDofsT; i++) { //temperature
      unsigned solTDof = msh->GetSolutionDof(i, iel, solTType);  //local to global solution dof
      solT[i] = (*sol->_Sol[solTIndex])(solTDof);  //global to local solution value
      solT0[i] = (*sol->_Sol[solT0Index])(solTDof);  //global to local solution value
      sysDof[i] = pdeSys->GetSystemDof(solTIndex, solTPdeIndex, i, iel);  //local to global system dof
      if(sysDof[i]==counter) {
	fT[i]=0.;
	std::cout<<counter<<" "<<"T"<<std::endl;
      }
    }

    for(unsigned i = 0; i < nDofsV; i++) { //velocity
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);  //local to global solution dof

      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);  //global to local solution value
	solV0[k][i] = (*sol->_Sol[solV0Index[k]])(solVDof);  //global to local solution value
        sysDof[i + nDofsT + k * nDofsV] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);  //local to global system dof
	if(sysDof[i + nDofsT + k * nDofsV]==counter) {
	  fV[k][i]=0.;
	  std::cout<<counter<<" V"<<k<<std::endl;
	}
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) { //pressure
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);  //local to global solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);  //global to local solution value
      sysDof[i + nDofsT + dim * nDofsV] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);  //local to global system dof
      if(sysDof[i + nDofsT + dim * nDofsV]==counter) {
	fP[i]=0.;
	std::cout<<counter<<" P"<<std::endl;
      }
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
      msh->_finiteElement[ielGeom][solTType]->Jacobian(coordX, ig, weight, phiT, phiT_x, phiT_xx);
      msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV, phiV_x, phiV_xx);
      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      double solT_gss = 0;
      vector < double > gradSolT_gss(dim, 0.);
      
      double solT0_gss = 0;
      vector < double > gradSolT0_gss(dim, 0.);

      for(unsigned i = 0; i < nDofsT; i++) {
        solT_gss += phiT[i] * solT[i];
	solT0_gss += phiT[i] * solT0[i];

        for(unsigned j = 0; j < dim; j++) {
          gradSolT_gss[j] += phiT_x[i * dim + j] * solT[i];
	  gradSolT0_gss[j] += phiT_x[i * dim + j] * solT0[i];
        }
      }

      vector < double > solV_gss(dim, 0);
      vector < vector < double > > gradSolV_gss(dim);
      
      vector < double > solV0_gss(dim, 0);
      vector < vector < double > > gradSolV0_gss(dim);

      for(unsigned  k = 0; k < dim; k++) {
        gradSolV_gss[k].resize(dim);
	gradSolV0_gss[k].resize(dim);
        std::fill(gradSolV_gss[k].begin(), gradSolV_gss[k].end(), 0);
	std::fill(gradSolV0_gss[k].begin(), gradSolV0_gss[k].end(), 0);
      }

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solV_gss[k] += phiV[i] * solV[k][i];
	  solV0_gss[k] += phiV[i] * solV0[k][i];
        }

        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            gradSolV_gss[k][j] += phiV_x[i * dim + j] * solV[k][i];
	    gradSolV0_gss[k][j] += phiV_x[i * dim + j] * solV0[k][i];
          }
        }
      }

      double solP_gss = 0;

      for(unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }


      double alpha = 1.;
      double beta = 1.;//40000.;


      double Pr = Prandtl;
      double Ra = Rayleigh;

      //BEGIN phiT_i loop: Energy balance
      for(unsigned i = 0; i < nDofsT; i++) {
        unsigned irow = i;
        
	Res[irow] +=  phiT[i] * fT[i] * weight;
	
        for(unsigned k = 0; k < dim; k++) {
          Res[irow] +=  -alpha / sqrt(Ra * Pr) * phiT_x[i * dim + k] * gradSolT_gss[k] * weight;
          Res[irow] +=  -  phiT[i] *(solV0_gss[k] * gradSolT_gss[k] + solV_gss[k] * gradSolT0_gss[k]) * weight; // why is gradsolT0_gass[k]

          if(assembleMatrix) {
            unsigned irowMat = irow * nDofsTVP;

            for(unsigned j = 0; j < nDofsT; j++) {
              Jac[ irowMat + j ] +=  alpha / sqrt(Ra * Pr) * phiT_x[i * dim + k] * phiT_x[j * dim + k] * weight;
              Jac[ irowMat + j ] +=  phiT[i] * solV0_gss[k] * phiT_x[j * dim + k] * weight; 
            }

            for(unsigned j = 0; j < nDofsV; j++) {
              unsigned jcol = nDofsT + k * nDofsV + j;
              Jac[ irowMat + jcol ] +=  phiT[i] * phiV[j] * gradSolT0_gss[k] * weight;
            }
          }

        }
      }

      //END phiT_i loop

      //BEGIN phiV_i loop: Momentum balance
      for(unsigned i = 0; i < nDofsV; i++) {

        for(unsigned k = 0; k < dim; k++) {
          unsigned irow = nDofsT + k * nDofsV + i;

	  Res[irow] +=  phiV[i] * fV[k][i] * weight;
	  
          for(unsigned l = 0; l < dim; l++) {
            Res[irow] +=  -sqrt(Pr / Ra) * phiV_x[i * dim + l] * (gradSolV_gss[k][l] + gradSolV_gss[l][k]) * weight;
            Res[irow] +=  - phiV[i] * ( solV0_gss[l] * gradSolV_gss[k][l] + solV_gss[l] * gradSolV0_gss[k][l] ) * weight;
          }

          Res[irow] += solP_gss * phiV_x[i * dim + k] * weight;

          if(k == 1) {
            Res[irow] += beta * solT_gss * phiV[i] * weight;
          }

          if(assembleMatrix) {
            unsigned irowMat = nDofsTVP * irow;

            for(unsigned l = 0; l < dim; l++) {
              for(unsigned j = 0; j < nDofsV; j++) {
                unsigned jcol1 = (nDofsT + k * nDofsV + j);
                unsigned jcol2 = (nDofsT + l * nDofsV + j);
                Jac[ irowMat + jcol1] += sqrt(Pr / Ra) * phiV_x[i * dim + l] * phiV_x[j * dim + l] * weight;
                Jac[ irowMat + jcol2] += sqrt(Pr / Ra) * phiV_x[i * dim + l] * phiV_x[j * dim + k] * weight;
//		Jac[ irowMat + jcol1] += sqrt(Pr / Ra) * phiV_x[i * dim + k] * phiV_x[j * dim + k] * weight;
//		Jac[ irowMat + jcol2] += sqrt(Pr / Ra) * phiV_x[i * dim + l] * phiV_x[j * dim + l] * weight;
                Jac[ irowMat + jcol1] +=  phiV[i] * solV0_gss[l] * phiV_x[j * dim + l] * weight;
//		Jac[ irowMat + jcol1] +=  phiV[i] * solV0_gss[l] * phiV_x[j * dim + k] * weight;
                Jac[ irowMat + jcol2] +=  phiV[i] * phiV[j] * gradSolV0_gss[k][l] * weight;
              }
            }

            for(unsigned j = 0; j < nDofsP; j++) {
              unsigned jcol = (nDofsT + dim * nDofsV) + j;
              Jac[ irowMat + jcol] += - phiV_x[i * dim + k] * phiP[j] * weight;
            }

            if(k == 1) {
              for(unsigned j = 0; j < nDofsT; j++) {
                Jac[ irowMat + j ] += -beta * phiV[i] * phiT[j] * weight;
              }
            }
          }

        }
      }

      //END phiV_i loop

      //BEGIN phiP_i loop: mass balance
      for(unsigned i = 0; i < nDofsP; i++) {
        unsigned irow = nDofsT + dim * nDofsV + i;
	Res[irow] += phiP[i] * fP[i]  * weight;
        for(int k = 0; k < dim; k++) {
          Res[irow] += -(gradSolV_gss[k][k]) * phiP[i]  * weight;

          if(assembleMatrix) {
            unsigned irowMat = nDofsTVP * irow;

            for(unsigned j = 0; j < nDofsV; j++) {
              unsigned jcol = (nDofsT + k * nDofsV + j);
              Jac[ irowMat + jcol ] += phiP[i] * phiV_x[j * dim + k] * weight;
            }
          }

        }
      }

      //END phiP_i loop

    }

    //END Gauss point loop

    //BEGIN local to global Matrix/Vector assembly
    RES->add_vector_blocked(Res, sysDof);

    if(assembleMatrix) {
      KK->add_matrix_blocked(Jac, sysDof, sysDof);
    }

    //END local to global Matrix/Vector assembly

  }

  //END element loop

  RES->close();

  if(assembleMatrix) {
    KK->close();
  }

  counter++;
  
  // ***************** END ASSEMBLY *******************
}
