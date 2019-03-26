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
#include "LinearImplicitSystem.hpp"
#include "adept.h"
#include "FieldSplitTree.hpp"
#include <stdlib.h>
#include "PetscMatrix.hpp"


using namespace femus;
unsigned counter = 0;
const double pi = 2.0 * acos(0.0); 

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;
  return dirichlet;
}

bool SetRefinementFlag(const std::vector < double >& x, const int& elemgroupnumber, const int& level) {

  bool refine = false;
  unsigned level0 = 0;

//   if (elemgroupnumber == 6 && level < 3) refine = 1;


  double a = static_cast<double>(rand())/RAND_MAX;
  if ( a < 0.5) refine	= true;
  return refine;

//  std::cout<<level<<std::endl;
// double radius = pi / 8.0 /(level - level0);

//  double radius = sqrt(2.0)/2.0/pow(2.0,level - level0);
  
//   unsigned powindex;
//   powindex = level -level0;
//   if (powindex % 2 == 0) powindex = powindex - 1;
//   double radius = sqrt(2.0)/2.0/pow(2.0,powindex);
  
 // double radius2 = radius * radius;
  
//   if ( (x[0]*x[0] + x[1] * x[1]) < radius2){
//     refine	= true;
//   }	 
//   return refine;  
  
  
//   if( fabs(x[0]) < 0.5/ pow(2,level) && fabs(x[1]) < 0.5/ pow(2,level) ){
//     refine = true;
//   }
//   return refine;
}


void AssembleBoussinesqAppoximation(MultiLevelProblem& ml_prob);

int main(int argc, char** args) {
  
  srand(1);
  
  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
 // mlMsh.ReadCoarseMesh("./input/square_quad.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/square_tri.neu","seventh",scalingFactor);
  mlMsh.ReadCoarseMesh("./input/cube_all_shapes.neu","seventh",scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 3;
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);
//   unsigned numberOfSelectiveLevels = 0;
//  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels, NULL);
//   erase all the coarse mesh levels
  //mlMsh.EraseCoarseLevels(1);
  //numberOfUniformLevels -= 1;		
  // print mesh info
  
  
  mlMsh.PrintInfo();
  MultiLevelSolution mlSol(&mlMsh);
  mlSol.AddSolution("U", LAGRANGE, SECOND);
  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);
  // add system Poisson in mlProb as a Linear Implicit System  
  LinearImplicitSystem& system = mlProb.add_system < LinearImplicitSystem > ("Poisson");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("U");

  //system.SetLinearEquationSolverType(FEMuS_DEFAULT);
  system.SetLinearEquationSolverType(FEMuS_ASM);
  // attach the assembling function to system
  system.SetAssembleFunction(AssembleBoussinesqAppoximation);
  
  //system.SetMaxNumberOfNonLinearIterations(1);
  //system.SetNonLinearConvergenceTolerance(1.e-8);
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
  system.SetTolerances(1.e-20, 1.e-20, 1.e+50, 1, 1); //GMRES tolerances
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");
  
  system.SetNumberOfSchurVariables(1);
  system.SetElementBlockNumber(2);
 
  //////////////////////////////////////////////////////////////////////
  //solution variable
  unsigned solUIndex;
  solUIndex = mlSol.GetIndex("U");	// get the position of "U" in the ml_sol object
  unsigned solUType = mlSol.GetSolutionType(solUIndex);		// get the finite element type for "U"
  Mesh* msh = mlMsh.GetLevel(numberOfUniformLevels+numberOfSelectiveLevels-1);
  
  unsigned nprocs = msh->n_processors();
  unsigned iproc = msh->processor_id();
  unsigned sizeU = msh->_dofOffset[solUType][nprocs];;  
  Solution* sol = mlSol.GetLevel(numberOfUniformLevels+numberOfSelectiveLevels-1);
  
  for(unsigned i = 0; i < sizeU; i++){
    system.SetOuterSolver(PREONLY);  
    system.MGsolve();
    mlSol.GenerateBdc("All");


    //BEGIN print solution vector
    std::vector<double> solLocal;
    sol->_Sol[solUIndex]->localize_to_one (solLocal, 0);
      
    if(iproc == 0){
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
      fout.precision(20);
    
      for(unsigned j = 0; j < sizeU; j++ ){
	fout << solLocal[j]<< " ";
	//fout << (*sol->_Sol[solUIndex])(j)<< " ";
      }  
      fout<<std::endl;
      fout.close();
    }
    //BEGIN print solution vector

    
    // print solutions
    std::vector < std::string > variablesToBePrinted;
    variablesToBePrinted.push_back("All");
  
    
    VTKWriter vtkIO(&mlSol);
    vtkIO.SetDebugOutput( true );
    vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, counter-1);
  }
  /////////////////////////////////////ultiLevelProb/////////////////////////
  
//   // print solutions
//   std::vector < std::string > variablesToBePrinted;
//   variablesToBePrinted.push_back("All");
// 
//   VTKWriter vtkIO(&mlSol);
//   vtkIO.SetDebugOutput( true );
//   vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);
  mlMsh.PrintInfo();
  
  return 0;
}

void AssembleBoussinesqAppoximation(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use
  LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<LinearImplicitSystem> ("Poisson");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*           msh         = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*           el          = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*   mlSol         = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*   sol         = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level];  // pointer to the equation (level) object
  bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  // call the adept stack object
  SparseMatrix*   KK          = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  
  SparseMatrix*   KKamr          = pdeSys->_KKamr;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*  RES         = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  
  //solution variable
  unsigned solUIndex;
  solUIndex = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  unsigned solUType = mlSol->GetSolutionType(solUIndex);    // get the finite element type for "U"
  
  unsigned solUPdeIndex;
  solUPdeIndex = mlPdeSys->GetSolPdeIndex("U");    // get the position of "T" in the pdeSys object
  
  vector < double >  solU; // local solution
  vector < vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
  
  solU.reserve(maxSize);
  for(unsigned  k = 0; k < dim; k++) {
    coordX[k].reserve(maxSize);
  }

  //RHS vectors
  vector < double >  fU; // local solution
  fU.reserve(maxSize);
  
  vector <double> phiU;  // local test function
  vector <double> phiU_x; // local test function first order partial derivatives
  vector <double> phiU_xx; // local test function second order partial derivatives

  phiU.reserve(maxSize);
  phiU_x.reserve(maxSize * dim);
  phiU_xx.reserve(maxSize * dim2);
  
  double weight; // gauss point weight

  vector< int > sysDof; // local to global pdeSys dofs
  sysDof.reserve((dim + 2) *maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve((dim + 2) *maxSize);

  vector < double > Jac;
  Jac.reserve((dim + 2) *maxSize * (dim + 2) *maxSize);

  if(assembleMatrix) KK->zero(); // Set to zero all the entries of the Global Matrix    
  sol->_Sol[solUIndex]->zero();  
  
  unsigned nprocs = msh->n_processors();  
  unsigned sizeU = msh->_dofOffset[solUType][nprocs];
//   if (counter < sizeU){
// //     if( (*sol->_Bdc[solUIndex])(counter) > 0.5){
//       sol->_Sol[solUIndex]->set(counter, 1.);
//       sol->_Sol[solUIndex]->close();
// //     }
//   }
  
//   (sol->_Eps[solUIndex])->matrix_mult( *sol->_Sol[solUIndex], *mlPdeSys->_PPamr[level]);
//   *(sol->_Sol[solUIndex]) = *(sol->_Eps[solUIndex]);
  
  
  //BEGIN element loop
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    //BEGIN local dof number extraction
    unsigned nDofsU = msh->GetElementDofNumber(iel, solUType);  //velocity
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType); // coordinates
    //END local dof number extraction

    //BEGIN memory allocation
    Res.resize(nDofsU);
    std::fill(Res.begin(), Res.end(), 0);
    Jac.resize(nDofsU * nDofsU);
    std::fill(Jac.begin(), Jac.end(), 0);
    sysDof.resize(nDofsU);
    solU.resize(nDofsU);

    //END memory allocation
    fU.assign(nDofsU,0.);
    
    //BEGIN global to local extraction
    for(unsigned i = 0; i < nDofsU; i++) { //velocity
      unsigned solUDof = msh->GetSolutionDof(i, iel, solUType);  //local to global solution dof
      solU[i] = (*sol->_Sol[solUIndex])(solUDof);  //global to local solution value
      sysDof[i] = pdeSys->GetSystemDof(solUIndex, solUPdeIndex, i, iel);  //local to global system dof
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
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solUType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solUType]->Jacobian(coordX, ig, weight, phiU, phiU_x, phiU_xx);
      
      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      double solU_gss = 0;
      vector < double > gradSolU_gss(dim, 0.);

      for(unsigned i = 0; i < nDofsU; i++) {
        solU_gss += phiU[i] * solU[i];
        for(unsigned j = 0; j < dim; j++) {
          gradSolU_gss[j] += phiU_x[i * dim + j] * solU[i];
        }
      }

      //BEGIN phiU_i loop
      for(unsigned i = 0; i < nDofsU; i++) {
        unsigned irow = i;       
	Res[irow] +=  phiU[i] * fU[i] * weight;	
        for(unsigned k = 0; k < dim; k++) {
          Res[irow] +=  - phiU_x[i * dim + k] * gradSolU_gss[k] * weight;
          if(assembleMatrix) {
            unsigned irowMat = irow * nDofsU;
            for(unsigned j = 0; j < nDofsU; j++) {
              Jac[ irowMat + j ] += phiU_x[i * dim + k] * phiU_x[j * dim + k] * weight;
            }
          }
        }
      }
      //END phiU_i loop
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
  RES->add(counter,1./nprocs);
  RES->close();

  if(assembleMatrix) {
    KK->close();
  }

  if(counter == 10){ 
    //KKamr->print_matlab("matrix.txt", "ascii");
    KK->print_matlab("matrix.txt", "ascii");
 //   Mat KKp = (static_cast< PetscMatrix* >(KK))->mat();  
//     PetscViewer    viewer;
//     PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,300,300,&viewer);
//     MatView(KKp,viewer);
  }   
  
  
  counter++;
  
  // ***************** END ASSEMBLY *******************
}
