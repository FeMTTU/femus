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
#include "Marker.hpp"
#include "MyVector.hpp"

//double Miu = 0.01;   int c0=-1; int cn=-1; //Re=100;
//double Miu = 0.002;  int c0=-1; int cn=-1; //Re=500;
double Miu = 0.001;  int c0=2; int cn=6; //Re=1000;

//unsigned numberOfUniformLevels = 7; unsigned numberOfSelectiveLevels = 0; //uniform
//unsigned numberOfUniformLevels = 8; unsigned numberOfSelectiveLevels = 0; //uniform
//unsigned numberOfUniformLevels = 9; unsigned numberOfSelectiveLevels = 0; //uniform
//<<<<<<< HEAD
//unsigned numberOfUniformLevels = 4; unsigned numberOfSelectiveLevels = 3; //non-uniform
//unsigned numberOfUniformLevels = 4; unsigned numberOfSelectiveLevels = 4; //non-uniform


//unsigned numberOfUniformLevels = 7; unsigned numberOfSelectiveLevels = 0; //non-uniform
// unsigned numberOfUniformLevels = 4; unsigned numberOfSelectiveLevels = 4; //non-uniform

//unsigned numberOfUniformLevels = 4; unsigned numberOfSelectiveLevels = 5; //non-uniform


int counter = 0 ;

using namespace femus;

double InitalValueU(const std::vector < double >& x)
{

  double value = 0.;
   if (x[1] > 0.492188 &&  x[0] > -0.5 + 1.0e-8 && x[0] < 0.5 - 1.0e-8) value = 1.*(x[1] - 0.492188) / (0.5 - 0.492188);


  return value;
}

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time)
{
  bool dirichlet = true; //dirichlet
  value = 0.;

  if (!strcmp(SolName, "U")) {
    if (facename == 4) {
      if (x[0] > -0.5 + 1.0e-8 && x[0] < 0.5 - 1.0e-8) value = 1.;
    }
  }
  else if (!strcmp(SolName, "P")) {
    dirichlet = false;
  }

  return dirichlet;
}

bool SetRefinementFlag(const std::vector < double >& x, const int& elemgroupnumber, const int& level)
{
  bool refine = false;
  unsigned level0 = 0;
//   double a = static_cast<double>(rand())/RAND_MAX;
//   if ( a < 0.25) refine	= true;

//   if( fabs(x[0] - 0.5) < 0.5/ pow(2,level) && fabs(x[1] - 0.5) < 0.5/ pow(2,level) ){
//     refine = true;
//   }

  double pi = acos(-1.);
  double radius = pi / 32.0 * (level - level0 - 2.0);
  double radius2 = radius * radius;

  if ( (x[0]*x[0] + x[1] * x[1]) > radius2) {
   std::cout << level << std::endl;
    refine	= true;
  }
  return refine;
}

void PrintConvergenceInfo(char *stdOutfile, char* outfile, const unsigned &numofrefinements);
void AssembleBoussinesqAppoximation(MultiLevelProblem& ml_prob);

int main(int argc, char** args)
{

  unsigned precType = 0;

  if (argc >= 2) {
    Miu = strtod(args[1], NULL);
    std::cout << Miu << std::endl;
  }

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

  unsigned numberOfUniformLevels = 8;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels, SetRefinementFlag);
 // mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(3);

  // print mesh info
  mlMsh.PrintInfo();
  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, SECOND);
  mlSol.AddSolution("V", LAGRANGE, SECOND);
  if (dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND);
  mlSol.AddSolution("P",  DISCONTINUOUS_POLYNOMIAL, FIRST);
  mlSol.AssociatePropertyToSolution("P", "Pressure");
  mlSol.Initialize("All");
//   mlSol.Initialize("U", InitalValueU);

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
  if (dim == 3) system.AddSolutionToSystemPDE("W");
  system.AddSolutionToSystemPDE("P");
  //system.SetLinearEquationSolverType(FEMuS_DEFAULT);
  system.SetLinearEquationSolverType(FEMuS_ASM);
  
  // attach the assembling function to system
  system.SetAssembleFunction(AssembleBoussinesqAppoximation);

  system.SetMaxNumberOfNonLinearIterations(20);
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
  system.SetPreconditionerFineGrids(MLU_PRECOND);

  system.SetTolerances(1.e-5, 1.e-8, 1.e+50, 30, 30); //GMRES tolerances

  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");
  system.SetNumberOfSchurVariables(1);
  system.SetElementBlockNumber("All");

  system.MGsolve();
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  mlMsh.PrintInfo();
  
  char *stdOutfile =  new char[100];
  char *outfile =  new char[100];
  sprintf(stdOutfile, "trueResidualMiu=%s.txt", args[1]);
  sprintf(outfile, "convergenceMiu=%s.txt", args[1]);
  std::cout << stdOutfile << std::endl;

  PrintConvergenceInfo(stdOutfile, outfile, numberOfUniformLevels);
  
  
  return 0;
}

void AssembleBoussinesqAppoximation(MultiLevelProblem& ml_prob)
{
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use

  double Mu;
  
  if(counter < c0 ) Mu = 2. * Miu;
  else if ( counter <= cn ) {
    Mu = 2*Miu*(cn-counter)/(cn-c0) + Miu*(counter-c0)/(cn-c0);
  }
  else{
    Mu = Miu;
  }
  std::cout << counter << " " << Mu <<std::endl;
  counter++;

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

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //solution variable
  vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object
  if (dim == 3) solVIndex[2] = mlSol->GetIndex("W");      // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"

  vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object
  if (dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W");

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object

  vector < vector < double > >  solV(dim);    // local solution
  vector < double >  solP; // local solution

  vector < vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned  k = 0; k < dim; k++) {
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

  double* phiP;
  double weight; // gauss point weight

  vector< int > sysDof; // local to global pdeSys dofs
  sysDof.reserve((dim + 2) *maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve((dim + 2) *maxSize);

  vector < double > Jac;
  Jac.reserve((dim + 2) *maxSize * (dim + 2) *maxSize);

  if (assembleMatrix)
    KK->zero(); // Set to zero all the entries of the Global Matrix

  //BEGIN element loop
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    //BEGIN local dof number extraction
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);  //velocity
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);  //pressure
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType); // coordinates
    unsigned nDofsVP = dim * nDofsV + nDofsP; // all solutions
    //END local dof number extraction

    //BEGIN memory allocation
    Res.resize(nDofsVP);
    std::fill(Res.begin(), Res.end(), 0);
    Jac.resize(nDofsVP * nDofsVP);
    std::fill(Jac.begin(), Jac.end(), 0);
    sysDof.resize(nDofsVP);

    for (unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsX);
    }
    
    std::vector <double> M(nDofsP);

    solP.resize(nDofsP);
    //END memory allocation

    //BEGIN global to local extraction
    for (unsigned i = 0; i < nDofsV; i++) { //velocity
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);  //local to global solution dof
      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);  //global to local solution value
        sysDof[i + k * nDofsV] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);  //local to global system dof
      }
    }

    for (unsigned i = 0; i < nDofsP; i++) { //pressure
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);  //local to global solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);  //global to local solution value
      sysDof[i + dim * nDofsV] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);  //local to global system dof
    }

    for (unsigned i = 0; i < nDofsX; i++) { //coordinates
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);  //local to global coordinate dof

      for (unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);  //global to local coordinate value
      }
    }

    //END global to local extraction

    //BEGIN Gauss point loop
    short unsigned ielGeom = msh->GetElementType(iel);

    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV, phiV_x, phiV_xx);
      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point

      vector < double > solV_gss(dim, 0);
      vector < vector < double > > gradSolV_gss(dim);

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

      double solP_gss = 0;

      for (unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }



//       //BEGIN phiT_i loop: Energy balance
//       for(unsigned i = 0; i < nDofsT; i++) {
//         unsigned irow = i;
//
//         for(unsigned k = 0; k < dim; k++) {
//           Res[irow] +=  -alpha / sqrt(Ra * Pr) * phiT_x[i * dim + k] * gradSolT_gss[k] * weight;
//           Res[irow] +=  -phiT[i] * solV_gss[k] * gradSolT_gss[k] * weight;
//
//           if(assembleMatrix) {
//             unsigned irowMat = irow * nDofsTVP;
//
//             for(unsigned j = 0; j < nDofsT; j++) {
//               Jac[ irowMat + j ] +=  alpha / sqrt(Ra * Pr) * phiT_x[i * dim + k] * phiT_x[j * dim + k] * weight;
//               Jac[ irowMat + j ] +=  phiT[i] * solV_gss[k] * phiT_x[j * dim + k] * weight;
//             }
//
//             for(unsigned j = 0; j < nDofsV; j++) {
//               unsigned jcol = nDofsT + k * nDofsV + j;
//               Jac[ irowMat + jcol ] += phiT[i] * phiV[j] * gradSolT_gss[k] * weight;
//             }
//           }
//
//         }
//       }
//
//       //END phiT_i loop

      //BEGIN phiV_i loop: Momentum balance
      for (unsigned i = 0; i < nDofsV; i++) {
        for (unsigned k = 0; k < dim; k++) {
          unsigned irow =  k * nDofsV + i;
          for (unsigned l = 0; l < dim; l++) {
            Res[irow] +=  -Mu * phiV_x[i * dim + l] * (gradSolV_gss[k][l] + gradSolV_gss[l][k]) * weight;
            Res[irow] +=  -phiV[i] * solV_gss[l] * gradSolV_gss[k][l] * weight;
          }
          Res[irow] += solP_gss * phiV_x[i * dim + k] * weight;

          if (assembleMatrix) {
            unsigned irowMat = nDofsVP * irow;

            for (unsigned l = 0; l < dim; l++) {
              for (unsigned j = 0; j < nDofsV; j++) {
                unsigned jcol1 = (k * nDofsV + j);
                unsigned jcol2 = (l * nDofsV + j);
                Jac[ irowMat + jcol1] += Mu * phiV_x[i * dim + l] * phiV_x[j * dim + l] * weight;
                Jac[ irowMat + jcol2] += Mu * phiV_x[i * dim + l] * phiV_x[j * dim + k] * weight;
                Jac[ irowMat + jcol1] += phiV[i] * solV_gss[l] * phiV_x[j * dim + l] * weight;
                Jac[ irowMat + jcol2] += phiV[i] * phiV[j] * gradSolV_gss[k][l] * weight;
              }
            }

            for (unsigned j = 0; j < nDofsP; j++) {
              unsigned jcol = (dim * nDofsV) + j;
              Jac[ irowMat + jcol] += - phiV_x[i * dim + k] * phiP[j] * weight;
            }
          }
        }
      }

      //END phiV_i loop

      //BEGIN phiP_i loop: mass balance
      for (unsigned i = 0; i < nDofsP; i++) {
        unsigned irow = dim * nDofsV + i;

        for (int k = 0; k < dim; k++) {
          Res[irow] += +(gradSolV_gss[k][k]) * phiP[i]  * weight;

          if (assembleMatrix) {
            unsigned irowMat = nDofsVP * irow;

            for (unsigned j = 0; j < nDofsV; j++) {
              unsigned jcol = ( k * nDofsV + j);
              Jac[ irowMat + jcol ] -= phiP[i] * phiV_x[j * dim + k] * weight;
            }
            
            M[i] += phiP[i] * phiP[i] * weight;
            
//             for (unsigned j = 0; j < nDofsP; j++) {
//               M[i] += phiP[i] * phiP[j] * weight;
//             }
            
          }

        }
      }

      //END phiP_i loop

    }
    //END Gauss point loop

    
    std::vector<std::vector<double> > BtMinv (dim * nDofsV);
    for(unsigned i=0; i< dim * nDofsV; i++){
      BtMinv[i].resize(nDofsP);
    }
    
    std::vector<std::vector<double> > B(nDofsP);
    for(unsigned i=0; i< nDofsP; i++){
      B[i].resize(dim * nDofsV);
    }
    
    for(unsigned i = 0; i < dim * nDofsV; i++){
      unsigned irow =  i * nDofsVP;
      for(unsigned j = 0; j < nDofsP; j++){
	unsigned jcol = (dim * nDofsV) + j;
	BtMinv[i][j] = .1 * Jac[ irow +  jcol] / M[j];
      }
    }
    
    for(unsigned i = 0; i < nDofsP; i++){
      unsigned irow = ( (dim * nDofsV) + i) * nDofsVP;
      for(unsigned j = 0; j < dim * nDofsV; j++){
	B[i][j] = Jac[irow + j];
      }
    }
    
    std::vector<std::vector<double> > Jg (dim * nDofsV);
    for(unsigned i=0; i< dim * nDofsV; i++){
      Jg[i].resize(dim * nDofsV);
    }
        
    for(unsigned i = 0; i < dim * nDofsV; i++){
      for(unsigned j = 0; j < dim * nDofsV; j++){
	Jg[i][j] = 0.;
	for(unsigned k = 0; k < nDofsP; k++){
	  Jg[i][j] += BtMinv[i][k] * B[k][j];
	}
      }
    }
   
    std::vector<double> fg (dim * nDofsV);
    for(unsigned i = 0; i < dim * nDofsV; i++){
      fg[i] = 0;
      for(unsigned j = 0; j < nDofsP; j++){
	fg[i] += BtMinv[i][j] * Res[dim * nDofsV + j];
      }
    }
    
    
     for(unsigned i = 0; i < dim * nDofsV; i++){
	unsigned irow =  i * nDofsVP;
	for(unsigned j = 0; j < dim * nDofsV; j++){
	  Jac[irow + j] += Jg[i][j];
	  //std::cout<< Jg[i][j]<<" ";
	  //std::cout<< Jac[irow + j]<<" ";
	}
	//std::cout<<"\n";
	Res[i] += fg[i];
     }
    //std::cout<<"\n";
    
    
    
    
    
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

void PrintConvergenceInfo(char *stdOutfile, char* outfile, const unsigned &numofrefinements) {

  std::cout << "END_COMPUTATION\n" << std::flush;

  std::ifstream inf;
  inf.open(stdOutfile);
  if(!inf) {
    std::cout << "Redirected standard output file not found\n";
    std::cout << "add option -std_output std_out_filename > std_out_filename\n";
    return;
  }

  unsigned counter1 = 0;
  std :: vector <unsigned> Level(numofrefinements, 0);
  std :: vector <unsigned> Num_Nonlinear(numofrefinements, 0);
  std :: vector <unsigned> Num_GMRES(numofrefinements, 0);
  std :: vector <unsigned> Ave_GMRES(numofrefinements, 0);

  std::ofstream outf;
  outf.open(outfile, std::ofstream::app);
  outf << std::endl << std::endl;
  outf << "Number_of_refinements=" << numofrefinements << std::endl;
  outf << "Nonlinear_Iteration,resid_norm0,resid_normN,N,convergence";

  std::string str1;
  inf >> str1;
  while(str1.compare("END_COMPUTATION") != 0) {

    if(str1.compare("Start") == 0) {
      inf >> str1;
      if(str1.compare("Level") == 0) {
        inf >> str1;
        if(str1.compare("Max") == 0) {
          int value;
          inf >> value;
          Level[counter1] = value;
        }
      }
    }

    if(str1.compare("Nonlinear") == 0) {
      inf >> str1;
      if(str1.compare("iteration") == 0) {
        inf >> str1;
        outf << std::endl << str1;
        Num_Nonlinear[counter1] += 1;
      }
    }
    else if(str1.compare("KSP") == 0) {
      inf >> str1;
      if(str1.compare("preconditioned") == 0) {
        inf >> str1;
        if(str1.compare("resid") == 0) {
          inf >> str1;
          if(str1.compare("norm") == 0) {
            double norm0 = 1.;
            double normN = 1.;
            unsigned counter = 0;
            inf >> norm0;
            outf << "," << norm0;
            for(unsigned i = 0; i < 11; i++) {
              inf >> str1;
            }
            while(str1.compare("norm") == 0) {
              inf >> normN;
              counter++;
              for(unsigned i = 0; i < 11; i++) {
                inf >> str1;
              }
              Num_GMRES[counter1] += 1;
            }
            outf << "," << normN;
            if(counter != 0) {
              outf << "," << counter << "," << pow(normN / norm0, 1. / counter);
            }
            else {
              outf << "Invalid solver, set -outer_ksp_solver \"gmres\"";
            }
          }
        }
      }
    }

    if(str1.compare("End") == 0) {
      inf >> str1;
      if(str1.compare("Level") == 0) {
        inf >> str1;
        if(str1.compare("Max") == 0) {
          counter1 ++ ;
        }
      }
    }
    inf >> str1;
  }


  outf << "\n" << "Level, Number of nonlinear,Number of GMRES, Average number of GMRES per nonlinear" << std::endl;
  for(unsigned i = 0; i < numofrefinements; i++) {
    outf << Level[i] << "," << Num_Nonlinear[i] << ","
         << Num_GMRES[i] << "," << double(Num_GMRES[i]) / double(Num_Nonlinear[i]) << std::endl;
  }

  outf.close();
  inf.close();
}
