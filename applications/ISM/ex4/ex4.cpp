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
#include "XDMFWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"
#include "Marker.hpp"
#include "Line.hpp"


using namespace femus;

unsigned dim;

bool SetBoundaryCondition(const std::vector < double >& x, const char name[],
                          double& value, const int FaceName, const double time)
{
  bool test = 1; //Dirichlet
  value = 0.;
  //   cout << "Time bdc : " <<  time << endl;
  if (!strcmp(name, "U")) {
    if (1 == FaceName) {  //inflow
      double um = 1;
      double r2 = x[1] * x[1] + x[2] * x[2];
      value = (1. - r2) * um;
    } else if (2 == FaceName) { //outflow
      test = 0;
      value = 0.;
    } else if (3 == FaceName) { // no-slip fluid wall
      test = 1;
      value = 0.;
    }
  } else if (!strcmp(name, "V")) {
    if (1 == FaceName) {        //inflow
      test = 1;
      value = 0.;
    } else if (2 == FaceName) { //outflow
      test = 0;
      value = 0.;
    } else if (3 == FaceName) { // no-slip fluid wall
      test = 1;
      value = 0.;
    }
  } else if (!strcmp(name, "W")) {
    if (1 == FaceName) {       //inflow
      test = 1;
      value = 0.;
    } else if (2 == FaceName) { //outflow
      test = 0;
      value = 0.;
    } else if (3 == FaceName) { // no-slip fluid wall
      test = 1;
      value = 0.;
    }
  } else if (!strcmp(name, "P")) {
    if (1 == FaceName) {
      test = 0;
      value = 0.;
    } else if (2 == FaceName) {
      test = 0;
      value = 0.;
    } else if (3 == FaceName) {
      test = 0;
      value = 0.;
    }
  }
  return test;
}

//------------------------------------------------------------------------------------------------------------
unsigned numberOfUniformLevels;

bool SetRefinementFlag(const std::vector < double >& x, const int& elemgroupnumber, const int& level)
{

  bool refine = 0;

  //------------------ 3D --------------------------//
  //if (elemgroupnumber == 6 && level < 2) refine = 1;

  if (elemgroupnumber == 7 && level < numberOfUniformLevels) refine = 1;

  if (elemgroupnumber == 8 && level < numberOfUniformLevels + 1) refine = 1;

  //------------------ 2D --------------------------//
  //if (elemgroupnumber == 7 && level < 2) refine = 1;

  //if (elemgroupnumber == 6 && level < 3) refine = 1;

  //if (elemgroupnumber == 5 && level < 4) refine = 1;

//   if (elemgroupnumber==6 && level<1) refine=1;
//   if (elemgroupnumber==7 && level<2) refine=1;
//   if (elemgroupnumber==8 && level<3) refine=1;

  return refine;

}

//------------------------------------------------------------------------------------------------------------


void AssembleIncompressibleNavierStokes(MultiLevelProblem& mlProb);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );


int main(int argc, char** args)
{

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  //mlMsh.ReadCoarseMesh("./input/cylinder2Dnew.neu", "seventh", scalingFactor);
  mlMsh.ReadCoarseMesh("./input/pipe.neu", "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  dim = mlMsh.GetDimension();

  numberOfUniformLevels = 3;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, SECOND);
  mlSol.AddSolution("V", LAGRANGE, SECOND);
  if (dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND);

  //mlSol.AddSolution("P", LAGRANGE, FIRST);
  mlSol.AddSolution("P",  DISCONTINUOUS_POLYNOMIAL, FIRST);

  mlSol.AssociatePropertyToSolution("P", "Pressure", false);
  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("NS");

  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  if (dim == 3) system.AddSolutionToSystemPDE("W");

  system.AddSolutionToSystemPDE("P");

  //system.SetLinearEquationSolverType(FEMuS_DEFAULT);
  system.SetLinearEquationSolverType(FEMuS_ASM);   // Additive Swartz Method
  // attach the assembling function to system
  system.SetAssembleFunction(AssembleIncompressibleNavierStokes);

  system.SetMaxNumberOfNonLinearIterations(10);
  system.SetNonLinearConvergenceTolerance(1.e-8);
  system.SetMaxNumberOfResidualUpdatesForNonlinearIteration(10);
  system.SetResidualUpdateConvergenceTolerance(1.e-15);

  system.SetMgType(F_CYCLE);

  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);
  // initialize and solve the system
  system.init();

  system.SetSolverFineGrids(GMRES);
  system.SetTolerances(1.e-20, 1.e-20, 1.e+50, 50, 10);

//   system.SetSolverFineGrids(RICHARDSON);
//   system.SetRichardsonScaleFactor(.5);
//   system.SetTolerances( 1.e-20, 1.e-20, 1.e+50, 50, 10 );

  system.SetPreconditionerFineGrids(ILU_PRECOND);

  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");
  system.SetNumberOfSchurVariables(1);
  system.SetElementBlockNumber(2);
  //system.UseSamePreconditioner();
  system.SetOuterSolver(PREONLY);
  system.MGsolve();

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  // print mesh info
  mlMsh.PrintInfo();


  //unsigned pSize = 7;
  unsigned theta_intervals = 10;
  unsigned radius_intervals = 9;
  unsigned size = radius_intervals * theta_intervals;

  std::vector < std::vector < double > > x; // marker
  std::vector < MarkerType > markerType;

  x.resize(size);
  markerType.resize(size);

  std::vector < std::vector < std::vector < double > > > streamline(size);

  for (unsigned j = 0; j < size; j++) {
    x[j].assign(dim, 0.);
    markerType[j] = VOLUME;
  }


//   for(unsigned j = 0; j < pSize; j++) {
//     std::vector < double > x(3);
//     x[0] = 0.;
//     x[1] = 0.;
//     x[2] = -0.75 + 0.25 * j;
//     particle[j] = new Marker(x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels-1), 2, true);
//   }


  clock_t start_time = clock();
  clock_t init_time = clock();
  //BEGIN INITIALIZE PARTICLES

  double pi = acos(-1.);
  unsigned counter = 0;
  for (unsigned k = 1; k < radius_intervals + 1 ; k++) {
    for (unsigned j = 0; j < theta_intervals; j++) {
      x[counter][0] = 0.;
      x[counter][1] = 0.1 * k * sin(2.*pi / theta_intervals * j);
      x[counter][2] = 0.1 * k * cos(2.*pi / theta_intervals * j);
      counter++;
    }
  }

  //Line linea0(x, markerType, mlMsh.GetLevel(numberOfUniformLevels - 1), 2);

  std::vector< Line* > linea(1);

  linea[0] =  new Line(x, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), 2);

  linea[0]->GetStreamLine(streamline, 0);
  linea[0]->GetStreamLine(streamline, 1);
  PrintLine(DEFAULT_OUTPUTDIR, streamline, true, 0);

  //END INITIALIZE PARTICLES

  double T = 240;
  unsigned n  = 160;

  std::cout << std::endl << " init in  " << std::setw(11) << std::setprecision(6) << std::fixed
            << static_cast<double>((clock() - init_time)) / CLOCKS_PER_SEC << " s" << std::endl;

  clock_t advection_time;
  for (unsigned k = 0; k < n; k++) {
    std::cout<< "Iteration = "<< k << std::endl;
    if(k == n/2) advection_time = clock();
    for(int i = linea.size() - 1; i>=0; i--){
      linea[i]->AdvectionParallel(40, T / n, 4);
      linea[i]->GetStreamLine(streamline, linea.size() - i );     
    }
    PrintLine(DEFAULT_OUTPUTDIR, streamline, true, k + 1);
    linea.resize(k+2);
    linea[k+1] =  new Line(x, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), 2);
    
  }

  std::cout << std::endl << " advection in: " << std::setw(11) << std::setprecision(6) << std::fixed
            << static_cast<double>((clock() - advection_time)) / CLOCKS_PER_SEC << " s" << std::endl;

  std::cout << std::endl << " RANNA in: " << std::setw(11) << std::setprecision(6) << std::fixed
            << static_cast<double>((clock() - start_time)) / CLOCKS_PER_SEC << " s" << std::endl;
  
  for (unsigned i = 0; i < linea.size(); i++) {
    delete linea[i];
  }

  return 0;
}


void AssembleIncompressibleNavierStokes(MultiLevelProblem& mlProb)
{
  //  mlProb is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem



  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &mlProb.get_system<NonLinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*           msh         = mlProb._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*           el          = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*   mlSol   = mlProb._ml_sol;  // pointer to the multilevel solution object
  Solution*   sol         = mlProb._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level];  // pointer to the equation (level) object


  bool assembleMatrix = mlPdeSys->GetAssembleMatrix();
  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;
  if (assembleMatrix) s.continue_recording();
  else s.pause_recording();


  SparseMatrix*   KK          = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*  RES         = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

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

  vector < vector < adept::adouble > >  solV(dim);    // local solution
  vector < adept::adouble >  solP; // local solution

  vector< vector < adept::adouble > > aResV(dim);    // local redidual vector
  vector< adept::adouble > aResP; // local redidual vector

  vector < vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(maxSize);
    aResV[k].reserve(maxSize);
    coordX[k].reserve(maxSize);
  }

  solP.reserve(maxSize);
  aResP.reserve(maxSize);


  vector <double> phiV;  // local test function
  vector <double> gradPhiV; // local test function first order partial derivatives
  vector <double> nablaPhiV; // local test function second order partial derivatives

  phiV.reserve(maxSize);
  gradPhiV.reserve(maxSize * dim);
  nablaPhiV.reserve(maxSize * dim2);

  double* phiP;
  double weight; // gauss point weight

  vector< int > sysDof; // local to global pdeSys dofs
  sysDof.reserve((dim + 1) *maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve((dim + 1) *maxSize);

  vector < double > Jac;
  Jac.reserve((dim + 1) *maxSize * (dim + 1) *maxSize);

  if (assembleMatrix) KK->zero();

  //BEGIN element loop for each process
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielType = msh->GetElementType(iel);

    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs

    unsigned nDofsTVP = dim * nDofsV + nDofsP;
    // resize local arrays
    sysDof.resize(nDofsTVP);

    for (unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsX);
    }
    solP.resize(nDofsP);

    for (unsigned  k = 0; k < dim; k++) {
      aResV[k].assign(nDofsV, 0.);   //resize
    }
    aResP.assign(nDofsP, 0.);   //resize

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

    //BEGIN a new recording of all the operations involving adept::adouble variables
    if (assembleMatrix) s.new_recording();

    //BEGIN Gauss point loop
    for (unsigned ig = 0; ig < msh->_finiteElement[ielType][solVType]->GetGaussPointNumber(); ig++) {
      msh->_finiteElement[ielType][solVType]->Jacobian(coordX, ig, weight, phiV, gradPhiV, nablaPhiV);
      phiP = msh->_finiteElement[ielType][solPType]->GetPhi(ig);

      // evaluate the solution, the solution derivatives and the coordinates in the Gauss point

      vector < adept::adouble > solVig(dim, 0);
      vector < vector < adept::adouble > > gradSolVig(dim);

      for (unsigned  k = 0; k < dim; k++) {
        gradSolVig[k].assign(dim, 0);
      }

      for (unsigned i = 0; i < nDofsV; i++) {
        for (unsigned  k = 0; k < dim; k++) {
          solVig[k] += phiV[i] * solV[k][i];
        }

        for (unsigned j = 0; j < dim; j++) {
          for (unsigned  k = 0; k < dim; k++) {
            gradSolVig[k][j] += gradPhiV[i * dim + j] * solV[k][i];
          }
        }
      }

      adept::adouble solPig = 0;

      for (unsigned i = 0; i < nDofsP; i++) {
        solPig += phiP[i] * solP[i];
      }

      //BEGIN phiV loop (momentum)
      double nu = 0.1;
      for (unsigned i = 0; i < nDofsV; i++) {
        vector < adept::adouble > NSV(dim, 0.);

        for (unsigned j = 0; j < dim; j++) {
          for (unsigned  k = 0; k < dim; k++) {
            NSV[k]   +=  nu * gradPhiV[i * dim + j] * (gradSolVig[k][j] + gradSolVig[j][k]);
            NSV[k]   +=  phiV[i] * (solVig[j] * gradSolVig[k][j]);
          }
        }

        for (unsigned  k = 0; k < dim; k++) {
          NSV[k] += -solPig * gradPhiV[i * dim + k];
        }

        for (unsigned  k = 0; k < dim; k++) {
          aResV[k][i] += - NSV[k] * weight;
        }
      }
      //END phiV loop

      //BEGIN phiP loop (continuity)
      for (unsigned i = 0; i < nDofsP; i++) {
        for (int k = 0; k < dim; k++) {
          aResP[i] += - (gradSolVig[k][k]) * phiP[i]  * weight;
        }
      }
      //END phiP loop

    }
    //END Gauss point loop

    //BEGIN Extract and store the residual
    Res.resize(nDofsTVP);    //resize
    for (int i = 0; i < nDofsV; i++) {
      for (unsigned  k = 0; k < dim; k++) {
        Res[ i + k * nDofsV ] = -aResV[k][i].value();
      }
    }
    for (int i = 0; i < nDofsP; i++) {
      Res[ i + dim * nDofsV ] = -aResP[i].value();
    }
    RES->add_vector_blocked(Res, sysDof);
    //END Extract and store the residual

    //BEGIN Extract and store the jacobian
    if (assembleMatrix) {
      Jac.resize(nDofsTVP * nDofsTVP);
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
    //END extract and store the jacobian
  }
  //END element loop for each process

  RES->close();
  if (assembleMatrix) KK->close();

  // ***************** END ASSEMBLY *******************
}
