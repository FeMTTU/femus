/** tutorial/Ex1
 * This example shows how to:
 * initialize a femus application;
 * define the multilevel-mesh object mlMsh;
 * read from the file ./input/square.neu the coarse-level mesh and associate it to mlMsh;
 * add in mlMsh uniform refined level-meshes;
 * define the multilevel-solution object mlSol associated to mlMsh;
 * add in mlSol different types of finite element solution variables;
 * initialize the solution varables;
 * define vtk and gmv writer objects associated to mlSol;
 * print vtk and gmv binary-format files in ./output directory.
 **/

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "PetscMatrix.hpp"
#include "slepceps.h"

#include "LinearImplicitSystem.hpp"
#include "Marker.hpp"
#include "Line.hpp"

#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"

using namespace femus;

Line* line1;
Line* line2;
Line* lineI;

unsigned DIM = 2;

void AssembleNitscheProblem_AD(MultiLevelProblem& mlProb);

void BuildFlag(MultiLevelSolution& mlSol);
void GetInterfaceElementEigenvalues(MultiLevelSolution& mlSol);

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = false; //dirichlet
  if(facename == 1)  dirichlet = true;
  value = 0.;
  return dirichlet;
}


int main(int argc, char** args) {

  if(DIM != 2 && DIM != 3) {
    std::cout << "Wrong Dimension!" << std::endl;
    return 0;
  }

  // init Petsc-MPI communicator

  SlepcInitialize(&argc, &args, PETSC_NULL, PETSC_NULL);
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;

  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;

  unsigned nx = 11; // this should always be a odd number
  unsigned ny = 11; // this should always be a odd number
  unsigned nz = 1;

  double length = 1.;
  double lengthx = 0.1;

  if(DIM == 2) {
    mlMsh.GenerateCoarseBoxMesh(nx, ny, 0, -lengthx / 2, lengthx / 2, -length / 2, length / 2, 0., 0., QUAD9, "seventh");
  }
  else if(DIM == 3) {
    nz = ny;
    mlMsh.GenerateCoarseBoxMesh(nx, ny, nz, -lengthx / 2, lengthx / 2, -length / 2, length / 2, -length / 2, length / 2, HEX27, "seventh");
  }

  double Lref = 1.;
  double Uref = 1.;
  double rhos1 = 7850;
  double rhos2 = 7850;
  double E1 = 2.e06;
  double E2 = 2.e07;
  double nu1 = 0.4;
  double nu2 = 0.4;

  Parameter par(Lref, Uref);

  // Generate Solid Object
  Solid solid1(par, E1, nu1, rhos1, "Neo-Hookean");
  Solid solid2(par, E2, nu2, rhos2, "Neo-Hookean");


  //mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);

  FEOrder femOrder = SECOND;

  mlSol.AddSolution("DX1", LAGRANGE, femOrder);
  mlSol.AddSolution("DY1", LAGRANGE, femOrder);
  if(DIM == 3) mlSol.AddSolution("DZ1", LAGRANGE, femOrder);

  mlSol.AddSolution("DX2", LAGRANGE, femOrder);
  mlSol.AddSolution("DY2", LAGRANGE, femOrder);
  if(DIM == 3) mlSol.AddSolution("DZ2", LAGRANGE, femOrder);

  mlSol.AddSolution("eflag", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("nflag", LAGRANGE, femOrder, 0, false);

  mlSol.AddSolution("CM1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("CM2", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

  mlSol.AddSolution("CL1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("CL2", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);


  mlSol.Initialize("All");

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");

  BuildFlag(mlSol);

  MultiLevelProblem ml_prob(&mlSol);

  ml_prob.parameters.set<Solid> ("Solid1") = solid1;
  ml_prob.parameters.set<Solid> ("Solid2") = solid2;


  // ******* Add FEM system to the MultiLevel problem *******
  LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ("Nitsche");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("DX1");
  system.AddSolutionToSystemPDE("DY1");
  if(DIM == 3) system.AddSolutionToSystemPDE("DZ1");
  system.AddSolutionToSystemPDE("DX2");
  system.AddSolutionToSystemPDE("DY2");
  if(DIM == 3) system.AddSolutionToSystemPDE("DZ2");


  // attach the assembling function to system
  system.SetAssembleFunction(AssembleNitscheProblem_AD);

  // time loop parameter
  system.SetMaxNumberOfLinearIterations(1);

  system.init();

  //init marker

  unsigned Ne = 4;

  double Lx = lengthx / nx;
  double Lx1 = 0.025 * Lx; //beam dimensions
  double Lx2 = Lx - Lx1; //beam dimensions

  double Ly = length;
  double Lz = (DIM == 3) ? length : 1.;

  unsigned rowy = ny * Ne;
  unsigned rowz = (DIM == 3) ?  nz * Ne : 1.;
  double dy = length / rowy; // size in y of each row
  double dz = Lz / rowz; // size in y of each row

  unsigned columns = Ne;
  double dx1 = Lx1 / Ne; // size in x of each column
  double dx2 = Lx2 / Ne; // size in x of each column
  unsigned size = rowy * rowz * columns;

  std::vector < std::vector < double > > x; // marker
  std::vector < MarkerType > markerType;

  x.resize(size);
  markerType.resize(size);

  for(unsigned j = 0; j < size; j++) {
    x[j].assign(DIM, 0.);
    markerType[j] = VOLUME;
  }



  //   double xc = 1.e-04 + 0.5 * H ; //should be this one to do COMSOL benchmark
  double x0 = -0.5 * Lx; //we are using this not to have markers on edges of elements from the beginning
  double y0 = -Ly / 2.;
  double z0 = -Lz / 2.;

  //BEGIN initialization
  for(unsigned i = 0; i < columns; i++) {
    for(unsigned j = 0; j < rowy; j++) {
      for(unsigned k = 0; k < rowz; k++) {
        x[ i * rowy * rowz + j * rowz + k][0] = x0 + 0.5 * dx1 + i * dx1;
        x[ i * rowy * rowz + j * rowz + k][1] = y0 + 0.5 * dy + j * dy;
        if(DIM == 3) x[i * rowy * rowz + j * rowz + k][2] = z0 + 0.5 * dz + k * dz;
      }
    }
  }

  //END

  unsigned solType = 2;
  std::vector < double > volume(x.size(), Lx1 * Ly * Lz / x.size());  // uniform marker volume
  line1 = new Line(x, volume, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), solType);

  std::vector < std::vector < std::vector < double > > >  line1Points(1);
  line1->GetLine(line1Points[0]);
  PrintLine(DEFAULT_OUTPUTDIR, "bulk1", line1Points, 0);
  //PrintLine("./output1", "bulk1", line1Points, 0);


  //BEGIN initialization bulk2 points
  x0 = -0.5 * Lx + Lx1; //we are using this not to have markers on edges of elements from the beginning
  for(unsigned i = 0; i < columns; i++) {
    for(unsigned j = 0; j < rowy; j++) {
      for(unsigned k = 0; k < rowz; k++) {
        x[ i * rowy * rowz + j * rowz + k][0] = x0 + 0.5 * dx2 + i * dx2;
      }
    }
  }
  volume.assign(x.size(), Lx2 * Ly * Lz / x.size());  // uniform marker volume
  line2 = new Line(x, volume, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), solType);
  //END initialization bulk2 points

  std::vector < std::vector < std::vector < double > > > line2Points(1);
  line2->GetLine(line2Points[0]);
  PrintLine(DEFAULT_OUTPUTDIR, "bulk2", line2Points, 0);
  //PrintLine("./output1", "bulk2", line2Points, 0);

  //interface marker initialization

  size = rowy * rowz;

  x.resize(size);
  markerType.resize(size);

  for(unsigned j = 0; j < size; j++) {
    x[j].assign(DIM, 0.);
    markerType[j] = INTERFACE;
  }

  std::vector < double > area(x.size(), Ly * Lz / x.size());  // uniform marker volume
  x0 = -0.5 * Lx + Lx1;; //we are using this not to have markers on edges of elements from the beginning
  std::vector < std::vector < std::vector < double > > > T;
  T.resize(x.size());
  for(unsigned i = 0; i < x.size(); i++) {
    T[i].resize(DIM - 1);
    for(unsigned k = 0; k < DIM - 1; k++) {
      T[i][k].resize(DIM, 0.);
    }
  }


  //BEGIN initialization
  for(unsigned i = 0; i < rowy; i++) {
    for(unsigned j = 0; j < rowz; j++) {
      x[ i * rowz + j][0] = x0;
      x[ i * rowz + j][1] = y0 + 0.5 * dy + i * dy;
      if(DIM == 3) x[ i * rowz + j][2] = z0 + 0.5 * dz + j * dz;

      T[i * rowz + j][0][0] = 0.;
      T[i * rowz + j][0][1] = -dy;
      if(DIM == 3) {
        T[i * rowz + j][0][2] = 0.;

        T[i * rowz + j][1][0] = 0.;
        T[i * rowz + j][1][1] = 0.;
        T[i * rowz + j][1][2] = dz;

      }
    }
  }
  //END

  lineI = new Line(x, T, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), solType);

  std::vector < std::vector < std::vector < double > > > lineIPoints(1);
  lineI->GetLine(lineIPoints[0]);
  PrintLine(DEFAULT_OUTPUTDIR, "interfaceLine", lineIPoints, 0);
  //PrintLine("./output1", "interfaceLine", lineIPoints, 0);
  //END interface markers


  // ******* Print solution *******
  mlSol.SetWriter(VTK);
  mlSol.GetWriter()->SetDebugOutput(true);

  std::vector<std::string> print_vars;
  print_vars.push_back("All");

  GetInterfaceElementEigenvalues(mlSol);

  system.MGsolve();
  
  
  mlSol.GetWriter()->Write("./output", "linear", print_vars, 0);

  std::vector<std::string> mov_vars1;
  mov_vars1.push_back("DX1");
  mov_vars1.push_back("DY1");
  if(DIM == 3) mov_vars1.push_back("DZ1");
  mlSol.GetWriter()->SetMovingMesh(mov_vars1);
  mlSol.GetWriter()->Write("./outputD1", "linear", print_vars, 0);


  std::vector<std::string> mov_vars2;
  mov_vars2.push_back("DX2");
  mov_vars2.push_back("DY2");
  if(DIM == 3) mov_vars2.push_back("DZ2");
  mlSol.GetWriter()->SetMovingMesh(mov_vars2);
  mlSol.GetWriter()->Write("./outputD2", "linear", print_vars, 0);


  ml_prob.clear();

  return 0;
}


void AssembleNitscheProblem_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object


  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("Nitsche");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble(); // We have different level of meshes. we assemble the problem on the specified one.

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  double rho1 = ml_prob.parameters.get<Solid> ("Solid1").get_density();
  double rho2 = ml_prob.parameters.get<Solid> ("Solid2").get_density();

  double mu1 = ml_prob.parameters.get<Solid> ("Solid1").get_lame_shear_modulus();
  double lambda1 = ml_prob.parameters.get<Solid> ("Solid1").get_lame_lambda();

  double mu2 = ml_prob.parameters.get<Solid> ("Solid2").get_lame_shear_modulus();
  double lambda2 = ml_prob.parameters.get<Solid> ("Solid2").get_lame_lambda();

  double g[DIM] = {1.};

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  std::vector< unsigned > solD1Index(dim);
  solD1Index[0] = mlSol->GetIndex("DX1");
  solD1Index[1] = mlSol->GetIndex("DY1");
  if(dim == 3)solD1Index[2] = mlSol->GetIndex("DZ1");

  std::vector< unsigned > solD2Index(dim);
  solD2Index[0] = mlSol->GetIndex("DX2");
  solD2Index[1] = mlSol->GetIndex("DY2");
  if(dim == 3)solD2Index[2] = mlSol->GetIndex("DZ2");

  unsigned solDType = mlSol->GetSolutionType(solD1Index[0]);

  std::vector< unsigned > solD1PdeIndex(dim);
  solD1PdeIndex[0] = mlPdeSys->GetSolPdeIndex("DX1");
  solD1PdeIndex[1] = mlPdeSys->GetSolPdeIndex("DY1");
  if(dim == 3) solD1PdeIndex[2] = mlPdeSys->GetSolPdeIndex("DZ1");

  std::vector< unsigned > solD2PdeIndex(dim);
  solD2PdeIndex[0] = mlPdeSys->GetSolPdeIndex("DX2");
  solD2PdeIndex[1] = mlPdeSys->GetSolPdeIndex("DY2");
  if(dim == 3) solD2PdeIndex[2] = mlPdeSys->GetSolPdeIndex("DZ2");

  std::vector < std::vector < adept::adouble > > solD1(dim); // local solution
  std::vector < std::vector < adept::adouble > > solD2(dim); // local solution

  unsigned CMIndex[2];
  unsigned CLIndex[2];

  CMIndex[0] = mlSol->GetIndex("CM1");
  CMIndex[1] = mlSol->GetIndex("CM2");
  
  CLIndex[0] = mlSol->GetIndex("CL1");
  CLIndex[1] = mlSol->GetIndex("CL2");

  unsigned eflagIndex = mlSol->GetIndex("eflag");
  unsigned nflagIndex = mlSol->GetIndex("nflag");

  vector < unsigned >  nodeFlag; // local solution

  vector < vector < double > > x(dim);    // local coordinates. x is now dim x m matrix.
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  std::vector< std::vector< adept::adouble > > aResD1(dim); // local redidual vector
  std::vector< std::vector< adept::adouble > > aResD2(dim); // local redidual vector

  vector< unsigned > l2GMap; // local to global mapping
  vector< double > Res; // local redidual vector
  vector < double > Jac;

  KK->zero(); // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual

  std::vector < std::vector < std::vector <double > > > aP(3);

  std::vector<Marker*> particle1 = line1->GetParticles();
  std::vector<unsigned> markerOffset1 = line1->GetMarkerOffset();
  unsigned imarker1 = markerOffset1[iproc];


  std::vector<Marker*> particle2 = line2->GetParticles();
  std::vector<unsigned> markerOffset2 = line2->GetMarkerOffset();
  unsigned imarker2 = markerOffset2[iproc];

  std::vector<Marker*> particleI = lineI->GetParticles();
  std::vector<unsigned> markerOffsetI = lineI->GetMarkerOffset();
  unsigned imarkerI = markerOffsetI[iproc];


  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofD  = msh->GetElementDofNumber(iel, solDType);  // number of solution element dofs

    unsigned eFlag = static_cast <unsigned>(floor((*sol->_Sol[eflagIndex])(iel) + 0.5));

    // resize local arrays
    l2GMap.resize(dim * 2 * nDofD);
    for(int k = 0; k < dim; k++) {
      solD1[k].resize(nDofD);
      solD2[k].resize(nDofD);
      aResD1[k].assign(nDofD, 0.);    //resize
      aResD2[k].assign(nDofD, 0.);    //resize
    }
    nodeFlag.resize(nDofD);

    for(int k = 0; k < dim; k++) {
      x[k].resize(nDofD);
    }



    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofD; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, solDType);
      nodeFlag[i] = (*sol->_Sol[nflagIndex])(iDof);
      for(unsigned k = 0; k < dim; k++) {
        solD1[k][i] = (*sol->_Sol[solD1Index[k]])(iDof);
        solD2[k][i] = (*sol->_Sol[solD2Index[k]])(iDof);

        l2GMap[k * nDofD + i] = pdeSys->GetSystemDof(solD1Index[k], solD1PdeIndex[k], i, iel);
        l2GMap[(dim + k) * nDofD + i] = pdeSys->GetSystemDof(solD2Index[k], solD2PdeIndex[k], i, iel);

      }
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofD; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, 2);    // global to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();
    if(eFlag == 0 || eFlag == 2) {
      // *** Element Gauss point loop ***
      for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solDType]->GetGaussPointNumber(); ig++) {
        // *** get gauss point weight, test function and test function partial derivatives ***
        msh->_finiteElement[ielGeom][solDType]->Jacobian(x, ig, weight, phi, phi_x);

        // evaluate the solution, the solution derivatives and the coordinates in the gauss point
        std::vector < std::vector < adept::adouble > > gradSolD1(dim);
        std::vector < std::vector < adept::adouble > > gradSolD2(dim);

        for(unsigned k = 0; k < dim; k++) {
          gradSolD1[k].assign(nDofD, 0.);
          gradSolD2[k].assign(nDofD, 0.);
        }

        for(unsigned i = 0; i < nDofD; i++) {
          for(unsigned k = 0; k < dim; k++) {
            for(unsigned j = 0; j < dim; j++) {
              gradSolD1[k][j] += phi_x[i * dim + j] * solD1[k][i];
              gradSolD2[k][j] += phi_x[i * dim + j] * solD2[k][i];
            }
          }
        }

        // *** phi_i loop ***
        for(unsigned i = 0; i < nDofD; i++) {

          adept::adouble divD1 = 0.;
          adept::adouble divD2 = 0.;
          for(unsigned k = 0; k < dim; k++) {
            divD1 += gradSolD1[k][k] * (eFlag == 0);
            divD2 += gradSolD2[k][k] * (eFlag == 2);
          }

          for(unsigned k = 0; k < dim; k++) {

            adept::adouble sigma1 = 0.;
            adept::adouble sigma2 = 0.;

            for(unsigned j = 0; j < dim; j++) {
              sigma1 += mu1 * (gradSolD1[k][j] + gradSolD1[j][k] * (eFlag == 0)) * phi_x[i * dim + j];
              sigma2 += mu2 * (gradSolD2[k][j] + gradSolD2[j][k] * (eFlag == 2)) * phi_x[i * dim + j];
            }
            sigma1 += lambda1 * divD1 * phi_x[i * dim + k];
            sigma2 += lambda2 * divD2 * phi_x[i * dim + k];

            if(eFlag == 0) {
              aResD1[k][i] += (- rho1 * g[k] * phi[i] + sigma1) * weight;
              if(nodeFlag[i] != 1) { // fake equation for sugular matrix
                aResD2[k][i] += sigma2 * weight;
              }
            }
            else if(eFlag == 2) {
              aResD2[k][i] += (- rho2 * g[k] * phi[i] + sigma2) * weight;
              if(nodeFlag[i] != 1) { // fake equation for sugular matrix
                aResD1[k][i] += sigma1 * weight;
              }
            }
          }

        } // end phi_i loop
      } // end gauss point loop
    }

    else {

      double iM1C1 = 1. / (mu1 * (*sol->_Sol[CMIndex[0]])(iel));
      double iM2C2 = 1. / (mu2 * (*sol->_Sol[CMIndex[1]])(iel));

      double denM = iM1C1 + iM2C2;

      double gammaM1 = iM1C1 / denM;
      double gammaM2 = iM2C2 / denM;

      double thetaM = 8. / denM;

      double iL1C1 = 1. / (lambda1 * (*sol->_Sol[CLIndex[0]])(iel));
      double iL2C2 = 1. / (lambda2 * (*sol->_Sol[CLIndex[1]])(iel));
      
      double denL = iL1C1 + iL2C2;
      
      double gammaL1 = iL1C1 / denL;
      double gammaL2 = iL2C2 / denL;
      
      double thetaL = 4. / denL;
      
      
//       double gammaL1 = 0.5;
//       double gammaL2 = 0.5;
// 
//       double gammaM1 = 0.5;
//       double gammaM2 = 0.5;
// 
// 
//       double thetaL = lambda1;
//       double thetaM = mu1;



      //bulk1
      while(imarker1 < markerOffset1[iproc + 1] && iel == particle1[imarker1]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particle1[imarker1]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][solDType]->Jacobian(x, xi, weight, phi, phi_x);
        double weight = particle1[imarker1]->GetMarkerMass();


        // evaluate the solution, the solution derivatives and the coordinates in the gauss point
        std::vector < std::vector < adept::adouble > > gradSolD1(dim);
        for(unsigned k = 0; k < dim; k++) {
          gradSolD1[k].assign(nDofD, 0.);
        }

        for(unsigned i = 0; i < nDofD; i++) {
          for(unsigned k = 0; k < dim; k++) {
            for(unsigned j = 0; j < dim; j++) {
              gradSolD1[k][j] += phi_x[i * dim + j] * solD1[k][i];
            }
          }
        }

        // *** phi_i loop ***
        for(unsigned i = 0; i < nDofD; i++) {

          adept::adouble divD1 = 0.;
          for(unsigned k = 0; k < dim; k++) {
            divD1 += gradSolD1[k][k];
          }

          for(unsigned k = 0; k < dim; k++) {
            adept::adouble sigma1 = 0.;
            for(unsigned j = 0; j < dim; j++) {
              sigma1 += mu1 * (gradSolD1[k][j] + gradSolD1[j][k]) * phi_x[i * dim + j];
            }
            sigma1 += lambda1 * divD1 * phi_x[i * dim + k];
            aResD1[k][i] += (- rho1 * g[k] * phi[i] + sigma1) * weight;
          }
        } // end phi_i loop
        imarker1++;
      }

      //bulk2
      while(imarker2 < markerOffset2[iproc + 1] && iel == particle2[imarker2]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particle2[imarker2]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][solDType]->Jacobian(x, xi, weight, phi, phi_x);
        double weight = particle2[imarker2]->GetMarkerMass();

        std::vector < std::vector < adept::adouble > > gradSolD2(dim);

        for(unsigned k = 0; k < dim; k++) {
          gradSolD2[k].assign(nDofD, 0.);
        }

        for(unsigned i = 0; i < nDofD; i++) {
          for(unsigned k = 0; k < dim; k++) {
            for(unsigned j = 0; j < dim; j++) {
              gradSolD2[k][j] += phi_x[i * dim + j] * solD2[k][i];
            }
          }
        }

        // *** phi_i loop ***
        for(unsigned i = 0; i < nDofD; i++) {
          adept::adouble divD2 = 0.;
          for(unsigned k = 0; k < dim; k++) {
            divD2 += gradSolD2[k][k];
          }
          for(unsigned k = 0; k < dim; k++) {
            adept::adouble sigma2 = 0.;
            for(unsigned j = 0; j < dim; j++) {
              sigma2 += mu2 * (gradSolD2[k][j] + gradSolD2[j][k]) * phi_x[i * dim + j];
            }
            sigma2 += lambda2 * divD2 * phi_x[i * dim + k];
            aResD2[k][i] += (- rho2 * g[k] * phi[i] + sigma2) * weight;
          }
        }
        imarker2++;
      }

      // interface
      while(imarkerI < markerOffsetI[iproc + 1] && iel == particleI[imarkerI]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particleI[imarkerI]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][solDType]->Jacobian(x, xi, weight, phi, phi_x);

        std::vector <std::vector < double > > T;
        particleI[imarkerI]->GetMarkerTangent(T);

        double weight;
        std::vector < double > N(dim);
        if(dim == 2) {
          N[0] =  T[0][1];
          N[1] = -T[0][0];
          weight = sqrt(N[0] * N[0] + N[1] * N[1]);
          N[0] /= weight;
          N[1] /= weight;
        }
        else {
          N[0] = T[0][1] * T[1][2] - T[0][2] * T[1][1];
          N[1] = T[0][2] * T[1][0] - T[0][0] * T[1][2];
          N[2] = T[0][0] * T[1][1] - T[0][1] * T[1][0];
          weight = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
          N[0] /= weight;
          N[1] /= weight;
          N[2] /= weight;
        }

        std::vector < adept::adouble > u1(dim, 0.);
        std::vector < adept::adouble > u2(dim, 0.);

        std::vector < adept::adouble > tau(dim, 0.);

        for(unsigned i = 0; i < nDofD; i++) {
          for(unsigned k = 0; k < dim; k++) {
            u1[k] += phi[i] * solD1[k][i];
            u2[k] += phi[i] * solD2[k][i];

            for(unsigned j = 0; j < dim; j++) {
              tau[k] += (gammaL1 * lambda1 * solD1[j][i] * phi_x[i * dim + j] +
                         gammaL2 * lambda2 * solD2[j][i] * phi_x[i * dim + j]) * N[k];

              tau[k] += (gammaM1 * mu1 * solD1[k][i] * phi_x[i * dim + j] +
                         gammaM2 * mu2 * solD2[k][i] * phi_x[i * dim + j]) * N[j];

              tau[k] += (gammaM1 * mu1 * solD1[j][i] * phi_x[i * dim + k] +
                         gammaM2 * mu2 * solD2[j][i] * phi_x[i * dim + k]) * N[j];

            }
          }
        }
        // *** phi_i loop ***
        for(unsigned i = 0; i < nDofD; i++) {
          for(unsigned k = 0; k < dim; k++) {
            aResD1[k][i] += tau[k] * phi[i] * weight;
            aResD1[k][i] += -thetaM * (u2[k] - u1[k]) * phi[i] * weight;

            aResD2[k][i] += -tau[k] * phi[i] * weight;
            aResD2[k][i] +=  thetaM * (u2[k] - u1[k]) * phi[i] * weight;

            for(unsigned j = 0; j < dim; j++) {
              aResD1[k][i] += -gammaL1 * (lambda1 * phi_x[i * dim + k] * N[j] * (u2[j] - u1[j])) * weight;
              aResD1[k][i] += -gammaM1 * (mu1 * phi_x[i * dim + j] * N[j] * (u2[k] - u1[k])) * weight;
              aResD1[k][i] += -gammaM1 * (mu1 * phi_x[i * dim + j] * N[k] * (u2[j] - u1[j])) * weight;

              aResD1[k][i] += -thetaL * (u2[j] - u1[j]) * N[j] * phi[i] * N[k] * weight;

              aResD2[k][i] += -gammaL2 * (lambda2 * phi_x[i * dim + k] * N[j] * (u2[j] - u1[j])) * weight;
              aResD2[k][i] += -gammaM2 * (mu2 * phi_x[i * dim + j] * N[j] * (u2[k] - u1[k])) * weight;
              aResD2[k][i] += -gammaM2 * (mu2 * phi_x[i * dim + j] * N[k] * (u2[j] - u1[j])) * weight;

              aResD2[k][i] +=  thetaL * (u2[j] - u1[j]) * N[j] * phi[i] * N[k] * weight;
            }


          }
        } // end phi_i loop
        imarkerI++;
      }

    }

    //copy the value of the adept::adoube aRes in double Res and store
    Res.resize(2 * dim * nDofD);    //resize

    for(unsigned k = 0; k < dim; k++) {
      for(int i = 0; i < nDofD; i++) {
        Res[k * nDofD + i] = - aResD1[k][i].value();
        Res[(k + dim) * nDofD + i] = - aResD2[k][i].value();
      }
    }

    RES->add_vector_blocked(Res, l2GMap);

    // define the dependent variables
    for(unsigned k = 0; k < dim; k++) {
      s.dependent(&aResD1[k][0], nDofD);
    }
    for(unsigned k = 0; k < dim; k++) {
      s.dependent(&aResD2[k][0], nDofD);
    }

    // define the independent variables
    for(unsigned k = 0; k < dim; k++) {
      s.independent(&solD1[k][0], nDofD);
    }
    for(unsigned k = 0; k < dim; k++) {
      s.independent(&solD2[k][0], nDofD);
    }

    // get the jacobian matrix (ordered by row major )
    Jac.resize(2 * dim * nDofD * 2 * dim * nDofD);   //resize
    s.jacobian(&Jac[0], true);

    //store K in the global matrix KK
    KK->add_matrix_blocked(Jac, l2GMap, l2GMap);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();

//   PetscViewer    viewer;
//   //PetscViewerDrawOpen (PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
//   //PetscObjectSetName ( (PetscObject) viewer, "FSI matrix");
//   //PetscViewerPushFormat (viewer, PETSC_VIEWER_DRAW_LG);
//   MatView ( (static_cast<PetscMatrix*> (KK))->mat(), viewer);
//   VecView ( (static_cast<PetscVector*> (RES))->vec(), viewer);

  //double a;
  //std::cin >> a;



  //***************** END ASSEMBLY *******************
}

void BuildFlag(MultiLevelSolution& mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol  = mlSol.GetSolutionLevel(level);
  Mesh     *msh   = mlSol._mlMesh->GetLevel(level);
  unsigned iproc  = msh->processor_id();

  unsigned eflagIndex = mlSol.GetIndex("eflag");
  unsigned nflagIndex = mlSol.GetIndex("nflag");

  unsigned nflagType = mlSol.GetSolutionType(nflagIndex);

  sol->_Sol[eflagIndex]->zero();
  sol->_Sol[nflagIndex]->zero();

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  std::vector < std::vector<double> >  x(dim);

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, nflagType);  // number of solution element dofs

    for(int k = 0; k < dim; k++) {
      x[k].resize(nDofu);  // Now we
    }

    for(unsigned i = 0; i < nDofu; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, 2);    // global to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

    bool interface = false;
    double signi = (x[0][0] < 0.) ? -1. : 1.;
    for(unsigned j = 1; j < nDofu; j++) {
      double signj = (x[0][j] < 0.) ? -1. : 1.;
      if(signi != signj) {
        interface = true;
        break;
      }
    }

    if(interface) {
      sol->_Sol[eflagIndex]->set(iel, 1.);
      for(unsigned i = 0; i < nDofu; i++) {
        unsigned iDof = msh->GetSolutionDof(i, iel, nflagType);
        sol->_Sol[nflagIndex]->set(iDof, 1.);
      }
    }
    else if(x[0][0] > 0. && x[0][1] > 0.) {
      sol->_Sol[eflagIndex]->set(iel, 2.);
    }
  }

  sol->_Sol[eflagIndex]->close();
  sol->_Sol[nflagIndex]->close();

}


void GetInterfaceElementEigenvalues(MultiLevelSolution& mlSol) {

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol  = mlSol.GetSolutionLevel(level);
  Mesh     *msh   = mlSol._mlMesh->GetLevel(level);
  unsigned iproc  = msh->processor_id();

  unsigned eflagIndex = mlSol.GetIndex("eflag");

  unsigned CMIndex[2];
  CMIndex[0] = mlSol.GetIndex("CM1");
  CMIndex[1] = mlSol.GetIndex("CM2");

  unsigned CLIndex[2];
  CLIndex[0] = mlSol.GetIndex("CL1");
  CLIndex[1] = mlSol.GetIndex("CL2");

  unsigned solIndex = mlSol.GetIndex("DX1");
  unsigned soluType = mlSol.GetSolutionType(solIndex);

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  std::vector < std::vector<double> >  x(dim);

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  std::vector < std::vector < std::vector <double > > > aP(3);

  std::vector<Marker*> particle1 = line1->GetParticles();
  std::vector<unsigned> markerOffset1 = line1->GetMarkerOffset();
  unsigned imarker1 = markerOffset1[iproc];

  std::vector<Marker*> particle2 = line2->GetParticles();
  std::vector<unsigned> markerOffset2 = line2->GetMarkerOffset();
  unsigned imarker2 = markerOffset2[iproc];

  std::vector<Marker*> particleI = lineI->GetParticles();
  std::vector<unsigned> markerOffsetI = lineI->GetMarkerOffset();
  unsigned imarkerI = markerOffsetI[iproc];

  sol->_Sol[CMIndex[0]]->zero();
  sol->_Sol[CMIndex[1]]->zero();

  sol->_Sol[CLIndex[0]]->zero();
  sol->_Sol[CLIndex[1]]->zero();

  std::vector < double > aM;
  std::vector < double > aM0;
  std::vector < std::vector < double > > bM(2);

  std::vector < double > aM1;
  std::vector < double > bM1;


  std::vector < double > aL;
  std::vector < double > aL0;
  std::vector < std::vector < double > > bL(2);

  std::vector < double > aL1;
  std::vector < double > bL1;


//   std::vector < double > a2;
//   std::vector < std::vector < double > > b2(2);

  Mat A, B;
  Vec v;
  EPS eps;

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned eFlag = static_cast <unsigned>(floor((*sol->_Sol[eflagIndex])(iel) + 0.5));
    if(eFlag == 1) {

      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);  // number of solution element dofs

      unsigned sizeAll = dim * nDofu;

      aM.assign(sizeAll * sizeAll, 0.);
      bM[0].assign(sizeAll * sizeAll, 0.);
      bM[1].assign(sizeAll * sizeAll, 0.);

      aL.assign(sizeAll * sizeAll, 0.);
      bL[0].assign(sizeAll * sizeAll, 0.);
      bL[1].assign(sizeAll * sizeAll, 0.);

      for(int k = 0; k < dim; k++) {
        x[k].resize(nDofu);  // Now we
      }

      for(unsigned i = 0; i < nDofu; i++) {
        unsigned xDof  = msh->GetSolutionDof(i, iel, 2);    // global to global mapping between coordinates node and coordinate dof
        for(unsigned k = 0; k < dim; k++) {
          x[k][i] = (*msh->_topology->_Sol[k])(xDof);      // global extraction and local storage for the element coordinates
        }
      }

      //bulk1
      while(imarker1 < markerOffset1[iproc + 1] && iel == particle1[imarker1]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particle1[imarker1]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][soluType]->Jacobian(x, xi, weight, phi, phi_x);
        double weight = particle1[imarker1]->GetMarkerMass();

        // *** phi_i loop ***

        for(unsigned k = 0; k < dim; k++) {
          for(unsigned i = 0; i < nDofu; i++) {
            for(unsigned l = 0; l < dim; l++) {
              for(unsigned j = 0; j < nDofu; j++) {
                bM[0][((nDofu * k) + i) * sizeAll + (k * nDofu + j)] += 0.5 * phi_x[i * dim + l] * phi_x[j * dim + l] * weight;
                bM[0][((nDofu * k) + i) * sizeAll + (l * nDofu + j)] += 0.5 * phi_x[i * dim + l] * phi_x[j * dim + k] * weight;

                bL[0][((nDofu * k) + i) * sizeAll + (l * nDofu + j)] += phi_x[i * dim + k] * phi_x[j * dim + l] * weight;

              }
            }
          }
        }
        imarker1++;
      }

      //bulk2
      while(imarker2 < markerOffset2[iproc + 1] && iel == particle2[imarker2]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particle2[imarker2]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][soluType]->Jacobian(x, xi, weight, phi, phi_x);
        double weight = particle2[imarker2]->GetMarkerMass();

        // *** phi_i loop ***

        for(unsigned k = 0; k < dim; k++) {
          for(unsigned i = 0; i < nDofu; i++) {
            for(unsigned l = 0; l < dim; l++) {
              for(unsigned j = 0; j < nDofu; j++) {
                bM[1][((nDofu * k) + i) * sizeAll + (k * nDofu + j)] += 0.5 * phi_x[i * dim + l] * phi_x[j * dim + l] * weight;
                bM[1][((nDofu * k) + i) * sizeAll + (l * nDofu + j)] += 0.5 * phi_x[i * dim + l] * phi_x[j * dim + k] * weight;

                bL[1][((nDofu * k) + i) * sizeAll + (l * nDofu + j)] += phi_x[i * dim + k] * phi_x[j * dim + l] * weight;
              }
            }
          }
        }
        imarker2++;
      }

      // interface
      while(imarkerI < markerOffsetI[iproc + 1] && iel == particleI[imarkerI]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particleI[imarkerI]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][soluType]->Jacobian(x, xi, weight, phi, phi_x);

        std::vector <std::vector < double > > T;
        particleI[imarkerI]->GetMarkerTangent(T);

        double weight;
        std::vector < double > N(dim);
        if(dim == 2) {
          N[0] =  T[0][1];
          N[1] = -T[0][0];
          weight = sqrt(N[0] * N[0] + N[1] * N[1]);
          N[0] /= weight;
          N[1] /= weight;
        }
        else {
          N[0] = T[0][1] * T[1][2] - T[0][2] * T[1][1];
          N[1] = T[0][2] * T[1][0] - T[0][0] * T[1][2];
          N[2] = T[0][0] * T[1][1] - T[0][1] * T[1][0];
          weight = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
          N[0] /= weight;
          N[1] /= weight;
          N[2] /= weight;
        }

        // *** phi_i loop ***

        for(unsigned k = 0; k < dim; k++) {
          for(int i = 0; i < nDofu; i++) {

            double gradPhiiDotN = 0.;
            for(unsigned l = 0; l < dim; l++) {
              gradPhiiDotN += phi_x[i * dim + l] * N[l];
            }
            for(int j = 0; j < nDofu; j++) {
              for(unsigned l = 0; l < dim; l++) {

                aM[((nDofu * k) + i) * sizeAll + (k * nDofu + j) ] += 0.5 * gradPhiiDotN * 0.5 * N[l] *  phi_x[j * dim + l]  * weight;
                aM[((nDofu * k) + i) * sizeAll + (l * nDofu + j) ] += 0.5 * gradPhiiDotN * 0.5 * N[l] *  phi_x[j * dim + k]  * weight;

                aL[((nDofu * k) + i) * sizeAll + (l * nDofu + j)] += phi_x[i * dim + k] * phi_x[j * dim + l] * weight;

              }
              for(unsigned l1 = 0; l1 < dim; l1++) {
                for(unsigned l2 = 0; l2 < dim; l2++) {
                  aM[((nDofu * k) + i) * sizeAll + (l1 * nDofu + j)] += 0.5 * N[k] * phi_x[i * dim + l1] * 0.5 * N[l2] *  phi_x[j * dim + l2]  * weight;
                  aM[((nDofu * k) + i) * sizeAll + (l2 * nDofu + j)] += 0.5 * N[k] * phi_x[i * dim + l1] * 0.5 * N[l2] *  phi_x[j * dim + l1]  * weight;
                }
              }

            }
          } // end phi_i loop
        }
        imarkerI++;
      }

      unsigned sizeAll0 = sizeAll;
      aM0 = aM;
      aL0 = aL;

      for(unsigned s = 0; s < 2; s++) {
        
        sizeAll = sizeAll0;

        //BEGIN DEFLATION

        unsigned sizeAll1 = dim * (nDofu - 1);
        aM1.resize(sizeAll1 * sizeAll1);
        bM1.resize(sizeAll1 * sizeAll1);

        MatCreateSeqDense(PETSC_COMM_SELF, sizeAll1, sizeAll1, NULL, &B);

        for(int k = 0; k < dim; k++) {
          for(int i = 0; i < nDofu - 1; i++) {

            int ip = i + 1;
            int i1 = (nDofu - 1) * k + i;
            for(int l = 0; l < dim; l++) {
              for(int j = 0; j < nDofu - 1; j++) {
                int jp = j + 1;
                int j1 = (nDofu - 1) * l + j;
                double value;
                value = aM0[((nDofu * k) + ip) * sizeAll0 + (nDofu * l + jp)] - aM0[(nDofu * k) * sizeAll0 + (nDofu * l + jp)];
                aM1[i1 * sizeAll1 + j1] = value;

                value = bM[s][((nDofu * k) + ip) * sizeAll0 + (nDofu * l + jp)] - bM[s][(nDofu * k) * sizeAll0 + (nDofu * l + jp)];
                bM1[i1 * sizeAll1 + j1] = value;
                MatSetValues(B, 1, &i1, 1, &j1, &value, INSERT_VALUES);

              }
            }
          }
        }

        MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);

        sizeAll = sizeAll1;
        aM.swap(aM1);
        bM[s].swap(bM1);

        double real = 0.;
        while(fabs(real) < 1.0e-12) {

          MatCreateVecs(B, &v, NULL);

          EPSCreate(PETSC_COMM_SELF, &eps);
          EPSSetOperators(eps, B, NULL);
          EPSSetFromOptions(eps);
          EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE);
          EPSSolve(eps);

          double imaginary;
          EPSGetEigenpair(eps, 0, &real, &imaginary, v, NULL);
          EPSDestroy(&eps);

          if(fabs(real) < 1.0e-12) {
            PetscScalar *pv;
            VecGetArray(v, &pv);
            unsigned ii = 0;
            for(unsigned i = 1; i < sizeAll; i++){
              if(fabs(pv[i]) > fabs(pv[ii])) ii = i;
            }

            unsigned sizeAll1 = sizeAll - 1;

            aM1.resize(sizeAll1 * sizeAll1);
            bM1.resize(sizeAll1 * sizeAll1);

            MatDestroy(&B);

            MatCreateSeqDense(PETSC_COMM_SELF, sizeAll1, sizeAll1, NULL, &B);

            for(unsigned i = 0; i < sizeAll; i++) {
              if(i != ii) {
                int i1 = (i < ii) ? i : i - 1;
                for(unsigned j = 0; j < sizeAll; j++) {
                  if(j != ii) {
                    int j1 = (j < ii) ? j : j - 1;
                    double value;
                    value = aM[i * sizeAll + j] - 1. / pv[ii] * pv[i] * aM[ii * sizeAll + j];
                    aM1[i1 * sizeAll1 + j1] = value;

                    value = bM[s][i * sizeAll + j] - 1. / pv[ii] * pv[i] * bM[s][ii * sizeAll + j];
                    bM1[i1 * sizeAll1 + j1] = value;
                    MatSetValues(B, 1, &i1, 1, &j1, &value, INSERT_VALUES);
                  }
                }
              }
            }
            MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
            VecRestoreArray(v, &pv);

            sizeAll = sizeAll1;
            aM.swap(aM1);
            bM[s].swap(bM1);
          }
          else {
            MatCreateSeqDense(PETSC_COMM_SELF, sizeAll, sizeAll, NULL, &A);
            for(int i = 0; i < sizeAll; i++) {
              for(int j = 0; j < sizeAll; j++) {
                MatSetValues(A, 1, &i, 1, &j, &aM[i * sizeAll + j], INSERT_VALUES);
              }
            }
            MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
          }

          VecDestroy(&v);
        }

        //END DEFLATION

        EPSCreate(PETSC_COMM_SELF, &eps);
        EPSSetOperators(eps, A, B);
        EPSSetFromOptions(eps);
        EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE);
        EPSSolve(eps);

        EPSGetEigenpair(eps, 0, &real, NULL, NULL, NULL);
        std::cout << real << " " << std::endl;

        sol->_Sol[CMIndex[s]]->set(iel, real);

        EPSDestroy(&eps);
        MatDestroy(&A);
        MatDestroy(&B);

      }

      for(unsigned s = 0; s < 2; s++) {
        
        sizeAll = sizeAll0;
   
        //BEGIN DEFLATION

        unsigned sizeAll1 = dim * (nDofu - 1);
        aL1.resize(sizeAll1 * sizeAll1);
        bL1.resize(sizeAll1 * sizeAll1);

        MatCreateSeqDense(PETSC_COMM_SELF, sizeAll1, sizeAll1, NULL, &B);

        for(int k = 0; k < dim; k++) {
          for(int i = 0; i < nDofu - 1; i++) {

            int ip = i + 1;
            int i1 = (nDofu - 1) * k + i;
            for(int l = 0; l < dim; l++) {
              for(int j = 0; j < nDofu - 1; j++) {
                int jp = j + 1;
                int j1 = (nDofu - 1) * l + j;
                double value;
                value = aL0[((nDofu * k) + ip) * sizeAll0 + (nDofu * l + jp)] - aL0[(nDofu * k) * sizeAll0 + (nDofu * l + jp)];
                aL1[i1 * sizeAll1 + j1] = value;

                value = bL[s][((nDofu * k) + ip) * sizeAll0 + (nDofu * l + jp)] - bL[s][(nDofu * k) * sizeAll0 + (nDofu * l + jp)];
                bL1[i1 * sizeAll1 + j1] = value;
                MatSetValues(B, 1, &i1, 1, &j1, &value, INSERT_VALUES);

              }
            }
          }
        }

        MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);

        sizeAll = sizeAll1;
        aL.swap(aL1);
        bL[s].swap(bL1);

        double real = 0.;
        while(fabs(real) < 1.0e-10) {

          MatCreateVecs(B, &v, NULL);

          EPSCreate(PETSC_COMM_SELF, &eps);
          EPSSetOperators(eps, B, NULL);
          EPSSetFromOptions(eps);
          EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE);
          EPSSolve(eps);

          double imaginary;
          EPSGetEigenpair(eps, 0, &real, &imaginary, v, NULL);
          EPSDestroy(&eps);

          if(fabs(real) < 1.0e-10) {
            PetscScalar *pv;
            VecGetArray(v, &pv);
            unsigned ii = 0;
            for(unsigned i = 1; i < sizeAll; i++){
              if(fabs(pv[i]) > fabs(pv[ii])) ii = i;
            }

            unsigned sizeAll1 = sizeAll - 1;

            aL1.resize(sizeAll1 * sizeAll1);
            bL1.resize(sizeAll1 * sizeAll1);

            MatDestroy(&B);

            MatCreateSeqDense(PETSC_COMM_SELF, sizeAll1, sizeAll1, NULL, &B);

            for(unsigned i = 0; i < sizeAll; i++) {
              if(i != ii) {
                int i1 = (i < ii) ? i : i - 1;
                for(unsigned j = 0; j < sizeAll; j++) {
                  if(j != ii) {
                    int j1 = (j < ii) ? j : j - 1;
                    double value;
                    value = aL[i * sizeAll + j] - 1. / pv[ii] * pv[i] * aL[ii * sizeAll + j];
                    aL1[i1 * sizeAll1 + j1] = value;

                    value = bL[s][i * sizeAll + j] - 1. / pv[ii] * pv[i] * bL[s][ii * sizeAll + j];
                    bL1[i1 * sizeAll1 + j1] = value;
                    MatSetValues(B, 1, &i1, 1, &j1, &value, INSERT_VALUES);
                  }
                }
              }
            }
            MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
            VecRestoreArray(v, &pv);

            sizeAll = sizeAll1;
            aL.swap(aL1);
            bL[s].swap(bL1);
          }
          else {
            MatCreateSeqDense(PETSC_COMM_SELF, sizeAll, sizeAll, NULL, &A);
            for(int i = 0; i < sizeAll; i++) {
              for(int j = 0; j < sizeAll; j++) {
                MatSetValues(A, 1, &i, 1, &j, &aL[i * sizeAll + j], INSERT_VALUES);
              }
            }
            MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
          }

          VecDestroy(&v);
        }

        //END DEFLATION

        EPSCreate(PETSC_COMM_SELF, &eps);
        EPSSetOperators(eps, A, B);
        EPSSetFromOptions(eps);
        EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE);
        EPSSolve(eps);

        EPSGetEigenpair(eps, 0, &real, NULL, NULL, NULL);
        std::cout << real << " " << std::endl;

        sol->_Sol[CLIndex[s]]->set(iel, real);

        EPSDestroy(&eps);
        MatDestroy(&A);
        MatDestroy(&B);

      }
    }
  }

  sol->_Sol[CMIndex[0]]->close();
  sol->_Sol[CMIndex[1]]->close();

  sol->_Sol[CLIndex[0]]->close();
  sol->_Sol[CLIndex[1]]->close();
}



