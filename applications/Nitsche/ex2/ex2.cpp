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

using namespace femus;

Line* line1;
Line* line2;
Line* lineI;

unsigned DIM = 3;

void AssembleNitscheProblem_AD(MultiLevelProblem& mlProb);

void BuildFlag(MultiLevelSolution& mlSol);
void GetInterfaceElementEigenvalues(MultiLevelSolution& mlSol);

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = false; //dirichlet
  if(DIM == 1) {
    if(facename == 1  && !strcmp(SolName, "u1")) dirichlet = true;
    if(facename == 2  && !strcmp(SolName, "u2")) dirichlet = true;
  }
  else if(DIM == 2) {
    if(facename == 4  && !strcmp(SolName, "u1")) dirichlet = true;
    if(facename == 2  && !strcmp(SolName, "u2")) dirichlet = true;
  }
  else if(DIM == 3) {
    if(facename == 5  && !strcmp(SolName, "u1")) dirichlet = true;
    if(facename == 3  && !strcmp(SolName, "u2")) dirichlet = true;
  }
  value = 0.;
  return dirichlet;
}


int main(int argc, char** args) {

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

  if(DIM == 1) {
    mlMsh.GenerateCoarseBoxMesh(nx, 0, 0, -length / 2, length / 2, 0., 0., 0., 0., EDGE3, "seventh");
  }
  else if(DIM == 2) {
    mlMsh.GenerateCoarseBoxMesh(nx, ny, 0, -length / 2, length / 2, -length / 2, length / 2, 0., 0., QUAD9, "seventh");
  }
  else if(DIM == 3) {
    nz = ny;
    mlMsh.GenerateCoarseBoxMesh(nx, ny, nz, -length / 2, length / 2, -length / 2, length / 2, -length / 2, length / 2, HEX27, "seventh");
  }

  //mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);

  FEOrder femOrder = SECOND;
  
  mlSol.AddSolution("u1", LAGRANGE, femOrder);
  mlSol.AddSolution("u2", LAGRANGE, femOrder);

  mlSol.AddSolution("eflag", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("nflag", LAGRANGE, femOrder, 0, false);

  mlSol.AddSolution("C1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("C2", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

  mlSol.Initialize("All");

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");

  BuildFlag(mlSol);

  MultiLevelProblem ml_prob(&mlSol);

  // ******* Add FEM system to the MultiLevel problem *******
  LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ("Nitsche");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("u1");
  system.AddSolutionToSystemPDE("u2");

  // attach the assembling function to system
  system.SetAssembleFunction(AssembleNitscheProblem_AD);

  // time loop parameter
  system.SetMaxNumberOfLinearIterations(1);

  system.init();

  //init marker

  unsigned Ne = 4;

  double Lx = length / nx;
  double Lx1 = 0.005 * Lx; //beam dimensions
  double Lx2 = Lx - Lx1; //beam dimensions

  double Ly = length;
  double Lz = length;

  unsigned rowy = ny * Ne;
  unsigned rowz = (DIM == 3) ?  nz * Ne : 1.;
  double dy = length / rowy; // size in y of each row
  double dz = (DIM == 3) ? length / rowz : 1.; // size in y of each row

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
  double y0 = -length / 2.;
  double z0 = -length / 2.;

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
  PrintLine("./output1", "bulk1", line1Points, 0);


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
  PrintLine("./output1", "bulk2", line2Points, 0);

  //interface marker initialization

  size = rowy * rowz;

  x.resize(size);
  markerType.resize(size);

  for(unsigned j = 0; j < size; j++) {
    x[j].assign(DIM, 0.);
    markerType[j] = INTERFACE;
  }

  std::vector < double > area(x.size(), length * length / x.size());  // uniform marker volume
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
  PrintLine("./output1", "interfaceLine", lineIPoints, 0);
  //END interface markers


  // ******* Print solution *******
  mlSol.SetWriter(VTK);
  mlSol.GetWriter()->SetDebugOutput(true);

  std::vector<std::string> print_vars;
  print_vars.push_back("All");

  GetInterfaceElementEigenvalues(mlSol);

  system.MGsolve();

  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "linear", print_vars, 0);

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

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned solu1Index = mlSol->GetIndex("u1");    // get the position of "u" in the ml_sol object
  unsigned solu2Index = mlSol->GetIndex("u2");    // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType(solu1Index);    // get the finite element type for "u"

  unsigned CIndex[2];
  
  CIndex[0] = mlSol->GetIndex("C1");
  CIndex[1] = mlSol->GetIndex("C2");
    
  
  unsigned solu1PdeIndex;
  solu1PdeIndex = mlPdeSys->GetSolPdeIndex("u1");    // get the position of "u" in the pdeSys object
  unsigned solu2PdeIndex;
  solu2PdeIndex = mlPdeSys->GetSolPdeIndex("u2");    // get the position of "u" in the pdeSys object

  vector < adept::adouble >  solu1; // local solution
  vector < adept::adouble >  solu2; // local solution

  unsigned eflagIndex = mlSol->GetIndex("eflag");
  unsigned nflagIndex = mlSol->GetIndex("nflag");

  vector < unsigned >  nodeFlag; // local solution

  vector < vector < double > > x(dim);    // local coordinates. x is now dim x m matrix.
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  vector< adept::adouble > aResu1; // local redidual vector
  vector< adept::adouble > aResu2; // local redidual vector

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
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);  // number of solution element dofs

    unsigned eFlag = static_cast <unsigned>(floor((*sol->_Sol[eflagIndex])(iel) + 0.5));

    // resize local arrays
    l2GMap.resize(2 * nDofu);
    solu1.resize(nDofu);
    solu2.resize(nDofu);
    nodeFlag.resize(nDofu);

    for(int k = 0; k < dim; k++) {
      x[k].resize(nDofu);
    }

    aResu1.assign(nDofu, 0.);    //resize
    aResu2.assign(nDofu, 0.);    //resize

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu1[i] = (*sol->_Sol[solu1Index])(solDof);                  // global extraction and local storage for the solution
      solu2[i] = (*sol->_Sol[solu2Index])(solDof);                  // global extraction and local storage for the solution
      nodeFlag[i] = (*sol->_Sol[nflagIndex])(solDof);
      l2GMap[i] = pdeSys->GetSystemDof(solu1Index, solu1PdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
      l2GMap[nDofu + i] = pdeSys->GetSystemDof(solu2Index, solu2PdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofu; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, 2);    // global to global mapping between coordinates node and coordinate dof
      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k])(xDof);      // global extraction and local storage for the element coordinates
      }
    }


    double alpha1 = .5;
    double alpha2 = 3.;

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();
    if(eFlag == 0 || eFlag == 2) {
      // *** Element Gauss point loop ***
      for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
        // *** get gauss point weight, test function and test function partial derivatives ***
        msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x);

        // evaluate the solution, the solution derivatives and the coordinates in the gauss point
        vector < adept::adouble > gradSolu1g(dim, 0.);
        vector < adept::adouble > gradSolu2g(dim, 0.);

        for(unsigned i = 0; i < nDofu; i++) {
          for(unsigned k = 0; k < dim; k++) {
            gradSolu1g[k] += phi_x[i * dim + k] * solu1[i];
            gradSolu2g[k] += phi_x[i * dim + k] * solu2[i];
          }
        }

        // *** phi_i loop ***
        for(unsigned i = 0; i < nDofu; i++) {

          adept::adouble graduGradphi1 = 0.;
          adept::adouble graduGradphi2 = 0.;

          for(unsigned k = 0; k < dim; k++) {
            graduGradphi1 += gradSolu1g[k] * phi_x[i * dim + k];
            graduGradphi2 += gradSolu2g[k] * phi_x[i * dim + k];
          }


          if(eFlag == 0) {
            aResu1[i] += (- phi[i] + alpha1 * graduGradphi1) * weight;
            if(nodeFlag[i] != 1) {
              aResu2[i] += (graduGradphi2) * weight;
            }
          }
          else if(eFlag == 2) {
            aResu2[i] += (- phi[i] + alpha2 * graduGradphi2) * weight;
            if(nodeFlag[i] != 1) {
              aResu1[i] += (graduGradphi1) * weight;
            }
          }

        } // end phi_i loop
      } // end gauss point loop
    }

    else {
      
      double ia1C1 = 1./ ( alpha1 * (*sol->_Sol[CIndex[0]])(iel) );
      double ia2C2 = 1./ ( alpha2 * (*sol->_Sol[CIndex[1]])(iel) );
      
      double den = ia1C1 + ia2C2;
      
      double gamma1 = ia1C1 / den;
      double gamma2 = ia2C2 / den;
      
      double theta = 2000. / den;
      
      

      //bulk1
      while(imarker1 < markerOffset1[iproc + 1] && iel == particle1[imarker1]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particle1[imarker1]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][soluType]->Jacobian(x, xi, weight, phi, phi_x);
        double weight = particle1[imarker1]->GetMarkerMass();

        // evaluate the solution, the solution derivatives and the coordinates in the gauss point
        vector < adept::adouble > gradSolu1g(dim, 0.);
        for(unsigned i = 0; i < nDofu; i++) {
          for(unsigned k = 0; k < dim; k++) {
            gradSolu1g[k] += phi_x[i * dim + k] * solu1[i];
          }
        }

        // *** phi_i loop ***
        for(unsigned i = 0; i < nDofu; i++) {
          adept::adouble graduGradphi1 = 0.;
          for(unsigned k = 0; k < dim; k++) {
            graduGradphi1 += gradSolu1g[k] * phi_x[i * dim + k];
          }
          aResu1[i] += (-phi[i] + alpha1 * graduGradphi1) * weight;
        }
        imarker1++;
      }

      //bulk2
      while(imarker2 < markerOffset2[iproc + 1] && iel == particle2[imarker2]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particle2[imarker2]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielGeom][soluType]->Jacobian(x, xi, weight, phi, phi_x);
        double weight = particle2[imarker2]->GetMarkerMass();

        // evaluate the solution, the solution derivatives and the coordinates in the gauss point
        vector < adept::adouble > gradSolu2g(dim, 0.);
        for(unsigned i = 0; i < nDofu; i++) {
          for(unsigned k = 0; k < dim; k++) {
            gradSolu2g[k] += phi_x[i * dim + k] * solu2[i];
          }
        }

        // *** phi_i loop ***
        for(unsigned i = 0; i < nDofu; i++) {
          adept::adouble graduGradphi2 = 0.;
          for(unsigned k = 0; k < dim; k++) {
            graduGradphi2 += gradSolu2g[k] * phi_x[i * dim + k];
          }
          aResu2[i] += (-phi[i] + alpha2 * graduGradphi2) * weight;
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

        adept::adouble solu1g  = 0.;
        adept::adouble solu2g  = 0.;
        adept::adouble alphaGradSoluDotN = 0.;

        for(unsigned i = 0; i < nDofu; i++) {
          solu1g += phi[i] * solu1[i];
          solu2g += phi[i] * solu2[i];
          for(unsigned k = 0; k < dim; k++) {
            alphaGradSoluDotN += (alpha1 * gamma1  * solu1[i] + alpha2 * gamma2 * solu2[i]) * phi_x[i * dim + k] * N[k];
          }
        }
        // *** phi_i loop ***
        for(unsigned i = 0; i < nDofu; i++) {

          adept::adouble gradPhiDotN = 0.;
          for(unsigned k = 0; k < dim; k++) {
            gradPhiDotN += phi_x[i * dim + k] * N[k];
          }

          aResu1[i] += (solu2g - solu1g) * (-phi[i] * theta - alpha1 * gamma1 * gradPhiDotN) * weight;
          aResu1[i] += alphaGradSoluDotN * (+phi[i]) * weight;
          aResu2[i] += (solu2g - solu1g) * (+phi[i] * theta - alpha2 * gamma2 * gradPhiDotN) * weight;
          aResu2[i] += alphaGradSoluDotN * (-phi[i]) * weight;
        } // end phi_i loop
        imarkerI++;
      }

    }

    //copy the value of the adept::adoube aRes in double Res and store
    Res.resize(2 * nDofu);    //resize

    for(int i = 0; i < nDofu; i++) {
      Res[i] = - aResu1[i].value();
      Res[nDofu + i] = - aResu2[i].value();
    }

    RES->add_vector_blocked(Res, l2GMap);

    // define the dependent variables
    s.dependent(&aResu1[0], nDofu);
    s.dependent(&aResu2[0], nDofu);

    // define the independent variables
    s.independent(&solu1[0], nDofu);
    s.independent(&solu2[0], nDofu);

    // get the jacobian matrix (ordered by row major )
    Jac.resize(2 * nDofu * 2 * nDofu);    //resize
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

  unsigned CIndex[2];
  CIndex[0] = mlSol.GetIndex("C1");
  CIndex[1] = mlSol.GetIndex("C2");

  unsigned solIndex = mlSol.GetIndex("u1");
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

  sol->_Sol[CIndex[0]]->zero();
  sol->_Sol[CIndex[1]]->zero();

  std::vector < double > a;
  std::vector < std::vector < double > > b(2);

  Mat A, B;
  EPS eps;

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    unsigned eFlag = static_cast <unsigned>(floor((*sol->_Sol[eflagIndex])(iel) + 0.5));
    if(eFlag == 1) {

      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);  // number of solution element dofs

      a.assign(nDofu * nDofu, 0.);
      b[0].assign(nDofu * nDofu, 0.);
      b[1].assign(nDofu * nDofu, 0.);

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
        for(int i = 0; i < nDofu; i++) {
          for(int j = 0; j < nDofu; j++) {
            for(unsigned k = 0; k < dim; k++) {
              b[0][ i * nDofu + j] += phi_x[i * dim + k] * phi_x[j * dim + k] * weight;
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
        for(int i = 0; i < nDofu; i++) {
          for(int j = 0; j < nDofu; j++) {
            for(unsigned k = 0; k < dim; k++) {
              b[1][ i * nDofu + j] += phi_x[i * dim + k] * phi_x[j * dim + k] * weight;
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
        for(int i = 0; i < nDofu; i++) {

          double gradPhiiDotN = 0.;
          for(unsigned k = 0; k < dim; k++) {
            gradPhiiDotN += phi_x[i * dim + k] * N[k];
          }
          for(int j = 0; j < nDofu; j++) {
            double gradPhijDotN = 0.;
            for(unsigned k = 0; k < dim; k++) {
              gradPhijDotN += phi_x[j * dim + k] * N[k];
            }
            a[ i * nDofu + j] += gradPhiiDotN * gradPhijDotN  * weight;
          }
        } // end phi_i loop
        imarkerI++;
      }

      
      /* Careful B has one zero eigenvalue (it is singular) with nullspace x = [1,1,1,...]^T, 
       * thus the generalized eigenvalue problem 
       * $$A u = \lambda B u$$ 
       * (or $B^{-1} A x = \lambda x$) has one indetermined eigenvalue, that makes the SLEPC solve very unstable, 
       * even using its built-in deflation method. 
       * Fortunately, A has at least one zero eigenvalue, with the same x = [1,1,1,...]^T being an element 
       * of its nullspace. Then, it is possible to deflate A and B simultaneously:
       * $Ad = A - x^T. a1$ and $Bd = B - x^T.b1$, where a1 and b1 are the first rows of A and B, respectively.
       * The generalized eigenvalue problem $Ab u = \lambda Bb u$, with matrices Ab and Bb, 
       * obtained as block matrices from Ad and Bd, removing the first row and the first column, 
       * has the same eigenvalues of the original one, except the indetermine one. 
       * Note that Bb is now invertible, and SLEPC has no problem in solving the deflated 
       * generalized eigenvalue problem.
       */
      
      for(unsigned k = 0; k < 2; k++) {
        MatCreateSeqDense(PETSC_COMM_SELF, nDofu - 1, nDofu - 1, NULL, &A);
        MatCreateSeqDense(PETSC_COMM_SELF, nDofu - 1, nDofu - 1, NULL, &B);

        for(int i = 0; i < nDofu - 1; i++) {
          for(int j = 0; j < nDofu - 1; j++) {
            double value;
            value = b[k][(i + 1) * nDofu + (j + 1)] - b[k][j + 1];
            MatSetValues(B, 1, &i, 1, &j, &value, INSERT_VALUES);
            value = a[(i + 1) * nDofu + (j + 1)] - a[j + 1];
            MatSetValues(A, 1, &i, 1, &j, &value, INSERT_VALUES);
          }
        }

        MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

        MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);

        EPSCreate(PETSC_COMM_SELF, &eps);
        EPSSetOperators(eps, A, B);
        EPSSetFromOptions(eps);
        EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE);
        EPSSolve(eps);

        double real, imaginary;
        EPSGetEigenpair(eps, 0, &real, &imaginary, NULL, NULL);
        //std::cout << real << " " << imaginary << std::endl;
        sol->_Sol[CIndex[k]]->set(iel, real);

        EPSDestroy(&eps);
        MatDestroy(&A);
        MatDestroy(&B);
        
      }
    }
  }

  sol->_Sol[CIndex[0]]->close();
  sol->_Sol[CIndex[1]]->close();
}


