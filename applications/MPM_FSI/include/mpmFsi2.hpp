#include "MultiLevelSolution.hpp"
#include "PetscMatrix.hpp"

using namespace femus;

double beta = 0.25;
double Gamma = 0.5;
//double gravity[3] = {9810, 0., 0.};
double gravity[3] = {0, 0., 0.};

Line* solidLine;
Line* interfaceLine;
Line* fluidLine;

const double faceNormal[6][6][3] = {
  {{0., -1., 0.}, {1., 0., 0.}, {0., 1., 0.}, { -1., 0., 0.}, {0., 0., -1.}, {0., 0., 1.}},
  {{0., 0., -1}, {0., -1., 0.}, {1. / sqrt (3.), 1. / sqrt (3.), 1. / sqrt (3.) }, { -1., 0., 0.}},
  {{0., -1., 0.}, {1. / sqrt (2.), 1. / sqrt (2.), 0.}, { -1., 0., 0.}, {0., 0., -1}, {0., 0., 1}},
  {{0., -1.}, {1., 0.}, {0., 1.}, { -1., 0.}},
  {{0., -1.}, {1. / sqrt (2.), 1. / sqrt (2.) }, { -1., 0.}},
  {{ -1}, {1}}
};

const unsigned faceNumber[6] = {6, 4, 5, 4, 3, 2};

const unsigned facePointNumber[6][3] = {{}, {}, {}, {5, 9, 9}, {4, 7, 7}, {}}; //facePointNumber[elemType][solType]

const unsigned facePoints[6][3][9] = { //facePointNumber[elemType][solType][localFaceIndex]
  {{ }},
  {{ }},
  {{ }},
  { { 0, 1, 2, 3, 0}, { 0, 4, 1, 5, 2, 6, 3, 7, 0 }, { 0, 4, 1, 5, 2, 6, 3, 7, 0 } },
  { {0, 1, 2, 0}, { 0, 3, 1, 4, 2, 5, 0 }, { 0, 3, 1, 4, 2, 5, 0 } },
  { }
};

const unsigned trianglesPerFace[6][3][6] = { //trainglesPerFace[elementType][solType][numberOfTrianglesPerFace]
  {{4, 4, 4, 4, 4, 4}, {8, 8, 8, 8, 8, 8}, {8, 8, 8, 8, 8, 8}},
  {{1, 1, 1, 1}, {6, 6, 6, 6}, {6, 6, 6, 6}},
  {{4, 4, 4, 1, 1}, {8, 8, 8, 6, 6}, {8, 8, 8, 6, 6}},
  {{4}, {8}, {8}},
  {{1}, {6}, {6}},
  {}
};

const unsigned faceTriangleNodes[6][3][6][8][4] = { // elementType - solType - number of faces - number of triangles per faces - 3 vertices of the triangle (the first is counted twice)

  { { {{0, 1, 20, 0}, {1, 5, 20, 1}, {5, 4, 20, 5}, {4, 0, 20, 4}},  //hex
      {{1, 2, 21, 1}, {2, 6, 21, 2}, {6, 5, 21, 6}, {5, 1, 21, 5}},
      {{2, 3, 22, 2}, {3, 7, 22, 3}, {7, 6, 22, 7}, {6, 2, 22, 6}},
      {{3, 0, 23, 3}, {0, 4, 23, 0}, {4, 7, 23, 4}, {7, 3, 23, 7}},
      {{0, 3, 24, 0}, {3, 2, 24, 3}, {2, 1, 24, 2}, {1, 0, 24, 1}},
      {{4, 5, 25, 4}, {5, 6, 25, 5}, {6, 7, 25, 6}, {7, 4, 25, 7}}
    },
    { {{0, 8, 20, 0}, {8, 1, 20, 8}, {1, 17, 20, 1}, {17, 5, 20, 17}, {5, 12, 20, 5}, {12, 4, 20, 12}, {4, 16, 20, 4 }, {16, 0, 20, 16}},
      {{1, 9, 21, 1}, {9, 2, 21, 9}, {2, 18, 21, 2}, {18, 6, 21, 18}, {6, 13, 21, 6}, {13, 5, 21, 13}, {5, 17, 21, 5}, {17, 1, 21, 17}},
      {{2, 10, 22, 2}, {10, 3, 22, 10}, {3, 19, 22, 3}, {19, 7, 22, 19}, {7, 14, 22, 7}, {14, 6, 22, 14}, {6, 18, 22, 6}, {18, 2, 22, 18}},
      {{3, 11, 23, 3}, {11, 0, 23, 11}, {0, 16, 23, 0}, {16, 4, 23, 16}, {4, 15, 23, 4}, {15, 7, 23, 15}, {7, 19, 23, 7}, {19, 3, 23, 19}},
      {{0, 11, 24, 0}, {11, 3, 24, 11}, {3, 10, 24, 3}, {10, 2, 24, 10}, {2, 9, 24, 2}, {9, 1, 24, 9}, {1, 8, 24, 1}, {8, 0, 24, 8}},
      {{4, 12, 25, 4}, {12, 5, 25, 12}, {5, 13, 25, 5}, {13, 6, 25, 13}, {6, 14, 25, 6}, {14, 7, 25, 14}, {7, 15, 25, 7}, {15, 4, 25, 15}}
    },
    { {{0, 8, 20, 0}, {8, 1, 20, 8}, {1, 17, 20, 1}, {17, 5, 20, 17}, {5, 12, 20, 5}, {12, 4, 20, 12}, {4, 16, 20, 4 }, {16, 0, 20, 16}},
      {{1, 9, 21, 1}, {9, 2, 21, 9}, {2, 18, 21, 2}, {18, 6, 21, 18}, {6, 13, 21, 6}, {13, 5, 21, 13}, {5, 17, 21, 5}, {17, 1, 21, 17}},
      {{2, 10, 22, 2}, {10, 3, 22, 10}, {3, 19, 22, 3}, {19, 7, 22, 19}, {7, 14, 22, 7}, {14, 6, 22, 14}, {6, 18, 22, 6}, {18, 2, 22, 18}},
      {{3, 11, 23, 3}, {11, 0, 23, 11}, {0, 16, 23, 0}, {16, 4, 23, 16}, {4, 15, 23, 4}, {15, 7, 23, 15}, {7, 19, 23, 7}, {19, 3, 23, 19}},
      {{0, 11, 24, 0}, {11, 3, 24, 11}, {3, 10, 24, 3}, {10, 2, 24, 10}, {2, 9, 24, 2}, {9, 1, 24, 9}, {1, 8, 24, 1}, {8, 0, 24, 8}},
      {{4, 12, 25, 4}, {12, 5, 25, 12}, {5, 13, 25, 5}, {13, 6, 25, 13}, {6, 14, 25, 6}, {14, 7, 25, 14}, {7, 15, 25, 7}, {15, 4, 25, 15}}
    }
  },
  { { {{0, 2, 1, 0}},   //tet
      {{0, 1, 3, 0}},
      {{1, 2, 3, 1}},
      {{2, 0, 3, 2}}
    },
    { {{0, 6, 10, 0}, {6, 2, 10, 6}, {2, 5, 10, 2}, {5, 1, 10, 5}, {1, 4, 10, 1}, {4, 0, 10, 4}},
      {{0, 4, 11, 0}, {4, 1, 11, 4}, {1, 8, 11, 1}, {8, 3, 11, 8}, {3, 7, 11, 3}, {7, 0, 11, 7}},
      {{1, 5, 12, 1}, {5, 2, 12, 5}, {2, 9, 12, 2}, {9, 3, 12, 9}, {3, 8, 12, 3}, {8, 1, 12, 8}},
      {{2, 6, 13, 2}, {6, 0, 13, 6}, {0, 7, 13, 0}, {7, 3, 13, 7}, {3, 9, 13, 3}, {9, 2, 13, 9}}
    },
    { {{0, 6, 10, 0}, {6, 2, 10, 6}, {2, 5, 10, 2}, {5, 1, 10, 5}, {1, 4, 10, 1}, {4, 0, 10, 4}},
      {{0, 4, 11, 0}, {4, 1, 11, 4}, {1, 8, 11, 1}, {8, 3, 11, 8}, {3, 7, 11, 3}, {7, 0, 11, 7}},
      {{1, 5, 12, 1}, {5, 2, 12, 5}, {2, 9, 12, 2}, {9, 3, 12, 9}, {3, 8, 12, 3}, {8, 1, 12, 8}},
      {{2, 6, 13, 2}, {6, 0, 13, 6}, {0, 7, 13, 0}, {7, 3, 13, 7}, {3, 9, 13, 3}, {9, 2, 13, 9}}
    }
  },
  { { {{0, 1, 15, 0}, {1, 4, 15, 1}, {4, 3, 15, 4}, {3, 0, 15, 3}},  //wedge
      {{1, 2, 16, 1}, {2, 5, 16, 2}, {5, 4, 16, 5}, {4, 1, 16, 4}},
      {{2, 0, 17, 2}, {0, 3, 17, 0}, {3, 5, 17, 3}, {5, 2, 17, 5}},
      {{0, 2, 1, 0}},
      {{3, 4, 5, 3}}
    },
    { {{0, 6, 15, 0 }, {6, 1, 15, 6}, {1, 13, 15, 1}, {13, 4, 15, 13}, {4, 9, 15, 4}, {9, 3, 15, 9}, {3, 12, 15, 3}, {12, 0, 15, 12}},
      {{1, 7, 16, 1}, {7, 2, 16, 7}, {2, 14, 16, 2}, {14, 5, 16, 14}, {5, 10, 16, 5}, {10, 4, 16, 10}, {4, 13, 16, 4}, {13, 1, 16, 13}},
      {{2, 8, 17, 2}, {8, 0, 17, 8}, {0, 12, 17, 0}, {12, 3, 17, 12}, {3, 11, 17, 3}, {11, 5, 17, 11}, {5, 14, 17, 5}, {14, 2, 17, 14}},
      {{0, 8, 18, 0}, {8, 2, 18, 8}, {2, 7, 18, 2}, {7, 1, 18, 7}, {1, 6, 18, 1}, {6, 0, 18, 6}},
      {{3, 9, 19, 3}, {9, 4, 19, 9}, {4, 10, 19, 4}, {10, 5, 19, 10}, {5, 11, 19, 5}, {11, 3, 19, 11}}
    },
    { {{0, 6, 15, 0 }, {6, 1, 15, 6}, {1, 13, 15, 1}, {13, 4, 15, 13}, {4, 9, 15, 4}, {9, 3, 15, 9}, {3, 12, 15, 3}, {12, 0, 15, 12}},
      {{1, 7, 16, 1}, {7, 2, 16, 7}, {2, 14, 16, 2}, {14, 5, 16, 14}, {5, 10, 16, 5}, {10, 4, 16, 10}, {4, 13, 16, 4}, {13, 1, 16, 13}},
      {{2, 8, 17, 2}, {8, 0, 17, 8}, {0, 12, 17, 0}, {12, 3, 17, 12}, {3, 11, 17, 3}, {11, 5, 17, 11}, {5, 14, 17, 5}, {14, 2, 17, 14}},
      {{0, 8, 18, 0}, {8, 2, 18, 8}, {2, 7, 18, 2}, {7, 1, 18, 7}, {1, 6, 18, 1}, {6, 0, 18, 6}},
      {{3, 9, 19, 3}, {9, 4, 19, 9}, {4, 10, 19, 4}, {10, 5, 19, 10}, {5, 11, 19, 5}, {11, 3, 19, 11}}
    }
  },
  { { {{0, 1, 8, 0}, {1, 2, 8, 1}, {2, 3, 8, 1}, {3, 0, 8, 3}}, //quad
      {{0, 4, 8, 0}, {4, 1, 8, 4}, {1, 5, 8, 1}, {5, 2, 8, 5}, {2, 6, 8, 2}, {6, 3, 8, 6}, {3, 7, 8, 3}, {7, 0, 8, 7}},
      {{0, 4, 8, 0}, {4, 1, 8, 4}, {1, 5, 8, 1}, {5, 2, 8, 5}, {2, 6, 8, 2}, {6, 3, 8, 6}, {3, 7, 8, 3}, {7, 0, 8, 7}}
    }
  },
  { { {{0, 1, 2, 0}},
      {{0, 3, 6, 0}, {3, 1, 6, 3}, {1, 4, 6, 1}, {4, 2, 6, 4}, {2, 5, 6, 2}, {5, 0, 6, 5}},
      {{0, 3, 6, 0}, {3, 1, 6, 3}, {1, 4, 6, 1}, {4, 2, 6, 4}, {2, 5, 6, 2}, {5, 0, 6, 5}}
    }
  },
  {{{{}}}},
};

void GetParticlesToNodeFlag (MultiLevelSolution &mlSol, Line & solidLine, Line & fluidLine);
void GetPressureNeighbor (MultiLevelSolution &mlSol, Line & solidLine, Line & fluidLine);

void ProjectGridVelocity (MultiLevelSolution &mlSol);
void ProjectGridVelocity2 (MultiLevelSolution &mlSol);

bool CheckInclusion2D (Solution* sol, const unsigned &elemToCheck, const unsigned &previousElem, std::vector <std::vector <double>> & xElement, std::vector <double> & xToCheck);
bool CheckInclusion3D (Solution* sol, const unsigned &elemToCheck, const unsigned &previousElem, const std::vector <std::vector <double>> & xElement, const std::vector <double> & xToCheck);
void FindLocalCoordinates (std::vector<double> & xi, std::vector < std::vector < std::vector < double > > >  & aX, const bool & pcElemUpdate, MultiLevelSolution &  mlSol, const unsigned &elemToCheck, const std::vector <double> &xToCheck, const std::vector<std::vector<double>> & xElement);

void AssembleMPMSys (MultiLevelProblem& ml_prob) {

  // ml_prob is the global object from/to where get/set all the data
  // level is the level of the PDE system to be assembled
  // levelMax is the Maximum level of the MultiLevelProblem
  // assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

  //pointers and references

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FSI");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = mlSol->GetSolutionLevel (level);    // pointer to the solution (level) object
  LinearEquationSolver* myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];  // pointer to the equation (level) object

  Mesh* msh = ml_prob._ml_msh->GetLevel (level);    // pointer to the mesh (level) object
  elem* el = msh->el;   // pointer to the elem object in msh (level)
  SparseMatrix* myKK = myLinEqSolver->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* myRES =  myLinEqSolver->_RES;  // pointer to the global residual vector object in pdeSys (level)

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  const unsigned dim = msh->GetDimension();

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned > (ceil (pow (3, dim)));       // conservative: based on line3, quad9, hex27

  // data
  unsigned iproc  = msh->processor_id();

  vector< vector< adept::adouble > > solD (dim);     // local solution (displacement)
  vector< vector< adept::adouble > > solV (dim);     // local solution (velocity)
  vector< adept::adouble > solP;     // local solution (velocity)
  vector< adept::adouble > solPOld;     // local solution (velocity)

  vector< vector< double > > solDOld (dim);     // local solution (displacement)
  vector< vector< double > > solVOld (dim);

  vector< double > rhs;    // local redidual vector
  vector< vector< adept::adouble > > aRhsD (dim);    // local redidual vector
  vector< vector< adept::adouble > > aRhsV (dim);    // local redidual vector
  vector< adept::adouble > aRhsP;    // local redidual vector
  vector < double > Jac;

  std::vector <unsigned> sysDofsAll;

  vector < bool > solidFlag;
  vector < bool > solidFlag1;
  vector < double > phi;
  vector < double > phiHat;
  vector < double > phiP;
  vector < adept::adouble> gradPhi;
  vector < double > gradPhiHat;

  vector <vector < adept::adouble> > vx (dim); //vx is coordX in assembly of ex30
  vector <vector < double> > vxHat (dim);

  adept::adouble weight;
  double weightHat;


  //reading parameters for MPM body
  double rhoMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_density();
  double EMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_young_module();
  double muMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_shear_modulus();
  double nuMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_poisson_coeff();
  double lambdaMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_lambda();
  double KMmp = EMpm / (3.* (1. - 2. * nuMpm)); //bulk modulus

  //reading parameters for fluid FEM domain
  double rhoFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_density();
  double muFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_viscosity();

  double dt =  my_nnlin_impl_sys.GetIntervalTime();

  std::cout.precision (10);

  //variable-name handling
  const char varname[10][5] = {"DX", "DY", "DZ", "VX", "VY", "VZ"};

  vector <unsigned> indexSolD (dim);
  vector <unsigned> indexSolV (dim);
  vector <unsigned> indexPdeD (dim);
  vector <unsigned> indexPdeV (dim);
  for (unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = mlSol->GetIndex (&varname[ivar][0]);
    indexSolV[ivar] = mlSol->GetIndex (&varname[ivar + 3][0]);
    indexPdeD[ivar] = my_nnlin_impl_sys.GetSolPdeIndex (&varname[ivar][0]);
    indexPdeV[ivar] = my_nnlin_impl_sys.GetSolPdeIndex (&varname[ivar + 3][0]);
  }
  unsigned solType = mlSol->GetSolutionType (&varname[0][0]);

  unsigned indexSolP = mlSol->GetIndex ("P");
  unsigned indexPdeP = my_nnlin_impl_sys.GetSolPdeIndex ("P");
  unsigned solTypeP = mlSol->GetSolutionType ("P");

  unsigned indexSolMat = mlSol->GetIndex ("Mat");
  unsigned solTypeMat = mlSol->GetSolutionType ("Mat");

  unsigned indexSolM = mlSol->GetIndex ("NodeFlag");
  unsigned indexSolM1 = mlSol->GetIndex ("M");

  //double  indexNodeDist = mlSol->GetIndex ("NodeDist");

  start_time = clock();

  myKK->zero();
  myRES->zero();

  std::vector<Marker*> particlesSolid = solidLine->GetParticles();
  std::vector<unsigned> markerOffsetSolid = solidLine->GetMarkerOffset();

  std::vector<unsigned> markerOffsetFluid = fluidLine->GetMarkerOffset();
  std::vector<Marker*> particlesFluid = fluidLine->GetParticles();

  std::vector<unsigned> markerOffsetInterface = interfaceLine->GetMarkerOffset();
  std::vector<Marker*> particlesInterface = interfaceLine->GetParticles();

  unsigned iSmarker = markerOffsetSolid[iproc];
  unsigned iFmarker = markerOffsetFluid[iproc];

  //BEGIN loop on elements (to initialize the "soft" stiffness matrix)
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielt = msh->GetElementType (iel);
    double  MPMmaterial = (*mysolution->_Sol[indexSolMat]) (iel);

    unsigned nDofs = msh->GetElementDofNumber (iel, solType);   // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber (iel, solTypeP);   // number of solution element dofs

    unsigned nDofsAll = 2 * dim * nDofs + nDofsP;
    // resize local arrays
    sysDofsAll.resize (nDofsAll);

    solidFlag.resize (nDofs);
    solidFlag1.resize (nDofs);

    for (unsigned  k = 0; k < dim; k++) {
      solD[k].resize (nDofs);
      solDOld[k].resize (nDofs);

      solV[k].resize (nDofs);
      solVOld[k].resize (nDofs);

      aRhsD[k].assign (nDofs, 0.);
      aRhsV[k].assign (nDofs, 0.);

      vx[k].resize (nDofs);
      vxHat[k].resize (nDofs);
    }
    solP.resize (nDofsP);
    solPOld.resize (nDofsP);
    aRhsP.assign (nDofsP, 0.);

    unsigned counter = 0;
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned idof = msh->GetSolutionDof (i, iel, solType);

      solidFlag[i] = ( (*mysolution->_Sol[indexSolM]) (idof) > 0.5) ? true : false;

      solidFlag1[i] = ( (*mysolution->_Sol[indexSolM1]) (idof) > 0.5) ? true : false;
      if (solidFlag1[i]) counter++;

      for (unsigned  k = 0; k < dim; k++) {
        solD[k][i] = (*mysolution->_Sol[indexSolD[k]]) (idof);
        solDOld[k][i] = (*mysolution->_SolOld[indexSolD[k]]) (idof);

        solV[k][i] = (*mysolution->_Sol[indexSolV[k]]) (idof);
        solVOld[k][i] = (*mysolution->_SolOld[indexSolV[k]]) (idof);

        sysDofsAll[i + k * nDofs] = myLinEqSolver->GetSystemDof (indexSolD[k], indexPdeD[k], i, iel);
        sysDofsAll[i + (k + dim) * nDofs] = myLinEqSolver->GetSystemDof (indexSolV[k], indexPdeV[k], i, iel);
      }
    }

    //bool test = (counter >= nDofs - 5) ? true : false;
    bool test = (counter > 0) ? true : false;

    for (unsigned i = 0; i < nDofsP; i++) {
      unsigned idof = msh->GetSolutionDof (i, iel, solTypeP);
      solP[i] = (*mysolution->_Sol[indexSolP]) (idof);
      solPOld[i] = (*mysolution->_SolOld[indexSolP]) (idof);
      sysDofsAll[i + (2 * dim) * nDofs] = myLinEqSolver->GetSystemDof (indexSolP, indexPdeP, i, iel);
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    for (unsigned i = 0; i < nDofs; i++) {
      unsigned idofX = msh->GetSolutionDof (i, iel, 2);
      for (unsigned  k = 0; k < dim; k++) {
        vxHat[k][i] = (*msh->_topology->_Sol[k]) (idofX) + solDOld[k][i];
        vx[k][i] = vxHat[k][i] + solD[k][i];
      }
    }

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {

      msh->_finiteElement[ielt][solType]->Jacobian (vxHat, ig, weightHat, phiHat, gradPhiHat);
      msh->_finiteElement[ielt][solType]->Jacobian (vx, ig, weight, phi, gradPhi);

      vector < adept::adouble > solVg (dim, 0.);
      vector < adept::adouble > solVgOld (dim, 0.);

      vector < adept::adouble > solDg (dim, 0.);
      vector < adept::adouble > solDgOld (dim, 0.);

      vector < vector < adept::adouble > > gradSolDgHat (dim);
      vector < vector < adept::adouble > > gradSolVg (dim);

      for (unsigned  k = 0; k < dim; k++) {
        gradSolDgHat[k].assign (dim, 0);
        gradSolVg[k].assign (dim, 0);
      }

      for (unsigned i = 0; i < nDofs; i++) {
        for (unsigned j = 0; j < dim; j++) {
          solVg[j] += phiHat[i] * solV[j][i];
          solDg[j] += phiHat[i] * solD[j][i];
          solVgOld[j] += phiHat[i] * solVOld[j][i];
          solDgOld[j] += phiHat[i] * solDOld[j][i];
          for (unsigned  k = 0; k < dim; k++) {
            gradSolDgHat[k][j] += gradPhiHat[i * dim + j] * solD[k][i];
            gradSolVg[k][j] += gradPhi[i * dim + j] * solV[k][i];
          }
        }
      }

      double *phiP = msh->_finiteElement[ielt][solTypeP]->GetPhi (ig);
      adept::adouble solPg = 0.;
      adept::adouble solPgOld = 0.;
      for (unsigned i = 0; i < nDofsP; i++) {
        solPg += phiP[i] * solP[i];
        solPgOld += phiP[i] * solPOld[i];
      }

      for (unsigned i = 0; i < nDofs; i++) {

        for (unsigned k = 0; k < dim; k++) {
          //adept::adouble  softStiffness  = 0.;
          adept::adouble  wlaplace1V  = 0.;
          adept::adouble  wlaplace1D  = 0.;
          for (unsigned  j = 0; j < dim; j++) {
            //softStiffness +=  muMpm * gradPhiHat[i * dim + j] * (gradSolDgHat[k][j] + gradSolDgHat[j][k]);
            wlaplace1V +=  gradPhi[i * dim + j] * (gradSolVg[k][j] + gradSolVg[j][k]);
            wlaplace1D +=  gradPhiHat[i * dim + j] * (gradSolDgHat[k][j] + gradSolDgHat[j][k]);
          }

//           if (!solidFlag[i]) { //kinematic equation in the fluid nodes
//             double stiffness = (MPMmaterial == 0) ? .00000000000000001 : 0.000001;
//             //aRhsD[k][i] += - stiffness * wlaplace1D * weightHat;
//             if (!solidFlag1[i]) {
//               aRhsD[k][i] += - stiffness * wlaplace1D * weightHat;
//             }
//             else {
//               //aRhsD[k][i] += phiHat[i] * (-0.00000000 * wlaplace1D + solVg[k] - (solDg[k] - solDgOld[k]) / dt) * weightHat;
//               //aRhsD[k][i] += phiHat[i] * (-muf * wlaplace1V + solP/*+ solV[k][i] - (solD[k][i] - solDOld[k][i]) / dt*/) * weightHat;
//
//               //std::cout << nodeDist[i]/dmax <<" ";
//
//               aRhsD[k][i] += ( (- muFluid * wlaplace1V + gradPhi[i * dim + k] * solPg) * (nodeDist[i]/dmax) +
//                               phiHat[i] * (solV[k][i] - (solD[k][i] - solDOld[k][i]) / dt)) * weight;
//             }
//           }
          if (!solidFlag1[i]) {
            double stiffness = (MPMmaterial == 0) ? 1 : 1;
            aRhsD[k][i] +=  - stiffness * wlaplace1D * weightHat;
          }
          else { //kinematic equation in the solid nodes
            aRhsV[k][i] += phiHat[i] * (solV[k][i] - (solD[k][i] - solDOld[k][i]) / dt) * weightHat;
          }
        }
      }

      if (MPMmaterial == 0) { //only cells that are completely fluid
        for (unsigned i = 0; i < nDofs; i++) {
          for (unsigned k = 0; k < dim; k++) {
            adept::adouble wlaplace = 0.;
            adept::adouble advection = 0.;
            for (unsigned j = 0; j < dim; j++) {
              wlaplace  +=  gradPhi[i * dim + j] * (gradSolVg[k][j] + gradSolVg[j][k]);
              advection  +=  phi[i] * (solVg[j] - (solDg[j] - solDgOld[j]) / dt) * gradSolVg[k][j];
            }
            if (!solidFlag1[i]) {
              aRhsV[k][i] += (- rhoFluid * phi[i] * (solVg[k] - solVgOld[k]) / dt - rhoFluid * advection - muFluid * wlaplace + gradPhi[i * dim + k] * solPg) * weight;
            }
            else {
              aRhsD[k][i] += (- rhoFluid * phi[i] * (solVg[k] - solVgOld[k]) / dt - rhoFluid * advection - muFluid * wlaplace + gradPhi[i * dim + k] * solPg) * weight;
            }
          }
        }
      }

      for (unsigned i = 0; i < nDofsP; i++) {
        for (unsigned  k = 0; k < dim; k++) {
//           if (MPMmaterial == nDofs) {  //all cells that are completely MPM solid
//             aRhsP[i] += phiP[i] * (solPg) * weight;
//           }
//           else if (MPMmaterial > 0 || test) {  //all cells that are completely MPM solid
//             aRhsP[i] += phiP[i] * (gradSolVg[k][k] + 0.01 * (solPg - solPgOld) / dt) * weight;
//           }
//           else { // if (MPMmaterial == 0) {
//             aRhsP[i] += phiP[i] *  gradSolVg[k][k] * weight;
//           }

          if (MPMmaterial == 0) {
            aRhsP[i] += phiP[i] *  gradSolVg[k][k] * weight;
          }
          else {  //all cells that are completely MPM solid
            aRhsP[i] += phiP[i] * (solPg) * weight;
          }

        }
      }
    } // end gauss point loop


    //BEGIN SOLID PARTICLE
    if (MPMmaterial > 0) { //solid markers
//       unsigned iSmarker = markerOffsetSolid[iproc];
//       while (iel != particlesSolid[imarker]->GetMarkerElement()) {
//         imarker++;
//       }
      while (iSmarker < markerOffsetSolid[iproc + 1] && iel == particlesSolid[iSmarker]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particlesSolid[iSmarker]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielt][solType]->Jacobian (vxHat, xi, weightHat, phiHat, gradPhiHat);

        std::vector <double> SolVpOld (dim);
        particlesSolid[iSmarker]->GetMarkerVelocity (SolVpOld);

        std::vector <double> SolApOld (dim);
        particlesSolid[iSmarker]->GetMarkerAcceleration (SolApOld);

        double mass = particlesSolid[iSmarker]->GetMarkerMass();

        msh->_finiteElement[ielt][solType]->Jacobian (vx, xi, weight, phi, gradPhi); //function to evaluate at the particles

        // displacement and velocity
        //BEGIN evaluates SolDp at the particle iSmarker
        vector<adept::adouble> SolDp (dim, 0.);
        vector<vector < adept::adouble > > gradSolDpHat (dim);
        for (int k = 0; k < dim; k++) {
          gradSolDpHat[k].assign (dim, 0.);
        }

        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nDofs; i++) {
            SolDp[j] += phi[i] * solD[j][i];
            for (int k = 0; k < dim; k++) {
              gradSolDpHat[j][k] +=  gradPhiHat[i * dim + k] * solD[j][i];
            }
          }
        }
        //END evaluates SolDp at the particle iSmarker

        //BEGIN computation of the Cauchy Stress
        std::vector < std::vector < double > > FpOld;
        FpOld = particlesSolid[iSmarker]->GetDeformationGradient(); //extraction of the deformation gradient

        adept::adouble FpNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
        adept::adouble F[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
        adept::adouble B[3][3];
        adept::adouble Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};
        adept::adouble Cauchy[3][3];

        for (unsigned j = 0; j < dim; j++) {
          for (unsigned k = 0; k < dim; k++) {
            FpNew[j][k] += gradSolDpHat[j][k];
          }
        }

        for (unsigned i = 0; i < dim; i++) {
          for (unsigned j = 0; j < dim; j++) {
            for (unsigned k = 0; k < dim; k++) {
              F[i][j] += FpNew[i][k] * FpOld[k][j];
            }
          }
        }

        if (dim == 2) F[2][2] = 1.;

        adept::adouble J_hat =  F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
                                - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];

        for (unsigned i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            B[i][j] = 0.;
            for (unsigned k = 0; k < 3; k++) {
              //left Cauchy-Green deformation tensor or Finger tensor (B = F*F^T)
              B[i][j] += F[i][k] * F[j][k];
            }
          }
        }

        adept::adouble I1_B = B[0][0] + B[1][1] + B[2][2];

        for (unsigned j = 0; j < 3; j++) {
          for (unsigned k = 0; k < 3; k++) {
            Cauchy[j][k] = lambdaMpm * log (J_hat) / J_hat * Id2th[j][k] + muMpm / J_hat * (B[j][k] - Id2th[j][k]); //alternative formulation
          }
        }
        //END computation of the Cauchy Stress

        //BEGIN redidual Solid Momentum in moving domain
        for (unsigned i = 0; i < nDofs; i++) {

          adept::adouble CauchyDIR[3] = {0., 0., 0.};

          for (unsigned j = 0.; j < dim; j++) {
            for (unsigned k = 0.; k < dim; k++) {
              CauchyDIR[j] += gradPhi[i * dim + k] * Cauchy[j][k];
            }
          }

          if (solidFlag1[i]) { // This is for diagonal dominance
            for (unsigned k = 0; k < dim; k++) {
              aRhsD[k][i] += (phi[i] * gravity[k] - J_hat * CauchyDIR[k] / rhoMpm
                              - phi[i] * (1. / (beta * dt * dt) * SolDp[k] - 1. / (beta * dt) * SolVpOld[k] - (1. - 2.* beta) / (2. * beta) * SolApOld[k])
                             ) * mass;
            }
          }
          else { // This is for the coupling with the fluid
            for (unsigned k = 0; k < dim; k++) {
              aRhsV[k][i] += (phi[i] * gravity[k] - J_hat * CauchyDIR[k] / rhoMpm
                              - phi[i] * (1. / (beta * dt * dt) * SolDp[k] - 1. / (beta * dt) * SolVpOld[k] - (1. - 2.* beta) / (2. * beta) * SolApOld[k])
                             ) * mass;
            }
          }
        }
        iSmarker++;
      }
    }
    //END SOLID PARTICLE

    //BEGIN FLUID PARTICLE
    if (MPMmaterial > 0 && MPMmaterial < nDofs) { //solid markers
      //unsigned iFmarker = markerOffsetFluid[iproc];
      //  while (iFmarker < markerOffsetFluid[iproc + 1] && iel != particlesFluid[iFmarker]->GetMarkerElement()) {
      //  iFmarker++;
      //}
      while (iFmarker < markerOffsetFluid[iproc + 1] && iel > particlesFluid[iFmarker]->GetMarkerElement()) {
        iFmarker++;
      }
      while (iFmarker < markerOffsetFluid[iproc + 1] && iel == particlesFluid[iFmarker]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particlesFluid[iFmarker]->GetMarkerLocalCoordinates();
        double mass = particlesFluid[iFmarker]->GetMarkerMass();

        msh->_finiteElement[ielt][solType]->Jacobian (vx, xi, weight, phi, gradPhi);
        msh->_finiteElement[ielt][solType]->Jacobian (vxHat, xi, weightHat, phiHat, gradPhiHat);

        //BEGIN evaluates SolDp at the particle iFmarker
        vector<adept::adouble> solDp (dim, 0.);
//         vector<double> solDpOld (dim, 0.);
//         vector<adept::adouble> solVp (dim, 0.);
        vector<double> solVpOld (dim, 0.);


//         vector<vector < adept::adouble > > gradSolVp (dim);
//         vector<vector < adept::adouble > > gradSolDp (dim);
        vector<vector < adept::adouble > > gradSolDpHat (dim);

        for (int j = 0; j < dim; j++) {
//           gradSolVp[j].assign (dim, 0.);
//           gradSolDp[j].assign (dim, 0.);
          gradSolDpHat[j].assign (dim, 0.);
        }

        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nDofs; i++) {
            solDp[j] += phi[i] * solD[j][i];
//             solDpOld[j] += phi[i] * solDOld[j][i];
//             solVp[j] += phi[i] * solV[j][i];
            solVpOld[j] += phi[i] * solVOld[j][i];
            for (int k = 0; k < dim; k++) {
//               gradSolVp[j][k] +=  gradPhi[i * dim + k] * solV[j][i];
//               gradSolDp[j][k] +=  gradPhi[i * dim + k] * solD[j][i];
              gradSolDpHat[j][k] += gradPhiHat[i * dim + k] * solD[j][i];
            }
          }
        }

//         adept::adouble solPp = 0.;
//         msh->_finiteElement[ielt][solTypeP]->GetPhi (phiP, xi);
//         for (unsigned i = 0; i < nDofsP; i++) {
//           solPp += phiP[i] * solP[i];
//         }

//         adept::adouble divV = 0.;
//         for (unsigned k = 0; k < dim; k++) {
//           divV +=  gradSolVp[k][k];
//         }

        adept::adouble F[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
        adept::adouble B[3][3];
        adept::adouble Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};
        adept::adouble Cauchy[3][3];

        for (unsigned j = 0; j < dim; j++) {
          for (unsigned k = 0; k < dim; k++) {
            F[j][k] += gradSolDpHat[j][k];
          }
        }
        adept::adouble J_hat =  F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
                                - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];

        for (unsigned i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            B[i][j] = 0.;
            for (unsigned k = 0; k < 3; k++) {
              //left Cauchy-Green deformation tensor or Finger tensor (B = F*F^T)
              B[i][j] += F[i][k] * F[j][k];
            }
          }
        }

        adept::adouble I1_B = B[0][0] + B[1][1] + B[2][2];

        for (unsigned j = 0; j < 3; j++) {
          for (unsigned k = 0; k < 3; k++) {
            Cauchy[j][k] = 100. * lambdaMpm * log (J_hat) / J_hat * Id2th[j][k] + muMpm / J_hat * (B[j][k] - Id2th[j][k]); //alternative formulation
          }
        }
        //END computation of the Cauchy Stress


        for (unsigned i = 0; i < nDofs; i++) {

          adept::adouble CauchyDIR[3] = {0., 0., 0.};
          for (unsigned j = 0.; j < dim; j++) {
            for (unsigned k = 0.; k < dim; k++) {
              CauchyDIR[j] += gradPhi[i * dim + k] * Cauchy[j][k];
            }
          }

          for (unsigned k = 0; k < dim; k++) {

//             adept::adouble Vlaplace = 0.;
//             adept::adouble Dlaplace = 0.;
//             adept::adouble advection = 0.;
//             for (unsigned j = 0; j < dim; j++) {
//               Vlaplace  +=  gradPhi[i * dim + j] * (gradSolVp[k][j] + gradSolVp[j][k]);
//               Dlaplace  +=  gradPhi[i * dim + j] * (gradSolDp[k][j] + gradSolDp[j][k]);
//               advection  +=  phi[i] * (solVp[j] - (solDp[j] - solDpOld[j]) / dt) * gradSolVp[k][j];
//             }
            if (!solidFlag1[i]) { // This is for diagonal dominance
              aRhsV[k][i] += (- phi[i] * (solDp[k] / dt - solVpOld[k]) / dt /*- advection +
                              muFluid / rhoFluid * (- Vlaplace + gradPhi[i * dim + k] * (2. / 3.) * divV)*/
                              //- 0.01 * muMpm / rhoFluid * Dlaplace
                              /*+gradPhi[i * dim + k] * solPp / rhoFluid*/
                              - J_hat * 0.01 * CauchyDIR[k] / rhoFluid) * mass;
            }
            else { // This is for the coupling with the solid
              aRhsD[k][i] += (- phi[i] * (solDp[k] / dt - solVpOld[k]) / dt /* - advection +
                              muFluid / rhoFluid * (- Vlaplace + gradPhi[i * dim + k] * (2. / 3.) * divV)*/
                              // - 0.01 * muMpm / rhoFluid * Dlaplace
                              /*+gradPhi[i * dim + k] * solPp / rhoFluid*/
                              - J_hat * 0.01 * CauchyDIR[k] / rhoFluid) * mass;
            }
          }
        }
//         if (!test) {
//           for (unsigned i = 0; i < nDofsP; i++) {
//             aRhsP[i] += phiP[i] * divV * mass / rhoFluid;
//           }
//         }

//
        iFmarker++;
      }
    }
    //END FLUID PARTICLE


    //BEGIN INTERFACE PARTICLE
//     if (MPMmaterial > 0 && MPMmaterial < nDofs) { //solid markers
//       unsigned imarker = markerOffsetInterface[iproc];
//       while (imarker < markerOffsetInterface[iproc + 1] && iel != particlesInterface[imarker]->GetMarkerElement()) {
//         imarker++;
//       }
//
//       while (imarker < markerOffsetInterface[iproc + 1] && iel == particlesInterface[imarker]->GetMarkerElement()) {
//
//         // the local coordinates of the particles are the Gauss points in this context
//         std::vector <double> xi = particlesFluid[imarker]->GetMarkerLocalCoordinates();
//
//         msh->_finiteElement[ielt][solType]->Jacobian (vxHat, xi, weightHat, phiHat, gradPhiHat);
//         msh->_finiteElement[ielt][solType]->Jacobian (vx, xi, weight, phi, gradPhi);
//         //BEGIN evaluates SolDp at the particle imarker
//
//
//
//         vector<vector < adept::adouble > > gradSolVp (dim);
//         vector<vector < adept::adouble > > gradSolDpHat (dim);
//
//         for (int j = 0; j < dim; j++) {
//           gradSolVp[j].assign (dim, 0.);
//           gradSolDpHat[j].assign (dim, 0.);
//         }
//
//         for (int j = 0; j < dim; j++) {
//           for (unsigned i = 0; i < nDofs; i++) {
//             for (int k = 0; k < dim; k++) {
//               gradSolVp[j][k] +=  gradPhi[i * dim + k] * solV[j][i];
//               gradSolDpHat[j][k] +=  gradPhiHat[i * dim + k] * solD[j][i];
//             }
//           }
//         }
//
//         adept::adouble solPp = 0.;
//         msh->_finiteElement[ielt][solTypeP]->GetPhi (phiP, xi);
//         for (unsigned i = 0; i < nDofsP; i++) {
//           solPp += phiP[i] * solP[i];
//         }
//
//         adept::adouble divV = 0.;
//         for (unsigned k = 0; k < dim; k++) {
//           divV +=  gradSolVp[k][k];
//         }
//
//         std::vector < std::vector < double > > FpOld;
//         FpOld = particlesSolid[imarker]->GetDeformationGradient(); //extraction of the deformation gradient
//
//         adept::adouble FpNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
//         adept::adouble F[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
//         adept::adouble B[3][3];
//         adept::adouble Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};
//         adept::adouble Cauchy[3][3];
//
//         for (unsigned j = 0; j < dim; j++) {
//           for (unsigned k = 0; k < dim; k++) {
//             FpNew[j][k] += gradSolDpHat[j][k];
//           }
//         }
//
//         for (unsigned i = 0; i < dim; i++) {
//           for (unsigned j = 0; j < dim; j++) {
//             for (unsigned k = 0; k < dim; k++) {
//               F[i][j] += FpNew[i][k] * FpOld[k][j];
//             }
//           }
//         }
//         if (dim == 2) F[2][2] = 1.;
//
//
//         std::vector <std::vector < double > > tangentHat;
//         particlesInterface[imarker]->GetMarkerTangent (tangentHat);
//
//         std::vector <std::vector < adept::adouble > > tangent;
//         tangent.resize (tangentHat.size());
//
//         for (unsigned k = 0; k < tangent.size(); k++) {
//           tangent[k].assign (dim, 0.);
//           for (unsigned i = 0; i < dim; i++) {
//             for (unsigned j = 0; j < dim; j++) {
//               tangent[k][i] += F[i][j] * tangentHat[k][j];
//             }
//           }
//         }
//
//         std::vector < adept::adouble > normal (dim);
//         if (dim == 2) {
//           normal[0] =  tangent[0][1];
//           normal[1] = -tangent[0][0];
//         }
//         else {
//           normal[0] = tangent[0][1] * tangent[1][2] - tangent[0][2] * tangent[1][1];
//           normal[1] = tangent[0][2] * tangent[1][0] - tangent[0][0] * tangent[1][2];
//           normal[2] = tangent[0][0] * tangent[1][1] - tangent[0][1] * tangent[1][0];
//         }
//
//         adept::adouble J_hat =  F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
//                                 - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];
//
//         for (unsigned i = 0; i < 3; i++) {
//           for (int j = 0; j < 3; j++) {
//             B[i][j] = 0.;
//             for (unsigned k = 0; k < 3; k++) {
//               //left Cauchy-Green deformation tensor or Finger tensor (B = F*F^T)
//               B[i][j] += F[i][k] * F[j][k];
//             }
//           }
//         }
//
//         adept::adouble I1_B = B[0][0] + B[1][1] + B[2][2];
//
//         for (unsigned j = 0; j < 3; j++) {
//           for (unsigned k = 0; k < 3; k++) {
//             Cauchy[j][k] = lambdaMpm * log (J_hat) / J_hat * Id2th[j][k] + muMpm / J_hat * (B[j][k] - Id2th[j][k]); //alternative formulation
//           }
//         }
//
//         adept::adouble sigmaF [3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
//         for (unsigned k = 0; k < dim; k++) {
//           sigmaF[k][k] = -solPp;
//           for (unsigned j = 0; j < dim; j++) {
//             sigmaF[k][j] += muFluid * (gradSolVp[k][j] + gradSolVp[j][k]);
//           }
//         }
//
//         adept::adouble tauSmF[3] = {0., 0., 0.};
//         for (unsigned k = 0; k < dim; k++) {
//           for (unsigned j = 0; j < dim; j++) {
//             tauSmF[k] += (Cauchy[k][j] - sigmaF[k][j]) * normal[j];
//           }
//         }
//
//         for (unsigned i = 0; i < nDofs; i++) {
//           if (!solidFlag[i]) {
//             for (unsigned k = 0; k < dim; k++) {
//
//               //aRhsD[k][i] += tauSmF[k] * phi[i];
//             }
//           }
//         }
//         imarker++;
//       }
//     }
    //END INTERFACE PARTICLES

    //copy the value of the adept::adoube aRes in double Res and store them in RES
    rhs.resize (nDofsAll); //resize

    for (int i = 0; i < nDofs; i++) {
      for (unsigned  k = 0; k < dim; k++) {
        rhs[ i +  k * nDofs ] = -aRhsD[k][i].value();
        rhs[ i + (k + dim) * nDofs ] = -aRhsV[k][i].value();
      }
    }
    for (int i = 0; i < nDofsP; i++) {
      rhs[ i + (2 * dim) * nDofs ] = -aRhsP[i].value();
    }

    myRES->add_vector_blocked (rhs, sysDofsAll);


    Jac.resize (nDofsAll * nDofsAll);
    // define the dependent variables

    for (unsigned  k = 0; k < dim; k++) {
      s.dependent (&aRhsD[k][0], nDofs);
    }
    for (unsigned  k = 0; k < dim; k++) {
      s.dependent (&aRhsV[k][0], nDofs);
    }
    s.dependent (&aRhsP[0], nDofsP);

    // define the independent variables
    for (unsigned  k = 0; k < dim; k++) {
      s.independent (&solD[k][0], nDofs);
    }
    for (unsigned  k = 0; k < dim; k++) {
      s.independent (&solV[k][0], nDofs);
    }
    s.independent (&solP[0], nDofsP);

    // get the and store jacobian matrix (row-major)
    s.jacobian (&Jac[0], true);
    myKK->add_matrix_blocked (Jac, sysDofsAll, sysDofsAll);

    s.clear_independents();
    s.clear_dependents();

  }
  //END building "soft" stiffness matrix

  myRES->close();
  myKK->close();


  //   PetscViewer    viewer1;
  //   PetscViewerDrawOpen (PETSC_COMM_WORLD, NULL, NULL, 0, 0, 1800, 1800, &viewer1);
  //   PetscObjectSetName ( (PetscObject) viewer1, "FSI matrix");
  //   PetscViewerPushFormat (viewer1, PETSC_VIEWER_DRAW_LG);
  //   MatView ( (static_cast< PetscMatrix* > (myKK))->mat(), viewer1);
  //
  //   double a;
  //   std::cin >> a;

  // *************************************
  end_time = clock();
  AssemblyTime += (end_time - start_time);
  // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

}



void GridToParticlesProjection (MultiLevelProblem & ml_prob, Line & solidLine, Line & fluidLine, Line & interfaceLine) {

  //pointers and references

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FSI");

  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = mlSol->GetSolutionLevel (level);    // pointer to the solution (level) object

  Mesh* msh = ml_prob._ml_msh->GetLevel (level);    // pointer to the mesh (level) object
  elem* el = msh->el;   // pointer to the elem object in msh (level)

  double dt =  my_nnlin_impl_sys.GetIntervalTime();
  const unsigned dim = msh->GetDimension();

  // data
  unsigned iproc  = msh->processor_id();

  // local objects
  vector< vector < double > > solD (dim);
  vector< vector < double > > solDOld (dim);
  vector< vector < double > > gradSolDHat (dim);

  for (int k = 0; k < dim; k++) {
    gradSolDHat[k].resize (dim);
  }

  vector < double > phiHat;
  vector < double > gradPhiHat;

  vector <vector < double> > vxHat (dim); //vx is coordX in assembly of ex30

  double weightHat;

  //variable-name handling
  const char varname[9][3] = {"DX", "DY", "DZ"};
  vector <unsigned> indexSolD (dim);
  unsigned solType = mlSol->GetSolutionType (&varname[0][0]);

  for (unsigned k = 0; k < dim; k++) {
    indexSolD[k] = mlSol->GetIndex (&varname[k][0]);
  }

  //line instances



  //BEGIN loop on solid particles
  std::vector<Marker*> particles = solidLine.GetParticles();
  std::vector<unsigned> markerOffset = solidLine.GetMarkerOffset();
  unsigned ielOld = UINT_MAX;
  for (unsigned iMarker = markerOffset[iproc]; iMarker < markerOffset[iproc + 1]; iMarker++) {
    unsigned iel = particles[iMarker]->GetMarkerElement();
    if (iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nDofs;
      //update element related quantities only if we are in a different element
      if (iel != ielOld) {
        ielt = msh->GetElementType (iel);
        nDofs = msh->GetElementDofNumber (iel, solType);
        for (int i = 0; i < dim; i++) {
          solD[i].resize (nDofs);
          solDOld[i].resize (nDofs);
          vxHat[i].resize (nDofs);
        }

        for (unsigned inode = 0; inode < nDofs; inode++) {
          unsigned idof = msh->GetSolutionDof (inode, iel, solType); //local 2 global solution
          unsigned idofX = msh->GetSolutionDof (inode, iel, 2); //local 2 global solution
          for (int i = 0; i < dim; i++) {
            solDOld[i][inode] = (*mysolution->_SolOld[indexSolD[i]]) (idof);
            solD[i][inode] = (*mysolution->_Sol[indexSolD[i]]) (idof) - solDOld[i][inode];
            //moving domain
            vxHat[i][inode] = (*msh->_topology->_Sol[i]) (idofX) + solDOld[i][inode];
          }
        }

      }
      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();

      msh->_finiteElement[ielt][solType]->Jacobian (vxHat, xi, weightHat, phiHat, gradPhiHat);

      std::vector <double> particleVelOld (dim);
      particles[iMarker]->GetMarkerVelocity (particleVelOld);

      std::vector <double> particleAccOld (dim);
      particles[iMarker]->GetMarkerAcceleration (particleAccOld);

      std::vector <double> particleDisp (dim, 0.);
      //update displacement and acceleration
      for (int i = 0; i < dim; i++) {
        for (unsigned inode = 0; inode < nDofs; inode++) {
          particleDisp[i] += phiHat[inode] * solD[i][inode];
        }
      }

      particles[iMarker]->SetMarkerDisplacement (particleDisp);
      particles[iMarker]->UpdateParticleCoordinates();

      std::vector <double> particleAcc (dim);
      std::vector <double> particleVel (dim);
      for (unsigned i = 0; i < dim; i++) {
        particleAcc[i] = 1. / (beta * dt * dt) * particleDisp[i] - 1. / (beta * dt) * particleVelOld[i] - (1. - 2.* beta) / (2. * beta) * particleAccOld[i];
        particleVel[i] = particleVelOld[i] + dt * ( (1. - Gamma) * particleAccOld[i] + Gamma * particleAcc[i]);
      }

      particles[iMarker]->SetMarkerVelocity (particleVel);
      particles[iMarker]->SetMarkerAcceleration (particleAcc);

      //   update the deformation gradient
      for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
          gradSolDHat[i][j] = 0.;
          for (unsigned inode = 0; inode < nDofs; inode++) {
            gradSolDHat[i][j] +=  gradPhiHat[inode * dim + j] * solD[i][inode];
          }
        }
      }
      std::vector < std::vector < double > > FpOld;
      FpOld = particles[iMarker]->GetDeformationGradient(); //extraction of the deformation gradient
      double FpNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
      std::vector < std::vector < double > > Fp (dim);
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = 0; j < dim; j++) {
          FpNew[i][j] += gradSolDHat[i][j];
        }
      }
      for (unsigned i = 0; i < dim; i++) {
        Fp[i].resize (dim);
        for (unsigned j = 0; j < dim; j++) {
          Fp[i][j] = 0.;
          for (unsigned k = 0; k < dim; k++) {
            Fp[i][j] += FpNew[i][k] * FpOld[k][j];
          }
        }
      }
      particles[iMarker]->SetDeformationGradient (Fp);
      ielOld = iel;
    }
    else {
      break;
    }
  }
  //END loop on solid particles


  //BEGIN loop on interface particles
  particles = interfaceLine.GetParticles();
  markerOffset = interfaceLine.GetMarkerOffset();
  ielOld = UINT_MAX;
  for (unsigned iMarker = markerOffset[iproc]; iMarker < markerOffset[iproc + 1]; iMarker++) {
    unsigned iel = particles[iMarker]->GetMarkerElement();
    if (iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nDofs;
      //update element related quantities only if we are in a different element
      if (iel != ielOld) {
        ielt = msh->GetElementType (iel);
        nDofs = msh->GetElementDofNumber (iel, solType);
        for (int i = 0; i < dim; i++) {
          solD[i].resize (nDofs);
          solDOld[i].resize (nDofs);
          vxHat[i].resize (nDofs);
        }

        for (unsigned inode = 0; inode < nDofs; inode++) {
          unsigned idof = msh->GetSolutionDof (inode, iel, solType); //local 2 global solution
          unsigned idofX = msh->GetSolutionDof (inode, iel, 2); //local 2 global solution
          for (int i = 0; i < dim; i++) {
            solDOld[i][inode] = (*mysolution->_SolOld[indexSolD[i]]) (idof);
            solD[i][inode] = (*mysolution->_Sol[indexSolD[i]]) (idof) - solDOld[i][inode];
            //moving domain
            vxHat[i][inode] = (*msh->_topology->_Sol[i]) (idofX) + solDOld[i][inode];
          }
        }

      }
      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();

      msh->_finiteElement[ielt][solType]->Jacobian (vxHat, xi, weightHat, phiHat, gradPhiHat);


      std::vector <double> particleDisp (dim, 0.);
      //update displacement and acceleration
      for (int i = 0; i < dim; i++) {
        for (unsigned inode = 0; inode < nDofs; inode++) {
          particleDisp[i] += phiHat[inode] * solD[i][inode];
        }
      }

      particles[iMarker]->SetMarkerDisplacement (particleDisp);
      particles[iMarker]->UpdateParticleCoordinates();

      //   update the deformation gradient
      for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
          gradSolDHat[i][j] = 0.;
          for (unsigned inode = 0; inode < nDofs; inode++) {
            gradSolDHat[i][j] +=  gradPhiHat[inode * dim + j] * solD[i][inode];
          }
        }
      }
      std::vector < std::vector < double > > FpOld;
      FpOld = particles[iMarker]->GetDeformationGradient(); //extraction of the deformation gradient
      double FpNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
      std::vector < std::vector < double > > Fp (dim);
      for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = 0; j < dim; j++) {
          FpNew[i][j] += gradSolDHat[i][j];
        }
      }
      for (unsigned i = 0; i < dim; i++) {
        Fp[i].resize (dim);
        for (unsigned j = 0; j < dim; j++) {
          Fp[i][j] = 0.;
          for (unsigned k = 0; k < dim; k++) {
            Fp[i][j] += FpNew[i][k] * FpOld[k][j];
          }
        }
      }
      particles[iMarker]->SetDeformationGradient (Fp);
      ielOld = iel;
    }
    else {
      break;
    }
  }
  //END loop on interface particles


  //BEGIN loop on fluid particles
  particles = fluidLine.GetParticles();
  markerOffset = fluidLine.GetMarkerOffset();
  ielOld = UINT_MAX;
  for (unsigned iMarker = markerOffset[iproc]; iMarker < markerOffset[iproc + 1]; iMarker++) {
    unsigned iel = particles[iMarker]->GetMarkerElement();
    if (iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nDofs;
      //update element related quantities only if we are in a different element
      if (iel != ielOld) {
        ielt = msh->GetElementType (iel);
        nDofs = msh->GetElementDofNumber (iel, solType);
        for (int i = 0; i < dim; i++) {
          solD[i].resize (nDofs);
          solDOld[i].resize (nDofs);
          vxHat[i].resize (nDofs);
        }

        for (unsigned inode = 0; inode < nDofs; inode++) {
          unsigned idof = msh->GetSolutionDof (inode, iel, solType); //local 2 global solution
          unsigned idofX = msh->GetSolutionDof (inode, iel, 2); //local 2 global solution
          for (int i = 0; i < dim; i++) {
            solDOld[i][inode] = (*mysolution->_SolOld[indexSolD[i]]) (idof);
            solD[i][inode] = (*mysolution->_Sol[indexSolD[i]]) (idof) - solDOld[i][inode];
            //moving domain
            vxHat[i][inode] = (*msh->_topology->_Sol[i]) (idofX) + solDOld[i][inode];
          }
        }
      }
      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();

      msh->_finiteElement[ielt][solType]->Jacobian (vxHat, xi, weightHat, phiHat, gradPhiHat);

      std::vector <double> particleDisp (dim, 0.);
      //update displacement and acceleration
      for (int i = 0; i < dim; i++) {
        for (unsigned inode = 0; inode < nDofs; inode++) {
          particleDisp[i] += phiHat[inode] * solD[i][inode];
        }
      }

      particles[iMarker]->SetMarkerDisplacement (particleDisp);
      particles[iMarker]->UpdateParticleCoordinates();

      ielOld = iel;
    }
    else {
      break;
    }
  }
  //END loop on fluid particles

  clock_t project_time = clock();
  ProjectGridVelocity (*mlSol);
//   ProjectGridVelocity2 (*mlSol);
  std::cout << std::endl << " total projection time : " << std::setw (11) << std::setprecision (6) << std::fixed
            << static_cast<double> ( (clock() - project_time)) / CLOCKS_PER_SEC << " s" << std::endl;

  //BEGIN loop on elements to update grid velocity and acceleration
  for (unsigned idof = msh->_dofOffset[solType][iproc]; idof < msh->_dofOffset[solType][iproc + 1]; idof++) {
    for (int i = 0; i < dim; i++) {
      mysolution->_Sol[indexSolD[i]]->set (idof, 0.);
    }
  }

  for (int i = 0; i < dim; i++) {
    mysolution->_Sol[indexSolD[i]]->close();
  }
  //END loop on elements to update grid velocity and acceleration

  fluidLine.UpdateLineMPM();
  interfaceLine.UpdateLineMPM();
  solidLine.UpdateLineMPM();

  bool updateMat = false;
  fluidLine.GetParticlesToGridMaterial (updateMat);
  interfaceLine.GetParticlesToGridMaterial (updateMat);
  updateMat = true;
  solidLine.GetParticlesToGridMaterial (updateMat);

  GetParticlesToNodeFlag (*mlSol, solidLine, fluidLine);

}


unsigned getNumberOfLayers (const double &a, const double &fac, const bool inverse = true) {

  double fac1  = (inverse) ? fac : 1. / fac;
  double da = 1. / fac1;
  double b =  da;
  unsigned n = 1;

  while (b < a) {
    da /= fac1;
    b += da;
    n++;
    if (n >= 100) {
      std::cout << "Error: number of layer is unbounded, try with a smaller factor\n";
      abort();
    }
  }
  return n;
}

void GetParticlesToNodeFlag (MultiLevelSolution &mlSol, Line & solidLine, Line & fluidLine) {

  const unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;
  Mesh* msh = mlSol._mlMesh->GetLevel (level);
  Solution* sol  = mlSol.GetSolutionLevel (level);

  unsigned solIndexNodeFlag = sol->GetIndex ("NodeFlag");
  unsigned solIndexNodeDist = sol->GetIndex ("NodeDist");
  unsigned solType = sol->GetSolutionType (solIndexNodeFlag);

  unsigned indexSolMat = sol->GetIndex ("Mat");

  sol->_Sol[solIndexNodeFlag]->zero(); // zero will mean fluid node

  const unsigned  dim = msh->GetDimension();
  unsigned    iproc = msh->processor_id();

  for (unsigned idof = msh->_dofOffset[solType][iproc]; idof < msh->_dofOffset[solType][iproc + 1]; idof++) {
    sol->_Sol[solIndexNodeDist]->set (idof, DBL_MAX);
  }
  sol->_Sol[solIndexNodeDist]->close();

  //BEGIN loop on solid particles

  vector <vector < double> > vxHat (dim);
  vector < unsigned > idof;

  std::vector<Marker*> particlesSolid = solidLine.GetParticles();
  std::vector<unsigned> markerOffset = solidLine.GetMarkerOffset();

  unsigned ielOld = UINT_MAX;
  for (unsigned iMarker = markerOffset[iproc]; iMarker < markerOffset[iproc + 1]; iMarker++) {

    unsigned iel = particlesSolid[iMarker]->GetMarkerElement();
    if (iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nDofs;
      if (iel != ielOld) {
        ielt = msh->GetElementType (iel);
        nDofs = msh->GetElementDofNumber (iel, solType);
        for (int k = 0; k < dim; k++) {
          vxHat[k].resize (nDofs);
        }
        idof.resize (nDofs);
        for (unsigned i = 0; i < nDofs; i++) {
          idof[i] = msh->GetSolutionDof (i, iel, solType);
          sol->_Sol[solIndexNodeFlag]->set (idof[i], 1.);
          unsigned idofX = msh->GetSolutionDof (i, iel, 2); //local 2 global solution
          for (int k = 0; k < dim; k++) {
            vxHat[k][i] = (*msh->_topology->_Sol[k]) (idofX);
          }
        }
      }
      std::vector<double> particleCoords (dim);
      particleCoords = particlesSolid[iMarker]->GetIprocMarkerCoordinates();

      for (unsigned i = 0; i < nDofs; i++) {
        double currentMinDist = (*sol->_Sol[solIndexNodeDist]) (idof[i]);
        double newDist = 0.;
        for (unsigned k = 0; k < dim; k++) {
          newDist  += pow ( (vxHat[k][i] - particleCoords[k]), 2.);
        }
        newDist  = sqrt (newDist);
//         if (newDist  < 0.125 * 1.562e-06) {
//           sol->_Sol[solIndexNodeFlag]->set (idof[i], 1.);
//         }
        if (newDist  < currentMinDist) {
          sol->_Sol[solIndexNodeDist]->set (idof[i], newDist);
        }


      }
      ielOld = iel;
    }
    else {
      break;
    }
  }
  sol->_Sol[solIndexNodeDist]->closeWithMinValues();
  //sol->_Sol[solIndexNodeDist]->close();
  sol->_Sol[solIndexNodeFlag]->close();
  //END

// BEGIN loop on the fluid particles
  markerOffset = fluidLine.GetMarkerOffset();
  std::vector<Marker*> particlesFluid = fluidLine.GetParticles();
  ielOld = UINT_MAX;

  for (unsigned iMarker = markerOffset[iproc]; iMarker < markerOffset[iproc + 1]; iMarker++) {
    unsigned iel = particlesFluid[iMarker]->GetMarkerElement();
    if (iel != UINT_MAX) {
      if ( (*sol->_Sol[indexSolMat]) (iel) != 0) { //only if it is an interface element
        short unsigned ielt;
        unsigned nDofs;
        if (iel != ielOld) {
          ielt = msh->GetElementType (iel);
          nDofs = msh->GetElementDofNumber (iel, solType);
          for (int k = 0; k < dim; k++) {
            vxHat[k].resize (nDofs);
          }
          idof.resize (nDofs);
          for (unsigned i = 0; i < nDofs; i++) {
            idof[i] = msh->GetSolutionDof (i, iel, solType);
            unsigned idofX = msh->GetSolutionDof (i, iel, 2); //local 2 global solution
            for (int k = 0; k < dim; k++) {
              vxHat[k][i] = (*msh->_topology->_Sol[k]) (idofX);
            }
          }
        }
        std::vector<double> particleCoords (dim);
        particleCoords = particlesFluid[iMarker]->GetIprocMarkerCoordinates();

        for (unsigned i = 0; i < nDofs; i++) {
          double currentMinDist = (*sol->_Sol[solIndexNodeDist]) (idof[i]);
          double newDist = 0.;
          for (unsigned k = 0; k < dim; k++) {
            newDist  += pow ( (vxHat[k][i] - particleCoords[k]), 2.);
          }
          newDist  = sqrt (newDist);
          if (newDist  < currentMinDist) {
            sol->_Sol[solIndexNodeFlag]->set (idof[i], 0.);
            sol->_Sol[solIndexNodeDist]->set (idof[i], newDist);
          }
        }
      }
      ielOld = iel;
    }
    else {
      break;
    }
  }
  sol->_Sol[solIndexNodeDist]->closeWithMinValues();
  sol->_Sol[solIndexNodeFlag]->closeWithMinValues();
  //END

}

void ProjectGridVelocity (MultiLevelSolution &mlSol) {

  const unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;
  Mesh* msh = mlSol._mlMesh->GetLevel (level);
  Solution* sol  = mlSol.GetSolutionLevel (level);
  const unsigned dim = msh->GetDimension();

  unsigned iproc  = msh->processor_id();
  unsigned nprocs = msh->n_processors();

  vector< vector< double > > solV (dim);     // local solution (velocity)

  vector< bool > nodeFlag;     // local solution (velocity)
  vector< unsigned > idof;     // local solution (velocity)

  vector < double > phi;
  vector <vector < double> > vx (dim); //vx is coordX in assembly of ex30
  vector <vector < double> > xp;


  //variable-name handling
  const char varname[10][5] = {"DX", "DY", "DZ", "VX", "VY", "VZ"};


  vector <unsigned> indexSolD (dim);
  vector <unsigned> indexSolV (dim);

  for (unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = mlSol.GetIndex (&varname[ivar][0]);
    indexSolV[ivar] = mlSol.GetIndex (&varname[ivar + 3][0]);
  }
  unsigned indexNodeFlag =  mlSol.GetIndex ("NodeFlag");

  unsigned solType = mlSol.GetSolutionType (&varname[0][0]);

  sol->_Sol[indexNodeFlag]->zero();

  for (unsigned k = 0; k < dim; k++) {
    (*sol->_SolOld[indexSolV[k]]) = (*sol->_Sol[indexSolV[k]]);
    sol->_Sol[indexSolV[k]]->zero();
  }

  unsigned counter = 0;

  std::vector < std::vector < std::vector <double > > > aP (3);

  //BEGIN loop on elements
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielType = msh->GetElementType (iel);

    unsigned nDofs = msh->GetElementDofNumber (iel, solType);   // number of solution element dofs

    for (unsigned  k = 0; k < dim; k++) {
      solV[k].resize (nDofs);
      vx[k].resize (nDofs);
    }
    nodeFlag.resize (nDofs);
    idof.resize (nDofs);
    xp.resize (nDofs);
    for (unsigned i = 0; i < nDofs; i++) {
      xp[i].resize (dim);
    }

    for (unsigned i = 0; i < nDofs; i++) {
      idof[i] = msh->GetSolutionDof (i, iel, solType);
      unsigned idofX = msh->GetSolutionDof (i, iel, 2);

      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_SolOld[indexSolV[k]]) (idof[i]);
        for (unsigned  k = 0; k < dim; k++) {
          xp[i][k] = (*msh->_topology->_Sol[k]) (idofX); // coordinates of the reference configuration;
          vx[k][i] = xp[i][k] + (*sol->_Sol[indexSolD[k]]) (idof[i]); // coordinates of the deformed configuration
        }
        nodeFlag[i] = ( (*sol->_Sol[indexNodeFlag]) (idof[i]) > 0.5) ? true : false;
      }
    }

    //BEGIN old
//     bool aPIsInitialized = false;
//
//     double r;
//     std::vector <double> xc;
//     GetConvexHullSphere (vx, xc, r, 0.0001); // get the ball that circumscribe the element
//     double r2 = r * r;
//
//     std::vector < std::vector< double > > xe; // get the box that encloses the element
//     GetBoundingBox (vx, xe, 0.0001);
//
//     for (unsigned i = 0; i < nDofs; i++) { // loop on the nodes of the reference elements now considered as independent points
//       if (!nodeFlag[i]) {
//
//         double d2 = 0.;
//         for (int k = 0; k < dim; k++) {
//           d2 += (xp[i][k] - xc[k]) * (xp[i][k] - xc[k]);
//         }
//         bool insideHull = true;
//         if (d2 > r2) {
//           insideHull = false;
//         }
//         for (unsigned k = 0; k < dim; k++) {
//           if (xp[i][k] < xe[k][0] || xp[i][k] > xe[k][1]) {
//             insideHull = false;
//           }
//         }
//
//         if (insideHull) { //rough test
//
//           if (!aPIsInitialized) {
//             aPIsInitialized = true;
//             std::vector < std::vector <double> > x1 (dim);
//             for (unsigned jtype = 0; jtype < solType + 1; jtype++) {
//               ProjectNodalToPolynomialCoefficients (aP[jtype], vx, ielType, jtype) ;
//             }
//           }
//           std::vector <double> xi;
//           GetClosestPointInReferenceElement (vx, xp[i], ielType, xi);
//
//           bool inverseMapping = GetInverseMapping (solType, ielType, aP, xp[i], xi, 100);
//           if (!inverseMapping) {
//             std::cout << "InverseMapping failed at " << iel << " " << idof[i] << std::endl;
//           }
//
//           bool insideDomain = CheckIfPointIsInsideReferenceDomain (xi, ielType, 1.e-3); // fine testing
//
// //           if( idof[i] == 939) {
// //             std::cout << iel << " " << xp[i][0] << " " << xp[i][1] << " "<< xi[0]<< " " <<xi[1]<<std::endl;
// //             for(unsigned jnode = 0; jnode < nDofs; jnode++){
// //               std::cout.precision(14);
// //               std::cout << idof[jnode] << " " << vx[0][jnode] <<" " << vx[1][jnode] << std::endl;
// //               std::cout.precision(7);
// //             }
// //           }
//
// //           std::cout << std::endl << " inclusion CPU time : " << std::setw (11) << std::setprecision (6) << std::fixed
// //                     << static_cast<double> ( (clock() - inclusion_time)) / CLOCKS_PER_SEC << " s" << std::endl;
//
//           if (inverseMapping && insideDomain) {
//             sol->_Sol[indexNodeFlag]->add (idof[i], 1.);
//             msh->_finiteElement[ielType][solType]->GetPhi (phi, xi);
//             //std::cout << iel << " " << i << "  ";
//             counter++;
//             for (unsigned k = 0; k < dim; k++) {
//               double solVk = 0.;
//               for (unsigned j = 0; j < nDofs; j++)    {
//                 solVk += phi[j] * solV[k][j];
//               }
//               sol->_Sol[indexSolV[k]]->add (idof[i], solVk);
//             }
//           }
//         }
//       }
//     }
    //END old

    //BEGIN new
    bool pcElemUpdate = true;
    std::vector < std::vector < std::vector < double > > > aX;
    for (unsigned i = 0; i < nDofs; i++) {
      if (!nodeFlag[i]) {
        bool pointIsInDeformedElement = (dim == 2) ? CheckInclusion2D (sol, iel, UINT_MAX, vx, xp[i]) : CheckInclusion3D (sol, iel, UINT_MAX, vx, xp[i]);
        if (pointIsInDeformedElement) {
          sol->_Sol[indexNodeFlag]->add (idof[i], 1.);
          std::vector <double> xi;
          clock_t local_coords_time = clock();
          FindLocalCoordinates (xi, aX, pcElemUpdate, mlSol, iel, xp[i], vx);
          pcElemUpdate = false;
          msh->_finiteElement[ielType][solType]->GetPhi (phi, xi);
          counter++;
          for (unsigned k = 0; k < dim; k++) {
            double solVk = 0.;
            for (unsigned j = 0; j < nDofs; j++)    {
              solVk += phi[j] * solV[k][j];
            }
            sol->_Sol[indexSolV[k]]->add (idof[i], solVk);
          }
        }
      }
    }
    //END new
  }

  sol->_Sol[indexNodeFlag]->close();

  for (unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexSolV[k]]->close();
  }

  unsigned c0 = 0;
  for (unsigned i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {
    unsigned cnt = static_cast < unsigned > (floor ( (*sol->_Sol[indexNodeFlag]) (i) + 0.5));
    if (cnt == 0) {
      c0++;
    }
    else if (cnt > 1) {
      counter -= (cnt - 1);
      for (unsigned k = 0; k < dim; k++) {
        double velk = (*sol->_Sol[indexSolV[k]]) (i) / cnt;
        sol->_Sol[indexSolV[k]]->set (i, velk);
      }
    }
  }

  unsigned counterAll;
  MPI_Reduce (&counter, &counterAll, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  std::cout << "COUNTER = " << counterAll << " " << msh->GetTotalNumberOfDofs (solType) << std::endl;

  idof.resize (c0);
  vector< double > xp0 (c0 * dim);

  unsigned c1 = 0;
  for (unsigned i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {
    if (static_cast < unsigned > (floor ( (*sol->_Sol[indexNodeFlag]) (i) + 0.5)) == 0) {
      idof[c1] = i;
      for (unsigned k = 0; k < dim; k++) {
        xp0[c1 * dim + k] = (*msh->_topology->_Sol[k]) (i);
      }
      c1++;
      if (c1 == c0) break;
    }
  }

//   std::vector < unsigned > c2(nprocs);
//   MPI_Allgather(&c0, 1, MPI_UNSIGNED, &c2[0], 1, MPI_UNSIGNED, MPI_COMM_WORLD);
//   std::vector < unsigned > c3(nprocs + 1);
//   c3[0] = 0;
//   unsigned c2Sum = 0.;
//   for (unsigned jproc = 0; jproc < nprocs; jproc++) {
//     c2Sum += c2[jproc];
//     c3[jproc + 1] = c3[jproc] + c2[jproc];
//   }
//   if(c2Sum > 0){
//     std::vector < double > xp2(c2Sum * dim);
//
//     for (unsigned i = 0; i < c1; i++) {
//       for (unsigned k = 0; k < dim; k++) {
//         xp2[ (c3[iproc] + i) * dim + k] = xp0[i * dim + k];
//       }
//     }
//     for (unsigned jproc = 0; jproc < nprocs; jproc++) {
//       MPI_Bcast (&xp2[c3[jproc] * dim], c2[jproc] * dim, MPI_DOUBLE, jproc, PETSC_COMM_WORLD);
//     }
//   }







  vector< double > xp1;
  for (unsigned jproc = 0; jproc < nprocs; jproc++) {
    c1 = c0;
    MPI_Bcast (&c1, 1, MPI_UNSIGNED, jproc, PETSC_COMM_WORLD);
    if (c1) {
      xp1.resize (c1 * dim);
      if (iproc == jproc) {
        for (unsigned i = 0; i < c1; i++) {
          for (unsigned k = 0; k < dim; k++) {
            xp1[i * dim + k] = xp0[i * dim + k];
          }
        }
      }
      MPI_Bcast (&xp1[0], c1 * dim, MPI_DOUBLE, jproc, PETSC_COMM_WORLD);

      for (unsigned i = 0; i < c1; i++) {
        std::vector < double > xp (dim);
        for (unsigned k = 0; k < dim; k++) {
          xp[k] = xp1[i * dim + k];
        }
        Marker p (xp, 1, VOLUME, sol, solType, true, 1.);
        unsigned mproc = p.GetMarkerProc (sol);
        if (iproc == mproc) {
          unsigned jel = p.GetMarkerElement();
          short unsigned jelType = msh->GetElementType (jel);
          unsigned nDofs = msh->GetElementDofNumber (jel, solType);   // number of solution element dofs

          for (unsigned  k = 0; k < dim; k++) {
            solV[k].resize (nDofs);
            vx[k].resize (nDofs);
          }

          for (unsigned j = 0; j < nDofs; j++) {
            unsigned jdof = msh->GetSolutionDof (j, jel, solType);
            unsigned jdofX = msh->GetSolutionDof (j, jel, 2);
            for (unsigned  k = 0; k < dim; k++) {
              solV[k][j] = (*sol->_SolOld[indexSolV[k]]) (jdof); //velocity to be projected
              vx[k][j] = (*msh->_topology->_Sol[k]) (jdofX) + (*sol->_Sol[indexSolD[k]]) (jdof); // coordinates of the deformed configuration
            }
          }
          std::vector <double> xi = p.GetMarkerLocalCoordinates();
          msh->_finiteElement[jelType][solType]->GetPhi (phi, xi);

          std::vector < double > vel (dim, 0.);
          for (unsigned k = 0; k < dim; k++) {
            for (unsigned j = 0; j < nDofs; j++)    {
              vel[k] += phi[j] * solV[k][j];
            }
          }
          MPI_Send (&vel[0], dim, MPI_DOUBLE, jproc, 1, PETSC_COMM_WORLD);

        }
        if (iproc == jproc) {
          std::vector < double > vel (dim);
          MPI_Recv (&vel[0], dim, MPI_DOUBLE, mproc, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
          for (unsigned k = 0; k < dim; k++) {
            sol->_Sol[indexSolV[k]]->set (idof[i], vel[k]);
          }
          counter++;
        }
      }
    }
  }

  for (unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexSolV[k]]->close();
  }

  MPI_Reduce (&counter, &counterAll, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  std::cout << "COUNTER = " << counterAll << " " << msh->GetTotalNumberOfDofs (solType) << std::endl;

}


void ProjectGridVelocity2 (MultiLevelSolution &mlSol) {

  const unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;
  Mesh* msh = mlSol._mlMesh->GetLevel (level);
  Solution* sol  = mlSol.GetSolutionLevel (level);
  const unsigned dim = msh->GetDimension();

  unsigned iproc  = msh->processor_id();
  unsigned nprocs = msh->n_processors();

  vector< vector< double > > solV (dim);     // local solution (velocity)
  vector< unsigned > idof;

  vector < double > phi;

  //variable-name handling
  const char varname[10][5] = {"DX", "DY", "DZ", "VX", "VY", "VZ"};

  vector <unsigned> indexSolV (dim);

  for (unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolV[ivar] = mlSol.GetIndex (&varname[ivar + 3][0]);
  }

  unsigned solType = mlSol.GetSolutionType (&varname[0][0]);
  unsigned solXType = 2;

  for (unsigned k = 0; k < dim; k++) {
    (*sol->_SolOld[indexSolV[k]]) = (*sol->_Sol[indexSolV[k]]);
    sol->_Sol[indexSolV[k]]->zero();
  }

  unsigned numberOfNodes = 0;
  for (unsigned jproc = 0; jproc < nprocs; jproc++) {
    numberOfNodes += (msh->_dofOffset[solXType][jproc + 1]) - (msh->_dofOffset[solXType][jproc]);
  }

  std::vector<std::vector<double>> particleCoordAndMore (numberOfNodes);
  for (unsigned i = 0; i < numberOfNodes; i++) {
    particleCoordAndMore[i].assign (dim + 1, 0.);
  }


  for (unsigned i = msh->_dofOffset[solXType][iproc]; i < msh->_dofOffset[solXType][iproc + 1]; i++) {
    for (unsigned k = 0; k < dim; k++) {
      particleCoordAndMore[i][k] = (*msh->_topology->_Sol[k]) (i);
    }
    particleCoordAndMore[i][dim] = iproc; //proc that contains dof i
  }

//   std::cout<<"iproc = " << iproc << std::endl;
//   for (unsigned isize = 0; isize < numberOfNodes; isize++) {
//       std::cout<<"particleCoordAndMore[ " << isize << "][" << 0 << " ]=" << particleCoordAndMore[isize][0] << std::endl;
//       std::cout<<"particleCoordAndMore[ " << isize << "][" << 1 << " ]=" << particleCoordAndMore[isize][1] << std::endl;
//       std::cout<<"particleCoordAndMore[ " << isize << "][" << 2 << " ]=" << particleCoordAndMore[isize][2] << std::endl;
//       std::cout<<"-------------------------------------------------------------------------------"<<std::endl;
// }

  for (unsigned isize = 0; isize < numberOfNodes; isize++) {
    for (unsigned k = 0; k < dim; k++) {
      double value = 0.;
      MPI_Allreduce (&particleCoordAndMore[isize][k], &value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      particleCoordAndMore[isize][k] = value;
    }
    double value = 0.;
    MPI_Allreduce (&particleCoordAndMore[isize][dim], &value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    particleCoordAndMore[isize][dim] = value;
  }

//     std::cout<<"iproc = " << iproc << std::endl;
//   for (unsigned isize = 0; isize < numberOfNodes; isize++) {
//       std::cout<<"particleCoordAndMore[ " << isize << "][" << 0 << " ]=" << particleCoordAndMore[isize][0] << std::endl;
//       std::cout<<"particleCoordAndMore[ " << isize << "][" << 1 << " ]=" << particleCoordAndMore[isize][1] << std::endl;
//       std::cout<<"particleCoordAndMore[ " << isize << "][" << 2 << " ]=" << particleCoordAndMore[isize][2] << std::endl;
//       std::cout<<"-------------------------------------------------------------------------------"<<std::endl;
// }

  for (unsigned isize = 0; isize < numberOfNodes; isize++) {
    std::vector<double> markerCoords (dim, 0);
    for (unsigned k = 0; k < dim; k++) {
      markerCoords[k] = particleCoordAndMore[isize][k];
    }
    Marker p (markerCoords, 1, VOLUME, sol, solType, true, 1.);
    unsigned mproc = p.GetMarkerProc (sol);
    unsigned elem = p.GetMarkerElement();

    unsigned dofproc = static_cast < unsigned > (particleCoordAndMore[isize][dim]);

    std::vector<double>VtoAdd (dim, 0.);

    if (iproc == mproc) {
      short unsigned ielType = msh->GetElementType (elem);
      unsigned nDofs = msh->GetElementDofNumber (elem, solType);   // number of solution element dofs
      for (unsigned  k = 0; k < dim; k++) {
        solV[k].resize (nDofs);
      }
      idof.resize (nDofs);
      for (unsigned i = 0; i < nDofs; i++) {
        idof[i] = msh->GetSolutionDof (i, elem, solType);
        for (unsigned  k = 0; k < dim; k++) {
          solV[k][i] = (*sol->_SolOld[indexSolV[k]]) (idof[i]);
        }
      }
      std::vector<double> xi;
      xi = p.GetMarkerLocalCoordinates();
      msh->_finiteElement[ielType][solType]->GetPhi (phi, xi);
      for (unsigned k = 0; k < dim; k++) {
        for (unsigned i = 0; i < nDofs; i++) {
          VtoAdd[k] += phi[i] * solV[k][i];
        }
      }

      MPI_Send (&VtoAdd[0], dim, MPI_DOUBLE, dofproc, 1, PETSC_COMM_WORLD);
    }

    if (iproc == dofproc) {
      MPI_Recv (&VtoAdd[0], dim, MPI_DOUBLE, mproc, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
      for (unsigned k = 0; k < dim; k++) {
        sol->_Sol[indexSolV[k]]->add (isize, VtoAdd[k]);
      }
    }
  }
}

bool CheckInclusion2D (Solution* sol, const unsigned &elemToCheck, const unsigned &previousElem,  std::vector <std::vector <double>> & xElement,  std::vector <double> &xToCheck) {

  unsigned solXType = 2;
  const unsigned dim = sol->GetMesh()->GetDimension();

  bool markerIsInElement = false;
  bool nextElementFound = false;
  short unsigned currentElementType = sol->GetMesh()->GetElementType (elemToCheck);
  double epsilon  = 10.e-10;
  double epsilon2  = epsilon * epsilon;
  double t;

  std::vector<double> xc (dim, 0); //stores the coordinates of the face node of elemToCheck
  unsigned faceNodeLocalIndex = (currentElementType == 3) ? 8 : 6;
  for (unsigned k = 0; k < dim; k++) {
    xc[k] = xElement[k][faceNodeLocalIndex] - xToCheck[k]; // coordinates are translated so that the point to check is the new origin
  }

  if (xc[0]*xc[0] < epsilon2 && xc[1]*xc[1] < epsilon2) {
    //   std::cout << "the marker is the central face node" << std::endl;
    markerIsInElement = true; //the marker is xc
  }
  else {
    unsigned faceNodeNumber = facePointNumber[currentElementType][solXType];
    std::vector< std::vector < double > > xv (dim);  //stores the coordinates of the vertices and midpoints of the element, the first and the last are the same
    for (unsigned k = 0; k < dim; k++) {
      xv[k].reserve (faceNodeNumber);
      xv[k].resize (faceNodeNumber - 1);
    }
    for (unsigned i = 0; i < faceNodeNumber - 1; i++) {
//       unsigned inodeDof  = sol->GetMesh()->GetSolutionDof (facePoints[currentElementType][solXType][i], elemToCheck, 2);
      for (unsigned k = 0; k < dim; k++) {
        xv[k][i] = xElement[k][facePoints[currentElementType][solXType][i]] - xToCheck[k];
      }
    }

    std::vector<double> r (dim, 0);  //coordinates of the intersection point between the line of the edges and the line that connects the marker and the face node
    for (unsigned k = 0; k < dim; k++) {
      xv[k].resize (faceNodeNumber);
      xv[k][faceNodeNumber - 1] = xv[k][0];
    }

    //BEGIN look for face intersection

    // rescaling coordinates to properly handle different scales of meshes

    double length = 0.;
    double sum = 0.;

    for (unsigned i = 0; i < faceNodeNumber - 1; i++) {
      for (unsigned k = 0; k < dim; k++) {
        sum += (xv[k][i + 1] - xv[k][i]) * (xv[k][i + 1] - xv[k][i]);
      }

      length += sqrt (sum);
    }

    length /= faceNodeNumber;

    for (unsigned k = 0; k < dim; k++) {
      xc[k] /= length;
      for (unsigned i = 0; i < faceNodeNumber; i++) {
        xv[k][i] /= length;
      }
    }

    for (unsigned i = 0 ; i < faceNodeNumber - 1; i++) {

      // let's find the plane passing through the points xv[][i], xv[][i+1] and xv[][2] = xv[][i] but with z = length .
      double A = (xv[1][i + 1] - xv[1][i]);
      double B = - (xv[0][i + 1] - xv[0][i]);

      //std::cout << "A= " << A << " , " <<"B= " << B <<std::endl;

      double tBottom = (A * xc[0] + B * xc[1]) ;
      double tTop = A * xv[0][i] + B * xv[1][i];
      //std::cout << "tBottom = " << tBottom << " , " << "A= " << A << " , " <<  "B= " << B << " , " << "xv[1][" << i << "] =" << xv[1][i] << " , " <<  "tTop = " <<   tTop << std::endl;

      if (fabs (tBottom) >= epsilon || fabs (tTop) < epsilon) {
        //now let's find the coordinates of the intersection point r
        t = tTop / tBottom ;
        //std::cout << "t = " << t << std::endl;

        for (unsigned k = 0; k < dim; k++) {
          r[k] = t * xc[k];
          //std::cout << "r[" << k << "] = " << r[k] <<std::endl;
        }

        if (t < 1.) {  //if not, it means the point r is far away from the marker, and we don't want to go in that direction

          std::vector< std::vector < double > > xvr (dim);

          for (unsigned k = 0; k < dim; k++) {
            xvr[k].reserve (9);
          }

          for (unsigned k = 0; k < dim; k++) {
            xvr[k].resize (faceNodeNumber);
          }

          //now we have to determine if r is inside edge i
          for (unsigned j = 0; j < faceNodeNumber; j++) {
            for (unsigned k = 0; k < dim; k++) {
              xvr[k][j] = xv[k][j] - r[k];     //transate again the reference frame so that the origin is r
            }
          }


          if ( (xvr[0][i] * xvr[0][i]  + xvr[1][i] * xvr[1][i]) < epsilon2 ||
               (xvr[0][i + 1]*xvr[0][i + 1] + xvr[1][i + 1]*xvr[1][i + 1]) < epsilon2) {
            // std::cout << "intersection on a vertex of the edge" << std::endl;

            if (fabs (t) < epsilon || t < 0) { //this means the marker is on one of the edges

              // if(fabs(t) < epsilon) std::cout << "setting markerIsInElement = true because the marker is one of the nodes" << std::endl;

              // if(t < 0) std::cout << "setting markerIsInElement = true because r is one of the nodes" << std::endl;

              markerIsInElement = true;
              break;
            }
            else {
              unsigned nodeIndex = (solXType == 0) ? i : i / 2;
              unsigned nextElem = (sol->GetMesh()->el->GetFaceElementIndex (elemToCheck, nodeIndex) - 1);
              if (nextElem != previousElem) {
                nextElementFound = true;
              }
              break;
            }

          }

          else if (xvr[0][i]*xvr[0][i + 1] < 0 || xvr[1][i]*xvr[1][i + 1] < 0) {
            // std::cout << "intersection on an edge" << std::endl;

            if (fabs (t) < epsilon || t < 0) { //this means the marker is on one of the edges

              //  if(fabs(t) < epsilon) std::cout << "setting markerIsInElement = true because the marker is on one of the edges " << std::endl;

              // if(t < 0) std::cout << "setting markerIsInElement = true because r is on one of the edges " << std::endl;

              markerIsInElement = true;
              break;
            }
            else {
              unsigned nodeIndex = (solXType == 0) ? i : i / 2;
              unsigned nextElem = (sol->GetMesh()->el->GetFaceElementIndex (elemToCheck, nodeIndex) - 1);
              if (nextElem != previousElem) {
                nextElementFound = true;
              }
              break;
            }
          }
        } // closes the " if t < 1 "
      } // closes if
    } //closes the for on the nodes
    //END look for face intersection
  }// closes the else before the for

  if (nextElementFound == true) markerIsInElement = false;

  return markerIsInElement;

}

bool CheckInclusion3D (Solution* sol, const unsigned &elemToCheck, const unsigned &previousElem, const std::vector <std::vector <double>> & xElement, const std::vector <double> &xToCheck) {

  const unsigned dim = sol->GetMesh()->GetDimension();
  unsigned solXType = 2;

  unsigned nDofs = sol->GetMesh()->GetElementDofNumber (elemToCheck, solXType);
  unsigned nFaceDofs = (solXType == 2) ? nDofs - 1 : nDofs;
  bool markerIsInElement = false;
  bool nextElementFound = false;
  short unsigned currentElementType = sol->GetMesh()->GetElementType (elemToCheck);
  double epsilon  = 10.e-10;
  double epsilon2  = epsilon * epsilon;
  double t;

  std::vector<double> xc (dim, 0); //stores the coordinates of the central node of elemToCheck
  unsigned centralNodeLocalIndex;

  if (currentElementType == 0) centralNodeLocalIndex = 26;
  else if (currentElementType == 1) centralNodeLocalIndex = 14;
  else if (currentElementType == 2) centralNodeLocalIndex = 20;

  for (unsigned k = 0; k < dim; k++) {
    xc[k] = xElement[k][centralNodeLocalIndex] - xToCheck[k]; // coordinates are translated so that the marker is the new origin
    //std::cout << " xc[" << k << "]= " <<xc[k] <<std::endl;
  }

  if (xc[0]*xc[0] < epsilon2 && xc[1]*xc[1] < epsilon2 && xc[2]*xc[2] < epsilon2) {
    //   std::cout << "the marker is the central element node" << std::endl;
    markerIsInElement = true; //the marker is xc
  }

  else {

    for (unsigned iface = 0; iface < sol->GetMesh()->GetElementFaceNumber (elemToCheck); iface++) {

      // std::cout << "iface = " << iface << std::endl;

      for (unsigned itri = 0; itri < trianglesPerFace[currentElementType][solXType][iface]; itri ++) {

        //  std::cout << "itri = " << itri << std::endl;

        std::vector<double> xcc (dim, 0); // will store the coordinates of the center scaled

        unsigned scalarCount = 0;
        std::vector<double> r (dim, 0);  //coordinates of the intersection point between the plane of itri and the line through the element center point and the marker
        std::vector< std::vector < double > > xv (dim);  //stores the coordinates of the nodes of the triangle itri

        // fill in the coordinates of the vertices of itri
        for (unsigned k = 0; k < dim; k++) {
          xv[k].reserve (4);
        }
        for (unsigned k = 0; k < dim; k++) {
          xv[k].resize (4);
        }

        for (unsigned i = 0; i < 4; i++) {
//           unsigned itriDof  = sol->GetMesh()->GetSolutionDof (faceTriangleNodes[currentElementType][solXType][iface][itri][i], elemToCheck, 2);
          // std::cout << "itriDof = " << itriDof << std::endl;
          for (unsigned k = 0; k < dim; k++) {
            xv[k][i] = xElement[k][faceTriangleNodes[currentElementType][solXType][iface][itri][i]] - xToCheck[k]; // coordinates are translated so that the marker is the new origin
          }
        }

        // rescaling coordinates to properly handle different scales of meshes
        double length = 0.;
        double sum = 0.;
        for (unsigned i = 0; i < 3; i++) {
          for (unsigned k = 0; k < dim; k++) {
            sum += (xv[k][i + 1] - xv[k][i]) * (xv[k][i + 1] - xv[k][i]);
          }
          length += sqrt (sum);
        }

        length /= 4;

        for (unsigned k = 0; k < dim; k++) {
          xcc[k] = xc[k] / length;
          //std::cout << " xcc[" << k << "]= " <<xcc[k] <<std::endl;
          for (unsigned i = 0; i < 4; i++) {
            xv[k][i] /= length;
          }
        }


        // let's find the plane passing through the vertices of the triangle itri
        double A = - (xv[1][2] - xv[1][0]) * (xv[2][1] - xv[2][0]) + (xv[2][2] - xv[2][0]) * (xv[1][1] - xv[1][0]);
        double B = - (xv[2][2] - xv[2][0]) * (xv[0][1] - xv[0][0]) + (xv[0][2] - xv[0][0]) * (xv[2][1] - xv[2][0]);
        double C = - (xv[0][2] - xv[0][0]) * (xv[1][1] - xv[1][0]) + (xv[1][2] - xv[1][0]) * (xv[0][1] - xv[0][0]);

        //std::cout << "A= " << A << " , " <<"B= " << B << " , " << "C = " << C << " , " <<std::endl;

        double tBottom = (A * xcc[0] + B * xcc[1] + C * xcc[2]);
        double tTop = A * xv[0][0] + B * xv[1][0] + C * xv[2][0];

        //std::cout << " tTop = " << tTop <<std::endl;
        //std::cout << " tBottom = " << tBottom <<std::endl;

        if (fabs (tBottom) < epsilon && fabs (tTop) >= epsilon) {
          // std::cout << "The plane of face" << itri << "does not intersect the line" <<std::endl;
          break; // must exit the loop on itri
        }

        else { //now let's find the coordinates of the intersection point r
          t = tTop / tBottom ;
          //std::cout << "t = " << t << std::endl;

          for (unsigned k = 0; k < dim; k++) {
            r[k] = t * xcc[k];
            // std::cout << "r[" << k << "] = " << r[k] <<std::endl;
          }

          if (t < 1) {  //if not, it means the point r is far away from the marker, and we don't want to go in that direction

            for (unsigned i = 0; i < 4; i++) { //now we have to determine if r is inside itri
              for (unsigned k = 0; k < dim; k++) {
                xv[k][i] = xv[k][i] - r[k];     //translate again the reference frame so that the origin is r
              }
            }

            for (unsigned i = 0; i < 3; i++) {
              double q0 = xv[1][i] * (xv[2][i] - xv[2][i + 1]) + xv[2][i] * (xv[1][i + 1] - xv[1][i]);
              double q1 = xv[2][i] * (xv[0][i] - xv[0][i + 1]) + xv[0][i] * (xv[2][i + 1] - xv[2][i]);
              double q2 = xv[0][i] * (xv[1][i] - xv[1][i + 1]) + xv[1][i] * (xv[0][i + 1] - xv[0][i]);

              // std::cout << "q0 = " << q0 << " , " << "q1 = " << q1 << " , " << " q2 = " << q2 <<  std::endl;

              double  scalarProduct = q0 * A + q1 * B + q2 * C;

              //   std::cout << "fabs(scalarProduct) = " << fabs(scalarProduct) << std::endl;

              if (scalarProduct > epsilon) {
                //   std::cout << "r is outside triangle " << itri <<  std::endl;
                break;

              }
              else if (fabs (scalarProduct) < epsilon) {  //scalarProduct == 0

                if ( (xv[0][i] * xv[0][i]  + xv[1][i] * xv[1][i] + xv[2][i] * xv[2][i]) < epsilon2 ||
                     (xv[0][i + 1]*xv[0][i + 1] + xv[1][i + 1]*xv[1][i + 1] + xv[2][i + 1]*xv[2][i + 1]) < epsilon2) {
                  //    std::cout << "intersection on a vertex of itri" << std::endl;
                  if (fabs (t) < epsilon || t < 0) { //this means the marker is on one of the faces

                    //     if(fabs(t) < epsilon) std::cout << "setting markerIsInElement = true because the marker is one vertex of triangle " << itri << std::endl;
                    //     if(t < 0) std::cout << "setting markerIsInElement = true because r is one vertex of triangle " << itri << std::endl;

                    markerIsInElement = true;
                    break;
                  }
                  else {
                    //     std::cout << "r is in triangle " << itri << std::endl;
                    unsigned nextElem = (sol->GetMesh()->el->GetFaceElementIndex (elemToCheck, iface) - 1);
                    if (nextElem != previousElem) {
                      nextElementFound = true;
                    }
                    break;
                  }
                }
                else if (xv[0][i]*xv[0][i + 1] < 0 || xv[1][i]*xv[1][i + 1] < 0 || xv[2][i]*xv[2][i + 1] < 0) {
                  //   std::cout << "intersection on an edge of itri" << std::endl;
                  if (fabs (t) < epsilon || t < 0) { //this means the marker is on one of the faces

                    //    if(fabs(t) < epsilon) std::cout << "setting markerIsInElement = true because the marker is on one of the edges of triangle " << itri << std::endl;
                    //    if(t < 0) std::cout << "setting markerIsInElement = true because r is on one of the edges of triangle " << itri << std::endl;

                    markerIsInElement = true;
                    break;
                  }
                  else {
                    //     std::cout << "r is in triangle " << itri << std::endl;
                    unsigned nextElem = (sol->GetMesh()->el->GetFaceElementIndex (elemToCheck, iface) - 1);
                    if (nextElem != previousElem) {
                      nextElementFound = true;
                    }
                    break;
                  }
                }
              }
              else if (scalarProduct < 0) {
                //    std::cout << " scalarProduct = " << scalarProduct << std::endl;
                scalarCount++;
              }
            } // closes the for loop
          } // closes " if t < 1 "
        } // closes the "else" on tBottom = 0


        if (scalarCount == 3) {
          if (fabs (t) < epsilon || t < 0) { //this means the marker is on one of the faces

            //   if(fabs(t) < epsilon) std::cout << "setting markerIsInElement = true because the marker is on one of the edges of triangle " << itri << std::endl;
            //  if(t < 0) std::cout << "setting markerIsInElement = true because r is on one of the edges of triangle " << itri << std::endl;

            markerIsInElement = true;
            break;
          }
          else {
            //    std::cout << "r is in triangle " << itri << std::endl;
            unsigned nextElem = (sol->GetMesh()->el->GetFaceElementIndex (elemToCheck, iface) - 1);
            if (nextElem != previousElem) {
              nextElementFound = true;
            }
            break;
          }
        }
        if (markerIsInElement == true) {
          break;
        }
        if (nextElementFound == true) {
          break;
        }
      } //end for on itri

      if (markerIsInElement == true) {
        break;
      }
      if (nextElementFound == true) {
        break;
      }
    } //end for on iface
  } // end of the first else

  if (nextElementFound == true) {
    markerIsInElement = false;
  }

  return markerIsInElement;

}

void FindLocalCoordinates (std::vector<double> & xi, std::vector < std::vector < std::vector < double > > >  & aX, const bool & pcElemUpdate, MultiLevelSolution & mlSol, const unsigned &elementToCheck, const std::vector <double> &xToCheck, const std::vector<std::vector<double>> & xElement) {

  const unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;
  Mesh* msh = mlSol._mlMesh->GetLevel (level);
  Solution* sol  = mlSol.GetSolutionLevel (level);
  const unsigned dim = msh->GetDimension();
  unsigned solXType = 2;

  short unsigned elemType = sol->GetMesh()->GetElementType (elementToCheck);
  unsigned nDofs = sol->GetMesh()->GetElementDofNumber (elementToCheck, solXType);

  if (pcElemUpdate) {

    //BEGIN projection nodal to polynomial coefficients
    aX.resize (solXType + 1);
    for (unsigned j = 0; j < solXType + 1; j++) {
      ProjectNodalToPolynomialCoefficients (aX[j], xElement, elemType, j);
    }
    //END projection nodal to polynomial coefficients
  }

  GetClosestPointInReferenceElement (xElement, xToCheck, elemType, xi);

  //BEGIN Inverse mapping loop
  for (unsigned j = 0; j < solXType; j++) {

    std::vector < double > phi;
    std::vector < std::vector < double > > gradPhi;
    bool convergence = false;
    while (!convergence) {
      GetPolynomialShapeFunctionGradient (phi, gradPhi, xi, elemType, solXType);
      convergence = GetNewLocalCoordinates (xi, xToCheck, phi, gradPhi, aX[solXType]);
    }
  }
  //END Inverse mapping loop

}

