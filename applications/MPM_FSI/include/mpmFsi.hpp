#include "MultiLevelSolution.hpp"
#include "PetscMatrix.hpp"

using namespace femus;

double beta = 0.25;
double Gamma = 0.5;
//double gravity[3] = {9810, 0., 0.};
double gravity[3] = {0, 0., 0.};
double NeumannFactor = .0;
Line* solidLine;
Line* fluidLine;

void GetParticlesToNodeFlag (MultiLevelSolution &mlSol, Line & solidLine, Line & fluidLine);
void GetParticlesToNodeFlag1 (MultiLevelSolution &mlSol, Line & solidLine, Line & fluidLine);
void GetPressureNeighbor (MultiLevelSolution &mlSol, Line & solidLine, Line & fluidLine);

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
  
  vector < double > phi;
  vector < double > phiHat;
  vector < double > phiP;
  vector < adept::adouble> gradPhi;
  vector < double > gradPhiHat;
  //vector < adept::adouble> gradPhiPres;
  
  phi.reserve (maxSize);
  phiHat.reserve (maxSize);
  //phiPres.reserve (maxSize);
  
  gradPhi.reserve (maxSize * dim);
  gradPhiHat.reserve (maxSize * dim);
  //gradPhiPres.reserve (maxSize * dim);
  
  vector <vector < adept::adouble> > vx (dim); //vx is coordX in assembly of ex30
  vector <vector < double> > vxHat (dim);
  
  vector< vector< adept::adouble > > solD (dim);     // local solution (displacement)
  vector< vector< adept::adouble > > solV (dim);     // local solution (velocity)
  vector< adept::adouble > solP;     // local solution (velocity)
  vector< double > solPOld;
  
  vector< vector< double > > solDOld (dim);     // local solution (displacement)
  vector< vector< double > > solVOld (dim);
  
  vector< double > rhs;    // local redidual vector
  vector< vector< adept::adouble > > aRhsD (dim);    // local redidual vector
  vector< vector< adept::adouble > > aRhsV (dim);    // local redidual vector
  vector< adept::adouble > aRhsP;    // local redidual vector
  
  std::vector <unsigned> sysDofsAll;
  
  vector < double > Jac;
  
  adept::adouble weight;
  double weightHat;
  //adept::adouble weightPres;
  
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
  //unsigned indexSolM = mlSol->GetIndex ("M");
  unsigned indexSolM = mlSol->GetIndex ("NodeFlag");
  
  //unsigned indexNodeFlag = mlSol->GetIndex ("NodeFlag");
  
  vector < bool > solidFlag;
  
  start_time = clock();
  
  myKK->zero();
  myRES->zero();
  
  //BEGIN loop on elements (to initialize the "soft" stiffness matrix)
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    
    short unsigned ielt = msh->GetElementType (iel);
    
    double  MPMmaterial = (*mysolution->_Sol[indexSolMat]) (iel);
    
    unsigned nDofsDV = msh->GetElementDofNumber (iel, solType);   // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber (iel, solTypeP);   // number of solution element dofs
    
    unsigned nDofsAll = 2 * dim * nDofsDV + nDofsP;
    // resize local arrays
    sysDofsAll.resize (nDofsAll);
    
    solidFlag.resize (nDofsDV);
    for (unsigned  k = 0; k < dim; k++) {
      solD[k].resize (nDofsDV);
      solDOld[k].resize (nDofsDV);
      
      solV[k].resize (nDofsDV);
      solVOld[k].resize (nDofsDV);
      
      aRhsD[k].assign (nDofsDV, 0.);
      aRhsV[k].assign (nDofsDV, 0.);
      
      vx[k].resize (nDofsDV);
      vxHat[k].resize (nDofsDV);
    }
    solP.resize (nDofsP);
    solPOld.resize (nDofsP);
    aRhsP.assign (nDofsP, 0.);
    
    
    unsigned counter = 0;
    for (unsigned i = 0; i < nDofsDV; i++) {
      unsigned idof = msh->GetSolutionDof (i, iel, solType);
      
      solidFlag[i] = ( (*mysolution->_Sol[indexSolM]) (idof) > 0.5) ? true : false;
      if (solidFlag[i]) counter++;
      
      for (unsigned  k = 0; k < dim; k++) {
        solD[k][i] = (*mysolution->_Sol[indexSolD[k]]) (idof);
        solDOld[k][i] = (*mysolution->_SolOld[indexSolD[k]]) (idof);
        
        solV[k][i] = (*mysolution->_Sol[indexSolV[k]]) (idof);
        solVOld[k][i] = (*mysolution->_SolOld[indexSolV[k]]) (idof);
        
        sysDofsAll[i + k * nDofsDV] = myLinEqSolver->GetSystemDof (indexSolD[k], indexPdeD[k], i, iel);
        sysDofsAll[i + (k + dim) * nDofsDV] = myLinEqSolver->GetSystemDof (indexSolV[k], indexPdeV[k], i, iel);
      }
    }
    
    bool test = (counter >= nDofsDV - 1) ? true : false;
    
    for (unsigned i = 0; i < nDofsP; i++) {
      unsigned idof = msh->GetSolutionDof (i, iel, solTypeP);
      solP[i] = (*mysolution->_Sol[indexSolP]) (idof);
      solPOld[i] = (*mysolution->_SolOld[indexSolP]) (idof);
      sysDofsAll[i + (2 * dim) * nDofsDV] = myLinEqSolver->GetSystemDof (indexSolP, indexPdeP, i, iel);
    }
    
    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();
    
    for (unsigned i = 0; i < nDofsDV; i++) {
      unsigned idofX = msh->GetSolutionDof (i, iel, 2);
      for (unsigned  k = 0; k < dim; k++) {
        vxHat[k][i] = (*msh->_topology->_Sol[k]) (idofX);
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
      
      for (unsigned i = 0; i < nDofsDV; i++) {
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
      
      for (unsigned i = 0; i < nDofsDV; i++) {
        
        for (unsigned k = 0; k < dim; k++) {
          //adept::adouble  softStiffness  = 0.;
          adept::adouble  wlaplace1V  = 0.;
          adept::adouble  wlaplace1D  = 0.;
          for (unsigned  j = 0; j < dim; j++) {
            //softStiffness +=  muMpm * gradPhiHat[i * dim + j] * (gradSolDgHat[k][j] + gradSolDgHat[j][k]);
            wlaplace1V +=  gradPhi[i * dim + j] * (gradSolVg[k][j] + gradSolVg[j][k]);
            wlaplace1D +=  gradPhiHat[i * dim + j] * (gradSolDgHat[k][j] + 0 * gradSolDgHat[j][k]);
          }
          if (!solidFlag[i]) { //kinematic equation in the fluid nodes
            double stiffness = (MPMmaterial == 0) ? .000001 : 1000;
            aRhsD[k][i] += - stiffness * wlaplace1D * weightHat;
          }
          else { //kinematic equation in the solid nodes
            aRhsV[k][i] += phiHat[i] * (solV[k][i] - (solD[k][i] - solDOld[k][i]) / dt) * weightHat;
          }
        }
      }
      
      if (MPMmaterial == 0) { //only cells that are completely fluid
        for (unsigned i = 0; i < nDofsDV; i++) {
          for (unsigned k = 0; k < dim; k++) {
            adept::adouble wlaplace = 0.;
            adept::adouble advection = 0.;
            for (unsigned j = 0; j < dim; j++) {
              wlaplace  +=  gradPhi[i * dim + j] * (gradSolVg[k][j] + gradSolVg[j][k]);
              advection  +=  phi[i] * (solVg[j] - (solDg[j] - solDgOld[j]) / dt) * gradSolVg[k][j];
            }
            if (!solidFlag[i]) {
              aRhsV[k][i] += (- rhoFluid * (solVg[k] - solVgOld[k]) / dt - rhoFluid * advection - muFluid * wlaplace + gradPhi[i * dim + k] * solPg) * weight;
            }
            else {
              aRhsD[k][i] += (- rhoFluid * (solVg[k] - solVgOld[k]) / dt - rhoFluid * advection - muFluid * wlaplace + gradPhi[i * dim + k] * solPg) * weight;
            }
          }
        }
      }
      
      for (unsigned i = 0; i < nDofsP; i++) {
        for (unsigned  k = 0; k < dim; k++) {
          if (MPMmaterial == nDofsDV || test) {  //all cells that are completely MPM solid
            aRhsP[i] += phiP[i] * (solPg) * weight;
          }
          else if (MPMmaterial == 0) {
            aRhsP[i] += phiP[i] *  gradSolVg[k][k] * weight;
          }
          //           else {
          //             aRhsP[i] += phiP[i] * (gradSolVg[k][k] + solPg * 0.01) * weight;
          //           }
        }
      }
    } // end gauss point loop
    
    
    //copy the value of the adept::adoube aRes in double Res and store them in RES
    rhs.resize (nDofsAll); //resize
    
    for (int i = 0; i < nDofsDV; i++) {
      for (unsigned  k = 0; k < dim; k++) {
        rhs[ i +  k * nDofsDV ] = -aRhsD[k][i].value();
        rhs[ i + (k + dim) * nDofsDV ] = -aRhsV[k][i].value();
      }
    }
    for (int i = 0; i < nDofsP; i++) {
      rhs[ i + (2 * dim) * nDofsDV ] = -aRhsP[i].value();
    }
    
    myRES->add_vector_blocked (rhs, sysDofsAll);
    
    
    Jac.resize (nDofsAll * nDofsAll);
    // define the dependent variables
    
    for (unsigned  k = 0; k < dim; k++) {
      s.dependent (&aRhsD[k][0], nDofsDV);
    }
    for (unsigned  k = 0; k < dim; k++) {
      s.dependent (&aRhsV[k][0], nDofsDV);
    }
    s.dependent (&aRhsP[0], nDofsP);
    
    // define the independent variables
    for (unsigned  k = 0; k < dim; k++) {
      s.independent (&solD[k][0], nDofsDV);
    }
    for (unsigned  k = 0; k < dim; k++) {
      s.independent (&solV[k][0], nDofsDV);
    }
    s.independent (&solP[0], nDofsP);
    
    // get the and store jacobian matrix (row-major)
    s.jacobian (&Jac[0] , true);
    myKK->add_matrix_blocked (Jac, sysDofsAll, sysDofsAll);
    
    s.clear_independents();
    s.clear_dependents();
    
  }
  //END building "soft" stiffness matrix
  
  
  //BEGIN loop on solid particles (used as Gauss points)
  
  std::vector<unsigned> markerOffsetSolid = solidLine->GetMarkerOffset();
  unsigned markerOffset1 = markerOffsetSolid[iproc];
  unsigned markerOffset2 = markerOffsetSolid[iproc + 1];
  std::vector<Marker*> particlesSolid = solidLine->GetParticles();
  
  std::vector < std::vector < adept::adouble > > solD1 (dim);
  
  unsigned ielOld = UINT_MAX;
  
  std::vector < unsigned > sysDofsAllD;
  std::vector < unsigned > sysDofsAllV;
  std::vector < unsigned > sysDofsP;
  
  for (unsigned iMarker = markerOffset1; iMarker < markerOffset2; iMarker++) {
    //element of particle iMarker
    unsigned iel = particlesSolid[iMarker]->GetMarkerElement();
    if (iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nDofs;
      unsigned nDofsP;
      double  MPMmaterial;
      if (iel != ielOld) { //update element related quantities only if we are in a different element
        
        MPMmaterial = (*mysolution->_Sol[indexSolMat]) (iel);
        
        ielt = msh->GetElementType (iel);
        nDofs = msh->GetElementDofNumber (iel, solType);
        nDofsP = msh->GetElementDofNumber (iel, solTypeP);
        
        solidFlag.resize (nDofs);
        for (unsigned k = 0; k < dim; k++) {
          solD[k].resize (nDofs);
          solDOld[k].resize (nDofs);
          aRhsD[k].assign (nDofs, 0.);
          aRhsV[k].assign (nDofs, 0.);
          vx[k].resize (nDofs);
          vxHat[k].resize (nDofs);
        }
        sysDofsAllD.resize (dim * nDofs);
        sysDofsAllV.resize (dim * nDofs);
        
        solP.resize (nDofsP);
        sysDofsP.resize (nDofsP);
        aRhsP.assign (nDofsP, 0.);
        
        //BEGIN copy of the value of Sol at the dofs idof of the element iel
        for (unsigned i = 0; i < nDofs; i++) {
          
          unsigned idof = msh->GetSolutionDof (i, iel, solType); //local 2 global solution
          
          solidFlag[i] = ( (*mysolution->_Sol[indexSolM]) (idof) > 0.5) ? true : false;
          
          unsigned idofX = msh->GetSolutionDof (i, iel, 2); //local 2 global coordinates
          for (unsigned k = 0; k < dim; k++) {
            solDOld[k][i] = (*mysolution->_SolOld[indexSolD[k]]) (idof);
            solD[k][i] = (*mysolution->_Sol[indexSolD[k]]) (idof) - solDOld[k][i];
            sysDofsAllD[k * nDofs + i ] = myLinEqSolver->GetSystemDof (indexSolD[k], indexPdeD[k], i, iel);
            sysDofsAllV[k * nDofs + i ] = myLinEqSolver->GetSystemDof (indexSolV[k], indexPdeV[k], i, iel);
            
            vxHat[k][i] = (*msh->_topology->_Sol[k]) (idofX) + solDOld[k][i];
          }
        }
        
        for (unsigned i = 0; i < nDofsP; i++) {
          unsigned idof = msh->GetSolutionDof (i, iel, solTypeP);
          solP[i] = (*mysolution->_Sol[indexSolP]) (idof);
          sysDofsP[i] = myLinEqSolver->GetSystemDof (indexSolP, indexPdeP, i, iel);
        }
        
        // start a new recording of all the operations involving adept::adouble variables
        s.new_recording();
      }
      
      // the local coordinates of the particles are the Gauss points in this context
      std::vector <double> xi = particlesSolid[iMarker]->GetMarkerLocalCoordinates();
      msh->_finiteElement[ielt][solType]->Jacobian (vxHat, xi, weightHat, phiHat, gradPhiHat);
      
      std::vector <double> SolVpOld (dim);
      particlesSolid[iMarker]->GetMarkerVelocity (SolVpOld);
      
      std::vector <double> SolApOld (dim);
      particlesSolid[iMarker]->GetMarkerAcceleration (SolApOld);
      
      double mass = particlesSolid[iMarker]->GetMarkerMass();
      
      unsigned ii[9][3][3] = {
        { {0, 3, 7}, {1, 2, 5}, {4, 6, 8}},
        { {1, 0, 4}, {2, 3, 6}, {5, 7, 8}},
        { {2, 1, 5}, {3, 0, 7}, {6, 4, 8}},
        { {3, 2, 6}, {0, 1, 4}, {7, 5, 8}}
      };
      
      
      for (int k = 0; k < dim; k++) {
        solD1[k].resize (nDofs);
        for (unsigned inode = 0; inode < nDofs; inode++) {
          solD1[k][inode] = solD[k][inode];
        }
      }
      
      for (unsigned iface = 0; iface < 4; iface++) {
        int faceIndex = el->GetBoundaryIndex (iel, iface);
        unsigned im = ii[iface][2][0];
        bool switchToNeumannBC = (faceIndex == 10  &&  SolVpOld[1] > 0) ? true : false;
        
        if (switchToNeumannBC) {
          
          for (unsigned inode = 0; inode < 3; inode++) {
            for (int k = 0; k < dim; k++) {
              unsigned i0 = ii[iface][inode][0];
              unsigned i1 = ii[iface][inode][1];
              unsigned i2 = ii[iface][inode][2];
              solD1[k][i0] = NeumannFactor * (- 1. / 3. * solD[k][i1] + 4. / 3 * solD[k][ i2 ]) + (1. - NeumannFactor) * solD[k][i0];
            }
          }
        }
      }
      
      for (int k = 0; k < dim; k++) {
        for (unsigned i = 0; i < nDofs; i++) {
          vx[k][i] = vxHat[k][i] + solD1[k][i];
        }
      }
      
      msh->_finiteElement[ielt][solType]->Jacobian (vx, xi, weight, phi, gradPhi); //function to evaluate at the particles
      
      
      // displacement and velocity
      //BEGIN evaluates SolDp at the particle iMarker
      vector<adept::adouble> SolDp (dim, 0.);
      vector<vector < adept::adouble > > gradSolDpHat (dim);
      for (int k = 0; k < dim; k++) {
        gradSolDpHat[k].assign (dim, 0.);
      }
      
      for (int j = 0; j < dim; j++) {
        for (unsigned i = 0; i < nDofs; i++) {
          SolDp[j] += phi[i] * solD[j][i];
          for (int k = 0; k < dim; k++) {
            gradSolDpHat[j][k] +=  gradPhiHat[i * dim + k] * solD1[j][i];
          }
        }
      }
      //END evaluates SolDp at the particle iMarker
      
      
      
      //BEGIN computation of the Cauchy Stress
      std::vector < std::vector < double > > FpOld;
      FpOld = particlesSolid[iMarker]->GetDeformationGradient(); //extraction of the deformation gradient
      
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
          //Cauchy[j][k] = mu_MPM * (B[j][k] - I1_B * Id2th[j][k] / 3.) / pow(J_hat, 5. / 3.)
          //             + K_MPM * (J_hat - 1.) * Id2th[j][k];  //Generalized Neo-Hookean solid, in Allan-Bower's book, for rubbers with very limited compressibility and K >> mu
          
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
        
        if (solidFlag[i]) { // This is for diagonal dominance
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
      
      adept::adouble solPp = 0.;
      msh->_finiteElement[ielt][solTypeP]->GetPhi (phiP, xi);
      for (unsigned i = 0; i < nDofsP; i++) {
        solPp += phiP[i] * solP[i];
      }
      
      for (unsigned i = 0; i < nDofsP; i++) {
        //aRhsP[i] += phiP[i] * solPp * mass / rhoMpm;
      }
      //END redidual Solid Momentum in moving domain
      
      
      if (iMarker == markerOffset2 - 1 || iel != particlesSolid[iMarker + 1]->GetMarkerElement()) {
        
        rhs.resize (nDofs * dim);
        //copy adouble aRhsD into double Rhs
        for (unsigned k = 0; k < dim; k++) {
          for (unsigned i = 0; i < nDofs; i++) {
            rhs[k * nDofs + i] = -aRhsD[k][i].value();
          }
        }
        myRES->add_vector_blocked (rhs, sysDofsAllD);
        
        
        //Store equations
        for (int k = 0; k < dim; k++) {
          s.dependent (&aRhsD[k][0], nDofs);
          s.independent (&solD[k][0], nDofs);
        }
        
        Jac.resize ( (dim * nDofs) * (dim * nDofs));
        
        s.jacobian (&Jac[0], true);
        
        myKK->add_matrix_blocked (Jac, sysDofsAllD, sysDofsAllD);
        s.clear_independents();
        s.clear_dependents();
        
        //BEGIN PRESSURE BLOCK
        
        rhs.resize (nDofsP);
        for (unsigned i = 0; i < nDofsP; i++) {
          rhs[i] = -aRhsP[i].value();
        }
        myRES->add_vector_blocked (rhs, sysDofsP);
        
        s.dependent (&aRhsP[0], nDofsP);
        s.independent (&solP[0], nDofsP);
        
        Jac.resize (nDofsP * nDofsP);
        
        s.jacobian (&Jac[0], true);
        
        myKK->add_matrix_blocked (Jac, sysDofsP, sysDofsP);
        s.clear_independents();
        s.clear_dependents();
        
        //END PRESSURE BLOCK
        
        if (MPMmaterial < 9) { // add this contribution only if this is an inteface element
          rhs.resize (nDofs * dim);
          for (unsigned k = 0; k < dim; k++) {
            for (unsigned i = 0; i < nDofs; i++) {
              rhs[k * nDofs + i] = -aRhsV[k][i].value();
            }
          }
          myRES->add_vector_blocked (rhs, sysDofsAllV);
          
          
          //Store equations
          for (int k = 0; k < dim; k++) {
            s.dependent (&aRhsV[k][0], nDofs);
            s.independent (&solD[k][0], nDofs);
          }
          
          Jac.resize ( (dim * nDofs) * (dim * nDofs));
          s.jacobian (&Jac[0], true);
          
          myKK->add_matrix_blocked (Jac, sysDofsAllV, sysDofsAllD);
          s.clear_independents();
          s.clear_dependents();
        }
      }
      //END local to global assembly
      
      ielOld = iel;
    }
    else {
      break;
    }
  }
  //END loop on solid particles
  
  //BEGIN loop on fluid particles (used as Gauss points)
  
  std::vector<unsigned> markerOffsetFluid = fluidLine->GetMarkerOffset();
  markerOffset1 = markerOffsetFluid[iproc];
  markerOffset2 = markerOffsetFluid[iproc + 1];
  std::vector<Marker*> particlesFluid = fluidLine->GetParticles();
  
  ielOld = UINT_MAX;
  
  for (unsigned iMarker = markerOffset1; iMarker < markerOffset2; iMarker++) {
    //element of particle iMarker
    unsigned iel = particlesFluid[iMarker]->GetMarkerElement();
    
    short unsigned ielt;
    unsigned nDofsDV;
    unsigned nDofsP;
    unsigned nDofsAll;
    
    unsigned counter = 0;
    bool test;
    
    if (iel != UINT_MAX) {
      
      if ( (*mysolution->_Sol[indexSolMat]) (iel) > 0) {
        
        //update element related quantities only if we are in a different element
        if (iel != ielOld) {
          
          ielt = msh->GetElementType (iel);
          nDofsDV = msh->GetElementDofNumber (iel, solType);   // number of solution element dofs
          nDofsP = msh->GetElementDofNumber (iel, solTypeP);
          nDofsAll = 2 * dim * nDofsDV + nDofsP;
          
          sysDofsAll.resize (nDofsAll);
          
          solidFlag.resize (nDofsDV);
          
          for (unsigned k = 0; k < dim; k++) {
            solD[k].resize (nDofsDV);
            solDOld[k].resize (nDofsDV);
            solV[k].resize (nDofsDV);
            solVOld[k].resize (nDofsDV);
            aRhsD[k].assign (nDofsDV, 0.);
            aRhsV[k].assign (nDofsDV, 0.);
            vx[k].resize (nDofsDV);
          }
          solP.resize (nDofsP);
          aRhsP.assign (nDofsP, 0.);
          
          
          counter = 0;
          for (unsigned i = 0; i < nDofsDV; i++) {
            unsigned idof = msh->GetSolutionDof (i, iel, solType);
            
            solidFlag[i] = ( (*mysolution->_Sol[indexSolM]) (idof) > 0.5) ? true : false;
            if (solidFlag[i]) counter++;
            
            for (unsigned  k = 0; k < dim; k++) {
              solD[k][i] = (*mysolution->_Sol[indexSolD[k]]) (idof);
              solDOld[k][i] = (*mysolution->_SolOld[indexSolD[k]]) (idof);
              
              solV[k][i] = (*mysolution->_Sol[indexSolV[k]]) (idof);
              solVOld[k][i] = (*mysolution->_SolOld[indexSolV[k]]) (idof);
              
              sysDofsAll[i + k * nDofsDV] = myLinEqSolver->GetSystemDof (indexSolD[k], indexPdeD[k], i, iel);
              sysDofsAll[i + (k + dim) * nDofsDV] = myLinEqSolver->GetSystemDof (indexSolV[k], indexPdeV[k], i, iel);
            }
            
            test = (counter >= nDofsDV - 1) ? true : false;
          }
          
          for (unsigned i = 0; i < nDofsP; i++) {
            unsigned idof = msh->GetSolutionDof (i, iel, solTypeP);
            solP[i] = (*mysolution->_Sol[indexSolP]) (idof);
            sysDofsAll[i + (2 * dim) * nDofsDV] = myLinEqSolver->GetSystemDof (indexSolP, indexPdeP, i, iel);
          }
          
          s.new_recording();
          
        }
        
        for (unsigned i = 0; i < nDofsDV; i++) {
          unsigned idofX = msh->GetSolutionDof (i, iel, 2);
          for (unsigned  k = 0; k < dim; k++) {
            vx[k][i] = (*msh->_topology->_Sol[k]) (idofX) + solD[k][i];
          }
        }
        
        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particlesFluid[iMarker]->GetMarkerLocalCoordinates();
        double mass = particlesFluid[iMarker]->GetMarkerMass();
        
        msh->_finiteElement[ielt][solType]->Jacobian (vx, xi, weight, phi, gradPhi);
        
        
        //std::vector <double> SolVpOld (dim);
        //     particlesFluid[iMarker]->GetMarkerVelocity (SolVpOld);//TODO
        
        //BEGIN evaluates SolDp at the particle iMarker
        vector<adept::adouble> solDp (dim, 0.);
        vector<double> solDpOld (dim, 0.);
        vector<adept::adouble> solVp (dim, 0.);
        vector<double> solVpOld (dim, 0.);
        
        
        vector<vector < adept::adouble > > gradSolVp (dim);
        vector<vector < adept::adouble > > gradSolDp (dim);
        
        for (int j = 0; j < dim; j++) {
          gradSolVp[j].assign (dim, 0.);
          gradSolDp[j].assign (dim, 0.);
        }
        
        for (int j = 0; j < dim; j++) {
          for (unsigned i = 0; i < nDofsDV; i++) {
            solDp[j] += phi[i] * solD[j][i];
            solDpOld[j] += phi[i] * solDOld[j][i];
            solVp[j] += phi[i] * solV[j][i];
            solVpOld[j] += phi[i] * solVOld[j][i];
            for (int k = 0; k < dim; k++) {
              gradSolVp[j][k] +=  gradPhi[i * dim + k] * solV[j][i];
              gradSolDp[j][k] +=  gradPhi[i * dim + k] * solD[j][i];
            }
          }
        }
        
        adept::adouble solPp = 0.;
        msh->_finiteElement[ielt][solTypeP]->GetPhi (phiP, xi);
        for (unsigned i = 0; i < nDofsP; i++) {
          solPp += phiP[i] * solP[i];
        }
        
        adept::adouble divV = 0.;
        for (unsigned k = 0; k < dim; k++) {
          divV +=  gradSolVp[k][k];
        }
        
        for (unsigned i = 0; i < nDofsDV; i++) {
          for (unsigned k = 0; k < dim; k++) {
            adept::adouble Vlaplace = 0.;
            adept::adouble Dlaplace = 0.;
            adept::adouble advection = 0.;
            for (unsigned j = 0; j < dim; j++) {
              Vlaplace  +=  gradPhi[i * dim + j] * (gradSolVp[k][j] + gradSolVp[j][k]);
              Dlaplace  +=  gradPhi[i * dim + j] * (gradSolDp[k][j] + gradSolDp[j][k]);
              advection  +=  phi[i] * (solVp[j] - (solDp[j] - solDpOld[j]) / dt) * gradSolVp[k][j];
            }
            if (!solidFlag[i]) { // This is for diagonal dominance
              aRhsV[k][i] += (- (solVp[k] - solVpOld[k]) / dt - advection +
              muFluid / rhoFluid * (- Vlaplace + gradPhi[i * dim + k] * (2. / 3.) * divV) +
              gradPhi[i * dim + k] * solPp / rhoFluid) * mass;
            }
            else { // This is for the coupling with the solid
              aRhsD[k][i] += (- (solVp[k] - solVpOld[k]) / dt - advection +
              muFluid / rhoFluid * (- Vlaplace + gradPhi[i * dim + k] * (2. / 3.) * divV) +
              gradPhi[i * dim + k] * solPp / rhoFluid) * mass;
            }
          }
        }
        
        if (!test) {
          for (unsigned i = 0; i < nDofsP; i++) {
            aRhsP[i] += phiP[i] * divV * mass / rhoFluid;
          }
          
          //           unsigned icase;
          //           std::vector < double > xi1 (dim);
          //
          //           if (xi[0] <= 0) {
          //             if (xi[1] <= 0) {
          //               icase = 0;
          //               xi1[0] = 2.* (xi[0] + 0.5);
          //               xi1[1] = 2.* (xi[1] + 0.5);
          //             }
          //             else {
          //               icase = 3;
          //               xi1[0] = 2.* (xi[0] + 0.5);
          //               xi1[1] = 2.* (xi[1] - 0.5);
          //             }
          //           }
          //           else {
          //             if (xi[1] <= 0) {
          //               icase = 1;
          //               xi1[0] = 2.* (xi[0] - 0.5);
          //               xi1[1] = 2.* (xi[1] + 0.5);
          //             }
          //             else {
          //               icase = 2;
          //               xi1[0] = 2.* (xi[0] - 0.5);
          //               xi1[1] = 2.* (xi[1] - 0.5);
          //             }
          //           }
          //           unsigned imap[4][4] = {
          //             {0, 4, 8, 7},
          //             {4, 1, 5, 8},
          //             {8, 5, 2, 6},
          //             {7, 8, 6, 3}
          //           };
          //           std::vector < std::vector < adept::adouble > > vx1 (dim);
          //           std::vector < std::vector < adept::adouble > > solV1 (dim);
          //           for (unsigned k = 0; k < dim; k++) {
          //             vx1[k].resize (4);
          //             solV1[k].resize (4);
          //             for (unsigned j = 0; j < 4; j++) {
          //               vx1[k][j] = vx[k][imap[icase][j]];
          //               solV1[k][j] = solV[k][imap[icase][j]];
          //             }
          //           }
          //
          //           msh->_finiteElement[ielt][0]->Jacobian (vx1, xi1, weight, phi, gradPhi);
          //
          //           vector<vector < adept::adouble > > gradSolVp1 (dim);
          //           for (int j = 0; j < dim; j++) {
          //             gradSolVp1[j].assign (dim, 0.);
          //           }
          //
          //           for (int j = 0; j < dim; j++) {
          //             for (unsigned i = 0; i < 4; i++) {
          //               for (int k = 0; k < dim; k++) {
          //                 //std::cout<< gradPhi[i * dim + k] << " " << solV1[j][i]<<"/t";
          //                 gradSolVp1[j][k] +=  gradPhi[i * dim + k] * solV1[j][i];
          //               }
          //             }
          //           }
          //
          //           adept::adouble divV1 = 0.;
          //           for (unsigned k = 0; k < dim; k++) {
          //             divV1 +=  gradSolVp1[k][k];
          //           }
          
          // std::cout << divV <<" "<< divV1<< " \t ";
          
          //           for (unsigned i = 0; i < nDofsP; i++) {
          //             aRhsP[i] += phiP[i] * divV1 * mass / rhoFluid;
          //           }
          
        }
        
        
        if (iMarker == markerOffset2 - 1 || iel != particlesFluid[iMarker + 1]->GetMarkerElement()) {
          
          rhs.resize (nDofsAll); //resize
          
          for (int i = 0; i < nDofsDV; i++) {
            for (unsigned  k = 0; k < dim; k++) {
              rhs[ i +  k * nDofsDV ] = -aRhsD[k][i].value();
              rhs[ i + (k + dim) * nDofsDV ] = -aRhsV[k][i].value();
            }
          }
          for (int i = 0; i < nDofsP; i++) {
            rhs[ i + (2 * dim) * nDofsDV ] = -aRhsP[i].value();
          }
          
          myRES->add_vector_blocked (rhs, sysDofsAll);
          
          
          Jac.resize (nDofsAll * nDofsAll);
          // define the dependent variables
          
          for (unsigned  k = 0; k < dim; k++) {
            s.dependent (&aRhsD[k][0], nDofsDV);
          }
          for (unsigned  k = 0; k < dim; k++) {
            s.dependent (&aRhsV[k][0], nDofsDV);
          }
          s.dependent (&aRhsP[0], nDofsP);
          
          // define the independent variables
          for (unsigned  k = 0; k < dim; k++) {
            s.independent (&solD[k][0], nDofsDV);
          }
          for (unsigned  k = 0; k < dim; k++) {
            s.independent (&solV[k][0], nDofsDV);
          }
          s.independent (&solP[0], nDofsP);
          
          // get the and store jacobian matrix (row-major)
          s.jacobian (&Jac[0] , true);
          myKK->add_matrix_blocked (Jac, sysDofsAll, sysDofsAll);
          
          s.clear_independents();
          s.clear_dependents();
          
        }
        //END local to global assembly
      }
      ielOld = iel;
    }
    else {
      break;
    }
  }
  //END loop on particles
  
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



void GridToParticlesProjection (MultiLevelProblem & ml_prob, Line & solidLine, Line & fluidLine) {
  
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
  std::vector < std::vector < double > > solD1 (dim);
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
      
      unsigned ii[9][3][3] = {
        { {0, 3, 7}, {1, 2, 5}, {4, 6, 8}},
        { {1, 0, 4}, {2, 3, 6}, {5, 7, 8}},
        { {2, 1, 5}, {3, 0, 7}, {6, 4, 8}},
        { {3, 2, 6}, {0, 1, 4}, {7, 5, 8}}
      };
      
      
      for (int k = 0; k < dim; k++) {
        solD1[k].resize (nDofs);
        for (unsigned inode = 0; inode < nDofs; inode++) {
          solD1[k][inode] = solD[k][inode];
        }
      }
      
      for (unsigned iface = 0; iface < 4; iface++) {
        int faceIndex = el->GetBoundaryIndex (iel, iface);
        unsigned im = ii[iface][2][0];
        bool switchToNeumannBC = (faceIndex == 10  &&  particleVelOld[1] > 0) ? true : false;
        if (switchToNeumannBC) {
          for (unsigned inode = 0; inode < 3; inode++) {
            for (int k = 0; k < dim; k++) {
              unsigned i0 = ii[iface][inode][0];
              unsigned i1 = ii[iface][inode][1];
              unsigned i2 = ii[iface][inode][2];
              solD1[k][i0] = NeumannFactor * (- 1. / 3. * solD[k][i1] + 4. / 3 * solD[k][ i2 ]) + (1. - NeumannFactor) * solD[k][i0];
            }
          }
        }
      }
      
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
            gradSolDHat[i][j] +=  gradPhiHat[inode * dim + j] * solD1[i][inode];
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
  solidLine.UpdateLineMPM();
  
  bool updateMat = false;
  fluidLine.GetParticlesToGridMaterial (updateMat);
  updateMat = true;
  solidLine.GetParticlesToGridMaterial (updateMat);
  
  GetParticlesToNodeFlag1 (*mlSol, solidLine, fluidLine);
  
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
  
  unsigned solIndexNodeFlag, solTypeNodeFlag, solIndexNodeDistF, solIndexNodeDistS;
  solIndexNodeFlag = sol->GetIndex ("NodeFlag");
  solTypeNodeFlag = sol->GetSolutionType (solIndexNodeFlag);
  
  solIndexNodeDistF = sol->GetIndex ("NodeDistF");
  solIndexNodeDistS = sol->GetIndex ("NodeDistS");
  
  sol->_Sol[solIndexNodeFlag]->zero(); // zero will mean fluid node
  
  const unsigned  dim = msh->GetDimension();
  unsigned    iproc = msh->processor_id();
  
  for (unsigned idof = msh->_dofOffset[solTypeNodeFlag][iproc]; idof < msh->_dofOffset[solTypeNodeFlag][iproc + 1]; idof++) {
    
    sol->_Sol[solIndexNodeDistF]->set (idof, 1.e8);
    sol->_Sol[solIndexNodeDistS]->set (idof, 1.e8);
    
  }
  
  sol->_Sol[solIndexNodeDistS]->close();
  sol->_Sol[solIndexNodeDistF]->close();
  
  //BEGIN loop on solid particles
  std::vector<unsigned> markerOffsetSolid = solidLine.GetMarkerOffset();
  unsigned markerOffset1 = markerOffsetSolid[iproc];
  unsigned markerOffset2 = markerOffsetSolid[iproc + 1];
  std::vector<Marker*> particlesSolid = solidLine.GetParticles();
  
  unsigned ielOld = UINT_MAX;
  
  for (unsigned iMarker = markerOffset1; iMarker < markerOffset2; iMarker++) {
    //element of particle iMarker
    unsigned iel = particlesSolid[iMarker]->GetMarkerElement();
    if (iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nDofs;
      vector <vector < double> > vxHat (dim);
      
      if (iel != ielOld) {
        ielt = msh->GetElementType (iel);
        nDofs = msh->GetElementDofNumber (iel, solTypeNodeFlag);
        for (int k = 0; k < dim; k++) {
          vxHat[k].resize (nDofs);
        }
        
        for (unsigned inode = 0; inode < nDofs; inode++) {
          unsigned idofX = msh->GetSolutionDof (inode, iel, 2); //local 2 global solution
          for (int k = 0; k < dim; k++) {
            vxHat[k][inode] = (*msh->_topology->_Sol[k]) (idofX);
          }
        }
        
      }
      
      std::vector<double> particleCoords (dim);
      particleCoords = particlesSolid[iMarker]->GetIprocMarkerCoordinates();
      
      for (unsigned inode = 0; inode < nDofs; inode++) {
        unsigned idof = msh->GetSolutionDof (inode, iel, solTypeNodeFlag);
        double currentMinDistS = (*sol->_Sol[solIndexNodeDistS]) (idof);
        double newDist = 0.;
        for (unsigned k = 0; k < dim; k++) {
          newDist  += (vxHat[k][inode] - particleCoords[k]) * (vxHat[k][inode] - particleCoords[k]);
        }
        newDist  = sqrt (newDist);
        std::cout.precision (16);
        //             std::cout<<"newDist = " << newDist << " , " << "currentMinDistS = " << currentMinDistS << std::endl;
        if (newDist  < currentMinDistS) {
          sol->_Sol[solIndexNodeDistS]->set (idof, newDist);
        }
      }
      ielOld = iel;
    }
    else {
      break;
    }
  }
  //   sol->_Sol[solIndexNodeDistS]->closeWithMinValues();
  sol->_Sol[solIndexNodeDistS]->close();
  //END
  
  //BEGIN loop on the fluid particles
  std::vector<unsigned> markerOffsetFluid = fluidLine.GetMarkerOffset();
  markerOffset1 = markerOffsetFluid[iproc];
  markerOffset2 = markerOffsetFluid[iproc + 1];
  std::vector<Marker*> particlesFluid = fluidLine.GetParticles();
  
  ielOld = UINT_MAX;
  
  for (unsigned iMarker = markerOffset1; iMarker < markerOffset2; iMarker++) {
    //element of particle iMarker
    unsigned iel = particlesFluid[iMarker]->GetMarkerElement();
    if (iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nDofs;
      vector <vector < double> > vxHat (dim);
      
      if (iel != ielOld) {
        ielt = msh->GetElementType (iel);
        nDofs = msh->GetElementDofNumber (iel, solTypeNodeFlag);
        for (int k = 0; k < dim; k++) {
          vxHat[k].resize (nDofs);
        }
        
        for (unsigned inode = 0; inode < nDofs; inode++) {
          unsigned idofX = msh->GetSolutionDof (inode, iel, 2); //local 2 global solution
          for (int k = 0; k < dim; k++) {
            vxHat[k][inode] = (*msh->_topology->_Sol[k]) (idofX);
          }
        }
        
      }
      
      std::vector<double> particleCoords (dim);
      particleCoords = particlesFluid[iMarker]->GetIprocMarkerCoordinates();
      
      for (unsigned inode = 0; inode < nDofs; inode++) {
        unsigned idof = msh->GetSolutionDof (inode, iel, solTypeNodeFlag);
        double currentMinDistF = (*sol->_Sol[solIndexNodeDistF]) (idof);
        double newDist = 0.;
        for (unsigned k = 0; k < dim; k++) {
          newDist  += (vxHat[k][inode] - particleCoords[k]) * (vxHat[k][inode] - particleCoords[k]);
        }
        newDist  = sqrt (newDist);
        if (newDist  < currentMinDistF) {
          sol->_Sol[solIndexNodeDistF]->set (idof, newDist);
        }
      }
      ielOld = iel;
    }
    else {
      break;
    }
  }
  //   sol->_Sol[solIndexNodeDistF]->closeWithMinValues();
  sol->_Sol[solIndexNodeDistF]->close();
  
  //END
  
  for (unsigned idof = msh->_dofOffset[solTypeNodeFlag][iproc]; idof < msh->_dofOffset[solTypeNodeFlag][iproc + 1]; idof++) {
    
    double minDistSolid = (*sol->_Sol[solIndexNodeDistS]) (idof);
    double minDistFluid = (*sol->_Sol[solIndexNodeDistF]) (idof);
    
    //         std::cout<<"minDistFluid = " << minDistFluid << " , " << "minDistSolid = " << minDistSolid << std::endl;
    
    if (minDistSolid < minDistFluid) sol->_Sol[solIndexNodeFlag]->set (idof, 1);
    
    
  }
  
  sol->_Sol[solIndexNodeFlag]->close();
  
  
}



void GetParticlesToNodeFlag1 (MultiLevelSolution &mlSol, Line & solidLine, Line & fluidLine) {
  
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
          //sol->_Sol[solIndexNodeFlag]->set (idof[i], 1.);
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
        if (newDist  < 0.75 * 1.562e-06) {
          sol->_Sol[solIndexNodeFlag]->set (idof[i], 1.);
        }
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
  //   sol->_Sol[solIndexNodeDist]->closeWithMinValues();
  sol->_Sol[solIndexNodeDist]->close();
  sol->_Sol[solIndexNodeFlag]->close();
  //END
  
  //BEGIN loop on the fluid particles
  //   markerOffset = fluidLine.GetMarkerOffset();
  //   std::vector<Marker*> particlesFluid = fluidLine.GetParticles();
  //   ielOld = UINT_MAX;
  //
  //   for (unsigned iMarker = markerOffset[iproc]; iMarker < markerOffset[iproc + 1]; iMarker++) {
  //     unsigned iel = particlesFluid[iMarker]->GetMarkerElement();
  //     if (iel != UINT_MAX) {
  //       if ( (*sol->_Sol[indexSolMat]) (iel) != 0) { //only if it is an interface element
  //         short unsigned ielt;
  //         unsigned nDofs;
  //         if (iel != ielOld) {
  //           ielt = msh->GetElementType (iel);
  //           nDofs = msh->GetElementDofNumber (iel, solType);
  //           for (int k = 0; k < dim; k++) {
  //             vxHat[k].resize (nDofs);
  //           }
  //           idof.resize (nDofs);
  //           for (unsigned i = 0; i < nDofs; i++) {
  //             idof[i] = msh->GetSolutionDof (i, iel, solType);
  //             unsigned idofX = msh->GetSolutionDof (i, iel, 2); //local 2 global solution
  //             for (int k = 0; k < dim; k++) {
  //               vxHat[k][i] = (*msh->_topology->_Sol[k]) (idofX);
  //             }
  //           }
  //         }
  //         std::vector<double> particleCoords (dim);
  //         particleCoords = particlesFluid[iMarker]->GetIprocMarkerCoordinates();
  //
  //         for (unsigned i = 0; i < nDofs; i++) {
  //           double currentMinDist = (*sol->_Sol[solIndexNodeDist]) (idof[i]);
  //           double newDist = 0.;
  //           for (unsigned k = 0; k < dim; k++) {
  //             newDist  += pow ( (vxHat[k][i] - particleCoords[k]), 2.);
  //           }
  //           newDist  = sqrt (newDist);
  //           if (newDist  < currentMinDist) {
  //             sol->_Sol[solIndexNodeFlag]->set (idof[i], 0.);
  //             sol->_Sol[solIndexNodeDist]->set (idof[i], newDist);
  //           }
  //         }
  //       }
  //       ielOld = iel;
  //     }
  //     else {
  //       break;
  //     }
  //   }
  //   sol->_Sol[solIndexNodeDist]->closeWithMinValues();
  //   sol->_Sol[solIndexNodeFlag]->closeWithMinValues();
  //END
  
  //   //BEGIN loop on elements (to initialize the "soft" stiffness matrix)
  //   for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
  //
  //     short unsigned ielt = msh->GetElementType (iel);
  //
  //     double  MPMmaterial = (*sol->_Sol[indexSolMat]) (iel);
  //
  //     if (MPMmaterial == 0) {
  //       unsigned nDofs = msh->GetElementDofNumber (iel, solType);   // number of solution element dofs
  //
  //       for (unsigned i = 0; i < nDofs; i++) {
  //         unsigned idof = msh->GetSolutionDof (i, iel, solType);
  //
  //         sol->_Sol[solIndexNodeFlag]->set (idof, 0.);
  //       }
  //     }
  //   }
  //   sol->_Sol[solIndexNodeFlag]->close();
  
  //  GetPressureNeighbor (mlSol, solidLine, fluidLine);
  
}



// void GetPressureNeighbor (MultiLevelSolution &mlSol, Line & solidLine, Line & fluidLine) {
//
//   const unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;
//   Mesh* msh = mlSol._mlMesh->GetLevel (level);
//   Solution* sol  = mlSol.GetSolutionLevel (level);
//
//   unsigned solIndexPNs = sol->GetIndex ("PNs");
//   unsigned solType = sol->GetSolutionType (solIndexPNs);
//
//   unsigned solIndexPNe = sol->GetIndex ("PNe");
//   unsigned indexSolMat = sol->GetIndex ("Mat");
//   unsigned solTypeMat = sol->GetSolutionType (indexSolMat);
//
//   (*sol->_Sol[solIndexPNs]) = -1.;
//   (*sol->_Sol[solIndexPNe]) = -1.;
//
//   const unsigned  dim = msh->GetDimension();
//   unsigned    iproc = msh->processor_id();
//
//
//   for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
//
//     short unsigned ielt = msh->GetElementType (iel);
//
//     double  MPMmaterial = (*sol->_Sol[indexSolMat]) (iel);
//
//     if (MPMmaterial == 0) {
//       unsigned nDofs = msh->GetElementDofNumber (iel, solType);   // number of solution element dofs
//
//       for (unsigned i = 0; i < nDofs; i++) {
//         unsigned idof = msh->GetSolutionDof (i, iel, solType);
//
//         sol->_Sol[solIndexPNs]->set (idof, iel);
//       }
//     }
//   }
//   sol->_Sol[solIndexPNs]->close();
//
// }




