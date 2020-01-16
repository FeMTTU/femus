#include "MultiLevelSolution.hpp"
#include "PetscMatrix.hpp"

using namespace femus;

double beta = 0.25;
double Gamma = 0.5;
//double gravity[3] = {9810, 0., 0.};
double gravity[3] = {0, 0., 0.};
double cTreshold = 0.01;

Line* solidLine;
Line* interfaceLine;
Line* fluidLine;

void GetParticlesToNodeFlag(MultiLevelSolution &mlSol, Line & solidLine, Line & fluidLine);
void GetPressureNeighbor(MultiLevelSolution &mlSol, Line & solidLine, Line & fluidLine);

void ProjectGridVelocity(MultiLevelSolution &mlSol);


double clamp(double x, double lowerlimit, double upperlimit) {
  if (x < lowerlimit)
    x = lowerlimit;
  if (x > upperlimit)
    x = upperlimit;
  return x;
}

double smoothstep(double edge0, double edge1, double x) {
  // Scale, bias and saturate x to 0..1 range
  x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0); 
  // Evaluate polynomial
  return x * x * (3 - 2 * x);
}


void AssembleMPMSys(MultiLevelProblem& ml_prob) {

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
  Solution* mysolution = mlSol->GetSolutionLevel(level);     // pointer to the solution (level) object
  LinearEquationSolver* myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];  // pointer to the equation (level) object

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);     // pointer to the mesh (level) object
  elem* el = msh->el;   // pointer to the elem object in msh (level)
  SparseMatrix* myKK = myLinEqSolver->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* myRES =  myLinEqSolver->_RES;  // pointer to the global residual vector object in pdeSys (level)

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  const unsigned dim = msh->GetDimension();

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  // data
  unsigned iproc  = msh->processor_id();

  vector< vector< adept::adouble > > solD(dim);      // local solution (displacement)
  vector< vector< adept::adouble > > solV(dim);      // local solution (velocity)
  vector< adept::adouble > solP;     // local solution (velocity)

  vector< vector< double > > solDOld(dim);      // local solution (displacement)
  vector< vector< double > > solVOld(dim);

  vector< double > rhs;    // local redidual vector
  vector< vector< adept::adouble > > aRhsD(dim);     // local redidual vector
  vector< vector< adept::adouble > > aRhsV(dim);     // local redidual vector
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

  vector <vector < adept::adouble> > vx(dim);   //vx is coordX in assembly of ex30
  vector <vector < double> > vxHat(dim);

  adept::adouble weight;
  double weightHat;


  //reading parameters for MPM body
  double rhoMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_density();
  double EMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_young_module();
  double muMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_shear_modulus();
  double nuMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_poisson_coeff();
  double lambdaMpm = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_lambda();
  double KMmp = EMpm / (3.* (1. - 2. * nuMpm));     //bulk modulus

  //reading parameters for fluid FEM domain
  double rhoFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_density();
  double muFluid = ml_prob.parameters.get<Fluid> ("FluidFEM").get_viscosity();

  double dt =  my_nnlin_impl_sys.GetIntervalTime();

  std::cout.precision(10);

  //variable-name handling
  const char varname[10][5] = {"DX", "DY", "DZ", "VX", "VY", "VZ"};

  vector <unsigned> indexSolD(dim);
  vector <unsigned> indexSolV(dim);
  vector <unsigned> indexPdeD(dim);
  vector <unsigned> indexPdeV(dim);
  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = mlSol->GetIndex(&varname[ivar][0]);
    indexSolV[ivar] = mlSol->GetIndex(&varname[ivar + 3][0]);
    indexPdeD[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);
    indexPdeV[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar + 3][0]);
  }
  unsigned solType = mlSol->GetSolutionType(&varname[0][0]);

  unsigned indexSolP = mlSol->GetIndex("P");
  unsigned indexPdeP = my_nnlin_impl_sys.GetSolPdeIndex("P");
  unsigned solTypeP = mlSol->GetSolutionType("P");

  unsigned indexSolMat = mlSol->GetIndex("Mat");
  unsigned indexC = mlSol->GetIndex("C");
  unsigned solTypeMat = mlSol->GetSolutionType("Mat");

  unsigned indexSolM = mlSol->GetIndex("NodeFlag");
  unsigned indexSolM1 = mlSol->GetIndex("M");

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
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielt = msh->GetElementType(iel);
    double  MPMmaterial = (*mysolution->_Sol[indexSolMat])(iel);

    double  C = (*mysolution->_Sol[indexC])(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solTypeP);    // number of solution element dofs

    unsigned nDofsAll = 2 * dim * nDofs + nDofsP;
    // resize local arrays
    sysDofsAll.resize(nDofsAll);

    solidFlag.resize(nDofs);
    solidFlag1.resize(nDofs);

    for(unsigned  k = 0; k < dim; k++) {
      solD[k].resize(nDofs);
      solDOld[k].resize(nDofs);

      solV[k].resize(nDofs);
      solVOld[k].resize(nDofs);

      aRhsD[k].assign(nDofs, 0.);
      aRhsV[k].assign(nDofs, 0.);

      vx[k].resize(nDofs);
      vxHat[k].resize(nDofs);
    }
    solP.resize(nDofsP);
    aRhsP.assign(nDofsP, 0.);


    for(unsigned i = 0; i < nDofs; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solType);

      solidFlag[i] = ((*mysolution->_Sol[indexSolM])(idof) > 0.5) ? true : false;       // sharp interface: divide fluid and solid dofs
      solidFlag1[i] = ((*mysolution->_Sol[indexSolM1])(idof) > 0.5) ? true : false;       // non-sharp interface: divide the dofs of cells that have at least one solid marker with the rest

      for(unsigned  k = 0; k < dim; k++) {
        solD[k][i] = (*mysolution->_Sol[indexSolD[k]])(idof);
        solDOld[k][i] = (*mysolution->_SolOld[indexSolD[k]])(idof);

        solV[k][i] = (*mysolution->_Sol[indexSolV[k]])(idof);
        solVOld[k][i] = (*mysolution->_SolOld[indexSolV[k]])(idof);

        sysDofsAll[i + k * nDofs] = myLinEqSolver->GetSystemDof(indexSolD[k], indexPdeD[k], i, iel);
        sysDofsAll[i + (k + dim) * nDofs] = myLinEqSolver->GetSystemDof(indexSolV[k], indexPdeV[k], i, iel);
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solTypeP);
      solP[i] = (*mysolution->_Sol[indexSolP])(idof);
      sysDofsAll[i + (2 * dim) * nDofs] = myLinEqSolver->GetSystemDof(indexSolP, indexPdeP, i, iel);
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    for(unsigned i = 0; i < nDofs; i++) {
      unsigned idofX = msh->GetSolutionDof(i, iel, 2);
      for(unsigned  k = 0; k < dim; k++) {
        vxHat[k][i] = (*msh->_topology->_Sol[k])(idofX) + solDOld[k][i];
        vx[k][i] = vxHat[k][i] + solD[k][i];
      }
    }

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {

      msh->_finiteElement[ielt][solType]->Jacobian(vxHat, ig, weightHat, phiHat, gradPhiHat);
      msh->_finiteElement[ielt][solType]->Jacobian(vx, ig, weight, phi, gradPhi);

      vector < adept::adouble > solVg(dim, 0.);
      vector < adept::adouble > solVgOld(dim, 0.);

      vector < adept::adouble > solDg(dim, 0.);
      vector < adept::adouble > solDgOld(dim, 0.);

      vector < vector < adept::adouble > > gradSolDgHat(dim);
      vector < vector < adept::adouble > > gradSolVg(dim);

      for(unsigned  k = 0; k < dim; k++) {
        gradSolDgHat[k].assign(dim, 0);
        gradSolVg[k].assign(dim, 0);
      }

      for(unsigned i = 0; i < nDofs; i++) {
        for(unsigned j = 0; j < dim; j++) {
          solVg[j] += phiHat[i] * solV[j][i];
          solDg[j] += phiHat[i] * solD[j][i];
          solVgOld[j] += phiHat[i] * solVOld[j][i];
          solDgOld[j] += phiHat[i] * solDOld[j][i];
          for(unsigned  k = 0; k < dim; k++) {
            gradSolDgHat[k][j] += gradPhiHat[i * dim + j] * solD[k][i];
            gradSolVg[k][j] += gradPhi[i * dim + j] * solV[k][i];
          }
        }
      }

      double *phiP = msh->_finiteElement[ielt][solTypeP]->GetPhi(ig);
      adept::adouble solPg = 0.;
      for(unsigned i = 0; i < nDofsP; i++) {
        solPg += phiP[i] * solP[i];
      }

      for(unsigned i = 0; i < nDofs; i++) {
        for(unsigned k = 0; k < dim; k++) {
          if(solidFlag[i]) {   //kinematic: v - dD/dt = 0 in the solid domain
            aRhsV[k][i] += phiHat[i] * (solV[k][i] - (solD[k][i] - solDOld[k][i]) / dt) * weightHat;
          }
          if(!solidFlag1[i]) {   //ALE equation: smooth displacement of the fluid mesh
            adept::adouble  wlaplaceD  = 0.;
            for(unsigned  j = 0; j < dim; j++) {
              wlaplaceD +=  gradPhiHat[i * dim + j] * (gradSolDgHat[k][j] + gradSolDgHat[j][k]);
            }
            aRhsD[k][i] +=  - wlaplaceD * weightHat;
          }
        }
      }

      if(MPMmaterial < nDofs) {   //only cells that are not completely solid
        for(unsigned i = 0; i < nDofs; i++) {
          for(unsigned k = 0; k < dim; k++) {
            adept::adouble wlaplace = 0.;
            adept::adouble advection = 0.;
            for(unsigned j = 0; j < dim; j++) {
              wlaplace  +=  gradPhi[i * dim + j] * (gradSolVg[k][j] + gradSolVg[j][k]);
              advection  +=  phi[i] * (solVg[j] - (solDg[j] - solDgOld[j]) / dt) * gradSolVg[k][j];
            }
            if(!solidFlag[i]) {   // contribution to a fluid momentum equation (solidFlag[i]) from fluid and interface cells
              aRhsV[k][i] += (- rhoFluid * phi[i] * (solVg[k] - solVgOld[k]) / dt
                              - rhoFluid * advection - muFluid * wlaplace
                              + gradPhi[i * dim + k] * solPg) * weight;
            }
            if(MPMmaterial == 0 && solidFlag1[i]) {   // contribution to a solid momentum equation (solidFlag1[i]) from a adjacent fluid cell
              aRhsD[k][i] += (- rhoFluid * phi[i] * (solVg[k] - solVgOld[k]) / dt
                              - rhoFluid * advection - muFluid * wlaplace
                              + gradPhi[i * dim + k] * solPg) * weight;
            }
          }
        }
      }



      for(unsigned i = 0; i < nDofsP; i++) {
        if(MPMmaterial == 0) {

          for(unsigned  k = 0; k < dim; k++) {
            aRhsP[i] += phiP[i] *  gradSolVg[k][k] * weight;
          }
        }
        else if (MPMmaterial < nDofs) {
          double nu = 0.49; //almost incompressible fluid in interface cell. This is because of the compressiblity of the solid and it relaxes a little bit the incompressibility of the fluid
          double lameFluidInverse = (1. - 2. * nu) / (2. * muFluid * nu);
          aRhsP[i] += phiP[i] * solPg * lameFluidInverse * weight;
          for (unsigned  k = 0; k < dim; k++) {
            aRhsP[i] += phiP[i] *  gradSolVg[k][k] * weight;
          }
        }
        else if(MPMmaterial == nDofs) {   //all cells that are completely MPM solid
          aRhsP[i] += phiP[i] * solPg * weight;
        }
      }
    } // end gauss point loop


    //BEGIN SOLID PARTICLE
    if(MPMmaterial > 0) {   //solid markers
      while(iSmarker < markerOffsetSolid[iproc + 1] && iel == particlesSolid[iSmarker]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particlesSolid[iSmarker]->GetMarkerLocalCoordinates();
        msh->_finiteElement[ielt][solType]->Jacobian(vxHat, xi, weightHat, phiHat, gradPhiHat);

        std::vector <double> SolVpOld(dim);
        particlesSolid[iSmarker]->GetMarkerVelocity(SolVpOld);

        std::vector <double> SolApOld(dim);
        particlesSolid[iSmarker]->GetMarkerAcceleration(SolApOld);

        double mass = particlesSolid[iSmarker]->GetMarkerMass();

        msh->_finiteElement[ielt][solType]->Jacobian(vx, xi, weight, phi, gradPhi);   //function to evaluate at the particles

        // displacement and velocity
        //BEGIN evaluates SolDp at the particle iSmarker
        vector<adept::adouble> SolDp(dim, 0.);
        vector<vector < adept::adouble > > gradSolDpHat(dim);
        for(int k = 0; k < dim; k++) {
          gradSolDpHat[k].assign(dim, 0.);
        }

        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nDofs; i++) {
            SolDp[j] += phi[i] * solD[j][i];
            for(int k = 0; k < dim; k++) {
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

//         if(C > cTreshold) {
          for(unsigned j = 0; j < dim; j++) {
            for(unsigned k = 0; k < dim; k++) {
              FpNew[j][k] += smoothstep(0., 0.5, C) * gradSolDpHat[j][k];
            }
          }
//         }

        for(unsigned i = 0; i < dim; i++) {
          for(unsigned j = 0; j < dim; j++) {
            for(unsigned k = 0; k < dim; k++) {
              F[i][j] += FpNew[i][k] * FpOld[k][j];
            }
          }
        }

        if(dim == 2) F[2][2] = 1.;

        adept::adouble J_hat =  F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
                                - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];

        for(unsigned i = 0; i < 3; i++) {
          for(int j = 0; j < 3; j++) {
            B[i][j] = 0.;
            for(unsigned k = 0; k < 3; k++) {
              //left Cauchy-Green deformation tensor or Finger tensor (B = F*F^T)
              B[i][j] += F[i][k] * F[j][k];
            }
          }
        }

        adept::adouble I1_B = B[0][0] + B[1][1] + B[2][2];

        for(unsigned j = 0; j < 3; j++) {
          for(unsigned k = 0; k < 3; k++) {
            Cauchy[j][k] = lambdaMpm * log(J_hat) / J_hat * Id2th[j][k] + muMpm / J_hat * (B[j][k] - Id2th[j][k]);     //alternative formulation
          }
        }
        //END computation of the Cauchy Stress

        //BEGIN redidual Solid Momentum in moving domain
        for(unsigned i = 0; i < nDofs; i++) {

          adept::adouble CauchyDIR[3] = {0., 0., 0.};

          for(unsigned j = 0.; j < dim; j++) {
            for(unsigned k = 0.; k < dim; k++) {
              CauchyDIR[j] += gradPhi[i * dim + k] * Cauchy[j][k];
            }
          }
          for(unsigned k = 0; k < dim; k++) {
            if(solidFlag1[i]) {
              aRhsD[k][i] += (phi[i] * gravity[k] - J_hat * CauchyDIR[k] / rhoMpm
                              - phi[i] * (1. / (beta * dt * dt) * SolDp[k] - 1. / (beta * dt) * SolVpOld[k] - (1. - 2.* beta) / (2. * beta) * SolApOld[k])
                             ) * mass;
            }
          }
        }

//         if(MPMmaterial < nDofs) {
// 
//           vector<vector < adept::adouble > > gradSolVp(dim);
//           for(int j = 0; j < dim; j++) {
//             gradSolVp[j].assign(dim, 0.);
//             for(unsigned i = 0; i < nDofs; i++) {
//               for(int k = 0; k < dim; k++) {
//                 gradSolVp[j][k] +=  gradPhi[i * dim + k] * solV[j][i];
//               }
//             }
//           }
// 
//           adept::adouble solPp = 0.;
//           msh->_finiteElement[ielt][solTypeP]->GetPhi(phiP, xi);
//           for(unsigned i = 0; i < nDofsP; i++) {
//             solPp += phiP[i] * solP[i];
//           }
//           double nu = 0.49;
//           double lameFluidInverse = (1. - 2. * nu) / (2. * muFluid * nu);
//           for(unsigned i = 0; i < nDofsP; i++) {
//             aRhsP[i] += phiP[i] * solPp * lameFluidInverse * mass / rhoMpm;
//             for(unsigned  k = 0; k < dim; k++) {
//               aRhsP[i] += phiP[i] *  gradSolVp[k][k] * mass / rhoMpm;
//             }
//           }
//         }

        iSmarker++;
      }
    }
    //END SOLID PARTICLE

    //BEGIN FLUID PARTICLE
    if(MPMmaterial > 0 && MPMmaterial < nDofs) {   //solid markers
      while(iFmarker < markerOffsetFluid[iproc + 1] && iel > particlesFluid[iFmarker]->GetMarkerElement()) {
        iFmarker++;
      }
      while(iFmarker < markerOffsetFluid[iproc + 1] && iel == particlesFluid[iFmarker]->GetMarkerElement()) {

        // the local coordinates of the particles are the Gauss points in this context
        std::vector <double> xi = particlesFluid[iFmarker]->GetMarkerLocalCoordinates();
        double mass = particlesFluid[iFmarker]->GetMarkerMass();

        msh->_finiteElement[ielt][solType]->Jacobian(vx, xi, weight, phi, gradPhi);
        //BEGIN evaluates SolDp at the particle iFmarker
        vector<adept::adouble> solDp(dim, 0.);
        vector<double> solDpOld(dim, 0.);
        vector<adept::adouble> solVp(dim, 0.);
        vector<double> solVpOld(dim, 0.);


        vector<vector < adept::adouble > > gradSolVp(dim);
        vector<vector < adept::adouble > > gradSolDp(dim);

        for(int j = 0; j < dim; j++) {
          gradSolVp[j].assign(dim, 0.);
          gradSolDp[j].assign(dim, 0.);
        }

        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nDofs; i++) {
            solDp[j] += phi[i] * solD[j][i];
            solDpOld[j] += phi[i] * solDOld[j][i];
            solVp[j] += phi[i] * solV[j][i];
            solVpOld[j] += phi[i] * solVOld[j][i];
            for(int k = 0; k < dim; k++) {
              gradSolDp[j][k] +=  gradPhi[i * dim + k] * solD[j][i];
              gradSolVp[j][k] +=  gradPhi[i * dim + k] * solV[j][i];
            }
          }
        }

        adept::adouble solPp = 0.;
        msh->_finiteElement[ielt][solTypeP]->GetPhi(phiP, xi);
        for(unsigned i = 0; i < nDofsP; i++) {
          solPp += phiP[i] * solP[i];
        }

        for(unsigned i = 0; i < nDofs; i++) {

          for(unsigned k = 0; k < dim; k++) {
            adept::adouble Vlaplace = 0.;
            adept::adouble Dlaplace = 0.;
            adept::adouble advection = 0.;
            for(unsigned j = 0; j < dim; j++) {
              Vlaplace  +=  gradPhi[i * dim + j] * (gradSolVp[k][j] + gradSolVp[j][k]);
              Dlaplace  +=  gradPhi[i * dim + j] * (gradSolDp[k][j] + gradSolDp[j][k]);
              advection  +=  phi[i] * (solVp[j] - (solDp[j] - solDpOld[j]) / dt) * gradSolVp[k][j];
            }
            if(solidFlag1[i]) {   // This is for the coupling with the solid
              aRhsD[k][i] += (- phi[i] * (solVp[k] - solVpOld[k]) / dt - advection +
                              muFluid / rhoFluid * (- Vlaplace) + gradPhi[i * dim + k] * solPp / rhoFluid
                              + 0.001 * (1. - smoothstep(0., 0.5, C) ) * muMpm / rhoFluid * (- Dlaplace)    //softStiffness (this is the only trick in this formulation)
                             ) * mass;
            }
          }
        }

//         for(unsigned i = 0; i < nDofsP; i++) {
//           for(unsigned  k = 0; k < dim; k++) {
//             aRhsP[i] += phiP[i] *  gradSolVp[k][k] * mass / rhoFluid;
//           }
//         }

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
//         std::vector <double> xi = particlesInterface[imarker]->GetMarkerLocalCoordinates();
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
    rhs.resize(nDofsAll);   //resize

    for(int i = 0; i < nDofs; i++) {
      for(unsigned  k = 0; k < dim; k++) {
        rhs[ i +  k * nDofs ] = -aRhsD[k][i].value();
        rhs[ i + (k + dim) * nDofs ] = -aRhsV[k][i].value();
      }
    }
    for(int i = 0; i < nDofsP; i++) {
      rhs[ i + (2 * dim) * nDofs ] = -aRhsP[i].value();
    }

    myRES->add_vector_blocked(rhs, sysDofsAll);


    Jac.resize(nDofsAll * nDofsAll);
    // define the dependent variables

    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&aRhsD[k][0], nDofs);
    }
    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&aRhsV[k][0], nDofs);
    }
    s.dependent(&aRhsP[0], nDofsP);

    // define the independent variables
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solD[k][0], nDofs);
    }
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solV[k][0], nDofs);
    }
    s.independent(&solP[0], nDofsP);

    // get the and store jacobian matrix (row-major)
    s.jacobian(&Jac[0], true);
    myKK->add_matrix_blocked(Jac, sysDofsAll, sysDofsAll);

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



void GridToParticlesProjection(MultiLevelProblem & ml_prob, Line & solidLine, Line & fluidLine, Line & interfaceLine) {

  //pointers and references

  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FSI");

  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = mlSol->GetSolutionLevel(level);     // pointer to the solution (level) object

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);     // pointer to the mesh (level) object
  elem* el = msh->el;   // pointer to the elem object in msh (level)

  double dt =  my_nnlin_impl_sys.GetIntervalTime();
  const unsigned dim = msh->GetDimension();

  // data
  unsigned iproc  = msh->processor_id();

  // local objects
  vector< vector < double > > solD(dim);
  vector< vector < double > > solDOld(dim);
  vector< vector < double > > gradSolDHat(dim);

  for(int k = 0; k < dim; k++) {
    gradSolDHat[k].resize(dim);
  }

  vector < double > phiHat;
  vector < double > gradPhiHat;

  vector <vector < double> > vxHat(dim);   //vx is coordX in assembly of ex30

  double weightHat;

  //variable-name handling
  const char varname[9][3] = {"DX", "DY", "DZ"};
  vector <unsigned> indexSolD(dim);
  unsigned solType = mlSol->GetSolutionType(&varname[0][0]);

  for(unsigned k = 0; k < dim; k++) {
    indexSolD[k] = mlSol->GetIndex(&varname[k][0]);
  }


  unsigned indexC = mlSol->GetIndex("C");
  //line instances



  //BEGIN loop on solid particles
  std::vector<Marker*> particles = solidLine.GetParticles();
  std::vector<unsigned> markerOffset = solidLine.GetMarkerOffset();
  unsigned ielOld = UINT_MAX;
  for(unsigned iMarker = markerOffset[iproc]; iMarker < markerOffset[iproc + 1]; iMarker++) {
    unsigned iel = particles[iMarker]->GetMarkerElement();
    if(iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nDofs;
      double  C;
      //update element related quantities only if we are in a different element
      if(iel != ielOld) {

        C = (*mysolution->_Sol[indexC])(iel);

        ielt = msh->GetElementType(iel);
        nDofs = msh->GetElementDofNumber(iel, solType);
        for(int i = 0; i < dim; i++) {
          solD[i].resize(nDofs);
          solDOld[i].resize(nDofs);
          vxHat[i].resize(nDofs);
        }

        for(unsigned inode = 0; inode < nDofs; inode++) {
          unsigned idof = msh->GetSolutionDof(inode, iel, solType);   //local 2 global solution
          unsigned idofX = msh->GetSolutionDof(inode, iel, 2);   //local 2 global solution
          for(int i = 0; i < dim; i++) {
            solDOld[i][inode] = (*mysolution->_SolOld[indexSolD[i]])(idof);
            solD[i][inode] = (*mysolution->_Sol[indexSolD[i]])(idof) - solDOld[i][inode];
            //moving domain
            vxHat[i][inode] = (*msh->_topology->_Sol[i])(idofX) + solDOld[i][inode];
          }
        }

      }
      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();

      msh->_finiteElement[ielt][solType]->Jacobian(vxHat, xi, weightHat, phiHat, gradPhiHat);

      std::vector <double> particleVelOld(dim);
      particles[iMarker]->GetMarkerVelocity(particleVelOld);

      std::vector <double> particleAccOld(dim);
      particles[iMarker]->GetMarkerAcceleration(particleAccOld);

      std::vector <double> particleDisp(dim, 0.);
      //update displacement and acceleration
      for(int i = 0; i < dim; i++) {
        for(unsigned inode = 0; inode < nDofs; inode++) {
          particleDisp[i] += phiHat[inode] * solD[i][inode];
        }
      }

      particles[iMarker]->SetMarkerDisplacement(particleDisp);
      particles[iMarker]->UpdateParticleCoordinates();

      std::vector <double> particleAcc(dim);
      std::vector <double> particleVel(dim);
      for(unsigned i = 0; i < dim; i++) {
        particleAcc[i] = 1. / (beta * dt * dt) * particleDisp[i] - 1. / (beta * dt) * particleVelOld[i] - (1. - 2.* beta) / (2. * beta) * particleAccOld[i];
        particleVel[i] = particleVelOld[i] + dt * ((1. - Gamma) * particleAccOld[i] + Gamma * particleAcc[i]);
      }

      particles[iMarker]->SetMarkerVelocity(particleVel);
      particles[iMarker]->SetMarkerAcceleration(particleAcc);

      //   update the deformation gradient
      for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          gradSolDHat[i][j] = 0.;
          for(unsigned inode = 0; inode < nDofs; inode++) {
            gradSolDHat[i][j] +=  gradPhiHat[inode * dim + j] * solD[i][inode];
          }
        }
      }
      //if(C > cTreshold) {
        std::vector < std::vector < double > > FpOld;
        FpOld = particles[iMarker]->GetDeformationGradient(); //extraction of the deformation gradient
        double FpNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
        std::vector < std::vector < double > > Fp(dim);
        for(unsigned i = 0; i < dim; i++) {
          for(unsigned j = 0; j < dim; j++) {
            FpNew[i][j] += smoothstep(0., 0.5, C) * gradSolDHat[i][j];
          }
        }
        for(unsigned i = 0; i < dim; i++) {
          Fp[i].resize(dim);
          for(unsigned j = 0; j < dim; j++) {
            Fp[i][j] = 0.;
            for(unsigned k = 0; k < dim; k++) {
              Fp[i][j] += FpNew[i][k] * FpOld[k][j];
            }
          }
        }
        particles[iMarker]->SetDeformationGradient(Fp);
      //}
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
  for(unsigned iMarker = markerOffset[iproc]; iMarker < markerOffset[iproc + 1]; iMarker++) {
    unsigned iel = particles[iMarker]->GetMarkerElement();
    if(iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nDofs;
      //update element related quantities only if we are in a different element
      if(iel != ielOld) {
        ielt = msh->GetElementType(iel);
        nDofs = msh->GetElementDofNumber(iel, solType);
        for(int i = 0; i < dim; i++) {
          solD[i].resize(nDofs);
          solDOld[i].resize(nDofs);
          vxHat[i].resize(nDofs);
        }

        for(unsigned inode = 0; inode < nDofs; inode++) {
          unsigned idof = msh->GetSolutionDof(inode, iel, solType);   //local 2 global solution
          unsigned idofX = msh->GetSolutionDof(inode, iel, 2);   //local 2 global solution
          for(int i = 0; i < dim; i++) {
            solDOld[i][inode] = (*mysolution->_SolOld[indexSolD[i]])(idof);
            solD[i][inode] = (*mysolution->_Sol[indexSolD[i]])(idof) - solDOld[i][inode];
            //moving domain
            vxHat[i][inode] = (*msh->_topology->_Sol[i])(idofX) + solDOld[i][inode];
          }
        }

      }
      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();

      msh->_finiteElement[ielt][solType]->Jacobian(vxHat, xi, weightHat, phiHat, gradPhiHat);


      std::vector <double> particleDisp(dim, 0.);
      //update displacement and acceleration
      for(int i = 0; i < dim; i++) {
        for(unsigned inode = 0; inode < nDofs; inode++) {
          particleDisp[i] += phiHat[inode] * solD[i][inode];
        }
      }

      particles[iMarker]->SetMarkerDisplacement(particleDisp);
      particles[iMarker]->UpdateParticleCoordinates();

      //   update the deformation gradient
      for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          gradSolDHat[i][j] = 0.;
          for(unsigned inode = 0; inode < nDofs; inode++) {
            gradSolDHat[i][j] +=  gradPhiHat[inode * dim + j] * solD[i][inode];
          }
        }
      }
      std::vector < std::vector < double > > FpOld;
      FpOld = particles[iMarker]->GetDeformationGradient(); //extraction of the deformation gradient
      double FpNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
      std::vector < std::vector < double > > Fp(dim);
      for(unsigned i = 0; i < dim; i++) {
        for(unsigned j = 0; j < dim; j++) {
          FpNew[i][j] += gradSolDHat[i][j];
        }
      }
      for(unsigned i = 0; i < dim; i++) {
        Fp[i].resize(dim);
        for(unsigned j = 0; j < dim; j++) {
          Fp[i][j] = 0.;
          for(unsigned k = 0; k < dim; k++) {
            Fp[i][j] += FpNew[i][k] * FpOld[k][j];
          }
        }
      }
      particles[iMarker]->SetDeformationGradient(Fp);
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
  for(unsigned iMarker = markerOffset[iproc]; iMarker < markerOffset[iproc + 1]; iMarker++) {
    unsigned iel = particles[iMarker]->GetMarkerElement();
    if(iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nDofs;
      //update element related quantities only if we are in a different element
      if(iel != ielOld) {
        ielt = msh->GetElementType(iel);
        nDofs = msh->GetElementDofNumber(iel, solType);
        for(int i = 0; i < dim; i++) {
          solD[i].resize(nDofs);
          solDOld[i].resize(nDofs);
          vxHat[i].resize(nDofs);
        }

        for(unsigned inode = 0; inode < nDofs; inode++) {
          unsigned idof = msh->GetSolutionDof(inode, iel, solType);   //local 2 global solution
          unsigned idofX = msh->GetSolutionDof(inode, iel, 2);   //local 2 global solution
          for(int i = 0; i < dim; i++) {
            solDOld[i][inode] = (*mysolution->_SolOld[indexSolD[i]])(idof);
            solD[i][inode] = (*mysolution->_Sol[indexSolD[i]])(idof) - solDOld[i][inode];
            //moving domain
            vxHat[i][inode] = (*msh->_topology->_Sol[i])(idofX) + solDOld[i][inode];
          }
        }
      }
      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();

      msh->_finiteElement[ielt][solType]->Jacobian(vxHat, xi, weightHat, phiHat, gradPhiHat);

      std::vector <double> particleDisp(dim, 0.);
      //update displacement and acceleration
      for(int i = 0; i < dim; i++) {
        for(unsigned inode = 0; inode < nDofs; inode++) {
          particleDisp[i] += phiHat[inode] * solD[i][inode];
        }
      }

      particles[iMarker]->SetMarkerDisplacement(particleDisp);
      particles[iMarker]->UpdateParticleCoordinates();

      ielOld = iel;
    }
    else {
      break;
    }
  }
  //END loop on fluid particles

  ProjectGridVelocity(*mlSol);
//   ProjectGridVelocity2 (*mlSol);

  //BEGIN loop on elements to update grid velocity and acceleration
  for(unsigned idof = msh->_dofOffset[solType][iproc]; idof < msh->_dofOffset[solType][iproc + 1]; idof++) {
    for(int i = 0; i < dim; i++) {
      mysolution->_Sol[indexSolD[i]]->set(idof, 0.);
    }
  }

  for(int i = 0; i < dim; i++) {
    mysolution->_Sol[indexSolD[i]]->close();
  }
  //END loop on elements to update grid velocity and acceleration

  fluidLine.UpdateLineMPM();
  interfaceLine.UpdateLineMPM();
  solidLine.UpdateLineMPM();

  bool updateMat = false;
  fluidLine.GetParticlesToGridMaterial(updateMat);
  interfaceLine.GetParticlesToGridMaterial(updateMat);
  updateMat = true;
  solidLine.GetParticlesToGridMaterial(updateMat);

  GetParticlesToNodeFlag(*mlSol, solidLine, fluidLine);

}


unsigned getNumberOfLayers(const double & a, const double & fac, const bool inverse = true) {

  double fac1  = (inverse) ? fac : 1. / fac;
  double da = 1. / fac1;
  double b =  da;
  unsigned n = 1;

  while(b < a) {
    da /= fac1;
    b += da;
    n++;
    if(n >= 100) {
      std::cout << "Error: number of layer is unbounded, try with a smaller factor\n";
      abort();
    }
  }
  return n;
}

void GetParticlesToNodeFlag(MultiLevelSolution & mlSol, Line & solidLine, Line & fluidLine) {

  const unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;
  Mesh* msh = mlSol._mlMesh->GetLevel(level);
  Solution* sol  = mlSol.GetSolutionLevel(level);

  unsigned solIndexNodeFlag = sol->GetIndex("NodeFlag");
  unsigned solIndexNodeDist = sol->GetIndex("NodeDist");
  unsigned solType = sol->GetSolutionType(solIndexNodeFlag);

  unsigned indexSolMat = sol->GetIndex("Mat");

  sol->_Sol[solIndexNodeFlag]->zero(); // zero will mean fluid node

  const unsigned  dim = msh->GetDimension();
  unsigned    iproc = msh->processor_id();

  for(unsigned idof = msh->_dofOffset[solType][iproc]; idof < msh->_dofOffset[solType][iproc + 1]; idof++) {
    sol->_Sol[solIndexNodeDist]->set(idof, DBL_MAX);
  }
  sol->_Sol[solIndexNodeDist]->close();

  //BEGIN loop on solid particles

  vector <vector < double> > vxHat(dim);
  vector < unsigned > idof;

  std::vector<Marker*> particlesSolid = solidLine.GetParticles();
  std::vector<unsigned> markerOffset = solidLine.GetMarkerOffset();

  unsigned ielOld = UINT_MAX;
  for(unsigned iMarker = markerOffset[iproc]; iMarker < markerOffset[iproc + 1]; iMarker++) {

    unsigned iel = particlesSolid[iMarker]->GetMarkerElement();
    if(iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nDofs;
      if(iel != ielOld) {
        ielt = msh->GetElementType(iel);
        nDofs = msh->GetElementDofNumber(iel, solType);
        for(int k = 0; k < dim; k++) {
          vxHat[k].resize(nDofs);
        }
        idof.resize(nDofs);
        for(unsigned i = 0; i < nDofs; i++) {
          idof[i] = msh->GetSolutionDof(i, iel, solType);
          sol->_Sol[solIndexNodeFlag]->set(idof[i], 1.);
          unsigned idofX = msh->GetSolutionDof(i, iel, 2);   //local 2 global solution
          for(int k = 0; k < dim; k++) {
            vxHat[k][i] = (*msh->_topology->_Sol[k])(idofX);
          }
        }
      }
      std::vector<double> particleCoords(dim);
      particleCoords = particlesSolid[iMarker]->GetIprocMarkerCoordinates();

      for(unsigned i = 0; i < nDofs; i++) {
        double currentMinDist = (*sol->_Sol[solIndexNodeDist])(idof[i]);
        double newDist = 0.;
        for(unsigned k = 0; k < dim; k++) {
          newDist  += pow((vxHat[k][i] - particleCoords[k]), 2.);
        }
        newDist  = sqrt(newDist);
//         if (newDist  < 0.125 * 1.562e-06) {
//           sol->_Sol[solIndexNodeFlag]->set (idof[i], 1.);
//         }
        if(newDist  < currentMinDist) {
          sol->_Sol[solIndexNodeDist]->set(idof[i], newDist);
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

  for(unsigned iMarker = markerOffset[iproc]; iMarker < markerOffset[iproc + 1]; iMarker++) {
    unsigned iel = particlesFluid[iMarker]->GetMarkerElement();
    if(iel != UINT_MAX) {
      if((*sol->_Sol[indexSolMat])(iel) != 0) {       //only if it is an interface element
        short unsigned ielt;
        unsigned nDofs;
        if(iel != ielOld) {
          ielt = msh->GetElementType(iel);
          nDofs = msh->GetElementDofNumber(iel, solType);
          for(int k = 0; k < dim; k++) {
            vxHat[k].resize(nDofs);
          }
          idof.resize(nDofs);
          for(unsigned i = 0; i < nDofs; i++) {
            idof[i] = msh->GetSolutionDof(i, iel, solType);
            unsigned idofX = msh->GetSolutionDof(i, iel, 2);   //local 2 global solution
            for(int k = 0; k < dim; k++) {
              vxHat[k][i] = (*msh->_topology->_Sol[k])(idofX);
            }
          }
        }
        std::vector<double> particleCoords(dim);
        particleCoords = particlesFluid[iMarker]->GetIprocMarkerCoordinates();

        for(unsigned i = 0; i < nDofs; i++) {
          double currentMinDist = (*sol->_Sol[solIndexNodeDist])(idof[i]);
          double newDist = 0.;
          for(unsigned k = 0; k < dim; k++) {
            newDist  += pow((vxHat[k][i] - particleCoords[k]), 2.);
          }
          newDist  = sqrt(newDist);
          if(newDist  < currentMinDist) {
            sol->_Sol[solIndexNodeFlag]->set(idof[i], 0.);
            sol->_Sol[solIndexNodeDist]->set(idof[i], newDist);
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


  unsigned indexSolC = sol->GetIndex("C");

  std::vector<unsigned> markerOffsetSolid = solidLine.GetMarkerOffset();
  std::vector<unsigned> markerOffsetFluid = fluidLine.GetMarkerOffset();


  unsigned iSmarker = markerOffsetSolid[iproc];
  unsigned iFmarker = markerOffsetFluid[iproc];
  unsigned solIndexNodeFlag1 = sol->GetIndex("M");


  double rhos = 7850;
  double rhof = 1000;


// sol->_Sol[solIndexNodeFlag1]->zero();

  //BEGIN loop on elements (to initialize the "soft" stiffness matrix)
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielt = msh->GetElementType(iel);
    double  MPMmaterial = (*sol->_Sol[indexSolMat])(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs

    double C;

    if(MPMmaterial == 0) {
      sol->_Sol[indexSolC]->set(iel, 0.);
      C = 0.;
    }
    else if(MPMmaterial == nDofs) {
      sol->_Sol[indexSolC]->set(iel, 1.);
      C = 1.;
    }
    else {
      double volumef = 0.;
      double volumes = 0.;

      //BEGIN SOLID PARTICLE
      while(iSmarker < markerOffsetSolid[iproc + 1] && iel > particlesSolid[iSmarker]->GetMarkerElement()) {
        iSmarker++;
      }
      while(iSmarker < markerOffsetSolid[iproc + 1] && iel == particlesSolid[iSmarker]->GetMarkerElement()) {
        volumes += particlesSolid[iSmarker]->GetMarkerMass() / rhos;
        iSmarker++;
      }
      //END SOLID PARTICLE

      //BEGIN FLUID PARTICLE
      while(iFmarker < markerOffsetFluid[iproc + 1] && iel > particlesFluid[iFmarker]->GetMarkerElement()) {
        iFmarker++;
      }
      while(iFmarker < markerOffsetFluid[iproc + 1] && iel == particlesFluid[iFmarker]->GetMarkerElement()) {
        volumef += particlesFluid[iFmarker]->GetMarkerMass() / rhof;
        iFmarker++;
      }
      //END FLUID PARTICLE
      C = volumes / (volumes + volumef);
      sol->_Sol[indexSolC]->set(iel, C);
    }
//     if(C >= cTreshold) {
//       for(unsigned i = 0; i < nDofs; i++) {
//         unsigned idof = msh->GetSolutionDof(i, iel, solType);
//         sol->_Sol[solIndexNodeFlag1]->set(idof, 1.);
//       }
//     }
  }
  sol->_Sol[indexSolC]->close();
//  sol->_Sol[solIndexNodeFlag1]->close();

}

void ProjectGridVelocity(MultiLevelSolution &mlSol) {

  const unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;
  Mesh* msh = mlSol._mlMesh->GetLevel(level);
  Solution* sol  = mlSol.GetSolutionLevel(level);
  const unsigned dim = msh->GetDimension();

  unsigned iproc  = msh->processor_id();
  unsigned nprocs = msh->n_processors();

  vector< vector< double > > solV(dim);      // local solution (velocity)

  vector< bool > nodeFlag;     // local solution (velocity)
  vector< unsigned > idof;     // local solution (velocity)

  vector < double > phi;
  vector <vector < double> > vx(dim);   //vx is coordX in assembly of ex30
  vector <vector < double> > xp;


  //variable-name handling
  const char varname[10][5] = {"DX", "DY", "DZ", "VX", "VY", "VZ"};


  vector <unsigned> indexSolD(dim);
  vector <unsigned> indexSolV(dim);

  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = mlSol.GetIndex(&varname[ivar][0]);
    indexSolV[ivar] = mlSol.GetIndex(&varname[ivar + 3][0]);
  }
  unsigned indexNodeFlag =  mlSol.GetIndex("NodeFlag");

  unsigned solType = mlSol.GetSolutionType(&varname[0][0]);

  sol->_Sol[indexNodeFlag]->zero();

  for(unsigned k = 0; k < dim; k++) {
    (*sol->_SolOld[indexSolV[k]]) = (*sol->_Sol[indexSolV[k]]);
    sol->_Sol[indexSolV[k]]->zero();
  }

  unsigned counter = 0;

  std::vector < std::vector < std::vector <double > > > aP(3);

  //BEGIN loop on elements
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielType = msh->GetElementType(iel);

    unsigned nDofs = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs

    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofs);
      vx[k].resize(nDofs);
    }
    nodeFlag.resize(nDofs);
    idof.resize(nDofs);
    xp.resize(nDofs);
    for(unsigned i = 0; i < nDofs; i++) {
      xp[i].resize(dim);
    }

    for(unsigned i = 0; i < nDofs; i++) {
      idof[i] = msh->GetSolutionDof(i, iel, solType);
      unsigned idofX = msh->GetSolutionDof(i, iel, 2);

      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_SolOld[indexSolV[k]])(idof[i]);
        for(unsigned  k = 0; k < dim; k++) {
          xp[i][k] = (*msh->_topology->_Sol[k])(idofX);     // coordinates of the reference configuration;
          vx[k][i] = xp[i][k] + (*sol->_Sol[indexSolD[k]])(idof[i]);     // coordinates of the deformed configuration
        }
        nodeFlag[i] = ((*sol->_Sol[indexNodeFlag])(idof[i]) > 0.5) ? true : false;
      }
    }

    bool aPIsInitialized = false;

    double r;
    std::vector <double> xc;
    GetConvexHullSphere(vx, xc, r, 0.0001);  // get the ball that circumscribe the element
    double r2 = r * r;

    std::vector < std::vector< double > > xe; // get the box that encloses the element
    GetBoundingBox(vx, xe, 0.0001);

    for(unsigned i = 0; i < nDofs; i++) {  // loop on the nodes of the reference elements now considered as independent points
      if(!nodeFlag[i]) {
        double d2 = 0.;
        for(int k = 0; k < dim; k++) {
          d2 += (xp[i][k] - xc[k]) * (xp[i][k] - xc[k]);
        }
        bool insideHull = true;
        if(d2 > r2) {
          insideHull = false;
        }
        for(unsigned k = 0; k < dim; k++) {
          if(xp[i][k] < xe[k][0] || xp[i][k] > xe[k][1]) {
            insideHull = false;
          }
        }

        if(insideHull) {  //rough test
          if(!aPIsInitialized) {
            aPIsInitialized = true;
            std::vector < std::vector <double> > x1(dim);
            for(unsigned jtype = 0; jtype < solType + 1; jtype++) {
              ProjectNodalToPolynomialCoefficients(aP[jtype], vx, ielType, jtype) ;
            }
          }
          std::vector <double> xi;
          GetClosestPointInReferenceElement(vx, xp[i], ielType, xi);

          bool inverseMapping = GetInverseMapping(solType, ielType, aP, xp[i], xi, 100);
          if(!inverseMapping) {
            std::cout << "InverseMapping failed at " << iel << " " << idof[i] << std::endl;
          }


          bool insideDomain = CheckIfPointIsInsideReferenceDomain(xi, ielType, 1.e-3);  // fine testing

          if(inverseMapping && insideDomain) {
            sol->_Sol[indexNodeFlag]->add(idof[i], 1.);
            msh->_finiteElement[ielType][solType]->GetPhi(phi, xi);
            //std::cout << iel << " " << i << "  ";
            counter++;
            for(unsigned k = 0; k < dim; k++) {
              double solVk = 0.;
              for(unsigned j = 0; j < nDofs; j++)    {
                solVk += phi[j] * solV[k][j];
              }
              sol->_Sol[indexSolV[k]]->add(idof[i], solVk);
            }
          }
        }
      }
    }
  }

  sol->_Sol[indexNodeFlag]->close();

  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexSolV[k]]->close();
  }

  unsigned c0 = 0;
  for(unsigned i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {
    unsigned cnt = static_cast < unsigned >(floor((*sol->_Sol[indexNodeFlag])(i) + 0.5));
    if(cnt == 0) {
      c0++;
    }
    else if(cnt > 1) {
      counter -= (cnt - 1);
      for(unsigned k = 0; k < dim; k++) {
        double velk = (*sol->_Sol[indexSolV[k]])(i) / cnt;
        sol->_Sol[indexSolV[k]]->set(i, velk);
      }
    }
  }

  unsigned counterAll;
  MPI_Reduce(&counter, &counterAll, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  std::cout << "COUNTER = " << counterAll << " " << msh->GetTotalNumberOfDofs(solType) << std::endl;

  idof.resize(c0);
  vector< double > xp0(c0 * dim);

  unsigned c1 = 0;
  for(unsigned i = msh->_dofOffset[solType][iproc]; i < msh->_dofOffset[solType][iproc + 1]; i++) {
    if(static_cast < unsigned >(floor((*sol->_Sol[indexNodeFlag])(i) + 0.5)) == 0) {
      idof[c1] = i;
      for(unsigned k = 0; k < dim; k++) {
        xp0[c1 * dim + k] = (*msh->_topology->_Sol[k])(i);
      }
      c1++;
      if(c1 == c0) break;
    }
  }

  vector< double > xp1;
  for(unsigned jproc = 0; jproc < nprocs; jproc++) {
    c1 = c0;
    MPI_Bcast(&c1, 1, MPI_UNSIGNED, jproc, PETSC_COMM_WORLD);
    if(c1) {
      xp1.resize(c1 * dim);
      if(iproc == jproc) {
        for(unsigned i = 0; i < c1; i++) {
          for(unsigned k = 0; k < dim; k++) {
            xp1[i * dim + k] = xp0[i * dim + k];
          }
        }
      }
      MPI_Bcast(&xp1[0], c1 * dim, MPI_DOUBLE, jproc, PETSC_COMM_WORLD);

      for(unsigned i = 0; i < c1; i++) {
        std::vector < double > xp(dim);
        for(unsigned k = 0; k < dim; k++) {
          xp[k] = xp1[i * dim + k];
        }
        Marker p(xp, 1, VOLUME, sol, solType, true, 1.);
        unsigned mproc = p.GetMarkerProc(sol);
        if(iproc == mproc) {
          unsigned jel = p.GetMarkerElement();
          short unsigned jelType = msh->GetElementType(jel);
          unsigned nDofs = msh->GetElementDofNumber(jel, solType);    // number of solution element dofs

          for(unsigned  k = 0; k < dim; k++) {
            solV[k].resize(nDofs);
            vx[k].resize(nDofs);
          }

          for(unsigned j = 0; j < nDofs; j++) {
            unsigned jdof = msh->GetSolutionDof(j, jel, solType);
            unsigned jdofX = msh->GetSolutionDof(j, jel, 2);
            for(unsigned  k = 0; k < dim; k++) {
              solV[k][j] = (*sol->_SolOld[indexSolV[k]])(jdof);     //velocity to be projected
              vx[k][j] = (*msh->_topology->_Sol[k])(jdofX) + (*sol->_Sol[indexSolD[k]])(jdof);         // coordinates of the deformed configuration
            }
          }
          std::vector <double> xi = p.GetMarkerLocalCoordinates();
          msh->_finiteElement[jelType][solType]->GetPhi(phi, xi);

          std::vector < double > vel(dim, 0.);
          for(unsigned k = 0; k < dim; k++) {
            for(unsigned j = 0; j < nDofs; j++)    {
              vel[k] += phi[j] * solV[k][j];
            }
          }
          MPI_Send(&vel[0], dim, MPI_DOUBLE, jproc, 1, PETSC_COMM_WORLD);

        }
        if(iproc == jproc) {
          std::vector < double > vel(dim);
          MPI_Recv(&vel[0], dim, MPI_DOUBLE, mproc, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
          for(unsigned k = 0; k < dim; k++) {
            sol->_Sol[indexSolV[k]]->set(idof[i], vel[k]);
          }
          counter++;
        }
      }
    }
  }

  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexSolV[k]]->close();
  }

  MPI_Reduce(&counter, &counterAll, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  std::cout << "COUNTER = " << counterAll << " " << msh->GetTotalNumberOfDofs(solType) << std::endl;

}



