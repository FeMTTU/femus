/* This example is for quasi-conformal minimization */
/* mu controls the beltrami coefficient */

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "TransientSystem.hpp"
#include "adept.h"
#include <cstdlib>
#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"

bool stopIterate = false;

unsigned conformalTriangleType = 2;
double eps = 1e-5;
const double normalSign = -1.;
bool O2conformal = false;
unsigned counter = 0;

using namespace femus;
void UpdateMu(MultiLevelSolution& mlSol);


#include "../include/supportFunctions.hpp"
#include "../include/assembleConformalMinimization.hpp"
// Comment back in for working code
//const double mu[2] = {0.8, 0.};



void AssembleConformalMinimizationOld(MultiLevelProblem& ml_prob);   //stable and not bad
void AssembleShearMinimization(MultiLevelProblem& ml_prob);

void AssembleConformalO1Minimization(MultiLevelProblem& ml_prob);

// IBVs.  No boundary, and IVs set to sphere (just need something).
bool SetBoundaryCondition(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {


  bool dirichlet = true;
  value = 0.;

  if(!strcmp(solName, "Dx1")) {
    if(1 == faceName || 3 == faceName) {
      dirichlet = false;
    }
    if(4 == faceName) {
      value = 0.5 * sin(x[1] / 0.5 * M_PI);
    }
  }
  else if(!strcmp(solName, "Dx2")) {
    if(2 == faceName) {
      dirichlet = false;
    }
  }


//   if (!strcmp (solName, "Dx1")) {
//     if (1 == faceName || 3 == faceName) {
//       dirichlet = false;
//     }
//   }
//   else if (!strcmp (solName, "Dx2")) {
//     if (2 == faceName || 4 == faceName) {
//       dirichlet = false;
//     }
//   }



  /*

  bool dirichlet = true;
  value = 0.;

  if(!strcmp(solName, "Dx1")) {
  if(1 == faceName) {
    //value = 0.04 * sin (4*(x[1] / 0.5 * acos (-1.)));
    value = 0. * sin(x[1] / 0.5 * M_PI);
    //dirichlet = false;
  }
  }*/


  return dirichlet;
}


// Main program starts here.
int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);


  // define multilevel mesh
  unsigned maxNumberOfMeshes;
  MultiLevelMesh mlMsh;

  // Read coarse level mesh and generate finer level meshes.
  double scalingFactor = 1.;

  //mlMsh.GenerateCoarseBoxMesh(32, 32, 0, -0.5, 0.5, -0.5, 0.5, 0., 0., QUAD9, "seventh");

  //mlMsh.ReadCoarseMesh("../input/squareReg3D.neu", "seventh", scalingFactor);
  mlMsh.ReadCoarseMesh("../input/square13D.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("../input/cylinder2.neu", "seventh", scalingFactor);


  unsigned numberOfUniformLevels = 4;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // Erase all the coarse mesh levels.
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  const unsigned  DIM = mlMsh.GetDimension();

  // Define the multilevel solution and attach the mlMsh object to it.
  MultiLevelSolution mlSol(&mlMsh);

  // Add variables X,Y,W to mlSol.

  FEOrder feOrder = FIRST;
  mlSol.AddSolution("Dx1", LAGRANGE, feOrder, 0);
  mlSol.AddSolution("Dx2", LAGRANGE, feOrder, 0);
  mlSol.AddSolution("Dx3", LAGRANGE, feOrder, 0);

  mlSol.AddSolution("Lambda1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0);

  mlSol.AddSolution("mu1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("mu2", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

  mlSol.AddSolution("muN1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("weight1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

  mlSol.AddSolution("muN2", LAGRANGE, feOrder, 0, false);
  mlSol.AddSolution("weight2", LAGRANGE, feOrder, 0, false);


  mlSol.AddSolution("theta1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("theta2", LAGRANGE, feOrder, 0, false);

  mlSol.AddSolution("phi1", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
  mlSol.AddSolution("phi2", LAGRANGE, feOrder, 0, false);

  // Initialize the variables and attach boundary conditions.
  mlSol.Initialize("All");

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");

  MultiLevelProblem mlProb(&mlSol);

  // Add system Conformal or Shear Minimization in mlProb.
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("conformal"); //for conformal

  // Add solutions newDX, Lambda1 to system.
  system.AddSolutionToSystemPDE("Dx1");
  system.AddSolutionToSystemPDE("Dx2");
  system.AddSolutionToSystemPDE("Dx3");
  system.AddSolutionToSystemPDE("Lambda1");

  // Parameters for convergence and # of iterations.
  system.SetMaxNumberOfNonLinearIterations(100);
  system.SetNonLinearConvergenceTolerance(1.e-10);

  system.init();

  mlSol.SetWriter(VTK);
  std::vector<std::string> mov_vars;
  mov_vars.push_back("Dx1");
  mov_vars.push_back("Dx2");
  mov_vars.push_back("Dx3");
  mlSol.GetWriter()->SetMovingMesh(mov_vars);

  // and this?
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");
  mlSol.GetWriter()->SetDebugOutput(true);
  //mlSol.GetWriter()->Write (Files::_application_output_directory, "linear", variablesToBePrinted, 0);
  mlSol.GetWriter()->Write(Files::_application_output_directory, "biquadratic", variablesToBePrinted, 0);

  // Attach the assembling function to system and initialize.
  //system.SetAssembleFunction(AssembleShearMinimization);
  //system.SetAssembleFunction (AssembleConformalMinimization);
  //system.MGsolve();
  mlSol.GetWriter()->Write(Files::_application_output_directory, "biquadratic", variablesToBePrinted, 1);

  system.SetAssembleFunction(AssembleConformalMinimization);
  //system.SetAssembleFunction(AssembleConformalO1Minimization);
  system.MGsolve();

  mlSol.GetWriter()->Write(Files::_application_output_directory, "biquadratic", variablesToBePrinted, 2);

  return 0;
}


void UpdateMu(MultiLevelSolution& mlSol) {

  //MultiLevelSolution*  mlSol = ml_prob._ml_sol;
  unsigned level = mlSol.GetMLMesh()->GetNumberOfLevels() - 1u;

  Solution* sol = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol.GetMLMesh()->GetLevel(level);
  elem* el = msh->GetMeshElements();

  //unsigned  dim = msh->GetDimension();
  unsigned dim = 2;
  unsigned DIM = 3;

  std::vector < unsigned > indexDx(DIM);
  indexDx[0] = mlSol.GetIndex("Dx1");
  indexDx[1] = mlSol.GetIndex("Dx2");
  indexDx[2] = mlSol.GetIndex("Dx3");
  unsigned solTypeDx = mlSol.GetSolutionType(indexDx[0]);

  std::vector < unsigned > indexMu(dim);
  indexMu[0] = mlSol.GetIndex("mu1");
  indexMu[1] = mlSol.GetIndex("mu2");

  unsigned indexMuN1 = mlSol.GetIndex("muN1"); //piecewice linear discontinuous

  unsigned indexW1 = mlSol.GetIndex("weight1");
  unsigned solType1 = mlSol.GetSolutionType(indexMu[0]);

  unsigned indexTheta1 = mlSol.GetIndex("theta1");
  unsigned indexTheta2 = mlSol.GetIndex("theta2");

  unsigned indexPhi1 = mlSol.GetIndex("phi1");
  unsigned indexPhi2 = mlSol.GetIndex("phi2");

  std::vector< double > dof1;

  std::vector < std::vector < double > > solx(DIM);
  std::vector < std::vector < double > > xHat(DIM);

  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexMu[k]]->zero();
  }
  sol->_Sol[indexW1]->zero();

  std::vector < double > phi;  // local test function for velocity
  std::vector < double > phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  unsigned iproc = msh->processor_id();
  unsigned nprocs = msh->n_processors();

// Setting the reference elements to be equilateral triangles.
  std::vector < std::vector < double > > xT(2);
  xT[0].resize(3);
  xT[0][0] = -0.5;
  xT[0][1] = 0.5;
  xT[0][2] = 0.;

  xT[1].resize(3);
  xT[1][0] = 0.;
  xT[1][1] = 0.;
  xT[1][2] = sqrt(3.) / 2.;

  std::vector<double> phi_uv0;
  std::vector<double> phi_uv1;

  std::vector< double > stdVectorPhi;
  std::vector< double > stdVectorPhi_uv;

  for(int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs1  = msh->GetElementDofNumber(iel, solType1);
    unsigned nDofsDx  = msh->GetElementDofNumber(iel, solTypeDx);

    dof1.resize(nDofs1);

    for(int K = 0; K < DIM; K++) {
      xHat[K].resize(nDofsDx);
      solx[K].resize(nDofsDx);
    }

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nDofs1; i++) {
      dof1[i] = msh->GetSolutionDof(i, iel, solType1);
    }
    // local storage of coordinates
    for(unsigned i = 0; i < nDofsDx; i++) {
      unsigned idof = msh->GetSolutionDof(i, iel, solTypeDx);
      unsigned xDof  = msh->GetSolutionDof(i, iel, 2);
      for(unsigned K = 0; K < DIM; K++) {
        xHat[K][i] = (*msh->GetTopology()->_Sol[K])(xDof);
        solx[K][i] = xHat[K][i] + (*sol->_Sol[indexDx[K]])(idof);
      }
    }

//####### DIDN'T CHANGE THIS
    // for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solTypeDx]->GetGaussPointNumber(); ig++) {

    // msh->_finiteElement[ielGeom][solTypeDx]->Jacobian(xHat, ig, weight, phi, phi_x);
    //
    // // std::vector < std::vector < double > > gradSolx(dim);
    // //
    // // for(unsigned  k = 0; k < dim; k++) {
    // //   gradSolx[k].assign(dim, 0.);
    // // }
    //
    // for(unsigned i = 0; i < nDofsDx; i++) {
    //   for(unsigned j = 0; j < dim; j++) {
    //     for(unsigned  k = 0; k < dim; k++) {
    //       gradSolx[k][j] += solx[k][i] * phi_x[i * dim + j];
    //     }
    //   }
    // }
// ##############

// *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solTypeDx]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function

      const double *phix_uv[dim]; // local test function first order partial derivatives

      double weight; // gauss point weight

      // Get Gauss point weight, test function, and first order derivatives.
      if(ielGeom == QUAD) {
        phix = msh->_finiteElement[ielGeom][solTypeDx]->GetPhi(ig);


        phix_uv[0] = msh->_finiteElement[ielGeom][solTypeDx]->GetDPhiDXi(ig);
        phix_uv[1] = msh->_finiteElement[ielGeom][solTypeDx]->GetDPhiDEta(ig);

        weight = msh->_finiteElement[ielGeom][solTypeDx]->GetGaussWeight(ig);
      }

      // Special adjustments for triangles.
      else {
        msh->_finiteElement[ielGeom][solTypeDx]->Jacobian(xT, ig, weight, stdVectorPhi, stdVectorPhi_uv);

        phix = &stdVectorPhi[0];
        phi_uv0.resize(nDofsDx);
        phi_uv1.resize(nDofsDx);


        for(unsigned i = 0; i < nDofsDx; i++) {
          phi_uv0[i] = stdVectorPhi_uv[i * dim];
          phi_uv1[i] = stdVectorPhi_uv[i * dim + 1];
        }

        phix_uv[0] = &phi_uv0[0];
        phix_uv[1] = &phi_uv1[0];

      }
      const double *phi1 = msh->_finiteElement[ielGeom][solType1]->GetPhi(ig);  // local test function

      // Initialize and compute values of x, Dx, NDx, x_uv at the Gauss points.
      //double solDxg[3] = {0., 0., 0.};
      double solxg[3] = {0., 0., 0.};
      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      for(unsigned K = 0; K < DIM; K++) {
        for(unsigned i = 0; i < nDofsDx; i++) {
          solxg[K] += phix[i] * solx[K][i];
        }
        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nDofsDx; i++) {
            solx_uv[K][j]    += phix_uv[j][i] * solx[K][i];
          }
        }
      }

      // Compute the metric, metric determinant, and area element.
      double g[2][2] = {{0., 0.}, {0., 0.}};
      for(unsigned i = 0; i < dim; i++) {
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }
      double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      //double Area = weight * sqrt(detg);
      //double Area2 = weight; // Trick to give equal weight to each element.

      // std::cout << detg << " ";

      double normal[DIM];
       normal[0] = 0; // (solx_uv[1][0] * solx_uv[2][1] - solx_uv[2][0] * solx_uv[1][1]) / sqrt(detg);
       normal[1] = 0; //(solx_uv[2][0] * solx_uv[0][1] - solx_uv[0][0] * solx_uv[2][1]) / sqrt(detg);
       normal[2] = 1; //(solx_uv[0][0] * solx_uv[1][1] - solx_uv[1][0] * solx_uv[0][1]) / sqrt(detg);
      //normal[0] = (solx_uv[1][0] * solx_uv[2][1] - solx_uv[2][0] * solx_uv[1][1]) * sqrt(detg) / detg;
      //normal[1] = (solx_uv[2][0] * solx_uv[0][1] - solx_uv[0][0] * solx_uv[2][1]) * sqrt(detg) / detg;
      //normal[2] = (solx_uv[0][0] * solx_uv[1][1] - solx_uv[1][0] * solx_uv[0][1]) * sqrt(detg) / detg;
      // normal[0] = 0;
      // normal[1] = solxg[1] / sqrt(solxg[1] * solxg[1] + solxg[2] * solxg[2]);
      // normal[2] = solxg[2] / sqrt(solxg[1] * solxg[1] + solxg[2] * solxg[2]);


      double dxPlus[DIM];
      dxPlus[0] = solx_uv[0][0] + solx_uv[1][1] * normal[2] - solx_uv[2][1] * normal[1];
      dxPlus[1] = solx_uv[1][0] + solx_uv[2][1] * normal[0] - solx_uv[0][1] * normal[2];
      dxPlus[2] = solx_uv[2][0] + solx_uv[0][1] * normal[1] - solx_uv[1][1] * normal[0];

      double sdxPlus[DIM];
      sdxPlus[0] = solx_uv[0][1] - solx_uv[1][0] * normal[2] + solx_uv[2][0] * normal[1];
      sdxPlus[1] = solx_uv[1][1] - solx_uv[2][0] * normal[0] + solx_uv[0][0] * normal[2];
      sdxPlus[2] = solx_uv[2][1] - solx_uv[0][0] * normal[1] + solx_uv[1][0] * normal[0];

      double dxMinus[DIM];
      dxMinus[0] = solx_uv[0][0] - solx_uv[1][1] * normal[2] + solx_uv[2][1] * normal[1];
      dxMinus[1] = solx_uv[1][0] - solx_uv[2][1] * normal[0] + solx_uv[0][1] * normal[2];
      dxMinus[2] = solx_uv[2][0] - solx_uv[0][1] * normal[1] + solx_uv[1][1] * normal[0];

      double norm2dxPlus = 0;
      double norm2sdxPlus = 0;
      double rhsmu1 = 0;
      double rhsmu2 = 0;

      //FAKE COMPUTATION NOTE THE SIGNS
      for(unsigned K = 0; K < DIM; K++) {
        norm2dxPlus += dxPlus[K] * dxPlus[K];
        norm2sdxPlus += sdxPlus[K] * sdxPlus[K];
        rhsmu1 += dxPlus[K] * dxMinus[K];
        rhsmu2 += sdxPlus[K] * dxMinus[K];

        //dxsdxp += dxMinus[K] * dxMinus[K];
      }

      // double norm2Xz = (1. / 4.) * (pow((gradSolx[0][0] + gradSolx[1][1]), 2) + pow((gradSolx[1][0] - gradSolx[0][1]), 2));
      // double XzBarXz_Bar[2];
      //
      // XzBarXz_Bar[0] = (1. / 4.) * (pow(gradSolx[0][0], 2) + pow(gradSolx[1][0], 2) - pow(gradSolx[0][1], 2) - pow(gradSolx[1][1], 2));
      // XzBarXz_Bar[1] = (1. / 2.) * (gradSolx[0][0] * gradSolx[0][1] + gradSolx[1][0] * gradSolx[1][1]);

      // Comment out for working code

      double mu[2] = {0., 0.};
      mu[0] = rhsmu1 / norm2dxPlus;
      mu[1] = rhsmu2 / norm2sdxPlus;

      if(iel == 4) {
        std::cout << mu[0] << " " << mu[1] << " " << norm2dxPlus << " " << norm2sdxPlus << "\n";
      }

      // for(unsigned k = 0; k < 2; k++) {
      //   if(norm2Xz > 0.) {
      //     mu[k] += (1. / norm2Xz) * XzBarXz_Bar[k];
      //   }
      // }

      for(unsigned i = 0; i < nDofs1; i++) {
        sol->_Sol[indexW1]->add(dof1[i], phi1[i] * weight);
        for(unsigned k = 0; k < dim; k++) {
          sol->_Sol[indexMu[k]]->add(dof1[i], mu[k] * phi1[i] * weight);
        }
      } // end phi_i loop

    } // end gauss point loop

  } //end element loop for each process*/

  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexMu[k]]->close();
  }
  sol->_Sol[indexW1]->close();

  sol->_Sol[indexTheta1]->zero();
  sol->_Sol[indexPhi1]->zero();

  for(unsigned i = msh->_dofOffset[solType1][iproc]; i < msh->_dofOffset[solType1][iproc + 1]; i++) {

    double weight = (*sol->_Sol[indexW1])(i);

    //std::cout << weight << " ";

    double mu[2];
    for(unsigned k = 0; k < dim; k++) {
      mu[k] = (*sol->_Sol[indexMu[k]])(i);
      sol->_Sol[indexMu[k]]->set(i, mu[k] / weight);
    }
    sol->_Sol[indexMuN1]->set(i, sqrt(mu[0] * mu[0] + mu[1] * mu[1]) / weight);
    sol->_Sol[indexTheta1]->set(i, atan2(mu[1] / weight, fabs(mu[0]) / weight));
    sol->_Sol[indexPhi1]->set(i, atan2(mu[0] / weight, fabs(mu[1]) / weight));
  }
  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexMu[k]]->close();
  }
  sol->_Sol[indexMuN1]->close();
  sol->_Sol[indexTheta1]->close();
  sol->_Sol[indexPhi1]->close();

  double norm = sol->_Sol[indexMuN1]->linfty_norm();
  std::cout << norm << std::endl;

  //BEGIN Iterative smoothing element -> nodes -> element

  for(unsigned smooth = 0; smooth < 0; smooth++) {
    unsigned indexW2 = mlSol.GetIndex("weight2");
    unsigned indexMuN2 = mlSol.GetIndex("muN2");  //smooth ni norm

    unsigned solType2 = mlSol.GetSolutionType(indexMuN2);

    std::vector< double > dof2;
    std::vector< double > sol1;
    std::vector< double > solTheta1;
    std::vector< double > solPhi1;

    sol->_Sol[indexMuN2]->zero();
    sol->_Sol[indexW2]->zero();
    sol->_Sol[indexTheta2]->zero();
    sol->_Sol[indexPhi2]->zero();

    //std::vector < double > phi2;  // local test function for velocity
    //std::vector < double > phi2_x; // local test function first order partial derivatives
    //double weight2; // gauss point weight

    for(int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned nDofs1  = msh->GetElementDofNumber(iel, solType1);
      unsigned nDofs2  = msh->GetElementDofNumber(iel, solType2);

      sol1.resize(nDofs1);
      solTheta1.resize(nDofs1);
      solPhi1.resize(nDofs1);

      dof2.resize(nDofs2);

      for(int K = 0; K < DIM; K++) {
        xHat[K].resize(nDofs2);
      }

      for(unsigned i = 0; i < nDofs1; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solType1);
        sol1[i] = (*sol->_Sol[indexMuN1])(idof);
        solTheta1[i] = (*sol->_Sol[indexTheta1])(idof);
        solPhi1[i] = (*sol->_Sol[indexPhi1])(idof);
      }

      // local storage of global mapping and solution
      for(unsigned i = 0; i < nDofs2; i++) {
        dof2[i] = msh->GetSolutionDof(i, iel, solType2);
      }
      // local storage of coordinates
      for(unsigned i = 0; i < nDofs2; i++) {
        unsigned xDof  = msh->GetSolutionDof(i, iel, 2);
        for(unsigned K = 0; K < DIM; K++) {
          xHat[K][i] = (*msh->GetTopology()->_Sol[K])(xDof);
        }
      }

      for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType2]->GetGaussPointNumber(); ig++) {
//         msh->_finiteElement[ielGeom][solType2]->Jacobian(xHat, ig, weight2, phi2, phi2_x);

        double *phi2 = msh->_finiteElement[ielGeom][solType2]->GetPhi(ig);
        double weight2 = msh->_finiteElement[ielGeom][solType2]->GetGaussWeight(ig);
        double *phi1 = msh->_finiteElement[ielGeom][solType1]->GetPhi(ig);
        double sol1g = 0.;
        double solTheta1g = 0.;
        double solPhi1g = 0.;
        for(unsigned i = 0; i < nDofs1; i++) {
          sol1g += phi1[i] * sol1[i];
          solTheta1g += phi1[i] * solTheta1[i];
          solPhi1g += phi1[i] * solPhi1[i];
        }

        // *** phi_i loop ***
        for(unsigned i = 0; i < nDofs2; i++) {
          sol->_Sol[indexW2]->add(dof2[i], phi2[i] * weight2);
          sol->_Sol[indexMuN2]->add(dof2[i], sol1g * phi2[i] * weight2);
          sol->_Sol[indexTheta2]->add(dof2[i], solTheta1g * phi2[i] * weight2);
          sol->_Sol[indexPhi2]->add(dof2[i], solPhi1g * phi2[i] * weight2);
        } // end phi_i loop
      } // end gauss point loop

    } //end element loop for each process*/
    sol->_Sol[indexW2]->close();
    sol->_Sol[indexMuN2]->close();
    sol->_Sol[indexTheta2]->close();
    sol->_Sol[indexPhi2]->close();

    for(unsigned i = msh->_dofOffset[solType2][iproc]; i < msh->_dofOffset[solType2][iproc + 1]; i++) {
      double weight = (*sol->_Sol[indexW2])(i);
      double value = (*sol->_Sol[indexMuN2])(i);
      sol->_Sol[indexMuN2]->set(i, value / weight);

      value = (*sol->_Sol[indexTheta2])(i);
      sol->_Sol[indexTheta2]->set(i, value / weight);

      value = (*sol->_Sol[indexPhi2])(i);
      sol->_Sol[indexPhi2]->set(i, value / weight);
    }
    sol->_Sol[indexMuN2]->close();
    sol->_Sol[indexTheta2]->close();
    sol->_Sol[indexPhi2]->close();


    sol->_Sol[indexMuN1]->zero();
    sol->_Sol[indexTheta1]->zero();
    sol->_Sol[indexPhi1]->zero();

    std::vector< double > sol2;
    std::vector< double > solTheta2;
    std::vector< double > solPhi2;

    for(int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

      short unsigned ielGeom = msh->GetElementType(iel);
      unsigned nDofs1  = msh->GetElementDofNumber(iel, solType1);
      unsigned nDofs2  = msh->GetElementDofNumber(iel, solType2);

      dof1.resize(nDofs1);

      for(int k = 0; k < dim; k++) {
        //xHat[k].resize(nDofs2);
        sol2.resize(nDofs2);
        solTheta2.resize(nDofs2);
        solPhi2.resize(nDofs2);
      }
      for(int K = 0; K < DIM; K++) {
        xHat[K].resize(nDofs2);
      }

      // local storage of global mapping and solution
      for(unsigned i = 0; i < nDofs1; i++) {
        dof1[i] = msh->GetSolutionDof(i, iel, solType1);
      }
      // local storage of coordinates
      for(unsigned i = 0; i < nDofs2; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solType2);
        unsigned xDof  = msh->GetSolutionDof(i, iel, 2);
        for(unsigned K = 0; K < DIM; K++) {
          xHat[K][i] = (*msh->GetTopology()->_Sol[K])(xDof);
        }
        sol2[i] = (*sol->_Sol[indexMuN2])(idof);
        solTheta2[i] = (*sol->_Sol[indexTheta2])(idof);
        solPhi2[i] = (*sol->_Sol[indexPhi2])(idof);

      }

      for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solTypeDx]->GetGaussPointNumber(); ig++) {

        //msh->_finiteElement[ielGeom][solTypeDx]->Jacobian(xHat, ig, weight, phi, phi_x);

        weight = msh->_finiteElement[ielGeom][solTypeDx]->GetGaussWeight(ig);
        double *phi = msh->_finiteElement[ielGeom][solTypeDx]->GetPhi(ig);

        double sol2g = 0;
        double solTheta2g = 0;
        double solPhi2g = 0;

        for(unsigned i = 0; i < nDofs2; i++) {
          sol2g += sol2[i] * phi[i];
          solTheta2g += solTheta2[i] * phi[i];
          solPhi2g += solPhi2[i] * phi[i];
        }

        double *phi1 = msh->_finiteElement[ielGeom][solType1]->GetPhi(ig);

        for(unsigned i = 0; i < nDofs1; i++) {
          sol->_Sol[indexMuN1]->add(dof1[i], sol2g * phi1[i] * weight);
          sol->_Sol[indexTheta1]->add(dof1[i], solTheta2g * phi1[i] * weight);
          sol->_Sol[indexPhi1]->add(dof1[i], solPhi2g * phi1[i] * weight);
        } // end phi_i loop

      } // end gauss point loop

    } //end element loop for each process*/

    sol->_Sol[indexMuN1]->close();
    sol->_Sol[indexTheta1]->close();
    sol->_Sol[indexPhi1]->close();

    for(unsigned i = msh->_dofOffset[solType1][iproc]; i < msh->_dofOffset[solType1][iproc + 1]; i++) {
      double weight = (*sol->_Sol[indexW1])(i);
      double radius = (*sol->_Sol[indexMuN1])(i) / weight;
      sol->_Sol[indexMuN1]->set(i,  radius);

      double theta = (*sol->_Sol[indexTheta1])(i) / weight;
      sol->_Sol[indexTheta1]->set(i,  theta);

      double phi = (*sol->_Sol[indexPhi1])(i) / weight;
      sol->_Sol[indexPhi1]->set(i,  phi);


      double mu[2];
      for(unsigned k = 0; k < dim; k++) {
        mu[k] = (*sol->_Sol[indexMu[k]])(i);
      }

//       if(fabs(theta) < fabs(phi)) {
//         if(mu[0] < 0) {
//           if(mu[1] < 0) theta = -theta - M_PI;
//           else theta = M_PI - theta;
//         }
//       }
//       else {
//         if(mu[1] < 0) {
//           theta = -M_PI / 2 + phi;
//         }
//         else {
//           theta =  M_PI / 2 - phi;
//         }
//       }

//       sol->_Sol[indexMu[0]]->set(i, radius * cos(theta));
//       sol->_Sol[indexMu[1]]->set(i, radius * sin(theta));

    }
    sol->_Sol[indexMuN1]->close();
    sol->_Sol[indexTheta1]->close();
    sol->_Sol[indexPhi1]->close();

//     for(unsigned k = 0; k < dim; k++) {
//       sol->_Sol[indexMu[k]]->close();
//     }
  }
  //END Iterative smoothing element -> nodes -> element


  //BEGIN mu update
  double MuNormLocalSum = 0.;
  for(unsigned i = msh->_dofOffset[solType1][iproc]; i < msh->_dofOffset[solType1][iproc + 1]; i++) {
    MuNormLocalSum += (*sol->_Sol[indexMuN1])(i);
  }

  double MuNormAverage;
  MPI_Allreduce(&MuNormLocalSum, &MuNormAverage, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MuNormAverage /= msh->_dofOffset[solType1][nprocs];

  for(unsigned i = msh->_dofOffset[solType1][iproc]; i < msh->_dofOffset[solType1][iproc + 1]; i++) {

    double theta = (*sol->_Sol[indexTheta1])(i);
    double phi = (*sol->_Sol[indexPhi1])(i);

    double mu[2];
    for(unsigned k = 0; k < dim; k++) {
      mu[k] = (*sol->_Sol[indexMu[k]])(i);
    }

    if(fabs(theta) < fabs(phi)) {
      if(mu[0] < 0) {
        if(mu[1] < 0) theta = -theta - M_PI;
        else theta = M_PI - theta;
      }
    }
    else {
      if(mu[1] < 0) {
        theta = -M_PI / 2 + phi;
      }
      else {
        theta =  M_PI / 2 - phi;
      }
    }

    sol->_Sol[indexMu[0]]->set(i, MuNormAverage * cos(theta));
    sol->_Sol[indexMu[1]]->set(i, MuNormAverage * sin(theta));

    //sol->_Sol[indexTheta1]->set(i, theta);

  }

  for(unsigned k = 0; k < dim; k++) {
    sol->_Sol[indexMu[k]]->close();
  }
  //sol->_Sol[indexTheta1]->close();
  //END mu update
}


void AssembleShearMinimization(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< LinearImplicitSystem> ("conformal");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh *msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem *el = msh->GetMeshElements();  // pointer to the elem object in msh (level)

  MultiLevelSolution *mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution *sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object
  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

  SparseMatrix *KK = pdeSys->_KK;  // pointer to the global stiffness matrix object in pdeSys (level)
  NumericVector *RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension();

  std::vector < double > phi;  // local test function for velocity
  std::vector <adept::adouble> phi_x; // local test function first order partial derivatives
  adept::adouble weight; // gauss point weight

  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  std::vector < unsigned > solDxIndex(dim);
  solDxIndex[0] = mlSol->GetIndex("Dx1");  // get the position of "DX" in the ml_sol object
  solDxIndex[1] = mlSol->GetIndex("Dx2");  // get the position of "DY" in the ml_sol object
  if(dim == 3) solDxIndex[2] = mlSol->GetIndex("Dx3");   // get the position of "DY" in the ml_sol object

  unsigned solType;
  solType = mlSol->GetSolutionType(solDxIndex[0]);   // get the finite element type for "U"

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector < unsigned > solDxPdeIndex(dim);
  solDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex("Dx1");    // get the position of "Dx1" in the pdeSys object
  solDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex("Dx2");    // get the position of "Dx2" in the pdeSys object
  if(dim == 3) solDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex("Dx3");    // get the position of "Dx3" in the pdeSys object

  std::vector < std::vector < adept::adouble > > solDx(dim);  // local Y solution
  std::vector < std::vector < adept::adouble > > x(dim);

  std::vector< int > SYSDOF; // local to global pdeSys dofs

  std::vector< double > Res; // local redidual vector
  std::vector< adept::adouble > aResDx[dim]; // local redidual vector

  std::vector < double > Jac; // local Jacobian matrix (ordered by column, adept)

  KK->zero();  // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual

  // element loop: each process loops only on the elements that owns
  for(int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nxDofs  = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs

    for(unsigned k = 0; k < dim; k++) {
      solDx[k].resize(nxDofs);
      x[k].resize(nxDofs);
    }

    // resize local arrays
    SYSDOF.resize(dim * nxDofs);
    Res.resize(dim * nxDofs);        //resize

    for(unsigned k = 0; k < dim; k++) {
      aResDx[k].assign(nxDofs, 0.);   //resize and zet to zero
    }


    // local storage of global mapping and solution
    for(unsigned i = 0; i < nxDofs; i++) {
      // Global-to-local mapping between X solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof(i, iel, solType);
      for(unsigned k = 0; k < dim; k++) {
        solDx[k][i] = (*sol->_Sol[solDxIndex[k]])(iDDof);
        // Global-to-global mapping between NDx solution node and pdeSys dof.
        SYSDOF[ k * nxDofs + i] = pdeSys->GetSystemDof(solDxIndex[k], solDxPdeIndex[k], i, iel);
      }
    }

    // start a new recording of all the operations involving adept variables.
    s.new_recording();
    for(unsigned i = 0; i < nxDofs; i++) {
      unsigned iXDof  = msh->GetSolutionDof(i, iel, xType);
      for(unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->GetTopology()->_Sol[k])(iXDof);
      }
    }

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      msh->_finiteElement[ielGeom][solType]->Jacobian(x, ig, weight, phi, phi_x);

      std::vector < std::vector < adept::adouble > > gradSolDx(dim);

      for(unsigned  k = 0; k < dim; k++) {
        gradSolDx[k].assign(dim, 0.);
      }

      for(unsigned i = 0; i < nxDofs; i++) {
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            gradSolDx[k][j] += (x[k][i] + solDx[k][i]) * phi_x[i * dim + j];
          }
        }
      }

      for(unsigned i = 0; i < nxDofs; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          adept::adouble term = 0.;
          term  +=  phi_x[i * dim + k] * (gradSolDx[k][k]);
          aResDx[k][i] += term * weight;
        }
      }
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store


    for(int k = 0; k < dim; k++) {
      for(int i = 0; i < nxDofs; i++) {
        Res[ k * nxDofs + i] = -aResDx[k][i].value();
      }
    }

    RES->add_vector_blocked(Res, SYSDOF);

    Jac.resize((dim * nxDofs) * (dim * nxDofs));

    // define the dependent variables

    for(int k = 0; k < dim; k++) {
      s.dependent(&aResDx[k][0], nxDofs);
    }

    // define the dependent variables

    for(int k = 0; k < dim; k++) {
      s.independent(&solDx[k][0], nxDofs);
    }

    // get the jacobian matrix (ordered by row)
    s.jacobian(&Jac[0], true);

    KK->add_matrix_blocked(Jac, SYSDOF, SYSDOF);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();

  counter++;

  // ***************** END ASSEMBLY *******************
}

// Building the Conformal Minimization system.
void AssembleConformalO1Minimization(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  Extract pointers to the several objects that we are going to use.
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem> ("conformal");   // pointer to the linear implicit system named "Poisson"

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  // Pointers to the mesh (level) object and elem object in mesh (level).
  Mesh *msh = ml_prob._ml_msh->GetLevel(level);
  elem *el = msh->GetMeshElements();

  // Pointers to the multilevel solution, solution (level) and equation (level).
  MultiLevelSolution *mlSol = ml_prob._ml_sol;

  if(counter > 0 && !stopIterate) {
    UpdateMu(*mlSol);
  }

  Solution *sol = ml_prob._ml_sol->GetSolutionLevel(level);
  LinearEquationSolver *pdeSys = mlPdeSys->_LinSolver[level];

  // Pointers to global stiffness matrix and residual vector in pdeSys (level).
  SparseMatrix *KK = pdeSys->_KK;
  NumericVector *RES = pdeSys->_RES;

  // Convenience variables to keep track of the dimension.
  const unsigned  dim = 2;
  const unsigned  DIM = 3;

  // Get the process_id (for parallel computation).
  unsigned iproc = msh->processor_id();

  // Setting the reference elements to be equilateral triangles.
  std::vector < std::vector < double > > xT(2);
  xT[0].resize(3);
  xT[0][0] = -0.5;
  xT[0][1] = 0.5;
  xT[0][2] = 0.;

  xT[1].resize(3);
  xT[1][0] = 0.;
  xT[1][1] = 0.;
  xT[1][2] = sqrt(3.) / 2.;

  std::vector<double> phi_uv0;
  std::vector<double> phi_uv1;

  std::vector< double > stdVectorPhi;
  std::vector< double > stdVectorPhi_uv;

  // Extract positions of Dx in ml_sol object.
//   unsigned solDxIndex[DIM];
//   solDxIndex[0] = mlSol->GetIndex ("Dx1");
//   solDxIndex[1] = mlSol->GetIndex ("Dx2");
//   solDxIndex[2] = mlSol->GetIndex ("Dx3");

  // Extract finite element type for the solution.


  // Local solution vectors for X, Dx, Xhat, XC.
  std::vector < double > solx[DIM];
  //std::vector < double > solDx[DIM];
  std::vector < double > xhat[DIM];
  std::vector < double > xc[DIM];

  // Get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC).
  unsigned xType = 2;

  // Get the poitions of Y in the ml_sol object.
  unsigned solNDxIndex[DIM];
  solNDxIndex[0] = mlSol->GetIndex("Dx1");
  solNDxIndex[1] = mlSol->GetIndex("Dx2");
  solNDxIndex[2] = mlSol->GetIndex("Dx3");

  unsigned solType;
  solType = mlSol->GetSolutionType(solNDxIndex[0]);

  // Get the positions of Y in the pdeSys object.
  unsigned solNDxPdeIndex[DIM];
  solNDxPdeIndex[0] = mlPdeSys->GetSolPdeIndex("Dx1");
  solNDxPdeIndex[1] = mlPdeSys->GetSolPdeIndex("Dx2");
  solNDxPdeIndex[2] = mlPdeSys->GetSolPdeIndex("Dx3");

  // Local solution vectors for Nx and NDx.
  std::vector < adept::adouble > solNDx[DIM];
  std::vector < adept::adouble > solNx[DIM];

  // Get the position of "Lambda1" in the ml_sol object.
  unsigned solLIndex;
  solLIndex = mlSol->GetIndex("Lambda1");

  // Get the finite element type for "Lambda1".
  unsigned solLType;
  solLType = mlSol->GetSolutionType(solLIndex);

  // Get the position of "Lambda1" in the pdeSys object.
  unsigned solLPdeIndex;
  solLPdeIndex = mlPdeSys->GetSolPdeIndex("Lambda1");

  // Local Lambda1 solution.
  std::vector < adept::adouble > solL;


  std::vector < unsigned > solMuIndex(dim);
  solMuIndex[0] = mlSol->GetIndex("mu1");
  solMuIndex[1] = mlSol->GetIndex("mu2");
  unsigned solType1 = mlSol->GetSolutionType(solMuIndex[0]);

  std::vector < std::vector < double > > solMu(dim);

  // Local-to-global pdeSys dofs.
  std::vector < int > SYSDOF;

  // Local residual vectors.
  std::vector< double > Res;
  std::vector< adept::adouble > aResNDx[3];
  std::vector< adept::adouble > aResL;

  // Local Jacobian matrix (ordered by column).
  std::vector < double > Jac;

  KK->zero();  // Zero all the entries of the Global Matrix
  RES->zero(); // Zero all the entries of the Global Residual

  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for(int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    // Numer of solution element dofs.
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nxDofs  = msh->GetElementDofNumber(iel, solType);
    unsigned nLDofs  = msh->GetElementDofNumber(iel, solLType);
    unsigned nDofs1  = msh->GetElementDofNumber(iel, solType1);

    // Resize local arrays.
    for(unsigned K = 0; K < DIM; K++) {

      xhat[K].resize(nxDofs);

      //solDx[K].resize (nxDofs);
      solx[K].resize(nxDofs);

      solNDx[K].resize(nxDofs);
      solNx[K].resize(nxDofs);

    }
    solL.resize(nLDofs);

    for(unsigned k = 0; k < dim; k++) {
      solMu[k].resize(nDofs1);
    }

    // Resize local arrays
    SYSDOF.resize(DIM * nxDofs + nLDofs);
    Res.resize(DIM * nxDofs + nLDofs);

    for(unsigned K = 0; K < DIM; K++) {
      aResNDx[K].assign(nxDofs, 0.);
    }
    aResL.assign(nLDofs, 0.);

    // local storage of global mapping and solution
    for(unsigned i = 0; i < nxDofs; i++) {

      // Global-to-local mapping between X solution node and solution dof.
      unsigned iDDof = msh->GetSolutionDof(i, iel, solType);
      unsigned iXDof  = msh->GetSolutionDof(i, iel, xType);

      for(unsigned K = 0; K < DIM; K++) {
        xhat[K][i] = (*msh->GetTopology()->_Sol[K])(iXDof);
        //solDx[K][i] = (*sol->_Sol[solDxIndex[K]]) (iDDof);
        solx[K][i] = xhat[K][i];// + solDx[K][i];
        solNDx[K][i] = (*sol->_Sol[solNDxIndex[K]])(iDDof);

        // Global-to-global mapping between NDx solution node and pdeSys dof.
        SYSDOF[ K * nxDofs + i] =
          pdeSys->GetSystemDof(solNDxIndex[K], solNDxPdeIndex[K], i, iel);
      }
    }

    // Local storage of global mapping and solution.
    for(unsigned i = 0; i < nLDofs; i++) {

      // Global-to-local mapping between Lambda solution node and solution dof.
      unsigned iLDof = msh->GetSolutionDof(i, iel, solLType);
      solL[i] = (*sol->_Sol[solLIndex])(iLDof);

      // Global-to-global mapping between Lambda solution node and pdeSys dof.
      SYSDOF[DIM * nxDofs + i] =
        pdeSys->GetSystemDof(solLIndex, solLPdeIndex, i, iel);
    }

    for(unsigned i = 0; i < nDofs1; i++) {
      unsigned iDof = msh->GetSolutionDof(i, iel, solType1);
      for(unsigned k = 0; k < dim; k++) {
        solMu[k][i] = (*sol->_Sol[solMuIndex[k]])(iDof);
      }
    }

    // start a new recording of all the operations involving adept variables.
    s.new_recording();

    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      const double *phix;  // local test function

      const double *phix_uv[dim]; // local test function first order partial derivatives

      double weight; // gauss point weight

      // Get Gauss point weight, test function, and first order derivatives.
      if(ielGeom == QUAD) {
        phix = msh->_finiteElement[ielGeom][solType]->GetPhi(ig);


        phix_uv[0] = msh->_finiteElement[ielGeom][solType]->GetDPhiDXi(ig);
        phix_uv[1] = msh->_finiteElement[ielGeom][solType]->GetDPhiDEta(ig);

        weight = msh->_finiteElement[ielGeom][solType]->GetGaussWeight(ig);
      }

      // Special adjustments for triangles.
      else {
        msh->_finiteElement[ielGeom][solType]->Jacobian(xT, ig, weight, stdVectorPhi, stdVectorPhi_uv);

        phix = &stdVectorPhi[0];
///////// WHAT ABOUT PHIL? ////////////////////

        phi_uv0.resize(nxDofs);
        phi_uv1.resize(nxDofs);


        for(unsigned i = 0; i < nxDofs; i++) {
          phi_uv0[i] = stdVectorPhi_uv[i * dim];
          phi_uv1[i] = stdVectorPhi_uv[i * dim + 1];
        }

        phix_uv[0] = &phi_uv0[0];
        phix_uv[1] = &phi_uv1[0];

      }

      const double *phiL = msh->_finiteElement[ielGeom][solLType]->GetPhi(ig);  // local test function
      const double *phi1 = msh->_finiteElement[ielGeom][solType1]->GetPhi(ig);  // local test function

      // Initialize and compute values of x, Dx, NDx, x_uv at the Gauss points.
      //double solDxg[3] = {0., 0., 0.};
      double solNxg[3] = {0., 0., 0.};
      adept::adouble solNDxg[3] = {0., 0., 0.};

      double solx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
      adept::adouble solNx_uv[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};

      for(unsigned K = 0; K < DIM; K++) {
        for(unsigned i = 0; i < nxDofs; i++) {
          //solDxg[K] += phix[i] * solDx[K][i];
          solNDxg[K] += phix[i] * solNDx[K][i];
          solNxg[K] += phix[i] * (xhat[K][i] + solNDx[K][i].value());
        }
        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nxDofs; i++) {
            solx_uv[K][j]    += phix_uv[j][i] * solx[K][i];
            solNx_uv[K][j]   += phix_uv[j][i] * (xhat[K][i] + solNDx[K][i]);
          }
        }
      }

      ///////// ADDED THIS /////////
      adept::adouble solLg = 0.;
      for(unsigned i = 0; i < nLDofs; i++) {
        solLg += phiL[i] * solL[i];
      }

      // Compute the metric, metric determinant, and area element.
      double g[dim][dim] = {{0., 0.}, {0., 0.}};
      for(unsigned i = 0; i < dim; i++) {
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned K = 0; K < DIM; K++) {
            g[i][j] += solx_uv[K][i] * solx_uv[K][j];
          }
        }
      }
      double detg = g[0][0] * g[1][1] - g[0][1] * g[1][0];
      double Area = weight * sqrt(detg);
      double Area2 = weight; // Trick to give equal weight to each element.

      // Compute the metric inverse.
      double gi[dim][dim];
      gi[0][0] =  g[1][1] / detg;
      gi[0][1] = -g[0][1] / detg;
      gi[1][0] = -g[1][0] / detg;
      gi[1][1] =  g[0][0] / detg;

      // Compute components of the unit normal N.
      double normal[DIM];
       normal[0] = 0; //(solx_uv[1][0] * solx_uv[2][1] - solx_uv[2][0] * solx_uv[1][1]) / sqrt(detg);
       normal[1] = 0; //(solx_uv[2][0] * solx_uv[0][1] - solx_uv[0][0] * solx_uv[2][1]) / sqrt(detg);
       normal[2] = 1; //(solx_uv[0][0] * solx_uv[1][1] - solx_uv[1][0] * solx_uv[0][1]) / sqrt(detg);
      // normal[0] = 0;
      // normal[1] = solNxg[1] / sqrt(solNxg[1] * solNxg[1] + solNxg[2] * solNxg[2]);
      // normal[2] = solNxg[2] / sqrt(solNxg[1] * solNxg[1] + solNxg[2] * solNxg[2]);


      //
      // // Discretize the equation \delta CD = 0 on the basis d/du, d/dv.
      // adept::adouble V[DIM];
      // V[0] = solNx_uv[0][1] - normal[1] * solNx_uv[2][0] + normal[2] * solNx_uv[1][0];
      // V[1] = solNx_uv[1][1] - normal[2] * solNx_uv[0][0] + normal[0] * solNx_uv[2][0];
      // V[2] = solNx_uv[2][1] - normal[0] * solNx_uv[1][0] + normal[1] * solNx_uv[0][0];
      //
      // adept::adouble W[DIM];
      // W[0] = solNx_uv[0][0] + normal[1] * solNx_uv[2][1] - normal[2] * solNx_uv[1][1];
      // W[1] = solNx_uv[1][0] + normal[2] * solNx_uv[0][1] - normal[0] * solNx_uv[2][1];
      // W[2] = solNx_uv[2][0] + normal[0] * solNx_uv[1][1] - normal[1] * solNx_uv[0][1];
      //
      // adept::adouble M[DIM][dim];
      // M[0][0] = W[0] - normal[2] * V[1] + normal[1] * V[2];
      // M[1][0] = W[1] - normal[0] * V[2] + normal[2] * V[0];
      // M[2][0] = W[2] - normal[1] * V[0] + normal[0] * V[1];
      //
      // M[0][1] = V[0] + normal[2] * W[1] - normal[1] * W[2];
      // M[1][1] = V[1] + normal[0] * W[2] - normal[2] * W[0];
      // M[2][1] = V[2] + normal[1] * W[0] - normal[0] * W[1];


      double mu[2] = {0., 0.};

      for(unsigned i = 0; i < nDofs1; i++) {
        for(unsigned k = 0; k < 2; k++) {
          mu[k] += phi1[i] * solMu[k][i];
        }
      }
      double muPlus = pow(1 + mu[0], 2) + mu[1] * mu[1];
      double muMinus = pow(1 - mu[0], 2) + mu[1] * mu[1];
      double mu1Sqr = mu[1] * mu[1];
      double sqrMin = ((1 - mu[0]) * (1 - mu[0]));
      double sqrPlus = ((1 + mu[0]) * (1 + mu[0]));
      double xMin = mu[1] * (1 - mu[0]);
      double xPlus = mu[1] * (1 + mu[0]);
      double muSqr = (mu[0] * mu[0] + mu[1] * mu[1] - 1.);

      //std::cout << mu[0] << " " << mu[1] << " ";


//       if (counter == 0) mu[0] = 0.8;
//
//       for (unsigned K = 0; K < DIM; K++) {
//         if (counter > 0 && norm2Xz.value() > 0.) {
//           mu[K] += (1. / norm2Xz.value()) * XzBarXz_Bar[K].value();
//         }
//         //if(counter % 2 == 0) mu[K]*=1.01;
//         //else mu[K]/=1.01;
//       }

      //std::cout << mu[0] <<" "<< mu[1]<<" ";

      adept::adouble M2[DIM][dim];
      double muSqrPlus = 1 + mu[0] * mu[0] + mu[1] * mu[1];
      double muSqrMinus = mu[0] * mu[0] + mu[1] * mu[1] - 1;

      M2[0][0] = muSqrPlus * solNx_uv[0][0] + muSqrMinus * solNx_uv[1][1];
      M2[1][0] = muSqrPlus * solNx_uv[1][0] - muSqrMinus * solNx_uv[0][1];
      M2[0][1] = muSqrPlus * solNx_uv[0][1] - muSqrMinus * solNx_uv[1][0];
      M2[1][1] = muSqrPlus * solNx_uv[1][1] + muSqrMinus * solNx_uv[0][0];


      adept::adouble M1[DIM][dim];

      M1[0][0] = (muMinus + muPlus * (normal[1] * normal[1] + normal[2] * normal[2])) * solNx_uv[0][0]
                 - muPlus * normal[0] * (solNx_uv[1][0] * normal[1] + solNx_uv[2][0] * normal[2])
                 - 2 * muSqr * (solNx_uv[2][1] * normal[1] - solNx_uv[1][1] * normal[2]);

      M1[1][0] = (muMinus + muPlus * (normal[2] * normal[2] + normal[0] * normal[0])) * solNx_uv[1][0]
                 - muPlus * normal[1] * (solNx_uv[2][0] * normal[2] + solNx_uv[0][0] * normal[0])
                 - 2 * muSqr * (solNx_uv[0][1] * normal[2] - solNx_uv[2][1] * normal[0]);

      M1[2][0] = (muMinus + muPlus * (normal[0] * normal[0] + normal[1] * normal[1])) * solNx_uv[2][0]
                 - muPlus * normal[2] * (solNx_uv[0][0] * normal[0] + solNx_uv[1][0] * normal[1])
                 - 2 * muSqr * (solNx_uv[1][1] * normal[0] - solNx_uv[0][1] * normal[1]);

      M1[0][1] = (muMinus + muPlus * (normal[1] * normal[1] + normal[2] * normal[2])) * solNx_uv[0][1]
                 - muPlus * normal[0] * (solNx_uv[2][1] * normal[2] + solNx_uv[1][1] * normal[1])
                 + 2 * muSqr * (solNx_uv[2][0] * normal[1] - solNx_uv[1][0] * normal[2]);

      M1[1][1] = (muMinus + muPlus * (normal[2] * normal[2] + normal[0] * normal[0])) * solNx_uv[1][1]
                 - muPlus * normal[1] * (solNx_uv[0][1] * normal[0] + solNx_uv[2][1] * normal[2])
                 + 2 * muSqr * (solNx_uv[0][0] * normal[2] - solNx_uv[2][0] * normal[0]);

      M1[2][1] = (muMinus + muPlus * (normal[0] * normal[0] + normal[1] * normal[1])) * solNx_uv[2][1]
                 - muPlus * normal[2] * (solNx_uv[1][1] * normal[1] + solNx_uv[0][1] * normal[0])
                 + 2 * muSqr * (solNx_uv[1][0] * normal[0] - solNx_uv[0][0] * normal[1]);


      adept::adouble N[DIM][dim];

      N[0][0] = (+ (sqrMin + mu1Sqr  * (normal[1] * normal[1] + normal[2] * normal[2])) * solNx_uv[0][0]
                 - (xMin   + xPlus   * (normal[1] * normal[1] + normal[2] * normal[2])) * solNx_uv[0][1])
                + (- (solNx_uv[1][0] * normal[1] + solNx_uv[2][0] * normal[2]) * mu1Sqr * normal[0]
                   + (solNx_uv[1][1] * normal[1] + solNx_uv[2][1] * normal[2]) * xPlus  * normal[0]
                   + (solNx_uv[1][1] * normal[2] - solNx_uv[2][1] * normal[1]) * muSqr);

      N[1][0] = (+ (sqrMin + mu1Sqr  * (normal[2] * normal[2] + normal[0] * normal[0])) * solNx_uv[1][0]
                 - (xMin   + xPlus   * (normal[2] * normal[2] + normal[0] * normal[0])) * solNx_uv[1][1])
                + (- (solNx_uv[2][0] * normal[2] + solNx_uv[0][0] * normal[0]) * mu1Sqr * normal[1]
                   + (solNx_uv[2][1] * normal[2] + solNx_uv[0][1] * normal[0]) * xPlus  * normal[1]
                   + (solNx_uv[2][1] * normal[0] - solNx_uv[0][1] * normal[2]) * muSqr);

      N[2][0] = (+ (sqrMin + mu1Sqr  * (normal[0] * normal[0] + normal[1] * normal[1])) * solNx_uv[2][0]
                 - (xMin   + xPlus   * (normal[0] * normal[0] + normal[1] * normal[1])) * solNx_uv[2][1])
                + (- (solNx_uv[0][0] * normal[0] + solNx_uv[1][0] * normal[1]) * mu1Sqr * normal[2]
                   + (solNx_uv[0][1] * normal[0] + solNx_uv[1][1] * normal[1]) * xPlus  * normal[2]
                   + (solNx_uv[0][1] * normal[1] - solNx_uv[1][1] * normal[0]) * muSqr);


      N[0][1] = (+ (mu1Sqr + sqrPlus * (normal[1] * normal[1] + normal[2] * normal[2])) * solNx_uv[0][1]
                 - (xMin   + xPlus   * (normal[1] * normal[1] + normal[2] * normal[2])) * solNx_uv[0][0])
                + (- (solNx_uv[1][1] * normal[1] + solNx_uv[2][1] * normal[2]) * sqrPlus * normal[0]
                   + (solNx_uv[1][0] * normal[1] + solNx_uv[2][0] * normal[2]) * xPlus   * normal[0]
                   - (solNx_uv[1][0] * normal[2] - solNx_uv[2][0] * normal[1]) * muSqr);

      N[1][1] = (+ (mu1Sqr + sqrPlus * (normal[2] * normal[2] + normal[0] * normal[0])) * solNx_uv[1][1]
                 - (xMin   + xPlus   * (normal[2] * normal[2] + normal[0] * normal[0])) * solNx_uv[1][0])
                + (- (solNx_uv[2][1] * normal[2] + solNx_uv[0][1] * normal[0]) * sqrPlus * normal[1]
                   + (solNx_uv[2][0] * normal[2] + solNx_uv[0][0] * normal[0]) * xPlus   * normal[1]
                   - (solNx_uv[2][0] * normal[0] - solNx_uv[0][0] * normal[2]) * muSqr);

      N[2][1] = (+ (mu1Sqr + sqrPlus * (normal[0] * normal[0] + normal[1] * normal[1])) * solNx_uv[2][1]
                 - (xMin   + xPlus   * (normal[0] * normal[0] + normal[1] * normal[1])) * solNx_uv[2][0])
                + (- (solNx_uv[0][1] * normal[0] + solNx_uv[1][1] * normal[1]) * sqrPlus * normal[2]
                   + (solNx_uv[0][0] * normal[0] + solNx_uv[1][0] * normal[1]) * xPlus   * normal[2]
                   - (solNx_uv[0][0] * normal[1] - solNx_uv[1][0] * normal[0]) * muSqr);


      adept::adouble V[DIM];
      V[0] = (1 - mu[0]) * solNx_uv[0][0] - (1 + mu[0]) * solNx_uv[1][1] + mu[1] * (solNx_uv[1][0] - solNx_uv[0][1]);
      V[1] = (1 - mu[0]) * solNx_uv[1][0] + (1 + mu[0]) * solNx_uv[0][1] - mu[1] * (solNx_uv[0][0] + solNx_uv[1][1]);
      V[2] = 0;


      adept::adouble M[DIM][dim];

      M[0][0] = (1 - mu[0]) * V[0] - mu[1] * V[1];
      M[1][0] = (1 - mu[0]) * V[1] + mu[1] * V[0];
      M[2][0] = 0;
      //M[0][0] = (1 - mu1) * V[0] - mu2 * V[1];
      //M[1][0] = (1 - mu1) * V[1] + mu2 * V[0];

      M[0][1] = (1 + mu[0]) * V[1] - mu[1] * V[0];
      M[1][1] = - (1 + mu[0]) * V[0] - mu[1] * V[1];
      M[2][1] = 0;

      // std::cout << 00 << " "<< M1[0][0]<<std::endl;

      if(iel == 4 && ig == 1) {

        std::cout << "derivatives " << std::endl;
        std::cout << 1 << " " << solNx_uv[0][0] << std::endl;
        std::cout << 2 << " " << solNx_uv[1][1] << std::endl;
        std::cout << 3 << " " << solNx_uv[1][0] << std::endl;
        std::cout << 4 << " " << solNx_uv[0][1] << std::endl;

        std::cout << "symmetric " << std::endl;
        std::cout << "00" << " " << M1[0][0] << " " << M2[0][0] << std::endl;
        std::cout << "01" << " " << M1[1][0] << " " << M2[1][0] << std::endl;
        std::cout << "01" << " " << M1[0][1] << " " << M2[0][1] << std::endl;
        std::cout << "11" << " " << M1[1][1] << " " << M2[1][1] << std::endl;

        std::cout << "asymmetric " << std::endl;
        std::cout << "00" << " " << N[0][0] << " " << M[0][0] << std::endl;
        std::cout << "01" << " " << N[1][0] << " " << M[1][0] << std::endl;
        std::cout << "01" << " " << N[0][1] << " " << M[0][1] << std::endl;
        std::cout << "11" << " " << N[1][1] << " " << M[1][1] << std::endl;

//         std::cout << "00" << " " << N[0][0] << " " << M[0][0] << std::endl;
//         std::cout << "10" << " " << N[1][0] << " " << M[1][0] << std::endl;
//         //std::cout << "20" << " "<< N[2][0]<<" "<< M[2][0]<<std::endl;
//         std::cout << "01" << " " << N[0][1] << " " << M[0][1] << std::endl;
//         std::cout << "11" << " " << N[1][1] << " " << M[1][1] << std::endl;
//         //std::cout << "21" << " "<< N[2][1]<<" "<< M[2][1]<<std::endl;

      }



//       adept::adouble Q[DIM][dim];
//       Q[0][0] = (gi[1][1] * W[0]
//                  + gi[0][0] * (normal[1] * V[2] - normal[2] * V[1])
//                  + gi[0][1] * (normal[2] * W[1] - normal[1] * W[2] - V[0]));
//
//       Q[1][0] = (gi[1][1] * W[1]
//                  + gi[0][0] * (normal[2] * V[0] - normal[0] * V[2])
//                  + gi[0][1] * (normal[0] * W[2] - normal[2] * W[0] - V[1]));
//
//       Q[2][0] = (gi[1][1] * W[2]
//                  + gi[0][0] * (normal[0] * V[1] - normal[1] * V[0])
//                  + gi[0][1] * (normal[1] * W[0] - normal[0] * W[1] - V[2]));
//
//       Q[0][1] = (gi[0][0] * V[0]
//                  + gi[1][1] * (normal[2] * W[1] - normal[1] * W[2])
//                  + gi[0][1] * (normal[1] * V[2] - normal[2] * V[1] - W[0]));
//
//       Q[1][1] = (gi[0][0] * V[1]
//                  + gi[1][1] * (normal[0] * W[2] - normal[2] * W[0])
//                  + gi[0][1] * (normal[2] * V[0] - normal[0] * V[2] - W[1]));
//
//       Q[2][1] = (gi[0][0] * V[2]
//                  + gi[1][1] * (normal[1] * W[0] - normal[0] * W[1])
//                  + gi[0][1] * (normal[0] * V[1] - normal[1] * V[0] - W[2]));


      // Compute new X minus old X dot N, for "reparametrization".
      adept::adouble DnXmDxdotN = 0.;
      for(unsigned K = 0; K < DIM; K++) {
        DnXmDxdotN += (/*solDxg[K]*/ - solNDxg[K]) * normal[K];
      }


      // Implement the Conformal Minimization equations.
      for(unsigned K = 0; K < DIM; K++) {
        for(unsigned i = 0; i < nxDofs; i++) {
          adept::adouble term1 = 0.;

          for(unsigned j = 0; j < dim; j++) {
            term1 += N[K][j] * phix_uv[j][i]; //asymmetric
            //term1 += M1[K][j] * phix_uv[j][i]; //symmetric
            //term1 += Q[K][j] * phix_uv[j][i];
          }

          // Conformal energy equation (with trick).
          aResNDx[K][i] += term1 * Area2
                           // + timederiv * (solNDxg[K] - solDxg[K]) * phix[i] * Area2
                           + solL[0] * phix[i] * normal[K] * Area; //2 occasionally better???
        }
      }

      // Lagrange multiplier equation (with trick).
      for(unsigned i = 0; i < nLDofs; i++) {
        aResL[i] += phiL[i] * (DnXmDxdotN + eps * solL[i]) * Area; // no2
      }
      //aResL[0] += (DnXmDxdotN + eps * solL[0]) * Area;

      if(iel == 4 && ig == 1) {
        for(unsigned i = 0; i < nxDofs; i++) {
          std::cout <<  mu[0] << " " << mu[1] << " " << aResNDx[0][i] << " " << aResNDx[1][i] << "\n";
        }
      }

    } // end GAUSS POINT LOOP

    //------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
    //copy the value of the adept::adoube aRes in double Res and store

    for(int K = 0; K < DIM; K++) {
      for(int i = 0; i < nxDofs; i++) {
        Res[ K * nxDofs + i] = -aResNDx[K][i].value();
      }
    }

    for(int i = 0; i < nLDofs; i++) {
      Res[DIM * nxDofs + i] = - aResL[i].value();
    }
    RES->add_vector_blocked(Res, SYSDOF);

    // Resize Jacobian.
    Jac.resize((DIM * nxDofs + nLDofs) * (DIM * nxDofs + nLDofs));

    // Define the dependent variables.
    for(int K = 0; K < DIM; K++) {
      s.dependent(&aResNDx[K][0], nxDofs);
    }
    s.dependent(&aResL[0], nLDofs);


    // Define the independent variables.
    for(int K = 0; K < DIM; K++) {
      s.independent(&solNDx[K][0], nxDofs);
    }
    s.independent(&solL[0], nLDofs);

    // Get the jacobian matrix (ordered by row).
    s.jacobian(&Jac[0], true);

    KK->add_matrix_blocked(Jac, SYSDOF, SYSDOF);

    s.clear_independents();
    s.clear_dependents();

  } //end ELEMENT LOOP for each process.

  RES->close();
  KK->close();

  counter++;

} // end AssembleConformalMinimization.
