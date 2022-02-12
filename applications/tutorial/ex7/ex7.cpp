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
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"


using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;

  if (!strcmp(SolName, "T")) {
    if (facename == 2) {
      value = 1.;
    } else if (facename == 3) {
      dirichlet = false; //Neumann
    }
  } else if (!strcmp(SolName, "P")) {
    dirichlet = false;
  }

  return dirichlet;
}


bool SetRefinementFlag(const std::vector < double >& x, const int& elemgroupnumber, const int& level) {

  bool refine = 0;

  if (elemgroupnumber == 6 && level < 4) refine = 1;

  if (elemgroupnumber == 7 && level < 5) refine = 1;

  if (elemgroupnumber == 8 && level < 6) refine = 1;

  return refine;

}



void AssembleBoussinesqAppoximation_AD(MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  //mlMsh.ReadCoarseMesh("./input/cube_hex.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/square_quad.neu", "seventh", scalingFactor);
  mlMsh.ReadCoarseMesh("./input/quadAMR.neu", "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

//   unsigned numberOfUniformLevels = 3;
//   unsigned numberOfSelectiveLevels = 0;
//   mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  unsigned numberOfUniformLevels = 4;
  unsigned numberOfSelectiveLevels = 3;
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);


  // erase all the coarse mesh levels
  //mlMsh.EraseCoarseLevels(numberOfUniformLevels - 3);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("T", LAGRANGE, FIRST);
  mlSol.AddSolution("U", LAGRANGE, SECOND);
  mlSol.AddSolution("V", LAGRANGE, SECOND);

  if (dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND);

  //mlSol.AddSolution("P", LAGRANGE, FIRST);
  mlSol.AddSolution("P",  DISCONTINUOUS_POLYNOMIAL, FIRST);

  mlSol.AssociatePropertyToSolution("P", "Pressure");
  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.FixSolutionAtOnePoint("P");
  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  // add system Poisson in mlProb as a Linear Implicit System
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("NS");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("T");
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");

  if (dim == 3) system.AddSolutionToSystemPDE("W");

  system.AddSolutionToSystemPDE("P");

  //system.SetLinearEquationSolverType(FEMuS_DEFAULT);
  system.SetLinearEquationSolverType(FEMuS_ASM); // Additive Swartz Method
  // attach the assembling function to system
  system.SetAssembleFunction(AssembleBoussinesqAppoximation_AD);

  system.SetMaxNumberOfNonLinearIterations(10);
  system.SetMaxNumberOfLinearIterations(3);
  system.SetAbsoluteLinearConvergenceTolerance(1.e-12);
  system.SetNonLinearConvergenceTolerance(1.e-8);
  system.SetMgType(F_CYCLE);

  system.SetNumberPreSmoothingStep(0);
  system.SetNumberPostSmoothingStep(2);
  // initilaize and solve the system
  system.init();

  system.SetSolverFineGrids(GMRES);
  system.SetPreconditionerFineGrids(ILU_PRECOND);
  //system.SetTolerances(1.e-20, 1.e-20, 1.e+50, 40);
  system.SetTolerances(1.e-3, 1.e-20, 1.e+50, 5);


  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");
  system.SetNumberOfSchurVariables(1);
  system.SetElementBlockNumber(4);
  //system.SetDirichletBCsHandling(ELIMINATION);
  //system.solve();
  system.MGsolve();

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  GMVWriter gmvIO(&mlSol);
  variablesToBePrinted.push_back("all");
  gmvIO.SetDebugOutput(true);
  gmvIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);


  return 0;
}


void AssembleBoussinesqAppoximation_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*           msh         = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*           el          = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*   mlSol         = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*   sol         = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level];  // pointer to the equation (level) object
  SparseMatrix*   KK          = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*  RES         = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //solution variable
  unsigned solTIndex;
  solTIndex = mlSol->GetIndex("T");    // get the position of "T" in the ml_sol object
  unsigned solTType = mlSol->GetSolutionType(solTIndex);    // get the finite element type for "T"

  vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object

  if (dim == 3) solVIndex[2] = mlSol->GetIndex("W");      // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"

  unsigned solTPdeIndex;
  solTPdeIndex = mlPdeSys->GetSolPdeIndex("T");    // get the position of "T" in the pdeSys object

  vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object

  if (dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W");

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object

  vector < adept::adouble >  solT; // local solution
  vector < vector < adept::adouble > >  solV(dim);    // local solution
  vector < adept::adouble >  solP; // local solution

  vector< adept::adouble > aResT; // local redidual vector
  vector< vector < adept::adouble > > aResV(dim);    // local redidual vector
  vector< adept::adouble > aResP; // local redidual vector

  vector < vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  solT.reserve(maxSize);
  aResT.reserve(maxSize);

  for (unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(maxSize);
    aResV[k].reserve(maxSize);
    coordX[k].reserve(maxSize);
  }

  solP.reserve(maxSize);
  aResP.reserve(maxSize);


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

  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsT = msh->GetElementDofNumber(iel, solTType);    // number of solution element dofs
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs

    unsigned nDofsTVP = nDofsT + dim * nDofsV + nDofsP;
    // resize local arrays
    sysDof.resize(nDofsTVP);

    solT.resize(nDofsV);

    for (unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsX);
    }

    solP.resize(nDofsP);

    aResT.resize(nDofsV);    //resize
    std::fill(aResT.begin(), aResT.end(), 0);    //set aRes to zero

    for (unsigned  k = 0; k < dim; k++) {
      aResV[k].resize(nDofsV);    //resize
      std::fill(aResV[k].begin(), aResV[k].end(), 0);    //set aRes to zero
    }

    aResP.resize(nDofsP);    //resize
    std::fill(aResP.begin(), aResP.end(), 0);    //set aRes to zero

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofsT; i++) {
      unsigned solTDof = msh->GetSolutionDof(i, iel, solTType);    // global to global mapping between solution node and solution dof
      solT[i] = (*sol->_Sol[solTIndex])(solTDof);      // global extraction and local storage for the solution
      sysDof[i] = pdeSys->GetSystemDof(solTIndex, solTPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dofs
    }

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
        sysDof[i + nDofsT + k * nDofsV] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    for (unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // global to global mapping between solution node and solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);      // global extraction and local storage for the solution
      sysDof[i + nDofsT + dim * nDofsV] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    // if (level == levelMax || !el->GetRefinedElementIndex(iel)) {      // do not care about this if now (it is used for the AMR)
    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solTType]->Jacobian(coordX, ig, weight, phiT, phiT_x, phiT_xx);
      msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV, phiV_x, phiV_xx);
      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      adept::adouble solT_gss = 0;
      vector < adept::adouble > gradSolT_gss(dim, 0.);

      for (unsigned i = 0; i < nDofsT; i++) {
        solT_gss += phiT[i] * solT[i];

        for (unsigned j = 0; j < dim; j++) {
          gradSolT_gss[j] += phiT_x[i * dim + j] * solT[i];
        }
      }

      vector < adept::adouble > solV_gss(dim, 0);
      vector < vector < adept::adouble > > gradSolV_gss(dim);

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

      adept::adouble solP_gss = 0;

      for (unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }

      double nu = 1.;
      double alpha = 1.;
      double beta = 2000.;

      // *** phiT_i loop ***
      for (unsigned i = 0; i < nDofsT; i++) {
        adept::adouble Temp = 0.;

        for (unsigned j = 0; j < dim; j++) {
          Temp +=  alpha * phiT_x[i * dim + j] * gradSolT_gss[j];
          Temp +=  phiT[i] * (solV_gss[j] * gradSolT_gss[j]);
        }

        aResT[i] += - Temp * weight;
      } // end phiT_i loop


      // *** phiV_i loop ***
      for (unsigned i = 0; i < nDofsV; i++) {
        vector < adept::adouble > NSV(dim, 0.);

        for (unsigned j = 0; j < dim; j++) {
          for (unsigned  k = 0; k < dim; k++) {
            NSV[k]   +=  nu * phiV_x[i * dim + j] * (gradSolV_gss[k][j] + gradSolV_gss[j][k]);
            NSV[k]   +=  phiV[i] * (solV_gss[j] * gradSolV_gss[k][j]);
          }
        }

        for (unsigned  k = 0; k < dim; k++) {
          NSV[k] += -solP_gss * phiV_x[i * dim + k];
        }

        NSV[1] += -beta * solT_gss * phiV[i];

        for (unsigned  k = 0; k < dim; k++) {
          aResV[k][i] += - NSV[k] * weight;
        }
      } // end phiV_i loop

      // *** phiP_i loop ***
      for (unsigned i = 0; i < nDofsP; i++) {
        for (int k = 0; k < dim; k++) {
          aResP[i] += - (gradSolV_gss[k][k]) * phiP[i]  * weight;
        }
      } // end phiP_i loop

    } // end gauss point loop

    // } // endif single element not refined or fine grid loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store them in RES
    Res.resize(nDofsTVP);    //resize

    for (int i = 0; i < nDofsT; i++) {
      Res[i] = -aResT[i].value();
    }

    for (int i = 0; i < nDofsV; i++) {
      for (unsigned  k = 0; k < dim; k++) {
        Res[ i + nDofsT + k * nDofsV ] = -aResV[k][i].value();
      }
    }

    for (int i = 0; i < nDofsP; i++) {
      Res[ i + nDofsT + dim * nDofsV ] = -aResP[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);

    //Extarct and store the Jacobian
    Jac.resize(nDofsTVP * nDofsTVP);
    // define the dependent variables
    s.dependent(&aResT[0], nDofsT);

    for (unsigned  k = 0; k < dim; k++) {
      s.dependent(&aResV[k][0], nDofsV);
    }

    s.dependent(&aResP[0], nDofsP);

    // define the independent variables
    s.independent(&solT[0], nDofsT);

    for (unsigned  k = 0; k < dim; k++) {
      s.independent(&solV[k][0], nDofsV);
    }

    s.independent(&solP[0], nDofsP);

    // get the and store jacobian matrix (row-major)
    s.jacobian(&Jac[0] , true);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();

  KK->close();

  // ***************** END ASSEMBLY *******************
}




