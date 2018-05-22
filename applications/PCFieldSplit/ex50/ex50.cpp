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
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"
#include "FieldSplitTree.hpp"
#include <stdlib.h>


double Re = 5000.0;
double Rm = 0.1;
double coefS = -1.0; 

using namespace femus;

double InitalValueT(const std::vector < double >& x){
  return (x[0]+0.5);
}

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
/*
  if(!strcmp(SolName, "T")) {
    if(facename == 2) {
      value = 1.;
    } else if(facename == 3) {
      dirichlet = false; //Neumann
    }
  } else if(!strcmp(SolName, "P")) {
    dirichlet = false;
  }*/
  
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
  else if {}
  
  
  return dirichlet;
}


void PrintConvergenceInfo(char *stdOutfile, char* infile, const unsigned &numofrefinements);
void AssembleBoussinesqAppoximation(MultiLevelProblem& ml_prob);

// enum PrecType {
//   FS_VTp = 1,
//   FS_TVp,
//   ASM_VTp,
//   ASM_TVp,
//   ILU_VTp,
//   ILU_TVp
// };

int main(int argc, char** args) {

//   unsigned precType = 0;
// 
//   if(argc >= 2) {
//     if(!strcmp("FS_VT", args[1])) precType = FS_VTp;
//     else if(!strcmp("FS_TV", args[1])) precType = FS_TVp;
//     else if(!strcmp("ASM_VT", args[1])) precType = ASM_VTp;
//     else if(!strcmp("ASM_TV", args[1])) precType = ASM_TVp;
//     else if(!strcmp("ILU_VT", args[1])) precType = ILU_VTp;
// 
//     if(!strcmp("ILU_TV", args[1])) precType = ILU_TVp;
// 
//     if(precType == 0) {
//       std::cout << "wrong input arguments!" << std::endl;
//       abort();
//     }
//   }
//   else {
//     std::cout << "No input argument set default preconditioner = NS+T" << std::endl;
//     precType = FS_VTp;
//   }
//   
//   if(argc >= 3) {
//     Prandtl = strtod(args[2], NULL);
//     std::cout << Prandtl<<std::endl;
//   }
//   
//   
//   if(argc >= 4) {
//     Rayleigh = strtod(args[3], NULL);
//     std::cout << Rayleigh <<std::endl;
//   }
    
  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  //mlMsh.ReadCoarseMesh("./input/cube_hex.neu","seventh",scalingFactor); 
  mlMsh.ReadCoarseMesh("./input/square_quad.neu", "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 6;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(3);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("B1",LAGRANGE, SECOND);
  mlSol.Addsoluiton{"B2",LAGRANGE, SECOND}; 
  mlSol.Addsolution{"R", DISCONTINOUS_POLYNOMIAL, FIRST);
  mlSol.AssociatePropertyToSolution("R", "Pressure");
  
  mlSol.AddSolution("U", LAGRANGE, SECOND);
  mlSol.AddSolution("V", LAGRANGE, SECOND);
  mlSol.AddSolution("P",  DISCONTINOUS_POLYNOMIAL, FIRST);
  mlSol.AssociatePropertyToSolution("P", "Pressure");
  mlSol.Initialize("All");
  
//   mlSol.Initialize("T",InitalValueT);
  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.FixSolutionAtOnePoint("R");
  mlSol.FixSolutionAtOnePoint("P");
  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);
  // add system Poisson in mlProb as a Linear Implicit System
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("NS");

  system.AddSolutionToSystemPDE("B1");
  system.AddSolutionToSystemPDE("B2");
  system.AddSolutionToSystemPDE("R");
  
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  system.AddSolutionToSystemPDE("P");

  //BEGIN buid fieldSplitTree (only for FieldSplitPreconditioner)
  std::vector <unsigned> fieldB1B2R(3);
  fieldB1B2R[0] = system.GetSolPdeIndex("B1");
  fieldB1B2R[1] = system.GetSolPdeIndex("B2");
  fieldB1B2R[2] = system.GetSolPdeIndex("R");
  
  std::vector <unsigned> solutionTypeB1B2R(3);
  solutionTypeB1B2R[0] = mlSol.GetSolutionType("B1");
  solutionTypeB1B2R[1] = mlSol.GetSolutionType("B2");
  solutionTypeB1B2R[2] = mlSol.GetSolutionType("R");
  
  FieldSplitTree FS_MF(PREONLY, ASM_PRECOND, fieldB1B2R, solutionTypeB1B2R, "Magnetic Field");
  FS_MF.SetASMBlockSize(4);
  FS_MF.SetAsmNumeberOfSchurVariables(1);
  
  std::vector < unsigned > fieldUVP(3);
  fieldUVP[0] = system.GetSolPdeIndex("U");
  fieldUVP[1] = system.GetSolPdeIndex("V");
  fieldUVP[2] = system.GetSolPdeIndex("P");

  std::vector < unsigned > solutionTypeUVP(3);
  solutionTypeUVP[0] = mlSol.GetSolutionType("U");
  solutionTypeUVP[1] = mlSol.GetSolutionType("V");
  solutionTypeUVP[2] = mlSol.GetSolutionType("P");

  FieldSplitTree FS_NS(PREONLY, ASM_PRECOND, fieldUVP, solutionTypeUVP, "Navier-Stokes");
  FS_NS.SetAsmBlockSize(4);
  FS_NS.SetAsmNumeberOfSchurVariables(1);

  std::vector < FieldSplitTree *> FS2;
  FS2.reserve(2);
  FS2.push_back(&FS_NS); // Navier-stokes block first
  FS2.push_back(&FS_MF); // Magnetic Field

  FieldSplitTree FS_MHD(RICHARDSON, FIELDSPLIT_PRECOND, FS2, "MHD");

  //END buid fieldSplitTree
  system.SetMgSmoother(FIELDSPLIT_SMOOTHER); // Field-Split preconditioner
  // system.SetMgSmoother(ASM_SMOOTHER);  // Additive Swartz preconditioner
  // ILU preconditioner 

  // attach the assembling function to system
  system.SetAssembleFunction(AssembleBoussinesqAppoximation);

  system.SetMaxNumberOfNonLinearIterations(10);
  system.SetNonLinearConvergenceTolerance(1.e-8);
  //system.SetMaxNumberOfResidualUpdatesForNonlinearIteration(2);
  //system.SetResidualUpdateConvergenceTolerance(1.e-12);

  //system.SetMaxNumberOfLinearIterations(10);
  //system.SetAbsoluteLinearConvergenceTolerance(1.e-15);
  
  system.SetMaxNumberOfLinearIterations(1);
  system.SetAbsoluteLinearConvergenceTolerance(1.e-15);
  

  system.SetMgType(F_CYCLE);

  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);
  // initilaize and solve the system
  system.init();

  system.SetSolverFineGrids(RICHARDSON);
  system.SetPreconditionerFineGrids(ILU_PRECOND);
  system.SetRichardsonScaleFactor(.6);
  
  system.SetFieldSplitTree(&FS_MHD);
  system.SetTolerances(1.e-5, 1.e-8, 1.e+50, 30, 30); //GMRES tolerances
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All")
  system.MGsolve();

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  mlMsh.PrintInfo();
  
  //char infile[256]="FSVTtrueResidual.txt";
  //char stdOutfile[256]="output.txt"; 

  return 0;
}

void AssembleBoussinesqAppoximation(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use
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
  
  unsigned solBIndex(dim);
  solBIndex[0] = mlSol->GetIndex("B1");
  solBIndex[1] = mlSol->GetIndex("B2");
  unsigned solBType = mlSol -> GetSolutionType(solBIndex[0]);  
  
  unsigned solRindex;
  solRindex = mlSol -> GetIndex{"R"};
  unsigned solRType = mlSol -> GetSolutionType(solRIndex);

  vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object
  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"

  vector < unsigned > solBPdeIndex(dim);
  solBPdeIndex[0] = mlPdeSys->GetSolPdeIndex("B1");    // get the position of "B1" in the pdeSys object
  solBPdeIndex[1] = mlPdeSys->GetSolPdeIndex("B2");    // get the position of "B2" in the pdeSys object
  unsigned solRPdeIndex;
  solRPdeIndex = mlPdeSys->GetSolPdeIndex("R");    // get the position of "R" in the pdeSys object

  vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object
  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object

  vector <vector<double>> solB(dim);
  vector <double> solR; 
  vector < vector < double > >  solV(dim);    // local solution
  vector < double >  solP; // local solution

  vector < vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)


  for(unsigned  k = 0; k < dim; k++) {
    solB[k].reserve(maxSize);
    coordX[k].reserve(maxSize);
  }
  solR.reserve(maxSize);
  
  for(unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(maxSize);
  }
  solP.reserve(maxSize);

  vector <double> phiB;  // local test function
  vector <double> phiB_x; // local test function first order partial derivatives
  vector <double> phiB_xx; // local test function second order partial derivatives

  phiB.reserve(maxSize);
  phiB_x.reserve(maxSize * dim);
  phiB_xx.reserve(maxSize * dim2);
  double* phiR;
  
  vector <double> phiV;  // local test function
  vector <double> phiV_x; // local test function first order partial derivatives
  vector <double> phiV_xx; // local test function second order partial derivatives

  phiV.reserve(maxSize);
  phiV_x.reserve(maxSize * dim);
  phiV_xx.reserve(maxSize * dim2);

//   vector <double> phiT;  // local test function
//   vector <double> phiT_x; // local test function first order partial derivatives
//   vector <double> phiT_xx; // local test function second order partial derivatives
// 
//   phiT.reserve(maxSize);
//   phiT_x.reserve(maxSize * dim);
//   phiT_xx.reserve(maxSize * dim2);

  double* phiP;
  double weight; // gauss point weight

  vector< int > sysDof; // local to global pdeSys dofs
  sysDof.reserve((dim + 2) *maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve((dim + 2) *maxSize);

  vector < double > Jac;
  Jac.reserve((dim + 2) *maxSize * (dim + 2) *maxSize);

  if(assembleMatrix)
    KK->zero(); // Set to zero all the entries of the Global Matrix

  //BEGIN element loop
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    //BEGIN local dof number extraction
    unsigned nDofsB = msh->GetElementDofNumber(iel, solBType);  //Magnetic 
    unsigned nDofsR = msh->GetElementDofNumber(iel, solRType);  //pressure
    
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);  //velocity
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);  //pressure
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType); // coordinates
    unsigned nDofsBRVP = dim * nDofsB + nDofsR + dim * nDofsV + nDofsP; // all solutions
    //END local dof number extraction

    //BEGIN memory allocation
    Res.resize(nDofsBRVP);
    std::fill(Res.begin(), Res.end(), 0);
    Jac.resize(nDofsBRVP * nDofsBRVP);
    std::fill(Jac.begin(), Jac.end(), 0);
    sysDof.resize(nDofsBRVP);

    for(unsigned  k = 0; k < dim; k++) {
      solB[k].resize(nDofsB);
      coordX[k].resize(nDofsX);
    }
    solR.resize(nDofsR);

    for(unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
    }
    solP.resize(nDofsP);
    //END memory allocation

    //BEGIN global to local extraction
     for(unsigned i = 0; i < nDofsB; i++) { //velocity
     unsigned solBDof = msh->GetSolutionDof(i, iel, solBType);  //local to global solution dof

     for(unsigned  k = 0; k < dim; k++) {
	solB[k][i] = (*sol->_Sol[solBIndex[k]])(solBDof);  //global to local solution value
	sysDof[i + k * nDofsB] = pdeSys->GetSystemDof(solBIndex[k], solBPdeIndex[k], i, iel);  //local to global system dof
     }
    }
    
    for(unsigned i = 0; i < nDofsR; i++) { //pressure
      unsigned solRDof = msh->GetSolutionDof(i, iel, solRType);  //local to global solution dof
      solR[i] = (*sol->_Sol[solRIndex])(solRDof);  //global to local solution value
      sysDof[i + dim * nDofsB] = pdeSys->GetSystemDof(solRIndex, solRPdeIndex, i, iel);  //local to global system dof
    }
    

    for(unsigned i = 0; i < nDofsV; i++) { //velocity
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);  //local to global solution dof

      for(unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);  //global to local solution value
        sysDof[i + dim * nDofsB + nDofsR + k * nDofsV] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);  //local to global system dof
      }
    }

    for(unsigned i = 0; i < nDofsP; i++) { //pressure
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);  //local to global solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);  //global to local solution value
      sysDof[i + dim * nDofsB + nDofsR + dim * nDofsV] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);  //local to global system dof
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

    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiB, phiB_x, phiB_xx);
      phiR = msh->_finiteElement[ielGeom][solRType]->GetPhi(ig);
      
      msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV, phiV_x, phiV_xx);
      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      vector < double > solB_gss(dim, 0);
      vector < vector < double > > gradSolB_gss(dim);

      for(unsigned  k = 0; k < dim; k++) {
        gradSolB_gss[k].resize(dim);
        std::fill(gradSolB_gss[k].begin(), gradSolB_gss[k].end(), 0);
      }

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solB_gss[k] += phiB[i] * solB[k][i];
        }

        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            gradSolB_gss[k][j] += phiB_x[i * dim + j] * solB[k][i];
          }
        }
      }

      double solR_gss = 0;
      for(unsigned i = 0; i < nDofsR; i++) {
        solR_gss += phiR[i] * solR[i];
      }
      
      vector < double > solV_gss(dim, 0);
      vector < vector < double > > gradSolV_gss(dim);

      for(unsigned  k = 0; k < dim; k++) {
        gradSolV_gss[k].resize(dim);
        std::fill(gradSolV_gss[k].begin(), gradSolV_gss[k].end(), 0);
      }

      for(unsigned i = 0; i < nDofsV; i++) {
        for(unsigned  k = 0; k < dim; k++) {
          solV_gss[k] += phiV[i] * solV[k][i];
        }

        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            gradSolV_gss[k][j] += phiV_x[i * dim + j] * solV[k][i];
          }
        }
      }

      double solP_gss = 0;

      for(unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }

      
      //BEGIN phiB_i loop: Momentum balance
      for(unsigned i = 0; i < nDofsV; i++) {

        for(unsigned k = 0; k < dim; k++) {
          unsigned irow = nDofsT + k * nDofsV + i;

          for(unsigned l = 0; l < dim; l++) {
            Res[irow] +=  -sqrt(Pr / Ra) * phiV_x[i * dim + l] * (gradSolV_gss[k][l] + gradSolV_gss[l][k]) * weight;
            Res[irow] +=  -phiV[i] * solV_gss[l] * gradSolV_gss[k][l] * weight;
          }

          Res[irow] += solP_gss * phiV_x[i * dim + k] * weight;

          if(k == 1) {
            Res[irow] += beta * solT_gss * phiV[i] * weight;
          }

          if(assembleMatrix) {
            unsigned irowMat = nDofsTVP * irow;

            for(unsigned l = 0; l < dim; l++) {
              for(unsigned j = 0; j < nDofsV; j++) {
                unsigned jcol1 = (nDofsT + k * nDofsV + j);
                unsigned jcol2 = (nDofsT + l * nDofsV + j);
                Jac[ irowMat + jcol1] += sqrt(Pr / Ra) * phiV_x[i * dim + l] * phiV_x[j * dim + l] * weight;
                Jac[ irowMat + jcol2] += sqrt(Pr / Ra) * phiV_x[i * dim + l] * phiV_x[j * dim + k] * weight;
                Jac[ irowMat + jcol1] += phiV[i] * solV_gss[l] * phiV_x[j * dim + l] * weight;
                Jac[ irowMat + jcol2] += phiV[i] * phiV[j] * gradSolV_gss[k][l] * weight;
              }
            }

            for(unsigned j = 0; j < nDofsP; j++) {
              unsigned jcol = (nDofsT + dim * nDofsV) + j;
              Jac[ irowMat + jcol] += - phiV_x[i * dim + k] * phiP[j] * weight;
            }

            if(k == 1) {
              for(unsigned j = 0; j < nDofsT; j++) {
                Jac[ irowMat + j ] += -beta * phiV[i] * phiT[j] * weight;
              }
            }
          }

        }
      }

      //END phiB_i loop
      
      
      
      
      
      for(unsigned i = 0; i < nDofsT; i++) {
        unsigned irow = i;

        for(unsigned k = 0; k < dim; k++) {
          Res[irow] +=  -alpha / sqrt(Ra * Pr) * phiT_x[i * dim + k] * gradSolT_gss[k] * weight;
          Res[irow] +=  -phiT[i] * solV_gss[k] * gradSolT_gss[k] * weight;

          if(assembleMatrix) {
            unsigned irowMat = irow * nDofsTVP;

            for(unsigned j = 0; j < nDofsT; j++) {
              Jac[ irowMat + j ] +=  alpha / sqrt(Ra * Pr) * phiT_x[i * dim + k] * phiT_x[j * dim + k] * weight;
              Jac[ irowMat + j ] +=  phiT[i] * solV_gss[k] * phiT_x[j * dim + k] * weight;
            }

            for(unsigned j = 0; j < nDofsV; j++) {
              unsigned jcol = nDofsT + k * nDofsV + j;
              Jac[ irowMat + jcol ] += phiT[i] * phiV[j] * gradSolT_gss[k] * weight;
            }
          }

        }
      }

      //END phiT_i loop

      //BEGIN phiV_i loop: Momentum balance
      for(unsigned i = 0; i < nDofsV; i++) {

        for(unsigned k = 0; k < dim; k++) {
          unsigned irow = nDofsT + k * nDofsV + i;

          for(unsigned l = 0; l < dim; l++) {
            Res[irow] +=  -sqrt(Pr / Ra) * phiV_x[i * dim + l] * (gradSolV_gss[k][l] + gradSolV_gss[l][k]) * weight;
            Res[irow] +=  -phiV[i] * solV_gss[l] * gradSolV_gss[k][l] * weight;
          }

          Res[irow] += solP_gss * phiV_x[i * dim + k] * weight;

          if(k == 1) {
            Res[irow] += beta * solT_gss * phiV[i] * weight;
          }

          if(assembleMatrix) {
            unsigned irowMat = nDofsTVP * irow;

            for(unsigned l = 0; l < dim; l++) {
              for(unsigned j = 0; j < nDofsV; j++) {
                unsigned jcol1 = (nDofsT + k * nDofsV + j);
                unsigned jcol2 = (nDofsT + l * nDofsV + j);
                Jac[ irowMat + jcol1] += sqrt(Pr / Ra) * phiV_x[i * dim + l] * phiV_x[j * dim + l] * weight;
                Jac[ irowMat + jcol2] += sqrt(Pr / Ra) * phiV_x[i * dim + l] * phiV_x[j * dim + k] * weight;
                Jac[ irowMat + jcol1] += phiV[i] * solV_gss[l] * phiV_x[j * dim + l] * weight;
                Jac[ irowMat + jcol2] += phiV[i] * phiV[j] * gradSolV_gss[k][l] * weight;
              }
            }

            for(unsigned j = 0; j < nDofsP; j++) {
              unsigned jcol = (nDofsT + dim * nDofsV) + j;
              Jac[ irowMat + jcol] += - phiV_x[i * dim + k] * phiP[j] * weight;
            }

            if(k == 1) {
              for(unsigned j = 0; j < nDofsT; j++) {
                Jac[ irowMat + j ] += -beta * phiV[i] * phiT[j] * weight;
              }
            }
          }

        }
      }

      //END phiV_i loop

      //BEGIN phiP_i loop: mass balance
      for(unsigned i = 0; i < nDofsP; i++) {
        unsigned irow = nDofsT + dim * nDofsV + i;

        for(int k = 0; k < dim; k++) {
          Res[irow] += -(gradSolV_gss[k][k]) * phiP[i]  * weight;

          if(assembleMatrix) {
            unsigned irowMat = nDofsTVP * irow;

            for(unsigned j = 0; j < nDofsV; j++) {
              unsigned jcol = (nDofsT + k * nDofsV + j);
              Jac[ irowMat + jcol ] += phiP[i] * phiV_x[j * dim + k] * weight;
            }
          }

        }
      }

      //END phiP_i loop

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

  RES->close();

  if(assembleMatrix) {
    KK->close();
  }

  // ***************** END ASSEMBLY *******************
}

// void PrintConvergenceInfo(char *stdOutfile, char* outfile, const unsigned &numofrefinements){
// 
//   std::cout<<"END_COMPUTATION\n"<<std::flush;
// 
//   std::ifstream inf;
//   inf.open(stdOutfile);
//   if (!inf) {
//     std::cout<<"Redirected standard output file not found\n";
//     std::cout<<"add option -std_output std_out_filename > std_out_filename\n";
//     return;
//   }
// 
//   std::ofstream outf;
// 
//   outf.open(outfile, std::ofstream::app);
//   outf << std::endl << std::endl;
//   outf << "Number_of_refinements="<<numofrefinements<<std::endl;
//   outf << "Nonlinear_Iteration,resid_norm0,resid_normN,N,convergence";
// 
//   std::string str1;
//   inf >> str1;
//   while (str1.compare("END_COMPUTATION") != 0) {
// 
//     if (str1.compare("Nonlinear") == 0) {
//       inf >> str1;
//       if (str1.compare("iteration") == 0) {
//         inf >> str1;
//         outf << std::endl << str1;
//       }
//     }
//     else if (str1.compare("KSP") == 0){
//       inf >> str1;
//       if (str1.compare("preconditioned") == 0){
//         inf >> str1;
//         if (str1.compare("resid") == 0){
//           inf >> str1;
//           if (str1.compare("norm") == 0){
//             double norm0 = 1.;
//             double normN = 1.;
//             unsigned counter = 0;
//             inf >> norm0;
//             outf <<","<< norm0;
//             for (unsigned i = 0; i < 11; i++){
//               inf >> str1;
//             }
//             while(str1.compare("norm") == 0){
//               inf >> normN;
//               counter++;
//               for (unsigned i = 0; i < 11; i++){
//                 inf >> str1;
//               }
//             }
//             outf <<","<< normN;
//             if(counter != 0){
//               outf << "," <<counter<< "," << pow(normN/norm0,1./counter);
//             }
//             else{
//               outf << "Invalid solver, set -outer_ksp_solver \"gmres\"";
//             }
//           }
//         }
//       }
//     }
//     inf >> str1;
//   }
// 
//   outf.close();
//   inf.close();
// 
// }


