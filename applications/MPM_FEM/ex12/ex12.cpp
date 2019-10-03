#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"

#include "PetscMatrix.hpp"

using namespace femus;

std::vector < SparseMatrix* > prol;

// bool SetBoundaryCondition (const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
//   bool dirichlet = true; //dirichlet
//   value = 0;
//
//   if (!strcmp (SolName, "weight")) {
//     dirichlet = false;
//   }
//   //if (facename == 2)
//   // dirichlet = false;
//
//   return dirichlet;
// }

bool SetBoundaryCondition (const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;

  if (!strcmp (SolName, "u")) { // strcmp compares two string in lexiographic sense. why !strcmp??? not false ==true?
    if (facename == 2) {
      value = - (x[1] - 0.5) * (x[1] + 0.5);
    }
    else if (facename == 3) {
      dirichlet = false;
    }
  }
  else if (!strcmp (SolName, "v")) {
    if (facename == 2) {
      value = 0.;
    }
    else if (facename == 3) {
      dirichlet = false;
    }
  }
  else if (!strcmp (SolName, "w")) {
    value = 0.;
  }
  else if (!strcmp (SolName, "P")) {
    dirichlet = false;
    value = 0.;
    if (facename == 3) {
      value = 0.;
    }
  }
  else {
    dirichlet = false;
    value = 0.;
  }

  return dirichlet;
}


double InitalValueU (const std::vector < double >& x) {
  return x[0] + x[1];
}

void BuidProjection (MultiLevelProblem& ml_prob, const unsigned &level);
void AssembleWithProjection (MultiLevelProblem& ml_prob);
void Assemble (MultiLevelProblem& ml_prob);
void ProjectSolutionIntoGradient (MultiLevelProblem& ml_prob, const unsigned &level);
std::pair < double, double > GetErrorNormWithProjection (MultiLevelSolution* mlSol);

int main (int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit (argc, args, MPI_COMM_WORLD);


  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;

  //mlMsh.GenerateCoarseBoxMesh (2., 0, 0, -0.5, 0.5, 0., 0., 0., 0., EDGE3, "seventh");

  mlMsh.ReadCoarseMesh ("./input/rectangle.neu", "seventh", scalingFactor);

  //mlMsh.ReadCoarseMesh ("./input/square_quad.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("./input/square_tri.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/square_mixed.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh ("./input/cube_hex.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cube_wedge.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cube_tet.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cube_mixed.neu","seventh",scalingFactor);

  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
   *    probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();
  unsigned dim2 = dim * dim;
  unsigned maxNumberOfMeshes;

  if (dim == 1) {
    maxNumberOfMeshes = 10;
  }
  else if (dim == 2) {
    maxNumberOfMeshes = 7;
  }
  else {
    maxNumberOfMeshes = 2;
  }

  vector < double > l2Norm (maxNumberOfMeshes);
  vector < double > semiNorm (maxNumberOfMeshes);

  for (unsigned i = maxNumberOfMeshes - 1; i < maxNumberOfMeshes; i++) {   // loop on the mesh level

    unsigned numberOfUniformLevels = i ;
    unsigned numberOfSelectiveLevels = 0;
    mlMsh.RefineMesh (numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

    // erase all the coarse mesh levels
    //mlMsh.EraseCoarseLevels (numberOfUniformLevels - 1);
    prol.resize (i);

    // print mesh info
    mlMsh.PrintInfo();

    // define the multilevel solution and attach the mlMsh object to it
    MultiLevelSolution mlSol (&mlMsh);


    FEOrder feOrder = SECOND;// FIRST; // SERENDIPITY; //SECOND;
    // add variables to mlSol

    std::string Uname[3] = {"u", "v", "w"};
    for (unsigned k = 0; k < dim; k++) {
      mlSol.AddSolution (Uname[k].c_str(), LAGRANGE, feOrder);
    }
    mlSol.AddSolution ("P", DISCONTINUOUS_POLYNOMIAL, FIRST);
    //mlSol.FixSolutionAtOnePoint ("P");

    std::string Uxname[3][3] = {{"ux", "uy", "uz"}, {"vx", "vy", "vz"}, {"wx", "wy", "wz"}};
    for (unsigned k = 0; k < dim; k++) {
      for (unsigned l = 0; l < dim; l++) {
        mlSol.AddSolution (Uxname[k][l].c_str(), LAGRANGE, feOrder);
      }
    }

    mlSol.AddSolution ("weight", LAGRANGE, feOrder, 0);
    mlSol.Initialize ("All");

//     mlSol.Initialize ("u", InitalValueU);
//     mlSol.Initialize ("v", InitalValueU);
//     mlSol.Initialize ("w", InitalValueU);

    // attach the boundary condition function and generate boundary data
    mlSol.AttachSetBoundaryConditionFunction (SetBoundaryCondition);
    mlSol.GenerateBdc ("All");

    // define the multilevel problem attach the mlSol object to it
    MultiLevelProblem mlProb (&mlSol);





    std::vector < LinearImplicitSystem* > systemP (dim2 + dim + 1);
    std::string Pxname[3][3] = { {"Pux", "Puy", "Puz"}, {"Pvx", "Pvy", "Pvz"}, {"Pwx", "Pwy", "Pwz"}};
    for (unsigned k = 0; k < dim; k++) {
      for (unsigned l = 0; l < dim; l++) {
        systemP[k * dim + l] = &mlProb.add_system < LinearImplicitSystem > (Pxname[k][l]);
        systemP[k * dim + l]->AddSolutionToSystemPDE (Uxname[k][l].c_str());
        systemP[k * dim + l]->init();
      }
    }
    std::string Pname[3] = {"Pu", "Pv", "Pw"};
    for (unsigned k = 0; k < dim; k++) {
      systemP[dim2 + k] = &mlProb.add_system < LinearImplicitSystem > (Pname[k]);
      systemP[dim2 + k]->AddSolutionToSystemPDE (Uname[k].c_str());
      systemP[dim2 + k]->init();
    }

    systemP[dim2 + dim] = &mlProb.add_system < LinearImplicitSystem > ("PP");
    systemP[dim2 + dim]->AddSolutionToSystemPDE ("P");
    systemP[dim2 + dim]->init();

    for (unsigned k = 0; k < i; k++) {
      BuidProjection (mlProb, k);
    }

    LinearImplicitSystem& system0 = mlProb.add_system < LinearImplicitSystem > ("auxiliary");

    for (unsigned k = 0; k < dim; k++) {
      for (unsigned l = 0; l < dim; l++) {
        system0.AddSolutionToSystemPDE (Uxname[k][l].c_str());
      }
    }
    for (unsigned k = 0; k < dim; k++) {
      system0.AddSolutionToSystemPDE (Uname[k].c_str());
    }
    system0.AddSolutionToSystemPDE ("P");
    system0.init();


    NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("NS");
    for (unsigned k = 0; k < dim; k++) {
      system.AddSolutionToSystemPDE (Uname[k].c_str());
    }
    system.AddSolutionToSystemPDE ("P");
    


    // attach the assembling function to system
    system.SetAssembleFunction (AssembleWithProjection);
    //system.SetAssembleFunction (Assemble);

    //system.SetLinearEquationSolverType(FEMuS_ASM);

    system.init();

    // ******* Set Smoother *******
    //system.SetSolverFineGrids(RICHARDSON);

    system.SetPreconditionerFineGrids(ILU_PRECOND);

    system.SetTolerances(1.e-12, 1.e-20, 1.e+50, 20, 10);
    
    // system.ClearVariablesToBeSolved();
    // system.AddVariableToBeSolved("All");

    // ******* Set the last (1) variables in system (i.e. P) to be a schur variable *******
    // system.SetNumberOfSchurVariables(1);

    // ******* Set block size for the ASM smoothers *******
    // system.SetElementBlockNumber(1);
    // initilaize and solve the system
    //system.SetOuterSolver (PREONLY);
    system.MGsolve();

    ProjectSolutionIntoGradient (mlProb, i - 1);

    std::pair< double , double > norm = GetErrorNormWithProjection (&mlSol);

    l2Norm[i]  = norm.first;
    semiNorm[i] = norm.second;
    // print solutions
    std::vector < std::string > variablesToBePrinted;
    variablesToBePrinted.push_back ("All");

    VTKWriter vtkIO (&mlSol);
    vtkIO.SetDebugOutput (true);
    // vtkIO.Write (DEFAULT_OUTPUTDIR, "linear", variablesToBePrinted, i);
    //vtkIO.Write (DEFAULT_OUTPUTDIR, "quadratic", variablesToBePrinted, i);
    vtkIO.Write (DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, i);
    for (unsigned k = 0; k < i; k++) {
      delete prol[k];
    }
  }


//   // print the seminorm of the error and the order of convergence between different levels
//   std::cout << std::endl;
//   std::cout << std::endl;
//   std::cout << "l2 ERROR and ORDER OF CONVERGENCE:\n\n";
//
//   for (unsigned i = 1; i < maxNumberOfMeshes; i++) {
//     std::cout << i + 1 << "\t";
//     std::cout.precision (14);
//     std::cout << l2Norm[i] << "\t";
//     std::cout << std::endl;
//     if (i < maxNumberOfMeshes - 1) {
//       std::cout.precision (3);
//       std::cout << "\t\t";
//       std::cout << log (l2Norm[i] / l2Norm[i + 1]) / log (2.) << "\t\t\t";
//       std::cout << std::endl;
//     }
//   }
//
//   std::cout << std::endl;
//   std::cout << std::endl;
//   std::cout << "SEMINORM ERROR and ORDER OF CONVERGENCE:\n\n";
//   //std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";
//
//   for (unsigned i = 1; i < maxNumberOfMeshes; i++) {
//     std::cout << i + 1 << "\t";
//     std::cout.precision (14);
//     std::cout << semiNorm[i] << "\t";
//     std::cout << std::endl;
//     if (i < maxNumberOfMeshes - 1) {
//       std::cout.precision (3);
//       std::cout << "\t\t";
//       std::cout << log (semiNorm[i] / semiNorm[i + 1]) / log (2.) << "\t\t\t";
//       std::cout << std::endl;
//     }
//   }
}


double GetExactSolutionValue (const std::vector < double >& x) {
  double pi = acos (-1.);
  if (x.size() == 1) {
    return cos (pi * x[0]);
  }
  else {
    return cos (pi * x[0]) * cos (pi * x[1]);
  }
};


void GetExactSolutionGradient (const std::vector < double >& x, vector < double >& solGrad) {
  double pi = acos (-1.);
  if (x.size() == 1) {
    solGrad[0]  = -pi * sin (pi * x[0]);
  }
  else {
    solGrad[0]  = -pi * sin (pi * x[0]) * cos (pi * x[1]);
    solGrad[1] = -pi * cos (pi * x[0]) * sin (pi * x[1]);
  }
};


double GetExactSolutionLaplace (const std::vector < double >& x) {
  double pi = acos (-1.);
  if (x.size() == 1) {
    return -pi * pi * cos (pi * x[0]);     // - pi*pi*cos(pi*x[0])*cos(pi*x[1]);
  }
  else {
    return -2.*pi * pi * cos (pi * x[0]) * cos (pi * x[1]);     // - pi*pi*cos(pi*x[0])*cos(pi*x[1]);
  }
};


std::pair < double, double > GetErrorNormWithProjection (MultiLevelSolution* mlSol) {
  unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1u;
  //  extract pointers to the several objects that we are going to use
  Mesh*          msh          = mlSol->_mlMesh->GetLevel (level);   // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)
  Solution*    sol        = mlSol->GetSolutionLevel (level);   // pointer to the solution (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

//solution variable


  std::string Uxname[3] = {"ux", "uy", "uz"};
  std::vector < unsigned > soluIndex (dim + 1);
  for (unsigned k = 0; k < dim; k++) {
    soluIndex[k] = mlSol->GetIndex (Uxname[k].c_str());
  }

  soluIndex[dim] = mlSol->GetIndex ("u");   // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType (soluIndex[dim]);   // get the finite element type for "u"

  std::vector < std::vector < double > >  solu (dim + 1); // local solution

  vector < vector < double > > x (dim);   // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  // reserve memory for the local standar vectors

  double seminorm = 0.;
  double l2norm = 0.;

// element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {


    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofs  = msh->GetElementDofNumber (iel, soluType);   // number of solution element dofs

    // resize local arrays
    for (int k = 0; k < dim; k++) {
      x[k].resize (nDofs);
      solu[k].resize (nDofs);
    }
    solu[dim].resize (nDofs);

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, soluType);   // global to global mapping between solution node and solution dof
      for (unsigned k = 0; k < dim + 1; k++)
        solu[k][i] = (*sol->_Sol[soluIndex[k]]) (solDof);     // global extraction and local storage for the solution
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);   // global to global mapping between coordinates node and coordinate dof

      for (unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k]) (xDof);     // global extraction and local storage for the element coordinates
      }
    }

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian (x, ig, weight, phi, phi_x);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      double soluGauss = 0;
      vector < double > soluGauss_x (dim, 0.);
      vector < double > xGauss (dim, 0.);

      for (unsigned i = 0; i < nDofs; i++) {
        soluGauss += phi[i] * solu[dim][i];

        for (unsigned k = 0; k < dim; k++) {
          soluGauss_x[k] += phi[i] * solu[k][i];
          xGauss[k] += x[k][i] * phi[i];
        }
      }

      vector <double> solGrad (dim);
      GetExactSolutionGradient (xGauss, solGrad);

      for (unsigned j = 0; j < dim ; j++) {
        seminorm   += ( (soluGauss_x[j] - solGrad[j]) * (soluGauss_x[j] - solGrad[j])) * weight;
      }

      double exactSol = GetExactSolutionValue (xGauss);
      l2norm += (exactSol - soluGauss) * (exactSol - soluGauss) * weight;
    } // end gauss point loop
  } //end element loop for each process

// add the norms of all processes
  NumericVector* norm_vec;
  norm_vec = NumericVector::build().release();
  norm_vec->init (msh->n_processors(), 1 , false, AUTOMATIC);

  norm_vec->set (iproc, l2norm);
  norm_vec->close();
  l2norm = norm_vec->l1_norm();

  norm_vec->set (iproc, seminorm);
  norm_vec->close();
  seminorm = norm_vec->l1_norm();

  delete norm_vec;

  std::pair < double, double > norm;
  norm.first  = sqrt (l2norm);
  norm.second = sqrt (seminorm);

  return norm;

}


void BuidProjection (MultiLevelProblem& ml_prob, const unsigned &level) {

  double scale = 4.;

  adept::Stack& s = FemusInit::_adeptStack;

  MultiLevelSolution*  mlSol = ml_prob._ml_sol;
  //unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1u;

  Solution* sol = ml_prob._ml_sol->GetSolutionLevel (level);
  Mesh* msh = ml_prob._ml_msh->GetLevel (level);
  elem* el = msh->el;

  unsigned  dim = msh->GetDimension();
  unsigned  dim2 = dim * dim;

  std::vector < LinearImplicitSystem* > mlSysP (dim2 + dim + 1);
  std::vector < LinearEquationSolver* > sysP (dim2 + dim + 1);
  std::vector < SparseMatrix*> P ( (dim2 + dim + 1) * (dim + 1));
  for (unsigned i = 0; i < P.size(); i++) P[i] = NULL;

  std::string Pxname[3][3] = {{"Pux", "Puy", "Puz"}, {"Pvx", "Pvy", "Pvz"}, {"Pwx", "Pwy", "Pwz"}};
  std::string Uxname[3] = {"ux", "uy", "uz"};
  std::vector < unsigned > soluxIndex (dim);
  for (unsigned k = 0; k < dim; k++) {
    soluxIndex[k] = mlSol->GetIndex (Uxname[k].c_str());
    for (unsigned l = 0; l < dim; l++) {
      mlSysP[k * dim + l] =  &ml_prob.get_system< LinearImplicitSystem > (Pxname[k][l]);
      sysP[k * dim + l] = mlSysP[k * dim + l]->_LinSolver[level];
      P[ (k * dim + l) * (dim + 1) + k ] = sysP[k * dim + l]->_KK;
      P[ (k * dim + l) * (dim + 1) + k ]->zero();
    }
  }

  std::string Pname[3] = {"Pu", "Pv", "Pw"};
  for (unsigned k = 0; k < dim; k++) {
    mlSysP[dim2 + k] =  &ml_prob.get_system< LinearImplicitSystem > (Pname[k]);
    sysP[dim2 + k] = mlSysP[dim2 + k]->_LinSolver[level];
    P[ (dim2 + k) * (dim + 1) + k] = sysP[dim2 + k]->_KK;
    P[ (dim2 + k) * (dim + 1) + k]->zero();
    NumericVector* D = sysP[dim2 + k]->_RES;
    *D = 1.;
    P[ (dim2 + k) * (dim + 1) + k]->matrix_set_diagonal_values (*D);
  }

  {
    mlSysP[dim2 + dim] =  &ml_prob.get_system< LinearImplicitSystem > ("PP");
    sysP[dim2 + dim] = mlSysP[dim2 + dim]->_LinSolver[level];
    P[ (dim2 + dim) * (dim + 1) + dim] = sysP[dim2 + dim]->_KK;
    P[ (dim2 + dim) * (dim + 1) + dim]->zero();
    NumericVector* D = sysP[dim2 + dim]->_RES;
    *D = 1.;
    P[ (dim2 + dim) * (dim + 1) + dim]->matrix_set_diagonal_values (*D);
  }

  vector < vector < double > > x (dim);
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector< int > sysDof;
  vector <double> phi;
  double* phi2;
  vector <double> phi_x;
  double weight;

  unsigned    iproc = msh->processor_id();

  unsigned solwIndex = mlSol->GetIndex ("weight");
  unsigned solType = mlSol->GetSolutionType (solwIndex);
  sol->_Sol[solwIndex]->zero();

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofs  = msh->GetElementDofNumber (iel, solType);
    sysDof.resize (nDofs);
    for (int k = 0; k < dim; k++) {
      x[k].resize (nDofs);
    }
    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, solType);
      sysDof[i] = msh->GetSolutionDof (i, iel, solType);
    }
    // local storage of coordinates
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);
      for (unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k]) (xDof);
      }
    }

    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
      msh->_finiteElement[ielGeom][solType]->Jacobian (x, ig, weight, phi, phi_x);
      phi2 = (solType != 1) ? msh->_finiteElement[ielGeom][solType]->GetPhi (ig) : msh->_finiteElement[ielGeom][2]->GetPhi (ig);

      vector < double > xGauss (dim, 0.);
      for (unsigned i = 0; i < nDofs; i++) {
        for (unsigned k = 0; k < dim; k++) {
          xGauss[k] += x[k][i] * phi[i];
        }
      }
      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {
        sol->_Sol[solwIndex]->add (sysDof[i], phi2[i] * weight);
      } // end phi_i loop
    } // end gauss point loop

  } //end element loop for each process*/
  sol->_Sol[solwIndex]->close();


  //solution variable

  unsigned soluIndex = mlSol->GetIndex ("u");
  {
    unsigned soluType = mlSol->GetSolutionType (soluIndex);
    if (soluType != solType) {
      std::cout << "error weight and u should be of the same type\n";
      abort();
    }
  }
  unsigned soluPdeIndex = mlSysP[0]->GetSolPdeIndex ("ux");
  std::vector < adept::adouble > solu;

  vector <double> Jac;
  std::vector < std::vector< adept::adouble > > aRes (dim); // local redidual vector

  std::vector < double > solw;

  //BEGIN element loop
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofs  = msh->GetElementDofNumber (iel, solType);
    solu.resize (nDofs);
    solw.resize (nDofs);
    sysDof.resize (nDofs);
    for (int k = 0; k < dim; k++) {
      x[k].resize (nDofs);
      aRes[k].assign (nDofs, 0.);
    }
    Jac.resize (nDofs * nDofs);
    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, solType);
      solu[i] = (*sol->_Sol[soluIndex]) (solDof);
      solw[i] = (*sol->_Sol[solwIndex]) (solDof);
      sysDof[i] = sysP[0]->GetSystemDof (soluIndex, soluPdeIndex, i, iel);
    }
    // local storage of coordinates
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);
      for (unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k]) (xDof);
      }
    }

    s.new_recording();
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {

      msh->_finiteElement[ielGeom][solType]->Jacobian (x, ig, weight, phi, phi_x);
      phi2 = (solType != 1) ? msh->_finiteElement[ielGeom][solType]->GetPhi (ig) : msh->_finiteElement[ielGeom][2]->GetPhi (ig);
      std::vector < adept::adouble > solux_g (dim, 0.);
      for (unsigned i = 0; i < nDofs; i++) {
        for (unsigned k = 0; k < dim; k++) {
          solux_g[k] += phi_x[i * dim + k] * solu[i];
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {
        for (unsigned k = 0; k < dim; k++) {
          aRes[k][i] += solux_g[k] * phi2[i] * weight / solw[i];
        }
      } // end phi_i loop
    } // end gauss point loop

    s.independent (&solu[0], nDofs);
    for (unsigned l = 0; l < dim; l++) {
      s.dependent (&aRes[l][0], nDofs);
      s.jacobian (&Jac[0], true);
      for (unsigned k = 0; k < dim; k++) {
        P[ (k * dim + l) * (dim + 1) + k]->add_matrix_blocked (Jac, sysDof, sysDof);
      }
      s.clear_dependents();
    }
    s.clear_independents();
  } //end element loop for each process*/

  for (unsigned k = 0; k < dim; k++) {
    for (unsigned l = 0; l < dim; l++) {
      P[ (k * dim + l) * (dim + 1) + k]->close();
    }
  }

  prol[level] = SparseMatrix::build().release();
  prol[level]->init (dim2 + dim + 1, dim + 1, P);

  double tolerance = 1.0e-14;
  prol[level]->RemoveZeroEntries (tolerance);

//   restr = SparseMatrix::build().release();
//   restr->init (dim2 + dim + 1, dim + 1, P);
//   restr->get_transpose(*restr);

  //MatView ( (static_cast< PetscMatrix* > (Proj))->mat(), PETSC_VIEWER_STDOUT_WORLD);

}

unsigned counter = 0;
unsigned levelOld = 0;

void AssembleWithProjection (MultiLevelProblem& ml_prob) {

  bool extended = true;

  adept::Stack& s = FemusInit::_adeptStack;
  
  NonLinearImplicitSystem* mlSys1 =  &ml_prob.get_system< NonLinearImplicitSystem > ("NS");
  const unsigned level = mlSys1->GetLevelToAssemble();
  
  if(level != levelOld){
    counter = 0;
    levelOld = level;    
  }

  LinearEquationSolver* sys1 = mlSys1->_LinSolver[level];
  SparseMatrix *K1 = sys1->_KK; ;
  NumericVector* RES1  = sys1->_RES;

  MultiLevelSolution*  mlSol = ml_prob._ml_sol;
  ProjectSolutionIntoGradient (ml_prob, level);

  Solution* solution = ml_prob._ml_sol->GetSolutionLevel (level);
  Mesh* msh = ml_prob._ml_msh->GetLevel (level);
  elem* el = msh->el;

  unsigned  dim = msh->GetDimension();
  unsigned  dim2 = dim * dim;

  LinearImplicitSystem* mlSys0 =  &ml_prob.get_system< LinearImplicitSystem > ("auxiliary");
  LinearEquationSolver* sys0 = mlSys0->_LinSolver[level];

  SparseMatrix *K = sys0->_KK; ;
  NumericVector* RES  = sys0->_RES;



  K->zero();
  RES->zero();

  K1->zero();
  RES1->zero();

  std::string Uxname[3][3] = {{"ux", "uy", "uz"}, {"vx", "vy", "vz"}, {"wx", "wy", "wz"}};
  std::string Uname[3] = {"u", "v", "w"};

  std::vector < unsigned > solIndex (dim2 + dim + 1);
  for (unsigned k = 0; k < dim; k++) {
    for (unsigned l = 0; l < dim; l++) {
      solIndex[k * dim + l] = mlSol->GetIndex (Uxname[k][l].c_str());
    }
  }
  for (unsigned k = 0; k < dim; k++) {
    solIndex[dim2 + k] = mlSol->GetIndex (Uname[k].c_str());
  }
  solIndex[dim2 + dim] = mlSol->GetIndex ("P");

  unsigned solTypeU = mlSol->GetSolutionType (solIndex[0]);
  unsigned solTypeP = mlSol->GetSolutionType (solIndex[dim2 + dim]);

  vector < vector < double > > x (dim);
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector< unsigned > sysDof;
  vector< unsigned > sysDof1;
  vector <double> phi;
  double* phiP;
  vector <double> phi_x;
  double weight;

  unsigned  iproc = msh->processor_id();

  vector < unsigned > solPdeIndex (dim2 + dim  + 1);
  for (unsigned k = 0; k < dim; k++) {
    for (unsigned l = 0; l < dim; l++) {
      solPdeIndex[ k * dim + l] = mlSys0->GetSolPdeIndex (Uxname[k][l].c_str());
    }
  }
  for (unsigned k = 0; k < dim; k++) {
    solPdeIndex[dim2 + k] = mlSys0->GetSolPdeIndex (Uname[k].c_str());
  }
  solPdeIndex[dim2 + dim] = mlSys0->GetSolPdeIndex ("P");


  vector < unsigned > solPde1Index (dim  + 1);
  for (unsigned k = 0; k < dim; k++) {
    solPde1Index[k] = mlSys1->GetSolPdeIndex (Uname[k].c_str());
  }
  solPde1Index[dim] = mlSys1->GetSolPdeIndex ("P");

  std::vector < std::vector < adept::adouble > > sol (dim2 + dim + 1);
  std::vector < std::vector< adept::adouble > > aRes (dim2 + dim + 1); // local residual vector

  vector <double> Jac;
  std::vector< double > res (dim2 + dim + 1);

  //BEGIN element loop
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofsU  = msh->GetElementDofNumber (iel, solTypeU);
    unsigned nDofsP  = msh->GetElementDofNumber (iel, solTypeP);

    for (int k = 0; k < dim2 + dim; k++) {
      sol[k].resize (nDofsU);
      aRes[k].assign (nDofsU, 0.);
    }
    sol[dim2 + dim].resize (nDofsP);
    aRes[dim2 + dim].assign (nDofsP, 0.);

    for (int k = 0; k < dim; k++) {
      x[k].resize (nDofsU);
    }
    sysDof.resize (nDofsU * (dim2 + dim) + nDofsP);
    sysDof1.resize (nDofsU * dim + nDofsP);

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofsU; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, solTypeU);
      for (int k = 0; k < dim2 + dim; k++) {
        sol[k][i] = (*solution->_Sol[solIndex[k]]) (solDof);
        sysDof[k * nDofsU + i] = sys0->GetSystemDof (solIndex[k], solPdeIndex[k], i, iel);
        if (k >= dim2) {
          sysDof1[ (k - dim2) * nDofsU + i] = sys1->GetSystemDof (solIndex[k], solPde1Index[k - dim2], i, iel);
        }
      }
    }

    for (unsigned i = 0; i < nDofsP; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, solTypeP);
      sol[dim2 + dim][i] = (*solution->_Sol[solIndex[dim2 + dim]]) (solDof);
      sysDof[ (dim2 + dim) * nDofsU + i] = sys0->GetSystemDof (solIndex[dim2 + dim], solPdeIndex[dim2 + dim], i, iel);
      sysDof1[dim * nDofsU + i] = sys1->GetSystemDof (solIndex[dim2 + dim], solPde1Index[dim], i, iel);
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofsU; i++) {
      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);
      for (unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k]) (xDof);
      }
    }

    s.new_recording();
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solTypeU]->GetGaussPointNumber(); ig++) {
      msh->_finiteElement[ielGeom][solTypeU]->Jacobian (x, ig, weight, phi, phi_x);
      phiP = msh->_finiteElement[ielGeom][solTypeP]->GetPhi (ig);

      std::vector < adept::adouble > solU (dim, 0.);
      std::vector < std::vector < adept::adouble > > solux (dim);
      std::vector < std::vector < adept::adouble > > solUx (dim);
      for (unsigned k = 0; k < dim; k++) {
        solux[k].assign (dim, 0.);
        solUx[k].assign (dim, 0.);
      }
      //vector < double > xGauss (dim, 0.);
      for (unsigned i = 0; i < nDofsU; i++) {
        for (unsigned k = 0; k < dim; k++) {
          for (unsigned l = 0; l < dim; l++) {
            solux[k][l] += phi[i] * sol[ k * dim + l][i]; //extended solution gradient
            solUx[k][l] += phi_x[i * dim + l] * sol[ dim2 + k][i]; //regular solution gradient
            //gradSolV_gss[k][j] += solV[k][i] * phiV_x[i * dim + j];
          }
          solU[k] += phi[i] * sol[dim2 + k][i];
        }
      }
      adept::adouble solP = 0.;
      for (unsigned i = 0; i < nDofsP; i++) {
        solP += phiP[i] * sol[dim2 + dim][i];
      }

      //double exactSolValue = GetExactSolutionValue (xGauss);
      //vector < double > exactSolGrad (dim);
      //GetExactSolutionGradient (xGauss , exactSolGrad);
      //double exactSolLaplace = GetExactSolutionLaplace (xGauss);

      double ReI = 100.;

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofsU; i++) {
        for (unsigned k = 0; k < dim; k++) {
          if (extended) {
            for (unsigned l = 0; l < dim; l++) {
              aRes[k * dim + l][i] += ReI * (solux[k][l] + solux[l][k]) * phi[i] * weight; //extended diffusivity
              //aRes[dim2 + k][i] += ReI * (solux[k][l] + solux[l][k]) * phi_x[ i * dim + l] * weight; //mixed 1 diffusivity
              //aRes[k * dim + l][i] += ReI * (solux[k][l] + solux[l][k]) * phi[i] * weight; //mixed 2 diffusivity
              aRes[dim2 + k][i] += phi[i] * solU[l] * solux[k][l] * weight; //extended convection
            }
            aRes[k * dim + k][i] -= solP * phi[i] * weight; //extended gradient
            //aRes[dim2 + k][i] -= solP * phi_x[ i * dim + k]  * weight; //regular gradient
          }
          else {
            for (unsigned l = 0; l < dim; l++) {
              aRes[dim2 + k][i] += ReI * (solUx[k][l] + solUx[l][k]) * phi_x[ i * dim + l] * weight; //regular diffusivity
              aRes[dim2 + k][i] += phi[i] * solU[l] * solUx[k][l] * weight; //regular convection
            }
            aRes[dim2 + k][i] -= solP * phi_x[ i * dim + k]  * weight; //regular gradient
          }
        }
      } // end phi_i loop

      for (unsigned i = 0; i < nDofsP; i++) {
        for (unsigned k = 0; k < dim; k++) {
          if (extended) {
            aRes[dim2 + dim][i] -= phiP[i] * solux[k][k] * weight; //extended divergence
          }
          else {
            aRes[dim2 + dim][i] -= phiP[i] * solUx[k][k] * weight; //regular divergence
          }
        }
      }
    } // end gauss point loop

    if (extended) {
      res.resize ( (dim2 + dim) * nDofsU + nDofsP);
      for (unsigned k = 0; k < dim2 + dim; k++) {
        for (int i = 0; i < nDofsU; i++) {
          res[k * nDofsU + i] = -aRes[k][i].value();
        }
      }
      for (int i = 0; i < nDofsP; i++) {
        res[ (dim2 + dim)  * nDofsU + i] = -aRes[dim2 + dim][i].value();
      }
      RES->add_vector_blocked (res, sysDof);


      Jac.resize ( ( (dim2 + dim) * nDofsU + nDofsP) * ( (dim2 + dim) * nDofsU + nDofsP));
      for (unsigned k = 0; k < dim2 + dim; k++) {
        s.independent (&sol[k][0], nDofsU);
      }
      s.independent (&sol[dim2 + dim][0], nDofsP);

      for (unsigned k = 0; k < dim2 + dim; k++) {
        s.dependent (&aRes[k][0], nDofsU);
      }
      s.dependent (&aRes[dim2 + dim][0], nDofsP);

      s.jacobian (&Jac[0], true);

      K->add_matrix_blocked (Jac, sysDof, sysDof);
      s.clear_dependents();
      s.clear_independents();
    }
    else {
      res.resize (dim * nDofsU + nDofsP);
      sysDof1.resize( dim * nDofsU + nDofsP );
      for (unsigned k = 0; k < dim; k++) {
        for (int i = 0; i < nDofsU; i++) {
          res[ k * nDofsU + i] = -aRes[dim2 + k][i].value();
          sysDof1[k * nDofsU + i] = sysDof[(dim2 + k) * nDofsU + i];
        }
      }
      for (int i = 0; i < nDofsP; i++) {
        res[ dim * nDofsU + i] = -aRes[dim2 + dim][i].value();
        sysDof1[dim * nDofsU + i] = sysDof[(dim2 + dim) * nDofsU + i];
      }
      RES->add_vector_blocked (res, sysDof1);


      Jac.resize ( (dim * nDofsU + nDofsP) * (dim * nDofsU + nDofsP));
      for (unsigned k = 0; k < dim; k++) {
        s.independent (&sol[dim2 + k][0], nDofsU);
      }
      s.independent (&sol[dim2 + dim][0], nDofsP);

      for (unsigned k = 0; k < dim; k++) {
        s.dependent (&aRes[dim2 + k][0], nDofsU);
      }
      s.dependent (&aRes[dim2 + dim][0], nDofsP);

      s.jacobian (&Jac[0], true);

      K->add_matrix_blocked (Jac, sysDof1, sysDof1);
      s.clear_dependents();
      s.clear_independents();
    }
  } //end element loop for each process*/

  K->close();
  RES->close();

  //K1->close();
  //RES1->close();

  ///if(extended){
    bool mathReuse = (counter == 0) ? false : true;

    //mathReuse = false;
    
    K1->matrix_PtAP (*prol[level], *K, mathReuse);
    RES1->matrix_mult_transpose (*RES, *prol[level]);
  //}

  counter++;


//   MatView((static_cast< PetscMatrix* > (K))->mat(),PETSC_VIEWER_STDOUT_WORLD);
//
//  VecView((static_cast< PetscVector* > (RES))->vec(),PETSC_VIEWER_STDOUT_WORLD);
//
//   VecView((static_cast< PetscVector* > (RES1))->vec(),PETSC_VIEWER_STDOUT_WORLD);


//   PetscViewer    viewer;
//   PetscViewerDrawOpen (PETSC_COMM_WORLD, NULL, NULL, 0, 0, 300, 700, &viewer);
//   PetscObjectSetName ( (PetscObject) viewer, "FSI matrix");
//   PetscViewerPushFormat (viewer, PETSC_VIEWER_DRAW_LG);
//   MatView ( (static_cast< PetscMatrix* > (prol[level]))->mat(), viewer);
// 
// 
// 
// 
//   PetscViewer    viewer0;
//   PetscViewerDrawOpen (PETSC_COMM_WORLD, NULL, NULL, 0, 0, 700, 700, &viewer0);
//   PetscObjectSetName ( (PetscObject) viewer0, "FSI matrix");
//   PetscViewerPushFormat (viewer0, PETSC_VIEWER_DRAW_LG);
//   MatView ( (static_cast< PetscMatrix* > (K))->mat(), viewer0);
// 
// 
// 
//   PetscViewer    viewer1;
//   PetscViewerDrawOpen (PETSC_COMM_WORLD, NULL, NULL, 0, 0, 1800, 1800, &viewer1);
//   PetscObjectSetName ( (PetscObject) viewer1, "FSI matrix");
//   PetscViewerPushFormat (viewer1, PETSC_VIEWER_DRAW_LG);
//   MatView ( (static_cast< PetscMatrix* > (K1))->mat(), viewer1);
// 
//   double a;
//   std::cin >> a;


}

void ProjectSolutionIntoGradient (MultiLevelProblem& ml_prob, const unsigned &level) {

  MultiLevelSolution*  mlSol = ml_prob._ml_sol;

  Solution* sol = ml_prob._ml_sol->GetSolutionLevel (level);
  Mesh* msh = ml_prob._ml_msh->GetLevel (level);
  elem* el = msh->el;

  unsigned  dim = msh->GetDimension();

  std::string Pxname[3][3] = {{"Pux", "Puy", "Puz"}, {"Pvx", "Pvy", "Pvz"}, {"Pwx", "Pwy", "Pwz"}};
  std::string Uxname[3][3] = {{"ux", "uy", "uz"}, {"vx", "vy", "vz"}, {"wx", "wy", "wz"}};
  std::string Uname[3] = {"u", "v", "w"};

  for (unsigned k = 0; k < dim; k++) {
    unsigned soluIndex = mlSol->GetIndex (Uname[k].c_str());
    for (unsigned l = 0; l < dim; l++) {
      LinearImplicitSystem* mlSysP =  &ml_prob.get_system< LinearImplicitSystem > (Pxname[k][l]);
      LinearEquationSolver* sysP = mlSysP->_LinSolver[level];
      unsigned solIndex = mlSol->GetIndex (Uxname[k][l].c_str());
      SparseMatrix* P = sysP->_KK;
      (*sol->_Sol[solIndex]).matrix_mult ( (*sol->_Sol[soluIndex]), (*P));
    }
  }

}

void Assemble (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned > (ceil (pow (3, dim)));       // conservative: based on line3, quad9, hex27

  //solution variable
  vector < unsigned > solVIndex (dim);
  solVIndex[0] = mlSol->GetIndex ("u");   // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex ("v");   // get the position of "V" in the ml_sol object

  if (dim == 3) solVIndex[2] = mlSol->GetIndex ("w");     // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType (solVIndex[0]);   // get the finite element type for "u"

  unsigned solPIndex;
  solPIndex = mlSol->GetIndex ("P");   // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType (solPIndex);   // get the finite element type for "u"

  vector < unsigned > solVPdeIndex (dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex ("u");   // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex ("v");   // get the position of "V" in the pdeSys object

  if (dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex ("w");

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex ("P");   // get the position of "P" in the pdeSys object

  vector < vector < adept::adouble > >  solV (dim);   // local solution
  vector < adept::adouble >  solP; // local solution

  vector< vector < adept::adouble > > aResV (dim);   // local redidual vector
  vector< adept::adouble > aResP; // local redidual vector

  vector < vector < double > > coordX (dim);   // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned  k = 0; k < dim; k++) {
    solV[k].reserve (maxSize);
    aResV[k].reserve (maxSize);
    coordX[k].reserve (maxSize);
  }

  solP.reserve (maxSize);
  aResP.reserve (maxSize);


  vector <double> phiV;  // local test function for velocity
  vector <double> phiV_x; // local test function first order partial derivatives

  phiV.reserve (maxSize);
  phiV_x.reserve (maxSize * dim);

  double* phiP; // local test function for the pressure
  double weight; // gauss point weight

  vector< int > sysDof; // local to global pdeSys dofs
  sysDof.reserve ( (dim + 1) * maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve ( (dim + 1) * maxSize);

  vector < double > Jac;
  Jac.reserve ( (dim + 1) * maxSize * (dim + 1) * maxSize);

  RES->zero(); // Set to zero all the entries of the Global Residual vector
  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);

    unsigned nDofsV = msh->GetElementDofNumber (iel, solVType);   // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber (iel, solPType);   // number of solution element dofs

    unsigned nDofsVP = dim * nDofsV + nDofsP;
    // resize local arrays
    sysDof.resize (nDofsVP);

    for (unsigned  k = 0; k < dim; k++) {
      solV[k].resize (nDofsV);
      coordX[k].resize (nDofsV);
    }
    solP.resize (nDofsP);

    for (unsigned  k = 0; k < dim; k++) {
      aResV[k].assign (nDofsV, 0.);
    }
    aResP.assign (nDofsP, 0.);

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof (i, iel, solVType);   // local to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]]) (solVDof);     // global extraction and local storage for the solution
        sysDof[k * nDofsV + i] = pdeSys->GetSystemDof (solVIndex[k], solVPdeIndex[k], i, iel);   // global to global mapping between solution node and pdeSys dof
      }
    }

    for (unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof (i, iel, solPType);   // local to global mapping between solution node and solution dof
      solP[i] = (*sol->_Sol[solPIndex]) (solPDof);     // global extraction and local storage for the solution
      sysDof[dim * nDofsV + i ] = pdeSys->GetSystemDof (solPIndex, solPPdeIndex, i, iel);   // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned coordXDof  = msh->GetSolutionDof (i, iel, coordXType);   // local to global mapping between coordinates node and coordinate dof
      for (unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k]) (coordXDof);     // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Face Gauss point loop (boundary Integral) ***
    for (unsigned jface = 0; jface < msh->GetElementFaceNumber (iel); jface++) {
      int faceIndex = el->GetBoundaryIndex (iel, jface);
      // look for boundary faces
      if (faceIndex == 3) {
        const unsigned faceGeom = msh->GetElementFaceType (iel, jface);
        unsigned faceDofs = msh->GetElementFaceDofNumber (iel, jface, solVType);

        vector  < vector  <  double> > faceCoordinates (dim);   // A matrix holding the face coordinates rowwise.
        for (int k = 0; k < dim; k++) {
          faceCoordinates[k].resize (faceDofs);
        }
        for (unsigned i = 0; i < faceDofs; i++) {
          unsigned inode = msh->GetLocalFaceVertexIndex (iel, jface, i);   // face-to-element local node mapping.
          for (unsigned k = 0; k < dim; k++) {
            faceCoordinates[k][i] =  coordX[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
          }
        }
        for (unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solVType]->GetGaussPointNumber(); ig++) {
          // We call the method GetGaussPointNumber from the object finiteElement in the mesh object msh.
          vector < double> normal;
          msh->_finiteElement[faceGeom][solVType]->JacobianSur (faceCoordinates, ig, weight, phiV, phiV_x, normal);

          vector< double > xg (dim, 0.);
          for (unsigned i = 0; i < faceDofs; i++) {
            for (unsigned k = 0; k < dim; k++) {
              xg[k] += phiV[i] * faceCoordinates[k][i]; // xg(ig)= \sum_{i=0}^faceDofs phi[i](xig) facecoordinates[i]
            }
          }
          double tau; // a(u)*grad_u\cdot normal
          SetBoundaryCondition (xg, "P", tau, faceIndex, 0.);  // return tau
          // *** phi_i loop ***
          for (unsigned i = 0; i < faceDofs; i++) {
            unsigned inode = msh->GetLocalFaceVertexIndex (iel, jface, i);
            for (unsigned k = 0; k < dim; k++) {
              aResV[k][inode] +=  -phiV[i] * tau * normal[k] * weight;
            }
          }
        }
      }
    }

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solVType]->Jacobian (coordX, ig, weight, phiV, phiV_x);
      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi (ig);

      vector < adept::adouble > solV_gss (dim, 0);
      vector < vector < adept::adouble > > gradSolV_gss (dim);

      for (unsigned  k = 0; k < dim; k++) {
        gradSolV_gss[k].assign (dim, 0.);
      }

      for (unsigned i = 0; i < nDofsV; i++) {
        for (unsigned  k = 0; k < dim; k++) {
          solV_gss[k] += solV[k][i] * phiV[i];
        }
        for (unsigned j = 0; j < dim; j++) {
          for (unsigned  k = 0; k < dim; k++) {
            gradSolV_gss[k][j] += solV[k][i] * phiV_x[i * dim + j];
          }
        }
      }

      adept::adouble solP_gss = 0;
      for (unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }

      double nu = 1;

      // *** phiV_i loop ***
      for (unsigned i = 0; i < nDofsV; i++) {
        vector < adept::adouble > NSV (dim, 0.);

        for (unsigned  k = 0; k < dim; k++) { //momentum equation in k
          for (unsigned j = 0; j < dim; j++) { // second index j in each equation
            NSV[k]   +=  nu * phiV_x[i * dim + j] * (gradSolV_gss[k][j] + gradSolV_gss[j][k]); // laplace
            NSV[k]   +=  phiV[i] * (solV_gss[j] * gradSolV_gss[k][j]); // non-linear term
          }
          NSV[k] += -solP_gss * phiV_x[i * dim + k]; // pressure gradient
        }
        for (unsigned  k = 0; k < dim; k++) {
          aResV[k][i] += NSV[k] * weight;
        }
      } // end phiV_i loop

      // *** phiP_i loop ***
      for (unsigned i = 0; i < nDofsP; i++) {
        for (int k = 0; k < dim; k++) {
          aResP[i] += - (gradSolV_gss[k][k]) * phiP[i]  * weight;
        }
      } // end phiP_i loop

    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store them in RES
    Res.resize (nDofsVP);   //resize

    for (int i = 0; i < nDofsV; i++) {
      for (unsigned  k = 0; k < dim; k++) {
        Res[ k * nDofsV + i ] = -aResV[k][i].value();
      }
    }

    for (int i = 0; i < nDofsP; i++) {
      Res[ dim * nDofsV + i ] = -aResP[i].value();
    }

    RES->add_vector_blocked (Res, sysDof);

    //Extract and store the Jacobian

    Jac.resize (nDofsVP * nDofsVP);
    // define the dependent variables

    for (unsigned  k = 0; k < dim; k++) {
      s.dependent (&aResV[k][0], nDofsV);
    }
    s.dependent (&aResP[0], nDofsP);

    // define the independent variables
    for (unsigned  k = 0; k < dim; k++) {
      s.independent (&solV[k][0], nDofsV);
    }
    s.independent (&solP[0], nDofsP);

    // get the and store jacobian matrix (row-major)
    s.jacobian (&Jac[0] , true); // This is rowwise order.
    KK->add_matrix_blocked (Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();

  PetscViewer    viewer;
  PetscViewerDrawOpen (PETSC_COMM_WORLD, NULL, NULL, 0, 0, 1800, 1800, &viewer);
  PetscObjectSetName ( (PetscObject) viewer, "NS matrix");
  PetscViewerPushFormat (viewer, PETSC_VIEWER_DRAW_LG);
  MatView ( (static_cast< PetscMatrix* > (KK))->mat(), viewer);

  double a;
  std::cin >> a;

  // ***************** END ASSEMBLY *******************
}


