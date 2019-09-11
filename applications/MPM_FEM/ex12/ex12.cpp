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

SparseMatrix *Proj;

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
  else{
    dirichlet = false;
    value = 0.;    
  }

  return dirichlet;
}


double InitalValueU (const std::vector < double >& x) {
  return x[0] + x[1];
}

void BuidProjection (MultiLevelProblem& ml_prob);
void AssembleWithProjection (MultiLevelProblem& ml_prob);
void ProjectSolutionIntoGradient (MultiLevelProblem& ml_prob);
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
    maxNumberOfMeshes = 2;
  }
  else {
    maxNumberOfMeshes = 2;
  }

 // std::cout<<"BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB";
  
  vector < double > l2Norm (maxNumberOfMeshes);
  vector < double > semiNorm (maxNumberOfMeshes);

  for (unsigned i = maxNumberOfMeshes - 1; i < maxNumberOfMeshes; i++) {   // loop on the mesh level

    unsigned numberOfUniformLevels = i ;
    unsigned numberOfSelectiveLevels = 0;
    mlMsh.RefineMesh (numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

    // erase all the coarse mesh levels
    mlMsh.EraseCoarseLevels (numberOfUniformLevels - 1);

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

    BuidProjection (mlProb);

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

    NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("Poisson");
    for (unsigned k = 0; k < dim; k++) {
      system.AddSolutionToSystemPDE (Uname[k].c_str());
    }
    system.AddSolutionToSystemPDE ("P");
    system.init();

    // attach the assembling function to system
    system.SetAssembleFunction (AssembleWithProjection);

    // initilaize and solve the system
    system.SetOuterSolver (PREONLY);
    system.MGsolve();

    ProjectSolutionIntoGradient (mlProb);

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
    delete Proj;
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


void BuidProjection (MultiLevelProblem& ml_prob) {

  double scale = 4.;

  adept::Stack& s = FemusInit::_adeptStack;

  MultiLevelSolution*  mlSol = ml_prob._ml_sol;
  unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1u;

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

  Proj = SparseMatrix::build().release();
  Proj->init (dim2 + dim + 1, dim + 1, P);
  
  double tolerance = 1.0e-14;
  Proj->RemoveZeroEntries (tolerance);

  MatView ( (static_cast< PetscMatrix* > (Proj))->mat(), PETSC_VIEWER_STDOUT_WORLD);

}

unsigned counter = 0;

void AssembleWithProjection (MultiLevelProblem& ml_prob) {

  ProjectSolutionIntoGradient (ml_prob);

  adept::Stack& s = FemusInit::_adeptStack;

  MultiLevelSolution*  mlSol = ml_prob._ml_sol;
  unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1u;

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

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofsU; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, solTypeU);
      for (int k = 0; k < dim2 + dim; k++) {
        sol[k][i] = (*solution->_Sol[solIndex[k]]) (solDof);
        sysDof[k * nDofsU + i] = sys0->GetSystemDof (solIndex[k], solPdeIndex[k], i, iel);
      }
    }

    for (unsigned i = 0; i < nDofsP; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, solTypeP);
      sol[dim2 + dim][i] = (*solution->_Sol[solIndex[dim2 + dim]]) (solDof);
      sysDof[ (dim2 + dim) * nDofsU + i] = sys0->GetSystemDof (solIndex[dim2 + dim], solPdeIndex[dim2 + dim], i, iel);
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
            solUx[k][l] += phi[i] * sol[ k * dim + l][i]; //extended solution gradient 
            solux[k][l] += phi_x[i * dim + l] * sol[ dim2 + k][i]; //regular solution gradient
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

      double ReI = 0.001;
      bool extended = true;
      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofsU; i++) {
        for (unsigned k = 0; k < dim; k++) {
          if(extended){
            for (unsigned l = 0; l < dim; l++) {  
              aRes[k * dim + l][i] += ReI * (solux[k][l] + solux[l][k]) * phi[i] * weight; //extended diffusivity 
              //aRes[dim2 + k][i] += ReI * (solUx[k][l] + solUx[l][k]) * phi_x[ i * dim + l] * weight; //mixed 1 diffusivity
              //aRes[k * dim + l][i] += ReI * (solux[k][l] + solux[l][k]) * phi[i] * weight; //mixed 2 diffusivity
              aRes[dim2 + k][i] += phi[i] * solU[l] * solux[k][l] * weight; //extended convection
            }
            aRes[k * dim + k][i] -= solP * phi[i] * weight; //extended gradient
          }
          else{
            for (unsigned l = 0; l < dim; l++) {  
              aRes[dim2 + k][i] += ReI * (solux[k][l] + solux[l][k]) * phi_x[ i * dim + l] * weight; //regular diffusivity
              aRes[dim2 + k][i] += phi[i] * solU[l] * solUx[k][l] * weight; //regular convection
            }
            aRes[dim2 + k][i] -= solP * phi_x[ i * dim + k]  * weight; //regular gradient
          }
        }
      } // end phi_i loop

      for (unsigned i = 0; i < nDofsP; i++) {
        for (unsigned k = 0; k < dim; k++) {
          if(extended){
            aRes[dim2 + dim][i] -= phiP[i] * solux[k][k] * weight; //extended divergence  
          }
          else{
            aRes[dim2 + dim][i] -= phiP[i] * solUx[k][k] * weight; //regular divergence  
          }
        }
      }
    } // end gauss point loop

    res.resize ( (dim2 + dim) * nDofsU + nDofsP);
    for (unsigned k = 0; k < dim2 + dim; k++) {
      for (int i = 0; i < nDofsU; i++) {
        res[k * nDofsU + i] = -aRes[k][i].value();
        //std::cout << k * nDofsU + i <<" "<< res[k * nDofsU + i]<<std::endl;
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

//     for(unsigned i = 0;i <( (dim2 + dim) * nDofsU + nDofsP);i++){
//       for(unsigned j = 0;j <( (dim2 + dim) * nDofsU + nDofsP);j++){
//          std::cout<<Jac[i*( (dim2 + dim) * nDofsU + nDofsP)+j]<<" "; 
//       }
//       std::cout<<std::endl;
//     }
    
    K->add_matrix_blocked (Jac, sysDof, sysDof);
    s.clear_dependents();
    s.clear_independents();
  } //end element loop for each process*/

  K->close();
  RES->close();

  NonLinearImplicitSystem* mlSys1 =  &ml_prob.get_system< NonLinearImplicitSystem > ("Poisson");
  LinearEquationSolver* sys1 = mlSys1->_LinSolver[level];

  SparseMatrix *K1 = sys1->_KK; ;
  NumericVector* RES1  = sys1->_RES;

  bool mathReuse = (counter == 0) ? false : true;
  K1->matrix_PtAP (*Proj, *K, mathReuse);
  RES1->matrix_mult_transpose (*RES, *Proj);
  counter++;
  

//  MatView((static_cast< PetscMatrix* > (K1))->mat(),PETSC_VIEWER_STDOUT_WORLD);
//
//  VecView((static_cast< PetscVector* > (RES))->vec(),PETSC_VIEWER_STDOUT_WORLD);
//
//   VecView((static_cast< PetscVector* > (RES1))->vec(),PETSC_VIEWER_STDOUT_WORLD);



}

void ProjectSolutionIntoGradient (MultiLevelProblem& ml_prob) {

  MultiLevelSolution*  mlSol = ml_prob._ml_sol;
  unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1u;

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



