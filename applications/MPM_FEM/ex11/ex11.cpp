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

bool SetBoundaryCondition (const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0;

  //if (facename == 2)
   // dirichlet = false;

  return dirichlet;
}


void AssembleStandardProblem (MultiLevelProblem& ml_prob);
void BuidProjection (MultiLevelProblem& ml_prob);
void AssembleWithProjection (MultiLevelProblem& ml_prob);

void ProjectSolutionIntoGradient (MultiLevelProblem& ml_prob);

std::pair < double, double > GetErrorNorm (MultiLevelSolution* mlSol);

int main (int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit (argc, args, MPI_COMM_WORLD);


  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  
  mlMsh.GenerateCoarseBoxMesh(2., 0, 0, -0.5, 0.5, 0., 0., 0., 0., EDGE3, "seventh");
  
  //mlMsh.ReadCoarseMesh ("./input/square_quad.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/square_tri.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/square_mixed.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cube_hex.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cube_wedge.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cube_tet.neu","seventh",scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cube_mixed.neu","seventh",scalingFactor);

  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
   *    probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();
  unsigned maxNumberOfMeshes;

  if (dim == 1) {
    maxNumberOfMeshes = 5;
  }
  else if (dim == 2) {
    maxNumberOfMeshes = 7;
  }
  else {
    maxNumberOfMeshes = 6;
  }

  vector < vector < double > > l2Norm;
  l2Norm.resize (maxNumberOfMeshes);

  vector < vector < double > > semiNorm;
  semiNorm.resize (maxNumberOfMeshes);

  for (unsigned i = 1; i < maxNumberOfMeshes; i++) {   // loop on the mesh level

    unsigned numberOfUniformLevels = i ;
    unsigned numberOfSelectiveLevels = 0;
    mlMsh.RefineMesh (numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

    // erase all the coarse mesh levels
    mlMsh.EraseCoarseLevels (numberOfUniformLevels - 1);

    // print mesh info
    mlMsh.PrintInfo();

    FEOrder feOrder[3] = {FIRST, SERENDIPITY, SECOND};
    l2Norm[i].resize (3);
    semiNorm[i].resize (3);

    for (unsigned j = 0; j < 3; j++) {   // loop on the FE Order
      // define the multilevel solution and attach the mlMsh object to it
      MultiLevelSolution mlSol (&mlMsh);



      // add variables to mlSol
      mlSol.AddSolution ("u", LAGRANGE, feOrder[j]);

      std::string Uxname[3] = {"ux", "uy", "uz"};
      for (unsigned k = 0; k < dim; k++) {
        mlSol.AddSolution (Uxname[k].c_str(), LAGRANGE, feOrder[j]);
      }

      mlSol.AddSolution ("weight", LAGRANGE, feOrder[j], 0);

      mlSol.Initialize ("All");

      // attach the boundary condition function and generate boundary data
      mlSol.AttachSetBoundaryConditionFunction (SetBoundaryCondition);
      mlSol.GenerateBdc ("All");

      // define the multilevel problem attach the mlSol object to it
      MultiLevelProblem mlProb (&mlSol);



      LinearImplicitSystem* systemP[4];
      std::string Pname[3] = {"Px", "Py", "Pz"};
      for (unsigned k = 0; k < dim; k++) {
        systemP[k] = &mlProb.add_system < LinearImplicitSystem > (Pname[k]);
        systemP[k]->AddSolutionToSystemPDE (Uxname[k].c_str());
        systemP[k]->init();
      }
      systemP[dim] = &mlProb.add_system < LinearImplicitSystem > ("Pu");
      systemP[dim]->AddSolutionToSystemPDE ("u");
      systemP[dim]->init();
      Proj = SparseMatrix::build().release();
       
      BuidProjection (mlProb);

      LinearImplicitSystem& system0 = mlProb.add_system < LinearImplicitSystem > ("auxiliary");

      for (unsigned k = 0; k < dim; k++) {
        system0.AddSolutionToSystemPDE (Uxname[k].c_str());
      }
      system0.AddSolutionToSystemPDE ("u");
      system0.init();

      NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("Poisson");
      system.AddSolutionToSystemPDE ("u");
      system.init();
      
      
      
      // attach the assembling function to system
      //system.SetAssembleFunction (AssembleWithProjection);
      system.SetAssembleFunction (AssembleStandardProblem);

      // initilaize and solve the system

      system.SetOuterSolver (PREONLY);
      system.MGsolve();
         
      ProjectSolutionIntoGradient (mlProb);  
      

      std::pair< double , double > norm = GetErrorNorm (&mlSol);
      l2Norm[i][j]  = norm.first;
      semiNorm[i][j] = norm.second;
      // print solutions
      std::vector < std::string > variablesToBePrinted;
      variablesToBePrinted.push_back ("All");

      VTKWriter vtkIO (&mlSol);
      vtkIO.SetDebugOutput (true);
      vtkIO.Write (DEFAULT_OUTPUTDIR, "linear", variablesToBePrinted, i + j * 10);
      vtkIO.Write (DEFAULT_OUTPUTDIR, "quadratic", variablesToBePrinted, i + j * 10);
      vtkIO.Write (DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, i + j * 10);

      delete Proj;
      
    }
  }

  // print the seminorm of the error and the order of convergence between different levels
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "l2 ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (unsigned i = 1; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision (14);

    for (unsigned j = 0; j < 3; j++) {
      std::cout << l2Norm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < maxNumberOfMeshes - 1) {
      std::cout.precision (3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < 3; j++) {
        std::cout << log (l2Norm[i][j] / l2Norm[i + 1][j]) / log (2.) << "\t\t\t";
      }

      std::cout << std::endl;
    }

  }

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "SEMINORM ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (unsigned i = 1; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision (14);

    for (unsigned j = 0; j < 3; j++) {
      std::cout << semiNorm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < maxNumberOfMeshes - 1) {
      std::cout.precision (3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < 3; j++) {
        std::cout << log (semiNorm[i][j] / semiNorm[i + 1][j]) / log (2.) << "\t\t\t";
      }
      std::cout << std::endl;
    }
  }
}


double GetExactSolutionValue (const std::vector < double >& x) {
  double pi = acos (-1.);
  if(x.size() == 1){  
    return cos (pi * x[0]);
  }
  else{
    return cos (pi * x[0]) * cos (pi * x[1]);    
  }
};


void GetExactSolutionGradient (const std::vector < double >& x, vector < double >& solGrad) {
  double pi = acos (-1.);
  if(x.size() == 1){
    solGrad[0]  = -pi * sin (pi * x[0]);
  }
  else{
    solGrad[0]  = -pi * sin (pi * x[0]) * cos (pi * x[1]);
    solGrad[1] = -pi * cos (pi * x[0]) * sin (pi * x[1]); 
  }
};


double GetExactSolutionLaplace (const std::vector < double >& x) {
  double pi = acos (-1.);
  if(x.size() == 1){
    return -pi * pi * cos (pi * x[0]);     // - pi*pi*cos(pi*x[0])*cos(pi*x[1]);
  }
  else{
    return -2.*pi * pi * cos (pi * x[0]) * cos (pi * x[1]);     // - pi*pi*cos(pi*x[0])*cos(pi*x[1]);  
  }
};

/**
 * Given the non linear problem
 *
 *      - \Delta u + < u, u, u > \cdot \nabla u = f(x)
 * in the unit box centered in the origin with
 *                      f(x) = - \Delta u_e + < u_e, u_e, u_e > \cdot \nabla u_e
 *                    u_e = \cos ( \pi * x ) * \cos( \pi * y ),
 * the following function assembles (explicitely) the Jacobian matrix J(u^i) and the residual vector Res(u^i)
 * for the Newton iteration, i.e.
 *
 *                  J(u^i) w = Res(u^i) = f(x) - ( - \Delta u^i + < u^i, u^i, u^i > \cdot \nabla u^i ),
 *        u^{i+1} = u^i + w
 *        where
 *        J(u^i) w = - \Delta w  + < w , w , w >  \cdot \nabla u^i + < u^i , u^i , u^i >  \cdot \nabla w
 *
 **/

void AssembleStandardProblem (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use

  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system< NonLinearImplicitSystem > ("Poisson");   // pointer to the linear implicit system named "Poisson"
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

//solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex ("u");   // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType (soluIndex);   // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex ("u");   // get the position of "u" in the pdeSys object

  vector < adept::adouble >  solu; // local solution

  vector < vector < double > > x (dim);   // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector< int > sysDof; // local to global pdeSys dofs
  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  double weight; // gauss point weight

  vector< double > Res; // local redidual vector
  vector< adept::adouble > aRes; // local redidual vector


// reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned > (ceil (pow (3, dim)));       // conservative: based on line3, quad9, hex27
  solu.reserve (maxSize);

  for (unsigned i = 0; i < dim; i++)
    x[i].reserve (maxSize);

  sysDof.reserve (maxSize);
  phi.reserve (maxSize);
  phi_x.reserve (maxSize * dim);
  Res.reserve (maxSize);
  aRes.reserve (maxSize);

  vector < double > Jac; // local Jacobian matrix (ordered by column, adept)
  Jac.reserve (maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Matrix

// element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofs  = msh->GetElementDofNumber (iel, soluType);   // number of solution element dofs
    // resize local arrays
    sysDof.resize (nDofs);
    solu.resize (nDofs);

    for (int i = 0; i < dim; i++) {
      x[i].resize (nDofs);
    }

    Res.assign (nDofs, 0.);  //resize
    aRes.assign (nDofs, 0.);  //resize
    Jac.resize (nDofs * nDofs);


    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, soluType);   // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex]) (solDof);     // global extraction and local storage for the solution
      sysDof[i] = pdeSys->GetSystemDof (soluIndex, soluPdeIndex, i, iel);   // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);   // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim]) (xDof);     // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian (x, ig, weight, phi, phi_x);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      adept::adouble soluGauss = 0;
      vector < adept::adouble > soluGauss_x (dim, 0.);
      vector < double > xGauss (dim, 0.);

      for (unsigned i = 0; i < nDofs; i++) {
        soluGauss += phi[i] * solu[i];

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          soluGauss_x[jdim] += phi_x[i * dim + jdim] * solu[i];
          xGauss[jdim] += x[jdim][i] * phi[i];
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {

        adept::adouble mLaplace = 0.;

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          mLaplace   +=  phi_x[i * dim + jdim] * soluGauss_x[jdim];
        }

        double exactSolValue = GetExactSolutionValue (xGauss);
        vector < double > exactSolGrad (dim);
        GetExactSolutionGradient (xGauss , exactSolGrad);
        double exactSolLaplace = GetExactSolutionLaplace (xGauss);

        double f = (- exactSolLaplace) * phi[i] ;
        //double f = - 1 * phi[i] ;
        aRes[i] += (f - (mLaplace)) * weight;

      } // end phi_i loop
    } // end gauss point loop


    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    for (int i = 0; i < nDofs; i++) {
      Res[i] = -aRes[i].value();
    }
    RES->add_vector_blocked (Res, sysDof);

    // define the dependent variables
    s.dependent (&aRes[0], nDofs);

    // define the independent variables
    s.independent (&solu[0], nDofs);

    // get the jacobian matrix (ordered by column)
    s.jacobian (&Jac[0], true);

    //store Jact in the global matrix KK
    KK->add_matrix_blocked (Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();
  KK->close();

// ***************** END ASSEMBLY *******************
}


std::pair < double, double > GetErrorNorm (MultiLevelSolution* mlSol) {
  unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1u;
  //  extract pointers to the several objects that we are going to use
  Mesh*          msh          = mlSol->_mlMesh->GetLevel (level);   // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)
  Solution*    sol        = mlSol->GetSolutionLevel (level);   // pointer to the solution (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

//solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex ("u");   // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType (soluIndex);   // get the finite element type for "u"

  vector < double >  solu; // local solution

  vector < vector < double > > x (dim);   // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

// reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned > (ceil (pow (3, dim)));       // conservative: based on line3, quad9, hex27
  solu.reserve (maxSize);

  for (unsigned i = 0; i < dim; i++)
    x[i].reserve (maxSize);

  phi.reserve (maxSize);
  phi_x.reserve (maxSize * dim);
  unsigned dim2 = (3 * (dim - 1) + ! (dim - 1));       // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  phi_xx.reserve (maxSize * dim2);

  double seminorm = 0.;
  double l2norm = 0.;

// element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {


    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofs  = msh->GetElementDofNumber (iel, soluType);   // number of solution element dofs
    unsigned nDofs2 = msh->GetElementDofNumber (iel, xType);   // number of coordinate element dofs

    // resize local arrays
    solu.resize (nDofs);

    for (int i = 0; i < dim; i++) {
      x[i].resize (nDofs2);
    }

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, soluType);   // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex]) (solDof);     // global extraction and local storage for the solution
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofs2; i++) {
      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);   // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim]) (xDof);     // global extraction and local storage for the element coordinates
      }
    }

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian (x, ig, weight, phi, phi_x, phi_xx);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      double soluGauss = 0;
      vector < double > soluGauss_x (dim, 0.);
      vector < double > xGauss (dim, 0.);

      for (unsigned i = 0; i < nDofs; i++) {
        soluGauss += phi[i] * solu[i];

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          soluGauss_x[jdim] += phi_x[i * dim + jdim] * solu[i];
          xGauss[jdim] += x[jdim][i] * phi[i];
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

  adept::Stack& s = FemusInit::_adeptStack;

  MultiLevelSolution*  mlSol = ml_prob._ml_sol;
  unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1u;

  Solution* sol = ml_prob._ml_sol->GetSolutionLevel (level);
  Mesh* msh = ml_prob._ml_msh->GetLevel (level);
  elem* el = msh->el;

  unsigned  dim = msh->GetDimension();

  std::vector < LinearImplicitSystem* > mlSysP (1 + dim);
  std::vector < LinearEquationSolver* > sysP (1 + dim);
  std::vector < SparseMatrix*> P (1 + dim);

  std::string Pname[3] = {"Px", "Py", "Pz"};
  std::string Uxname[3] = {"ux", "uy", "uz"};
  std::vector < unsigned > soluxIndex (dim);
  for (unsigned k = 0; k < dim; k++) {
    mlSysP[k] =  &ml_prob.get_system< LinearImplicitSystem > (Pname[k]);
    sysP[k] = mlSysP[k]->_LinSolver[level];
    soluxIndex[k] = mlSol->GetIndex (Uxname[k].c_str());
    P[k] = sysP[k]->_KK;
    P[k]->zero();
  }

  mlSysP[dim] =  &ml_prob.get_system< LinearImplicitSystem > ("Pu");
  sysP[dim] = mlSysP[dim]->_LinSolver[level];
  P[dim] = sysP[dim]->_KK;
  P[dim]->zero();
  NumericVector* D  = sysP[dim]->_RES;
  *D = 1.;
  P[dim]->matrix_set_diagonal_values (*D);

  unsigned solwIndex = mlSol->GetIndex ("weight");
  unsigned solType = mlSol->GetSolutionType (solwIndex);

  sol->_Sol[solwIndex]->zero();

  vector < vector < double > > x (dim);
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector< int > sysDof;
  vector <double> phi;
  vector <double> phi_x;
  double weight;

  unsigned    iproc = msh->processor_id();

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofs  = msh->GetElementDofNumber (iel, solType);
    sysDof.resize (nDofs);
    for (int k = 0; k < dim; k++) {
      x[k].resize (nDofs);
    }
    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
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
      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {
        sol->_Sol[solwIndex]->add (sysDof[i], phi[i] * weight);
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
  std::vector < adept::adouble > solw;
  vector <double> Jac;
  std::vector < std::vector< adept::adouble > > aRes (dim); // local redidual vector


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
      std::vector < adept::adouble > solux_g (dim, 0.);
      for (unsigned i = 0; i < nDofs; i++) {
        for (unsigned k = 0; k < dim; k++) {
          solux_g[k] += phi_x[i * dim + k] * solu[i];
        }
      }
      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {
        for (unsigned k = 0; k < dim; k++) {
          aRes[k][i] += solux_g[k] * phi[i] * weight / solw[i];
        }
      } // end phi_i loop
    } // end gauss point loop

    s.independent (&solu[0], nDofs);
    for (unsigned k = 0; k < dim; k++) {
      s.dependent (&aRes[k][0], nDofs);
      s.jacobian (&Jac[0], true);
      P[k]->add_matrix_blocked (Jac, sysDof, sysDof);
      s.clear_dependents();
    }
    s.clear_independents();
  } //end element loop for each process*/

  for (unsigned k = 0; k < dim; k++) {
    P[k]->close();
  }

  Proj->init (dim + 1, 1, P);
  
  //SparseMatrix* nProj;
  //nProj = SparseMatrix::build().release();
  //nProj->init (dim + 1, 1, P);
   
  //Mat pnProj = (static_cast< PetscMatrix* > (nProj))->mat();
  Mat psProj = (static_cast< PetscMatrix* > (Proj))->mat();
//   
  
  MatConvert( psProj, MATSEQAIJ, MAT_INPLACE_MATRIX, &psProj);
  
  //std::cout << "AAAAAAAAAAAAAAAAAAAA"<<std::endl;
    
  //delete nProj;

//
  //PetscViewer viewer;
  //PetscViewerDrawOpen (PETSC_COMM_WORLD, NULL, NULL, 0, 0, dim * 300, 300, &viewer);
  //PetscObjectSetName ( (PetscObject) viewer, "P matrix");
  //PetscViewerPushFormat (viewer, PETSC_VIEWER_DRAW_LG);
  //MatView ((static_cast< PetscMatrix* > (Proj))->mat(), viewer);
//   MatView((static_cast< PetscMatrix* > (Proj))->mat(),PETSC_VIEWER_STDOUT_WORLD);

  //double a;
  //std::cin>>a;




}

void AssembleWithProjection (MultiLevelProblem& ml_prob) {

  ProjectSolutionIntoGradient (ml_prob);  
    
  adept::Stack& s = FemusInit::_adeptStack;

  MultiLevelSolution*  mlSol = ml_prob._ml_sol;
  unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1u;

  Solution* solution = ml_prob._ml_sol->GetSolutionLevel (level);
  Mesh* msh = ml_prob._ml_msh->GetLevel (level);
  elem* el = msh->el;

  unsigned  dim = msh->GetDimension();

  LinearImplicitSystem* mlSys0 =  &ml_prob.get_system< LinearImplicitSystem > ("auxiliary");
  LinearEquationSolver* sys0 = mlSys0->_LinSolver[level];

  SparseMatrix *K = sys0->_KK; ;
  NumericVector* RES  = sys0->_RES;

  K->zero();
  RES->zero();


  std::string Uxname[3] = {"ux", "uy", "uz"};
  std::vector < unsigned > solIndex (dim + 1);
  for (unsigned k = 0; k < dim; k++) {
    solIndex[k] = mlSol->GetIndex (Uxname[k].c_str());
  }
  solIndex[dim] = mlSol->GetIndex ("u");
  
  unsigned solType = mlSol->GetSolutionType (solIndex[0]);

  vector < vector < double > > x (dim);
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector< unsigned > sysDof;
  vector <double> phi;
  vector <double> phi_x;
  double weight;

  unsigned  iproc = msh->processor_id();

  vector < unsigned > solPdeIndex (dim + 1);
  for (unsigned k = 0; k < dim; k++) {
    solPdeIndex[k] = mlSys0->GetSolPdeIndex (Uxname[k].c_str());
  }
  solPdeIndex[dim] = mlSys0->GetSolPdeIndex ("u");

  std::vector < std::vector < adept::adouble > > sol (dim + 1);
  std::vector < std::vector< adept::adouble > > aRes (dim + 1); // local residual vector

  vector <double> Jac;
  std::vector< double > res (dim + 1);

  //BEGIN element loop
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofs  = msh->GetElementDofNumber (iel, solType);
   
    for (int k = 0; k < dim + 1; k++) {
      sol[k].resize (nDofs);
      aRes[k].assign (nDofs, 0.);
    }

    for (int k = 0; k < dim; k++) {
      x[k].resize (nDofs);
    }
    sysDof.resize (nDofs * (dim + 1));
    
    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, solType);
      for (int k = 0; k < dim + 1; k++) {
        sol[k][i] = (*solution->_Sol[solIndex[k]]) (solDof);
        sysDof[k * nDofs + i] = sys0->GetSystemDof (solIndex[k], solPdeIndex[k], i, iel);
        //std::cout << k * nDofs + i << " "<< sysDof[k * nDofs + i] <<std::endl;
      }
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
      adept::adouble solu = 0.;
      std::vector < adept::adouble > solux (dim,0.);
      vector < double > xGauss (dim, 0.);
      for (unsigned i = 0; i < nDofs; i++) {
        solu += phi[i] * sol[dim][i];
        for (unsigned k = 0; k < dim; k++) {
          solux[k] += phi[i] * sol[k][i];
          xGauss[k] += x[k][i] * phi[i];
        }
      }
      
      double exactSolValue = GetExactSolutionValue (xGauss);
      vector < double > exactSolGrad (dim);
      GetExactSolutionGradient (xGauss , exactSolGrad);
      double exactSolLaplace = GetExactSolutionLaplace (xGauss);
      
      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {
        for (unsigned k = 0; k < dim; k++) {
          aRes[k][i] += -solux[k] * phi[i] * weight;
        }
        aRes[dim][i] += -exactSolLaplace * phi[i] * weight;
        //aRes[dim][i] += -1 * phi[i] * weight;
      } // end phi_i loop
    } // end gauss point loop

    res.resize (nDofs * (dim + 1));
    for (unsigned k = 0; k < dim + 1; k++) {
      for (int i = 0; i < nDofs; i++) {
        res[k * nDofs + i] = -aRes[k][i].value();
      }
    }
    RES->add_vector_blocked (res, sysDof);

    
    Jac.resize (nDofs * (dim + 1) * nDofs * (dim + 1));
    for (unsigned k = 0; k < dim + 1; k++) {
      s.independent (&sol[k][0], nDofs);
    }
    for (unsigned k = 0; k < dim + 1; k++) {
      s.dependent (&aRes[k][0], nDofs);
    }
    s.jacobian (&Jac[0], true);
    
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
  
  K1->matrix_PtAP (*Proj, *K, false);
  RES1->matrix_mult_transpose ( *RES, *Proj);
  
 // MatView((static_cast< PetscMatrix* > (K1))->mat(),PETSC_VIEWER_STDOUT_WORLD);
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

  std::vector < LinearImplicitSystem* > mlSysP (dim + 1);
  std::vector < LinearEquationSolver* > sysP (dim + 1);
  std::vector < SparseMatrix*> P (dim + 1);

  std::string Pname[3] = {"Px", "Py", "Pz"};
  std::string Uxname[3] = {"ux", "uy", "uz"};
  std::vector < unsigned > solIndex (dim);
  for (unsigned k = 0; k < dim; k++) {
    mlSysP[k] =  &ml_prob.get_system< LinearImplicitSystem > (Pname[k]);
    sysP[k] = mlSysP[k]->_LinSolver[level];
    solIndex[k] = mlSol->GetIndex (Uxname[k].c_str());
    P[k] = sysP[k]->_KK;
  }

  unsigned soluIndex = mlSol->GetIndex ("u");

//   mlSysP[dim] =  &ml_prob.get_system< LinearImplicitSystem > ("Pu");
//   sysP[dim] = mlSysP[dim]->_LinSolver[level];
//   solIndex[dim] = mlSol->GetIndex ("u");
//   P[dim] = sysP[dim]->_KK;
//   (*sol->_Sol[solIndex[0]]).matrix_mult ( (*sol->_Sol[soluIndex]), (*P[dim]));
//   (*sol->_Sol[solIndex[dim]])=(*sol->_Sol[solIndex[0]]);

  for (unsigned k = 0; k < dim; k++) {
    (*sol->_Sol[solIndex[k]]).matrix_mult ( (*sol->_Sol[soluIndex]), (*P[k]));
  }



}























//#include "./include/gmpm.hpp"
// int main (int argc, char** args) {
//
//   FemusInit (argc, args, MPI_COMM_WORLD);
//
//   //std::vector < double > x {0., 0.43, 0.81, 1.25, 1.6};
//   unsigned nve = 5;
//   double L = 2.;
//   std::vector < double > x (nve);
//
//   double H = L / (nve - 1.);
//
//   x[0] = 0.;
//   for (unsigned i = 1; i < nve; i++) {
//     x[i] = x[i - 1] + H;
//   }
//
//   std::vector < std::vector < double > > P (nve);
//   std::vector < std::vector < double > > M (nve);
//   std::vector < std::vector < double > > K (nve);
//   std::vector < std::vector < double > > MP (nve);
//   std::vector < std::vector < double > > Kt (nve);
//
//
//   for (unsigned i = 0; i < nve; i++) {
//     P[i].assign (nve, 0.);
//     M[i].assign (nve, 0.);
//     K[i].assign (nve + 1, 0.);
//     MP[i].assign (nve, 0.);
//     Kt[i].assign (nve + 1, 0.);
//   }
//
//   for (unsigned iel = 0; iel < nve - 1; iel++) {
//
//     unsigned i0 = iel;
//     unsigned i1 = iel + 1;
//     double x0 = x[i0];
//     double x1 = x[i1];
//     double h = x1 - x0;
//
//     M[i0][i0] += 1. / 6. * h * 2;
//     M[i0][i1] += 1. / 6. * h;
//     M[i1][i0] += 1. / 6. * h;
//     M[i1][i1] += 1. / 6. * h * 2;
//
//     K[i0][i0] +=  1. / h;
//     K[i0][i1] += -1. / h;
//     K[i1][i0] += -1. / h;
//     K[i1][i1] +=  1. / h;
//
//   }
//   P[0][0] = -1. / (x[1] - x[0]);
//   P[0][1] =  1. / (x[1] - x[0]);
//   for (unsigned i = 1; i < nve - 1; i++) {
//     double h = x[i + 1] - x[i - 1];
//     P[i][i - 1] = -1. / h;
//     P[i][i + 1] =  1. / h;
//   }
//   P[nve - 1][nve - 2] = -1. / (x[nve - 1] - x[nve - 2]);
//   P[nve - 1][nve - 1] =  1. / (x[nve - 1] - x[nve - 2]);
//
//   for (unsigned i = 0; i < nve; i++) {
//     for (unsigned j = 0; j < nve; j++) {
//       for (unsigned k = 0; k < nve; k++) {
//         MP[i][j] += M[i][k] * P[k][j];
//       }
//     }
//   }
//
//   for (unsigned i = 0; i < nve; i++) {
//     for (unsigned j = 0; j < nve; j++) {
//       for (unsigned k = 0; k < nve; k++) {
//         Kt[i][j] += P[k][i] * MP[k][j];
//       }
//     }
//   }
//
//   /*//point load
//   for(unsigned j = 1; j< nve + 1; j++){
//     K[0][j] = 0.;
//     Kt[0][j] = 0.;
//   }
//   K[nve - 1][nve] = 1.;
//   Kt[nve - 1][nve] = 1.;
//   */
//
//
//   //parabola
//   for (unsigned j = 0; j < nve + 1; j++) {
//     K[0][j] = 0.;
//     K[nve - 1][j] = 0.;
//
//     Kt[0][j] = 0.;
//     Kt[nve - 1][j] = 0.;
//   }
//   K[0][0] = 1.;
//   K[nve - 1][nve - 1] = 1.;
//   Kt[0][0] = 1.;
//   Kt[nve - 1][nve - 1] = 1.;
//
//   for (unsigned i = 1; i < nve - 1; i++) {
//     K[i][nve] = 2. * H;
//     Kt[i][nve] = 2. * H;
//   }
//
//   std::vector < double> u (nve);
//   GaussianEleminationWithPivoting (K, u, true);
//   std::cout << u[ (nve - 1) / 2] << std::endl;
//
//   GaussianEleminationWithPivoting (Kt, u, true);
//   std::cout << u[ (nve - 1) / 2] << std::endl;
//
// //   for (unsigned i = 0; i < nve; i++) {
// //     for (unsigned j = 0; j < nve + 1; j++) {
// //       std::cout << Kt[i][j] <<" ";
// //     }
// //     std::cout << std::endl;
// //   }
// //   std::cout << std::endl;
// //
//
//
// }
//
//
//
//
