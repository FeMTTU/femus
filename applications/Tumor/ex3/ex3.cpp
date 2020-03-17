/** 
 * This example shows how to set and solve the weak form of the Poisson problem
 *          $$ \dfrac{\partial u}{ \partial t}\=\nabla \cdot (a(u)\nabla u)   $$
 *          $$ \nabla u.n=-\epsilon \text{ on} \partial B(0,1) $$
 *          $$ u= V0*exp(1-1/(1-r^2)) for r< 1, u=0 for r>=1 where r=||x|| $$
 * on a sphere domain B(0,1) with boundary $\Gamma$;
 * all the coarse-level meshes are removed;
 * a multilevel problem and an equation system are initialized;
 * a direct solver is used to solve the problem.
 **/

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"
#include "Marker.hpp"

using namespace femus;

bool CheckIfPositiveDefinite (double K[6]);

void GetKFromFileISO (MultiLevelSolution &mlSol, const unsigned & split, const unsigned &patient);
void GetKFromFileANISO (MultiLevelSolution &mlSol);

double GetTimeStep (const double time) {
  double dt = .02;
  return dt;
}

bool SetBoundaryCondition (const std::vector < double >& x, const char solName[], double& value, const int faceIndex, const double time) {
  bool dirichlet = false; //Neumann
  value = 0.;

  return dirichlet;
}

double Keps = 1.0e-08;
double V0;



unsigned jInitialPosition = 0;
//std::vector <double>  xc = {0,0,0,0,0,0,0,0,0,0,-0.9,-0.4,0.4,-0.9,-0.4,0.4,-0.9,-0.4,0.4,-0.9};
//std::vector <double>  yc = {-0.9,-0.5,0.3,0.9,-0.9,-0.9,-0.9,-0.5,-0.5,-0.5,0,0,0,0,0,0,0.4,0.4,0.4,-0.9};
//std::vector <double>  zc = {0,0,0,0,0.4,-0.6,-1.2,0.4,-0.6,-1.2,-0.4,-0.4,-0.4,-1.2,-1.2,-1.2,-0.4,-0.4,-0.4,-0.4};


double InitalValueU3D (const std::vector < double >& x) {
    
  unsigned j = jInitialPosition;
  //std::cout<< "J value is: " << j << std::endl;

  //original=(0.5,0,0), xcentered = (0.5,-0.3,-0.5), ycentered = (-0.3,0,-0.6), zcentered = (-0.3,-0.3,-0.6), badcentered = (0,0.7,0.6)
  double xc = 0.5;
  double yc = 0.0;
  double zc = 0.0;
  double r = sqrt ( (x[0] - xc) * (x[0] - xc) + (x[1] - yc) * (x[1] - yc) + (x[2] - zc) * (x[2] - zc));
  //double r = sqrt ( (x[0] - xc[j]) * (x[0] - xc[j]) + (x[1] - yc[j]) * (x[1] - yc[j]) + (x[2] - zc[j]) * (x[2] - zc[j]));
  double r2 = r * r;
  double R = 1.; //radius of the tumor
  double R2 = R * R;
  double R3 = R2 * R;
  double Vb;

  if (R == 1.) {
    Vb = 1.1990039070212866;
  }
  else if (R == 0.5) {
    Vb = 0.149875;
  }
  else if (R == 0.35) {
    Vb = 0.0514073;
  }
  else {
    std::cout << "wrong R, evaluate the integral \n $Vb = 4 \\pi \\int_0^{R} ( 1 - \\frac{R^2}{ R^2 - \\rho^2} \\rho^2 d\\rho$" << std::endl;
    exit (0);
  }
  if (r2 > R2) return 0.;
  return (V0 * M_PI * 4. / 3.* R3) / Vb * exp ( (1. - R2 / (R2 - r2)));
  
}

double InitalValueD (const std::vector < double >& x) {
  return 100.;
}

double GetSmootK (const double & kmin, const double & kmax, const double & h, const double & r0, const std::vector < double >& x) {
  double value = kmin;
  double r = sqrt (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  if (r <= r0 - h) {
    value = kmax;
  }
  else if (r <= r0 + h) {
    value = kmin + (kmax - kmin) * 0.5 * (1. - atan ( (r - r0) / h) / atan (1.));
  }
  return value;
}

bool GetDeadCells (const double &time, MultiLevelSolution &mlSol, const bool & last);

void AssemblePoissonProblem_AD (MultiLevelProblem& ml_prob);

std::pair < double, double > GetErrorNorm (MultiLevelSolution* mlSol);

int main (int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit (argc, args, MPI_COMM_WORLD);

  //const unsigned split[4] = {8, 4, 2, 1};

  double scalingFactor = 1.;
  unsigned numberOfUniformLevels = 4; //We apply uniform refinement.
  unsigned numberOfSelectiveLevels = 0; // We may want to see the solution on some levels.

  MultiLevelMesh mlMsh (numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels,
                        "./input/cube8x8x3scaled.neu", "fifth", 1., NULL);

  unsigned dim = mlMsh.GetDimension();
  // erase all the coarse mesh levels
  // mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1); // We check the solution on the finest mesh.
 //std::vector <double> doses = {0.0100000,0.01170000,0.0128074,0.01375,0.0163183,0.0207917,0.0264914,
  //   0.0337537,0.0430067,0.05,0.2,0.3,0.43,0.47,0.57,0.68,0.83,0.91,1.1,1.5,2.3,2.8,3.6,3.775, 3.95,4.125,4.3,4.475,4.65,4.825,5};
     
  for (unsigned simulation = 0; simulation < 1 ; simulation++) {
    // unsigned l_base = xc.size();
    //for (unsigned simulation = 0; simulation < l_base ; simulation++) {
        
    jInitialPosition = simulation; 
    //V0 = 0.05 * (simulation + 1) ;   
    //V0 = doses[simulation];
    V0 = 1.5;
    // define the multilevel solution and attach the mlMsh object to it
    MultiLevelSolution mlSol (&mlMsh); // Here we provide the mesh info to the problem.

    // add variables to mlSol
    mlSol.AddSolution ("u", LAGRANGE, SECOND, 2);
    mlSol.AddSolution ("d", LAGRANGE, SECOND,  0, false);

    mlSol.AddSolution ("K11", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
    mlSol.AddSolution ("K12", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
    mlSol.AddSolution ("K13", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
    mlSol.AddSolution ("K22", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
    mlSol.AddSolution ("K23", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
    mlSol.AddSolution ("K33", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
    //mlSol.AddSolution ("AD", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

    mlSol.Initialize ("All");
    mlSol.Initialize ("u", InitalValueU3D);
    mlSol.Initialize ("d", InitalValueD);
    GetKFromFileANISO (mlSol);

    // attach the boundary condition function and generate boundary data
    mlSol.AttachSetBoundaryConditionFunction (SetBoundaryCondition);
    mlSol.GenerateBdc ("u", "Steady");

    // define the multilevel problem attach the mlSol object to it
    MultiLevelProblem mlProb (&mlSol); //

    // add system Poisson in mlProb as a Non Linear Implicit System
    TransientNonlinearImplicitSystem & system = mlProb.add_system < TransientNonlinearImplicitSystem > ("Poisson");

    // add solution "u" to system
    system.AddSolutionToSystemPDE ("u");

    // attach the assembling function to system
    system.SetAssembleFunction (AssemblePoissonProblem_AD);

    // time loop parameter
    system.AttachGetTimeIntervalFunction (GetTimeStep);
    const unsigned int n_timesteps = 60;

    system.SetMaxNumberOfNonLinearIterations (1);
    system.SetMaxNumberOfLinearIterations (1);

    system.init();
    system.SetTolerances (1.e-10, 1.e-20, 1.e+50, 200, 40);

    system.SetSolverFineGrids (RICHARDSON);
    system.SetPreconditionerFineGrids (JACOBI_PRECOND);
    system.SetMgType (V_CYCLE);

    // ******* Print solution *******
    mlSol.SetWriter (VTK);
    mlSol.GetWriter()->SetDebugOutput (true);

    std::vector<std::string> print_vars;
    print_vars.push_back ("All");

    double time = system.GetTime();
    GetDeadCells (time, mlSol, false);

    mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);

    for (unsigned time_step = 0; time_step < n_timesteps; time_step++) {

      system.CopySolutionToOldSolution();

      system.MGsolve();

      bool last = (time_step == n_timesteps - 1) ? true : false;

      double time = system.GetTime();
      bool stop = GetDeadCells (time, mlSol, last);
      mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "biquadratic", print_vars, time_step + 1);
      if (stop) break;
    }
    // mlProb.clear();

  }

  return 0;
}



/**
 * This function assemble the stiffnes matrix KK and the residual vector Res
 * Using automatic differentiation for Newton iterative scheme
 *                  J(u0) w =  - F(u0)  ,
 *                  with u = u0 + w
 *                  - F = f(x) - J u = Res
 *                  J = \grad_u F
 *
 * thus
 *                  J w = f(x) - J u0
 **/
void AssemblePoissonProblem_AD (MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object


  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use

  TransientNonlinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<TransientNonlinearImplicitSystem> ("Poisson");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble(); // We have different level of meshes. we assemble the problem on the specified one.

  Mesh*                    msh = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + ! (dim - 1));       // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned > (ceil (pow (3, dim))); // Return a value of unsigned // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex ("u");   // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType (soluIndex);   // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex ("u");   // get the position of "u" in the pdeSys object

  vector < adept::adouble >  solu; // local solution
  solu.reserve (maxSize);

  vector < double >  soluOld; // local solution
  soluOld.reserve (maxSize);


  std::string kname[6] = {"K11", "K12", "K13", "K22", "K23", "K33"};

  unsigned kIndex[6];
  for (unsigned i = 0; i < 6; i++) {
    kIndex[i] = mlSol->GetIndex (kname[i].c_str());
  }
  unsigned kType = mlSol->GetSolutionType (kIndex[0]);
  double kT[3][3];

  vector < vector < double > > x (dim);   // local coordinates. x is now dim x m matrix.
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned k = 0; k < dim; k++) {
    x[k].reserve (maxSize); // dim x maxsize is reserved for x.
  }

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives

  double weight; // gauss point weight
  phi.reserve (maxSize);
  phi_x.reserve (maxSize * dim); // This is probably gradient but he is doing the life difficult for me!

  vector< adept::adouble > aRes; // local redidual vector
  aRes.reserve (maxSize);

  vector< int > l2GMap; // local to global mapping
  l2GMap.reserve (maxSize);
  vector< double > Res; // local redidual vector
  Res.reserve (maxSize);
  vector < double > Jac;
  Jac.reserve (maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual

  // element loop: each process loops only on the elements that owns
  // Adventure starts here!

  double dt = GetTimeStep (0.);

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofu  = msh->GetElementDofNumber (iel, soluType); // number of solution element dofs

    double K = (*sol->_Sol[kIndex[0]]) (iel);

    if (K <= Keps) {
      // local storage of global mapping and solution
      for (unsigned i = 0; i < nDofu; i++) {
        unsigned solDof = msh->GetSolutionDof (i, iel, soluType);   // global to global mapping between solution node and solution dof
        sol->_Sol[soluIndex]->set (solDof, 0.);
        sol->_SolOld[soluIndex]->set (solDof, 0.);
      }
    }
  }
  sol->_Sol[soluIndex]->close();


  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofu  = msh->GetElementDofNumber (iel, soluType); // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber (iel, xType);   // number of coordinate element dofs

    // resize local arrays
    l2GMap.resize (nDofu);
    solu.resize (nDofu);
    soluOld.resize (nDofu);

    for (int k = 0; k < dim; k++) {
      x[k].resize (nDofx); // Now we
    }

    aRes.resize (nDofu);   //resize
    std::fill (aRes.begin(), aRes.end(), 0);   //set aRes to zero

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, soluType);   // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex]) (solDof);                 // global extraction and local storage for the solution
      soluOld[i] = (*sol->_SolOld[soluIndex]) (solDof);           // global extraction and local storage for the solution
      l2GMap[i] = pdeSys->GetSystemDof (soluIndex, soluPdeIndex, i, iel);   // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);   // global to global mapping between coordinates node and coordinate dof

      for (unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k]) (xDof);     // global extraction and local storage for the element coordinates
      }
    }

    kT[0][0] = (*sol->_Sol[kIndex[0]]) (iel);
    kT[0][1] = (*sol->_Sol[kIndex[1]]) (iel);
    kT[0][2] = (*sol->_Sol[kIndex[2]]) (iel);
    kT[1][0] = kT[0][1];
    kT[1][1] = (*sol->_Sol[kIndex[3]]) (iel);
    kT[1][2] = (*sol->_Sol[kIndex[4]]) (iel);
    kT[2][0] = kT[0][2];
    kT[2][1] = kT[1][2];
    kT[2][2] = (*sol->_Sol[kIndex[5]]) (iel);


    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Face Gauss point loop (boundary Integral) ***
    for (unsigned jface = 0; jface < msh->GetElementFaceNumber (iel); jface++) {
      int faceIndex = el->GetBoundaryIndex (iel, jface);
      // look for boundary faces
      if (faceIndex == 1) {
        const unsigned faceGeom = msh->GetElementFaceType (iel, jface);
        unsigned faceDofs = msh->GetElementFaceDofNumber (iel, jface, soluType);

        vector  < vector  <  double> > faceCoordinates (dim);   // A matrix holding the face coordinates rowwise.
        for (int k = 0; k < dim; k++) {
          faceCoordinates[k].resize (faceDofs);
        }
        for (unsigned i = 0; i < faceDofs; i++) {
          unsigned inode = msh->GetLocalFaceVertexIndex (iel, jface, i);   // face-to-element local node mapping.
          for (unsigned k = 0; k < dim; k++) {
            faceCoordinates[k][i] =  x[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
          }
        }
        for (unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][soluType]->GetGaussPointNumber(); ig++) {
          // We call the method GetGaussPointNumber from the object finiteElement in the mesh object msh.
          vector < double> normal;
          msh->_finiteElement[faceGeom][soluType]->JacobianSur (faceCoordinates, ig, weight, phi, phi_x, normal);

          adept::adouble solu_gss = 0;
          double soluOld_gss = 0;



          for (unsigned i = 0; i < faceDofs; i++) {
            unsigned inode = msh->GetLocalFaceVertexIndex (iel, jface, i);   // face-to-element local node mapping.
            solu_gss += phi[i] * solu[inode];
            soluOld_gss += phi[i] * soluOld[inode];
          }

          // *** phi_i loop ***
          double eps = 0.0002;
          for (unsigned i = 0; i < faceDofs; i++) {
            unsigned inode = msh->GetLocalFaceVertexIndex (iel, jface, i);
            aRes[inode] +=  phi[i] * eps * (1.0 * solu_gss + 0. * soluOld_gss) * weight;
          }
        }
      }
    }

    // *** Element Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian (x, ig, weight, phi, phi_x);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      adept::adouble solu_gss = 0;
      vector < adept::adouble > gradSolu_gss (dim, 0.);

      double soluOld_gss = 0;
      vector < double > gradSoluOld_gss (dim, 0.);


      for (unsigned i = 0; i < nDofu; i++) {
        solu_gss += phi[i] * solu[i];
        soluOld_gss += phi[i] * soluOld[i];

        for (unsigned k = 0; k < dim; k++) {
          gradSolu_gss[k] += phi_x[i * dim + k] * solu[i];
          gradSoluOld_gss[k] += phi_x[i * dim + k] * soluOld[i];
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofu; i++) {

        adept::adouble graduGradphi = 0.;
        double graduOldGradphi = 0.;

        for (unsigned k = 0; k < dim; k++) {
          for (unsigned l = 0; l < dim; l++) {
            graduGradphi     +=   kT[k][l] * gradSolu_gss[l] * phi_x[i * dim + k] ;
            //graduOldGradphi  +=   kT[k][l] * gradSoluOld_gss[l] * phi_x[i * dim + k];
          }
        }

        aRes[i] += ( (solu_gss - soluOld_gss) * phi[i] / dt + (1. * graduGradphi + 0. * graduOldGradphi)) * weight;

      } // end phi_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    Res.resize (nDofu);   //resize

    for (int i = 0; i < nDofu; i++) {
      Res[i] = - aRes[i].value();
    }

    RES->add_vector_blocked (Res, l2GMap);

    // define the dependent variables
    s.dependent (&aRes[0], nDofu);

    // define the independent variables
    s.independent (&solu[0], nDofu);

    // get the jacobian matrix (ordered by row major )
    Jac.resize (nDofu * nDofu);   //resize
    s.jacobian (&Jac[0], true);

    //store K in the global matrix KK
    KK->add_matrix_blocked (Jac, l2GMap, l2GMap);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();

  KK->close();

  // ***************** END ASSEMBLY *******************
}


bool GetDeadCells (const double &time, MultiLevelSolution &mlSol, const bool & last) {

  double a = 0.0471;
  double b = -1.4629;
  double c = 0.0866;

  std::pair<double, double> uT[3];

  for (unsigned i = 0; i < 3; i++) {
    uT[i].first = (i + 1) * 24;
    uT[i].second =  a - b * exp (-c * uT[i].first); //0.1 * exp( - uT[i].first / 24);
  }

  double treshold = a - b * exp (-c * time);
  std::cout << "time = " << time << " treshold = " << treshold << std::endl;

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol  = mlSol.GetSolutionLevel (level);
  Mesh     *msh   = mlSol._mlMesh->GetLevel (level);
  unsigned iproc  = msh->processor_id();

  unsigned soluIndex = mlSol.GetIndex ("u");
  unsigned soldIndex = mlSol.GetIndex ("d");
  unsigned solType = mlSol.GetSolutionType (soluIndex);

  for (int inode = msh->_dofOffset[solType][iproc]; inode < msh->_dofOffset[solType][iproc + 1]; inode++) {
    double d = (*sol->_Sol[soldIndex]) (inode);
    double u = (*sol->_Sol[soluIndex]) (inode);

    for (unsigned i = 0; i < 3; i++) {
      if (d > uT[i].first && u > uT[i].second) {
        sol->_Sol[soldIndex]->set (inode, uT[i].first);
        break;
      }
    }
  }
  sol->_Sol[soldIndex]->close();

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  std::vector<double> sold;
  std::vector < std::vector<double> >  x (dim);

  unsigned soldType = mlSol.GetSolutionType (soldIndex);   // get the finite element type for "u"
  unsigned xType = 2;

  double weight;
  std::vector<double> phi;
  std::vector<double> phi_x;

  double volume = 0;
  double volumeUT[3] = { 0., 0., 0.};

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofd  = msh->GetElementDofNumber (iel, soldType); // number of solution element dofs


    sold.resize (nDofd);

    for (int k = 0; k < dim; k++) {
      x[k].resize (nDofd); // Now we
    }


    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofd; i++) {
      unsigned solDof = msh->GetSolutionDof (i, iel, soldType);   // global to global mapping between solution node and solution dof
      sold[i] = (*sol->_Sol[soldIndex]) (solDof);                 // global extraction and local storage for the solution
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofd; i++) {
      unsigned xDof  = msh->GetSolutionDof (i, iel, xType);   // global to global mapping between coordinates node and coordinate dof

      for (unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->_topology->_Sol[k]) (xDof);     // global extraction and local storage for the element coordinates
      }
    }

    double K11 = (*sol->_Sol[ mlSol.GetIndex ("K11")]) (iel);
    double K22 = (*sol->_Sol[ mlSol.GetIndex ("K22")]) (iel);
    double K33 = (*sol->_Sol[ mlSol.GetIndex ("K33")]) (iel);

    if (x[0][nDofd - 1] > -0.6 && x[0][nDofd - 1] < 1.7 &&
        x[1][nDofd - 1] > -1.5 && x[1][nDofd - 1] < 1.5 &&
        x[2][nDofd - 1] > -1.1 && x[0][nDofd - 1] < 2.2 && (K11+K22+K33)/3. > 0.6) {
      // *** Element Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soldType]->GetGaussPointNumber(); ig++) {
        // *** get gauss point weight, test function and test function partial derivatives ***
        msh->_finiteElement[ielGeom][soldType]->Jacobian (x, ig, weight, phi, phi_x);

        // evaluate the solution, the solution derivatives and the coordinates in the gauss point
        double sold_gss = 0;

        for (unsigned i = 0; i < nDofd; i++) {
          sold_gss += phi[i] * sold[i];
        }

        volume += weight;

        if (sold_gss <= 24) volumeUT[0] += weight;
        if (sold_gss <= 48) volumeUT[1] += weight;
        if (sold_gss <= 72) volumeUT[2] += weight;

      } // end gauss point loop
    }

  }

  double volumeAll = 0.;
  MPI_Reduce (&volume, &volumeAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  double volumeUTAll[3] = { 0., 0., 0.};
  MPI_Reduce (volumeUT, volumeUTAll, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  std::cout << "VOLUME FRACTIONS = \n";
  std::cout << volumeUTAll[0] / volumeAll << " " << volumeUTAll[1] / volumeAll << " " << volumeUTAll[2] / volumeAll << " " << std::endl;


  double lInfinityNorm    = sol->_Sol[soluIndex]->linfty_norm();

  std::cout << "Max = " << lInfinityNorm << " Treshold = " << uT[2].second << std::endl;
  


  bool stop = (lInfinityNorm < uT[2].second) ? true : false;

  if ( (stop || last) && iproc == 0) {
    std::ofstream fout;
    unsigned j = jInitialPosition;
    
    fout.open ("DoseResponseCurve.csv", std::ofstream::app);
    //fout << xc[j] <<","<<  yc[j] <<","<< zc[j] << ","<<  V0 << "," << volumeUTAll[0] / volumeAll << "," << volumeUTAll[1] / volumeAll << "," << volumeUTAll[2] / volumeAll << "," << std::endl;
    fout <<  V0 << "," << volumeUTAll[0] / volumeAll << "," << volumeUTAll[1] / volumeAll << "," << volumeUTAll[2] / volumeAll << "," << std::endl;

    fout.close();

    fout.open ("DoseResponseCurve.txt", std::ofstream::app);
    //fout << xc[j] <<","<<  yc[j] <<","<< zc[j] << ","<<  V0 << "," << volumeUTAll[0] / volumeAll << "," << volumeUTAll[1] / volumeAll << "," << volumeUTAll[2] / volumeAll << "," << std::endl;
    fout << V0 << "," << volumeUTAll[0] / volumeAll << "," << volumeUTAll[1] / volumeAll << "," << volumeUTAll[2] / volumeAll << "," << std::endl;
    fout.close();
  }
  return stop;
}


bool CheckIfPositiveDefinite (double K[6]) {
  bool pDefinite = true;
  double det1 = K[0];
  double det2 = K[0] * K[3] - K[1] * K[1];
  double det3 =   K[0] * (K[3] * K[5] - K[4] * K[4])
                  - K[1] * (K[1] * K[5] - K[4] * K[2])
                  + K[2] * (K[1] * K[4] - K[3] * K[2]);
  if (det1 <= 0. || det2 <= 0. || det3 <= 0.) {
    pDefinite = false;
  }
  return pDefinite;
};


void GetKFromFileISO (MultiLevelSolution &mlSol, const unsigned & split, const unsigned &patient) {

  unsigned Level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol  = mlSol.GetSolutionLevel (Level);
  Mesh     *msh   = mlSol._mlMesh->GetLevel (Level);

  std::ostringstream filename;
  filename << "./input/CroppedMD_Data0" << patient << ".txt";

  std::ifstream fin;

  //fin.open ("./input/MeanDiffData.txt");
  fin.open (filename.str().c_str());
  if (!fin.is_open()) {
    std::cout << std::endl << " The output file " << "./input/MeanDiffData.txt" << " cannot be opened.\n";
    abort();
  }

  unsigned n1, n2, n3;
  fin >> n1;
  fin >> n2;
  fin >> n3;

  double h1, h2, h3;

  h1 = 2. / (n1 / split);
  h2 = 2. / (n2 / split);
  h3 = 2. / (n3 / split);

  std::string kname[6] = {"K11", "K12", "K13", "K22", "K23", "K33"};

  unsigned kIndex[6];
  for (unsigned i = 0; i < 6; i++) {
    kIndex[i] = mlSol.GetIndex (kname[i].c_str());
  }

  unsigned kType = mlSol.GetSolutionType (kIndex[0]);

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  std::vector<double> x (dim);

  unsigned iproc  = msh->processor_id();
  unsigned nprocs  = msh->n_processors();



  for (unsigned i = 0; i < n1; i++) {
    for (unsigned j = 0; j < n2; j++) {
      for (unsigned k = 0; k < n3; k++) {

        if (i < n1 / split && j < n2 / split & k < n3 / split) {
          x[0] = -1. + h1 * i + 0.5 * h1;
          x[1] = -1. + h2 * j + 0.5 * h2;
          x[2] = -1. + h3 * k + 0.5 * h3;

          double K[6];
          //for (unsigned l = 0; l < 1; l++) {
          fin >> K[0];
          if (K[0] < Keps) K[0] = Keps;
          K[1] = K[2] = K[4] = 0.;
          K[3] = K[5] = K[0];
          //}

          Marker center = Marker (x, 1., VOLUME , sol, 0);
          unsigned mproc = center.GetMarkerProc (sol);



          if (mproc == iproc) {
            unsigned iel = center.GetMarkerElement();
//             while (!CheckIfPositiveDefinite (K)) {
//               std::cout << " warning k[" << iel << "] is not positive definite\n";
//
//               K[0] += (K[0] < 0.) ? fabs (K[0]) + fabs (K[1]) + fabs (K[2]) : fabs (K[1]) + fabs (K[2]);
//               K[3] += (K[3] < 0.) ? fabs (K[1]) + fabs (K[3]) + fabs (K[4]) : fabs (K[1]) + fabs (K[4]);
//               K[5] += (K[5] < 0.) ? fabs (K[2]) + fabs (K[4]) + fabs (K[5]) : fabs (K[2]) + fabs (K[4]);
//               if (K[0] < eps) K[0] = eps;
//               if (K[3] < eps) K[3] = eps;
//               if (K[5] < eps) K[5] = eps;
//
//             }
            for (unsigned l = 0; l < 6; l++) {
              sol->_Sol[kIndex[l]]->set (iel, K[l]);
            }
          }
        }
      }
    }
  }

  fin.close();
  for (unsigned l = 0; l < 6; l++) {
    sol->_Sol[kIndex[l]]->close();
  }

//   double trace = 0.;
//   for (int iel = msh->_dofOffset[kType][iproc]; iel < msh->_dofOffset[kType][iproc + 1]; iel++) {
//     trace += ( (*sol->_Sol[kIndex[0]]) (iel) + (*sol->_Sol[kIndex[3]]) (iel)
//                + (*sol->_Sol[kIndex[5]]) (iel)) / 3.;
//   }
//   double traceAll = 0.;
//   MPI_Allreduce (&trace, &traceAll, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//
//   traceAll /= msh->_dofOffset[kType][nprocs];
//
//
//   for (int iel = msh->_dofOffset[kType][iproc]; iel < msh->_dofOffset[kType][iproc + 1]; iel++) {
//
//     for (unsigned j = 0; j < 6; j++) {
//       double value = (*sol->_Sol[kIndex[j]]) (iel) / traceAll;
//       sol->_Sol[kIndex[j]]->set (iel, value);
//     }
//   }
//   for (unsigned l = 0; l < 6; l++) {
//     sol->_Sol[kIndex[l]]->close();
//   }



}


void GetKFromFileANISO (MultiLevelSolution &mlSol) {

  unsigned Level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution *sol  = mlSol.GetSolutionLevel (Level);
  Mesh     *msh   = mlSol._mlMesh->GetLevel (Level);

  std::ostringstream filename;
  std::ostringstream fileAD;
  
  double treshold = 0.;
  filename << "./input/NewCorrectedTensorSPD1.txt";
  treshold = 0.002;
  
//   filename << "./input/NewCorrectedTensorSPD2.txt";
//   treshold = 0.0002;
 
//   filename << "./input/NewCorrectedTensorSPD3.txt";
//   treshold = 0.00002;
  
  
  //filename << "./input/NewCorrectedTensorSPD_0_1.txt";
  //fileAD << "/home/erdi/FEMuS/MyFEMuS/applications/Tumor/ex3/input/AxialDiffusivity.txt";

  std::ifstream fin;

  fin.open (filename.str().c_str());
  
  if (!fin.is_open()) {
    std::cout << " The input file " << filename.str().c_str() << " cannot be opened.\n" << std::endl;
    abort();
  }

  unsigned n1, n2, n3;
  fin >> n3;
  fin >> n2;
  fin >> n1;

  double h1, h2, h3;

  h1 = 5. / n1;
  h2 = 5. / n2;
  h3 = 5. / n3;

  std::string kname[6] = {"K11", "K12", "K13", "K22", "K23", "K33"};
  //std::string ADName = "AD";

  unsigned kIndex[6];
  for (unsigned i = 0; i < 6; i++) {
    kIndex[i] = mlSol.GetIndex (kname[i].c_str());
  }
  //unsigned ADIndex;
  //ADIndex = mlSol.GetIndex(ADName.c_str());
  

  unsigned kType = mlSol.GetSolutionType (kIndex[0]);
  //unsigned ADType = mlSol.GetSolutionType (ADIndex);

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  std::vector<double> x (dim);

  unsigned iproc  = msh->processor_id();
  unsigned nprocs  = msh->n_processors();


  for (unsigned k = 0; k < n3; k++) {
    for (unsigned j = 0; j < n2; j++) {
      for (unsigned i = 0; i < n1; i++) {

        double K[6];
        for (unsigned l = 0; l < 6; l++) {
          fin >> K[l];
        }
        //double AD;
        //fAD >> AD;
        if (K[0] < Keps) K[0] = Keps;
        if (K[3] < Keps) K[3] = Keps;
        if (K[5] < Keps) K[5] = Keps;

        x[0] = -2.5 + h1 * i + 0.5 * h1;
        x[1] = -2.5 + h2 * j + 0.5 * h2;
        x[2] = -2.5 + h3 * k + 0.5 * h3;

        Marker center = Marker (x, 1., VOLUME , sol, 0);
        unsigned mproc = center.GetMarkerProc (sol);

        if (mproc == iproc) {
          unsigned iel = center.GetMarkerElement();
          while (!CheckIfPositiveDefinite (K)) {
            std::cout << " warning k[" << iel << "] is not positive definite\n";

            K[0] += (K[0] < 0.) ? fabs (K[0]) + fabs (K[1]) + fabs (K[2]) : fabs (K[1]) + fabs (K[2]);
            K[3] += (K[3] < 0.) ? fabs (K[1]) + fabs (K[3]) + fabs (K[4]) : fabs (K[1]) + fabs (K[4]);
            K[5] += (K[5] < 0.) ? fabs (K[2]) + fabs (K[4]) + fabs (K[5]) : fabs (K[2]) + fabs (K[4]);
            if (K[0] < Keps) K[0] = Keps;
            if (K[3] < Keps) K[3] = Keps;
            if (K[5] < Keps) K[5] = Keps;
          }
          //sol->_Sol[ADIndex]->set (iel, AD);
          for (unsigned l = 0; l < 6; l++) {
            sol->_Sol[kIndex[l]]->set (iel, K[l]);
          }
        }
      }
    }

  }

  fin.close();
  for (unsigned l = 0; l < 6; l++) {
    sol->_Sol[kIndex[l]]->close();
  }
  //fAD.close();
  //sol->_Sol[ADIndex]->close();

  //rescale so that the average of all inSkull DTI traces is 1

  double trace = 0.;
  unsigned counter = 0;
  for (int iel = msh->_dofOffset[kType][iproc]; iel < msh->_dofOffset[kType][iproc + 1]; iel++) {

    unsigned xType = 2;
    unsigned n  = msh->GetElementDofNumber (iel, xType); // number of solution element dofs

    unsigned xDof  = msh->GetSolutionDof (n - 1, iel, xType);   // local to global mapping between coordinates node and coordinate dof

    std::vector<double> x (3);
    for (unsigned k = 0; k < dim; k++) {
      x[k] = (*msh->_topology->_Sol[k]) (xDof);     // global extraction and local storage for the element coordinates
    }



    if (x[0] > -0.6 && x[0] < 1.7 &&
        x[1] > -1.5 && x[1] < 1.5 &&
        x[2] > -1.1 && x[0] < 2.2)  {

      //if (r < 1.3) {
      double traceIel = ( (*sol->_Sol[kIndex[0]]) (iel) + (*sol->_Sol[kIndex[3]]) (iel)
                          + (*sol->_Sol[kIndex[5]]) (iel)) / 3.;
      if (traceIel > treshold) {
        trace += traceIel;
        counter++;
      }
    }
  }
  std::cout <<"Trace is: " << trace <<"counter is :"<< counter <<std::endl;

  double traceAll = 0.;
  unsigned counterAll;
  MPI_Allreduce (&trace, &traceAll, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&counter, &counterAll, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  
  std::cout << "TraceAll is : " << traceAll << "--" << " CounterAll = "<< counterAll <<std::endl;
  
  traceAll *= 1. / counterAll;
  

  std::cout << "TraceAll after is : " << traceAll << "--" << " CounterAll after = "<< counterAll <<std::endl;


  for (int iel = msh->_dofOffset[kType][iproc]; iel < msh->_dofOffset[kType][iproc + 1]; iel++) {
    for (unsigned j = 0; j < 6; j++) {
      double value = (*sol->_Sol[kIndex[j]]) (iel) / traceAll;
      sol->_Sol[kIndex[j]]->set (iel, value);
    }
  }
  for (unsigned l = 0; l < 6; l++) {
    sol->_Sol[kIndex[l]]->close();
  }

}
