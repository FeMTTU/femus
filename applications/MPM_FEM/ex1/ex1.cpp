// tests the inclusion algorithm for all 2D and 3D elements

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "XDMFWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "LinearImplicitSystem.hpp"
#include "adept.h"
#include "Marker.hpp"
#include "Line.hpp"
#include "TransientSystem.hpp"
#include "OprtrTypeEnum.hpp"

using namespace femus;


bool SetRefinementFlag(const std::vector < double >& x, const int& elemgroupnumber, const int& level) {

  bool refine = 0;

  if(elemgroupnumber == 6 && level < 4) refine = 1;
  if(elemgroupnumber == 7 && level < 5) refine = 1;
  if(elemgroupnumber == 8 && level < 6) refine = 1;

  return refine;

}


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  unsigned numberOfUniformLevels = 3; //for refinement in 3D
  //unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;


  /* element types
  0 = HEX
  1 = TET
  2 = WEDGE
  3 = QUAD
  4 = TRI
   */

  unsigned solType = 2;

  std::cout << " --------------------------------------------------     TEST     --------------------------------------------------" << std::endl;

  mlMsh.ReadCoarseMesh("./input/adaptiveRef6.neu", "seventh", scalingFactor);
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);

  unsigned dim = mlMsh.GetDimension();

//   PrintLine(DEFAULT_OUTPUTDIR, linea.GetLine(), false, 0);

  MultiLevelSolution mlSol(&mlMsh);
  // add variables to mlSol
  mlSol.AddSolution("DX", LAGRANGE, SECOND, 2);
  if(dim > 2) mlSol.AddSolution("DY", LAGRANGE, SECOND, 2);
  if(dim > 3) mlSol.AddSolution("DZ", LAGRANGE, SECOND, 2);

  mlSol.AddSolution("U", LAGRANGE, SECOND, 2);
  if(dim > 2) mlSol.AddSolution("V", LAGRANGE, SECOND, 2);
  if(dim > 3) mlSol.AddSolution("W", LAGRANGE, SECOND, 2);

  mlSol.AddSolution("AU", LAGRANGE, SECOND, 2);
  if(dim > 2) mlSol.AddSolution("AV", LAGRANGE, SECOND, 2);
  if(dim > 3) mlSol.AddSolution("AW", LAGRANGE, SECOND, 2);

  mlSol.Initialize("All");

  std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;


  unsigned rows = 20;
  unsigned columns = 50;
  unsigned size = rows * columns;

  std::vector < std::vector < double > > x; // marker
  std::vector < MarkerType > markerType;

  x.resize(size);
  markerType.resize(size);

  std::vector < std::vector < std::vector < double > > > line(1);
  std::vector < std::vector < std::vector < double > > > line0(1);

  for(unsigned j = 0; j < size; j++) {
    x[j].assign(dim, 0.);
    markerType[j] = VOLUME;
  }

  // double pi = acos(-1.);

  //BEGIN initialization
  for(unsigned i = 0; i < rows; i++) {
    for(unsigned j = 0; j < columns; j++) {


      x[i * columns + j][0] = -0.5 + (0.75 / (columns - 1)) * j;
      x[i * columns + j][1] = -0.1 + (0.2 / (rows - 1)) * i;
      if(dim == 3) {
        x[j][2] = 0.;
      }
    }
  }

  //END


  Line linea(x, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), solType);

  //linea.GetPointsToGridProjections();

  linea.GetLine(line0[0]);
  PrintLine(DEFAULT_OUTPUTDIR, line0, false, 0);

  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);


//BEGIN stuff for the advection that we will have to remove
//   double T = 2 * acos(-1.);
//
//   unsigned n = 4;

//   std::cout << std::endl << " init in  " << std::setw(11) << std::setprecision(6) << std::fixed
//             << static_cast<double>((clock() - init_time)) / CLOCKS_PER_SEC << " s" << std::endl;
//
//
//   clock_t advection_time = clock();

//   for (unsigned k = 1; k <= n; k++) {
//     std::cout << "Iteration = " << k << std::endl;
//     //uncomment for  vortex test
//     mlSol.CopySolutionToOldSolution();
//     mlSol.UpdateSolution("U" , InitalValueU, pi * k / n);
//     mlSol.UpdateSolution("V" , InitalValueV, pi * k / n);
//     if (dim == 3) mlSol.UpdateSolution("W" , InitalValueW, pi * k / n);
//     linea.AdvectionParallel(40, T / n, 4);
//     linea.GetLine(line[0]);
//     PrintLine(DEFAULT_OUTPUTDIR, line, false, k);
//   }


//   std::cout << std::endl << " advection in: " << std::setw(11) << std::setprecision(6) << std::fixed
//             << static_cast<double>((clock() - advection_time)) / CLOCKS_PER_SEC << " s" << std::endl;
//
//   std::cout << std::endl << " RANNA in: " << std::setw(11) << std::setprecision(6) << std::fixed
//             << static_cast<double>((clock() - start_time)) / CLOCKS_PER_SEC << " s" << std::endl;

  //computing the geometric error
//   double error = 0.;
//   for (unsigned j = 0; j < size + 1; j++) {
//     double tempError = 0.;
//     for (unsigned i = 0; i < dim; i++) {
//       tempError += (line0[0][j][i] - line[0][j][i]) * (line0[0][j][i] - line[0][j][i]);
//     }
//     error += sqrt(tempError);
//   }
  /*
    error = error / size;*/
  /*
    std::cout << " ERROR = " << std::setprecision(15) << error << std::endl;*/

//   for(unsigned j = 0; j < size; j++) {
//     std::vector <double> trial(dim);
//     trial = linea._particles[linea._printList[j]]->GetIprocMarkerCoordinates();
//     for(unsigned i=0; i<dim; i++){
//       std::cout << " x[" << j << "][" << i << "]=" << trial[i] << std::endl;
//     }
//   }
//END stuff for the advection that we will have to remove


  return 0;

} //end main



void AssembleMPMSys(MultiLevelProblem & ml_prob, Line &linea, const unsigned &numberOfUniformLevels) {

  // ml_prob is the global object from/to where get/set all the data
  // level is the level of the PDE system to be assembled
  // levelMax is the Maximum level of the MultiLevelProblem
  // assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

  //pointers and references

  //TODO fix the one below for our case
  TransientNonlinearImplicitSystem & my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("Fluid-Structure-Interaction");
  const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution * ml_sol = ml_prob._ml_sol; // pointer to the multilevel solution object
  Solution * mysolution = ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object
  LinearEquationSolver * myLinEqSolver = my_nnlin_impl_sys._LinSolver[level]; // pointer to the equation (level) object

  Mesh * mymsh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem * myel = mymsh->el;  // pointer to the elem object in msh (level)
  SparseMatrix * myKK = myLinEqSolver->_KK; // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector * myRES =  myLinEqSolver->_RES; // pointer to the global residual vector object in pdeSys (level)
  bool assembleMatrix = my_nnlin_impl_sys.GetAssembleMatrix();

// call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;
  if(assembleMatrix) s.continue_recording();
  else s.pause_recording();

  const unsigned dim = mymsh->GetDimension();
  const unsigned nabla_dim = 3 * (dim - 1);   //TODO don't know if we need it

  // reserve memory for the local standar vectors
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  // data
  unsigned nel    = mymsh->GetNumberOfElements();
  unsigned iproc  = mymsh->processor_id();
  double dt =  my_nnlin_impl_sys.GetIntervalTime();
  double time =  my_nnlin_impl_sys.GetTime();

  // local objects
  vector<adept::adouble> SolVAR(3 * dim);
  vector<double> SolVAR_old(3 * dim);

  vector<vector < adept::adouble > > GradSolVAR(3 * dim);
  vector<vector < double > > GradSolVAR_old(3 * dim);

  vector<vector < adept::adouble > > GradSolhatVAR(3 * dim);    //TODO what is this?
  vector<vector < double > > GradSolhatVAR_old(3 * dim);

  vector<vector<adept::adouble> > NablaSolVAR(3 * dim);    //TODO what is this?
  vector<vector < double > > NablaSolVAR_old(3 * dim);

  for(int i = 0; i < 3 * dim; i++) {
    GradSolVAR[i].resize(dim);
    GradSolVAR_old[i].resize(dim);

    GradSolhatVAR[i].resize(dim);
    GradSolhatVAR_old[i].resize(dim);

    NablaSolVAR[i].resize(nabla_dim);
    NablaSolVAR_old[i].resize(nabla_dim);
  }

  //BEGIN do we need these since our phi's are evaluated at the particle local coordinates?
  vector < double > phi;
  vector < double > phi_hat;
  vector < double > phi_old; //TODO I don't think we need this
  vector < adept::adouble> gradphi;
  vector < double > gradphi_hat;
  vector < double > gradphi_old; //TODO I don't think we need this
  vector < adept::adouble> nablaphi;
  vector < double > nablaphi_hat;
  vector < double > nablaphi_old; //TODO I don't think we need this

  phi.reserve(max_size);
  phi_hat.reserve(max_size);
  phi_old.reserve(max_size);

  gradphi.reserve(max_size * dim);
  gradphi_hat.reserve(max_size * dim);
  gradphi_old.reserve(max_size * dim);

  nablaphi.reserve(max_size * 3 * (dim - 1));
  nablaphi_hat.reserve(max_size * 3 * (dim - 1));
  nablaphi_old.reserve(max_size * 3 * (dim - 1));
  //END

  vector <vector < adept::adouble> > vx(dim); //vx is coordX in assembly of ex30
  vector <vector < double> > vx_hat(dim);
  vector <vector < double> > vx_old(dim);

  //TODO what is this?
  vector <vector < adept::adouble > > vx_face(dim);
  vector <vector < double > > vx_face_old(dim);

  for(int i = 0; i < dim; i++) {
    vx[i].reserve(max_size);
    vx_hat[i].reserve(max_size);
    vx_old[i].reserve(max_size);

    vx_face[i].resize(9);
    vx_face_old[i].resize(9);
  }

  vector< vector< adept::adouble > > Soli(3 * dim);      // local solution
  vector< vector< double > > Soli_old(3 * dim);      // local solution old
  vector< vector< int > > dofsVAR(3 * dim);

  for(int i = 0; i < 3 * dim ; i++) {
    Soli[i].reserve(max_size);
    Soli_old[i].reserve(max_size);
    dofsVAR[i].reserve(max_size);
  }

  vector< vector< double > > Rhs(3 * dim);     // local redidual vector
  vector< vector< adept::adouble > > aRhs(3 * dim);     // local redidual vector

  for(int i = 0; i < 3 * dim; i++) {
    aRhs[i].reserve(max_size);
    Rhs[i].reserve(max_size);
  }

  vector < int > dofsAll;
  dofsAll.reserve(max_size * (3 * dim));

  vector < double > Jac;
  Jac.reserve(dim * max_size * (3 * dim) *dim * max_size * (3 * dim));


  // gravity
  double _gravity[3] = {0., 0., 0.};


  // space discretization parameters
  unsigned SolType2 = ml_sol->GetSolutionType(ml_sol->GetIndex("DX"));


  //variable-name handling
  const char varname[9][3] = {"DX", "DY", "DZ", "U", "V", "W", "AU", "AV", "AW"};
  vector <unsigned> indexVAR(3 * dim);
  vector <unsigned> indVAR(3 * dim);
  vector <unsigned> SolType(3 * dim);

  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indVAR[ivar] = ml_sol->GetIndex(&varname[ivar][0]);
    indVAR[ivar + dim] = ml_sol->GetIndex(&varname[ivar + 3][0]);
    indVAR[ivar + 2 * dim] = ml_sol->GetIndex(&varname[ivar + 6][0]);
    SolType[ivar] = ml_sol->GetSolutionType(&varname[ivar][0]);
    SolType[ivar + dim] = ml_sol->GetSolutionType(&varname[ivar + 3][0]);
    SolType[ivar + 2 * dim] = ml_sol->GetSolutionType(&varname[ivar + 6][0]);
    indexVAR[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);
    indexVAR[ivar + dim] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar + 3][0]);
    indexVAR[ivar + 2 * dim] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar + 6][0]);
  }

  start_time = clock();

  if(assembleMatrix) myKK->zero();

  //line instances
  std::vector<unsigned> markerOffset = linea.GetMarkerOffset();
  unsigned markerOffset1 = markerOffset[iproc];
  unsigned markerOffset2 = markerOffset[iproc + 1];
  std::vector<Marker*> particles = linea.GetParticles();
  std::map<unsigned, std::vector < std::vector < std::vector < std::vector < double > > > > > aX;

  //initialization of iel
  unsigned ielOld = particles[markerOffset2]->GetMarkerElement();
  //declaration of element instances
  short unsigned ielt;
  unsigned nve;

//BEGIN loop on particles (used as Gauss points)
  for(unsigned iMarker = markerOffset1; iMarker < markerOffset2; iMarker++) {

    //element of particle iMarker
    unsigned iel = particles[iMarker]->GetMarkerElement();

    //BEGIN buliding Cauchy stress (using Saint Venant)
    //physical quantity
    adept::adouble J_hat;
    double J_hat_old;
    adept::adouble I_e;
    double I_e_old;
    adept::adouble Cauchy[3][3];
    adept::adouble Cauchy_old[3][3]; //TODO do we need this
    double Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};
    double mus = 1.; //TODO a real initialization

    adept::adouble e[3][3];
    double e_old[3][3];

    //computation of the stress tensor
    for(int i = 0; i < dim; i++) {
      for(int j = 0; j < dim; j++) {
        e[i][j]    = 0.5 * (GradSolhatVAR[i][j]    + GradSolhatVAR[j][i]);
        e_old[i][j] = 0.5 * (GradSolhatVAR_old[i][j] + GradSolhatVAR_old[j][i]);
      }
    }

    I_e = 0;
    I_e_old = 0;

    for(int i = 0; i < dim; i++) {
      I_e     += e[i][i];
      I_e_old += e_old[i][i];
    }

    for(int i = 0; i < dim; i++) {
      for(int j = 0; j < dim; j++) {
        //incompressible
        Cauchy[i][j]     = 2 * mus * e[i][j] - 2 * mus * I_e * SolVAR[2 * dim] * Id2th[i][j];
        Cauchy_old[i][j] = 2 * mus * e_old[i][j] - 2 * mus * I_e_old * SolVAR[2 * dim] * Id2th[i][j];
        //+(penalty)*lambda*I_e*Id2th[i][j];
      }
    }
//END building Cauchy stress

    //update element related quantities only if we are in a different element
    if(iel != ielOld) {

      ielt = mymsh->GetElementType(iel);
      nve = mymsh->GetElementDofNumber(iel, SolType2);

      //initialization of everything is in common fluid and solid

      //Rhs
      for(int i = 0; i < 3 * dim; i++) {
        dofsVAR[i].resize(nve);
        Soli[indexVAR[i]].resize(nve);
        Soli_old[indexVAR[i]].resize(nve);

        aRhs[indexVAR[i]].resize(nve);
        //  Rhs[indexVAR[i]].resize ( nve );
      }

      dofsAll.resize(0);

      Jac.resize((3 * dim * nve) * (3 * dim * nve));

      // ----------------------------------------------------------------------------------------
      // coordinates, solutions, displacement, velocity dofs

      for(int i = 0; i < dim; i++) {
        vx[i].resize(nve);
        vx_hat[i].resize(nve);
        vx_old[i].resize(nve);
      }

      //BEGIN copy of the value of Sol at idof in the local vector
      for(unsigned i = 0; i < nve; i++) {
        unsigned idof = mymsh->GetSolutionDof(i, iel, SolType2);
        for(int j = 0; j < dim; j++) {
          Soli[indexVAR[j]][i]     = (*mysolution->_Sol[indVAR[j]])(idof);
          Soli[indexVAR[j + dim]][i] = (*mysolution->_Sol[indVAR[j + dim]])(idof);
          Soli[indexVAR[j + 2 * dim]][i] = (*mysolution->_Sol[indVAR[j + 2 * dim]])(idof);

          Soli_old[indexVAR[j]][i]     = (*mysolution->_SolOld[indVAR[j]])(idof);
          Soli_old[indexVAR[j + dim]][i] = (*mysolution->_SolOld[indVAR[j + dim]])(idof);
          Soli_old[indexVAR[j + 2 * dim]][i] = (*mysolution->_SolOld[indVAR[j + 2 * dim]])(idof);

          aRhs[indexVAR[j]][i]     = 0.;
          aRhs[indexVAR[j + dim]][i] = 0.; //TODO maybe we don't need this
          aRhs[indexVAR[j + 2 * dim]][i] = 0.;  //TODO maybe we don't need this

          //Fixed coordinates (Reference frame) //TODO do we need this?
          vx_hat[j][i] = (*mymsh->_topology->_Sol[j])(idof);
          // displacement dofs
          dofsVAR[j][i] = myLinEqSolver->GetSystemDof(indVAR[j], indexVAR[j], i, iel);
          // velocity dofs
          dofsVAR[j + dim][i] = myLinEqSolver->GetSystemDof(indVAR[j + dim], indexVAR[j + dim], i, iel);
          // acceleration dofs
          dofsVAR[j + 2 * dim][i] = myLinEqSolver->GetSystemDof(indVAR[j + 2 * dim], indexVAR[j + 2 * dim], i, iel);
        }
      }
      //END

      // build dof composition
      for(int idim = 0; idim < 3 * dim; idim++) {
        dofsAll.insert(dofsAll.end(), dofsVAR[idim].begin(), dofsVAR[idim].end());
      }

      dofsAll.insert(dofsAll.end(), dofsVAR[2 * dim].begin(), dofsVAR[2 * dim].end());

      if(assembleMatrix) s.new_recording();

      //TODO do we need this?
      for(unsigned idim = 0; idim < dim; idim++) {
        for(int j = 0; j < nve; j++) {
          vx[idim][j]    = vx_hat[idim][j] + Soli[indexVAR[idim]][j];
          vx_old[idim][j] = vx_hat[idim][j] + Soli_old[indexVAR[idim]][j];
        }
      }
    }

    // start a new recording of all the operations involving adept::adouble variables //TODO not sure if this should be here
    if(assembleMatrix) s.new_recording();

    bool elementUpdate = (aX.find(iel) != aX.end()) ? false : true;     //update if iel was never updated
    particles[iMarker]->FindLocalCoordinates(SolType2, aX[iel], elementUpdate, ml_sol->GetLevel(numberOfUniformLevels - 1), 0);


    // the local coordinates of the particles are the Gauss points in this context
    std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();
    // the mass of the particles acts as weight
    double mass = particles[iMarker]->GetMarkerMass();
    adept::adouble massAdept = mass;
    mymsh->_finiteElement[ielt][SolType2]->Jacobian(vx, iMarker, massAdept, phi, gradphi, nablaphi);
    mymsh->_finiteElement[ielt][SolType2]->Jacobian(vx_hat, iMarker, mass, phi_hat, gradphi_hat, nablaphi_hat);

    // displacement and velocity
    //BEGIN evaluates SolVar at idof of element iel
    //TODO maybe we don't need this for the displacement
    for(int i = 0; i < 3 * dim; i++) {
      SolVAR[i] = 0.;
      SolVAR_old[i] = 0.;

      for(int j = 0; j < dim; j++) {
        GradSolVAR[i][j] = 0.;
        GradSolVAR_old[i][j] = 0.;

        GradSolhatVAR[i][j] = 0.;
        GradSolhatVAR_old[i][j] = 0.;
      }

      for(int j = 0; j < nabla_dim; j++) {
        NablaSolVAR[i][j] = 0.;
        NablaSolVAR_old[i][j] = 0.;
      }

      for(unsigned inode = 0; inode < nve; inode++) {
        SolVAR[i] += phi[inode] * Soli[indexVAR[i]][inode];
        SolVAR_old[i] += phi_old[inode] * Soli_old[indexVAR[i]][inode];

        for(int j = 0; j < dim; j++) {
          GradSolVAR[i][j]     += gradphi[inode * dim + j]    * Soli[indexVAR[i]][inode];
          GradSolVAR_old[i][j] += gradphi_old[inode * dim + j] * Soli_old[indexVAR[i]][inode];

          GradSolhatVAR[i][j]     += gradphi_hat[inode * dim + j] * Soli[indexVAR[i]][inode];
          GradSolhatVAR_old[i][j] += gradphi_hat[inode * dim + j] * Soli_old[indexVAR[i]][inode];
        }
        for(int j = 0; j < nabla_dim; j++) {
          NablaSolVAR[i][j]     += nablaphi[inode * nabla_dim + j] * Soli[indexVAR[i]][inode];
          NablaSolVAR_old[i][j] += nablaphi_old[inode * nabla_dim + j] * Soli_old[indexVAR[i]][inode];
        }
      }
    }
    //END evaluates SolVar at idof of element iel


    //TODO i only goes up to dim because we only care about displacement here

    //contribution from the mass of particles
    basis* base = mymsh->GetBasis(ielt, SolType2);
    for(unsigned i = 0; i <  dim; i++) {
      for(int j = 0; j < nve; j++) {
        for(int k = 0; k < nve; k++) {
          double value = base->eval_phi(k, xi);
          aRhs[indexVAR[i]][j] -= mass * value;
        }
        aRhs[indexVAR[i]][j] *= (4 / (dt * dt)) * SolVAR[i] - (4 / dt) * SolVAR[i + dim] - SolVAR[i + 2 * dim] ;
      }
    }

    //contirbution from f_int //TODO
    //need to evaluate the Cauchy stress at the particles

    //contribution from f_ext //TODO
    double body_prtcl = 0.; // value of the body force at the particle
    for(int j = 0; j < nve; j++) {
      body_prtcl += phi[j] ;
    }
    for(unsigned i = 0; i <  dim; i++) {
      for(int j = 0; j < nve; j++) {
	  double value = base->eval_phi(j, xi);
        aRhs[indexVAR[i]][j] -= mass * value * body_prtcl;
      }
    }

    if(iel != ielOld) {

      //copy adouble aRhs into double Rhs
      for(unsigned i = 0; i < 3 * dim; i++) {
        Rhs[indexVAR[i]].resize(nve);

        for(int j = 0; j < nve; j++) {
          Rhs[indexVAR[i]][j] = -aRhs[indexVAR[i]][j].value();
        }
      }

      for(int i = 0; i < 3 * dim + 1; i++) {
        myRES->add_vector_blocked(Rhs[indexVAR[i]], dofsVAR[i]);
      }

      if(assembleMatrix) {
        //Store equations
        for(int i = 0; i < 3 * dim; i++) {
          s.dependent(&aRhs[indexVAR[i]][0], nve);
          s.independent(&Soli[indexVAR[i]][0], nve);
        }

        Jac.resize((3 * dim * nve) * (3 * dim * nve));

        s.jacobian(&Jac[0], true);

        myKK->add_matrix_blocked(Jac, dofsAll, dofsAll);
        s.clear_independents();
        s.clear_dependents();

        //END local to global assembly
      }

      ielOld = iel;
    }

  }
  //END loop on particles

  myRES->close();

  if(assembleMatrix) {
    myKK->close();
  }

  // *************************************
  end_time = clock();
  AssemblyTime += (end_time - start_time);
  // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

}



