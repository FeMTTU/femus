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
  NonLinearImplicitSystem & my_nnlin_impl_sys = ml_prob.get_system<NonLinearImplicitSystem> ("MPM_FEM");
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

  // reserve memory for the local standar vectors
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  // data
  unsigned nel    = mymsh->GetNumberOfElements();
  unsigned iproc  = mymsh->processor_id();

  // local objects
  vector<adept::adouble> SolDp(dim);
  vector<double> SolVp(dim);
  vector<double> SolAp(dim);


  vector<vector < adept::adouble > > GradSolDp(dim);

  for(int i = 0; i < dim; i++) {
    GradSolDp[i].resize(dim);
  }

  //BEGIN do we need these since our phi's are evaluated at the particle local coordinates?
  vector < double > phi;
  vector < double > phi_hat;
  vector < adept::adouble> gradphi;
  vector < double > gradphi_hat;
  vector < adept::adouble> nablaphi;
  vector < double > nablaphi_hat;

  phi.reserve(max_size);
  phi_hat.reserve(max_size);

  gradphi.reserve(max_size * dim);
  gradphi_hat.reserve(max_size * dim);

  //END

  vector <vector < adept::adouble> > vx(dim); //vx is coordX in assembly of ex30
  vector <vector < double> > vx_hat(dim);

  vector< vector< adept::adouble > > SolDd(dim);      // local solution
  vector< vector< double > > SolVd(dim);      // local solution
  vector< vector< double > > SolAd(dim);      // local solution

  vector< vector< int > > dofsVAR(dim);

  vector< vector< double > > Rhs(dim);     // local redidual vector
  vector< vector< adept::adouble > > aRhs(dim);     // local redidual vector

  vector < int > dofsAll;

  vector < double > Jac;

  // gravity
  double _gravity[3] = {0., -1., 0.};

  //variable-name handling
  const char varname[9][3] = {"DX", "DY", "DZ", "U", "V", "W", "AU", "AV", "AW"}; //TODO is there a reason why varname[9][3] and not just varname[9] ?
  vector <unsigned> indexSolD(dim);
  vector <unsigned> indexSolV(dim);
  vector <unsigned> indexSolA(dim);
  vector <unsigned> indexPdeD(dim);
  unsigned solType = ml_sol->GetSolutionType(&varname[0][0]);

  for(unsigned ivar = 0; ivar < dim; ivar++) {
    indexSolD[ivar] = ml_sol->GetIndex(&varname[ivar][0]);
    indexPdeD[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);
    indexSolV[ivar] = ml_sol->GetIndex(&varname[3 + ivar][0]);
    indexSolA[ivar] = ml_sol->GetIndex(&varname[6 + ivar][0]);
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
  unsigned ielOld = UINT_MAX;
  //declaration of element instances


//BEGIN loop on particles (used as Gauss points)
  for(unsigned iMarker = markerOffset1; iMarker < markerOffset2; iMarker++) {
    //element of particle iMarker
    unsigned iel = particles[iMarker]->GetMarkerElement();
    if(iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nve;
      //update element related quantities only if we are in a different element
      if(iel != ielOld) {

        ielt = mymsh->GetElementType(iel);
        nve = mymsh->GetElementDofNumber(iel, solType);
        //initialization of everything is in common fluid and solid

        //Rhs
        for(int i = 0; i < dim; i++) {
          dofsVAR[i].resize(nve);
          SolDd[indexPdeD[i]].resize(nve);
          SolVd[i].resize(nve);
          SolAd[i].resize(nve);
          aRhs[indexPdeD[i]].resize(nve);
        }

        dofsAll.resize(0);

        // ----------------------------------------------------------------------------------------
        // coordinates, solutions, displacement, velocity dofs

        for(int i = 0; i < dim; i++) {
          vx[i].resize(nve);
          vx_hat[i].resize(nve);
        }

        //BEGIN copy of the value of Sol at the dofs idof of the element iel
        for(unsigned i = 0; i < nve; i++) {
          unsigned idof = mymsh->GetSolutionDof(i, iel, solType); //local 2 global solution
          for(int j = 0; j < dim; j++) {
            SolDd[indexPdeD[j]][i] = (*mysolution->_Sol[indexSolD[j]])(idof);
            dofsVAR[j][i] = myLinEqSolver->GetSystemDof(indexSolD[j], indexPdeD[j], i, iel); //local 2 global Pde
            aRhs[indexPdeD[j]][i] = 0.;

            //Fixed coordinates (Reference frame) //TODO do we need this?
            vx_hat[j][i] = (*mymsh->_topology->_Sol[j])(idof);

            SolVd[j][i] = (*mysolution->_Sol[indexSolV[j]])(idof);
            SolAd[j][i] = (*mysolution->_Sol[indexSolA[j]])(idof);
          }
        }
        //END

        // build dof composition
        for(int idim = 0; idim < dim; idim++) {
          dofsAll.insert(dofsAll.end(), dofsVAR[idim].begin(), dofsVAR[idim].end());
        }

        if(assembleMatrix) s.new_recording();

        //TODO do we need this?
        for(unsigned idim = 0; idim < dim; idim++) {
          for(int j = 0; j < nve; j++) {
            vx[idim][j]    = vx_hat[idim][j] + SolDd[indexPdeD[idim]][j];
          }
        }
        // start a new recording of all the operations involving adept::adouble variables
        if(assembleMatrix) s.new_recording();
      }


      bool elementUpdate = (aX.find(iel) != aX.end()) ? false : true;     //update if iel was never updated
      particles[iMarker]->FindLocalCoordinates(solType, aX[iel], elementUpdate, ml_sol->GetLevel(numberOfUniformLevels - 1), 0);

      // the local coordinates of the particles are the Gauss points in this context
      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();
      // the mass of the particle acts as weight
      double mass = particles[iMarker]->GetMarkerMass();
      double density = particles[iMarker]->GetMarkerDensity();
      adept::adouble massAdept;
      mymsh->_finiteElement[ielt][solType]->Jacobian(vx, 0, massAdept, phi, gradphi, nablaphi);
      mymsh->_finiteElement[ielt][solType]->Jacobian(vx_hat, 0, mass, phi_hat, gradphi_hat, nablaphi_hat); //TODO do we need this?

      // displacement and velocity
      //BEGIN evaluates SolDp at the particle iMarker
      //TODO maybe we don't need this for the displacement
      for(int i = 0; i < dim; i++) {
        SolDp[i] = 0.;

        for(int j = 0; j < dim; j++) {
          GradSolDp[i][j] = 0.;
        }

        for(unsigned inode = 0; inode < nve; inode++) {
          SolDp[i] += phi[inode] * SolDd[indexPdeD[i]][inode];

          for(int j = 0; j < dim; j++) {
            GradSolDp[i][j] += gradphi[inode * dim + j] * SolDd[indexPdeD[i]][inode];
          }
        }
      }
      //END evaluates SolVar at idof of element iel

      //BEGIN computation of local residual
      //TODO to double check
      double mu = 500.; //dynamic viscosity
      for(unsigned k = 0; k < nve; k++) {
        for(unsigned i = 0; i < dim; i++) {
          adept::adouble Laplacian = 0.;
          for(unsigned j = 0; j < dim; j++) {
            Laplacian += mu * gradphi[k * dim + j] * (GradSolDp[i][j] + GradSolDp[j][i]);
          }
          aRhs[indexPdeD[i]][k] += - (Laplacian * mass) / density;
        }
      }

      //END

      if(iMarker == markerOffset2 - 1 || iel != particles[iMarker + 1]->GetMarkerElement()) {

        //copy adouble aRhs into double Rhs
        for(unsigned i = 0; i < dim; i++) {
          Rhs[indexPdeD[i]].resize(nve);

          for(int j = 0; j < nve; j++) {
            Rhs[indexPdeD[i]][j] = -aRhs[indexPdeD[i]][j].value();
          }
        }

        for(int i = 0; i < dim; i++) {
          myRES->add_vector_blocked(Rhs[indexPdeD[i]], dofsVAR[i]);
        }

        if(assembleMatrix) {
          //Store equations
          for(int i = 0; i < dim; i++) {
            s.dependent(&aRhs[indexPdeD[i]][0], nve);
            s.independent(&SolDd[indexPdeD[i]][0], nve);
          }

          Jac.resize((dim * nve) * (dim * nve));

          s.jacobian(&Jac[0], true);

          myKK->add_matrix_blocked(Jac, dofsAll, dofsAll);
          s.clear_independents();
          s.clear_dependents();

          //END local to global assembly
        }
      }
      ielOld = iel;
    }
    else {
      break;
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



