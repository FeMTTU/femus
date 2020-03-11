
#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "TransientSystem.hpp"
#include "LinearImplicitSystem.hpp"

#include "NumericVector.hpp"
#include "adept.h"

#include "CurrentElem.hpp"
#include "ElemType_template.hpp"

#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"

#include "slepceps.h"

#include "PetscMatrix.hpp"

using namespace femus;

#define N_UNIFORM_LEVELS  3
#define N_ERASED_LEVELS   2
#define S_FRAC 0.75

#define OP_L2       0
#define OP_H1       0
#define OP_Hhalf    1
#define RHS_ONE     1

#define USE_Cns     1

#define Nsplit      5

#define EX_1       -1.
#define EX_2        1.
#define EY_1       -1.
#define EY_2        1.


double Antiderivative1(const double &theta, const double &s, const double &y);
double Antiderivative2(const double &theta, const double &s, const double &x);
void GetElementPartition1D(const std::vector <double >  & xg1, const std::vector < std::vector <double > > & x1, const unsigned &split,  std::vector < std::vector < std::vector<double>>> &x);
void GetElementPartition2D(const std::vector <double >  & xg1, const std::vector < std::vector <double > > & x1, const unsigned &split,  std::vector < std::vector < std::vector<double>>> &x);
void GetElementPartitionQuad(const std::vector <double >  & xg1, const std::vector < std::vector <double > > & xNodes, const unsigned & split, const unsigned & totalNumberofSplits,  std::vector < std::vector < std::vector<double>>> &x);

double InitialValueU(const std::vector < double >& x)
{
  return 0. * x[0] * x[0];
}

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time)
{
  bool dirichlet = true; //dirichlet
  value = 0.;

//   if(facename == 1) {
//     dirichlet = true; //dirichlet
//     value = 0.;
//   }
//   else if(facename == 2) {
//     dirichlet = true; //dirichlet
//     value = 0.;
//   }

  return dirichlet;
}

double hypergeometric(double a, double b, double c, double x)
{
  const double TOLERANCE = 1.0e-10;
  double term = a * b * x / c;
  double value = 1.0 + term;
  int n = 1;

  while(abs(term) > TOLERANCE) {
    a++, b++, c++, n++;
    term *= a * b * x / c / n;
    value += term;
  }

  return value;
}

void GetHsNorm(const unsigned level, MultiLevelProblem& ml_prob);

void AssembleFracProblem(MultiLevelProblem& ml_prob);


int main(int argc, char** argv)
{


  //quadr rule order
//   const std::string fe_quad_rule_1 = "fifth";
//   const std::string fe_quad_rule_2 = "sixth";
  const std::string fe_quad_rule_1 = "seventh";
  const std::string fe_quad_rule_2 = "eighth";


  //BEGIN deterministic FEM instances

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);


  unsigned numberOfUniformLevels = N_UNIFORM_LEVELS;


  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  unsigned numberOfSelectiveLevels = 0;
  //const std::string mesh_file = "./input/Mesh_1_x.med";
//   const std::string mesh_file = "./input/Mesh_1_x_dir_neu_200_elem.med";
//const std::string mesh_file = "./input/Mesh_1_x_dir_neu.med";
//   const std::string mesh_file = "./input/disk.neu";
//   mlMsh.ReadCoarseMesh(mesh_file.c_str(), fe_quad_rule_1.c_str(), scalingFactor);

//   mlMsh.GenerateCoarseBoxMesh(2, 0, 0, EX_1, EX_2, 0., 0., 0., 0., EDGE3, fe_quad_rule_1.c_str());
//   mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels, NULL);

  mlMsh.GenerateCoarseBoxMesh(2, 2, 0, EX_1, EX_2, EY_1, EY_2, 0., 0., QUAD9, fe_quad_rule_1.c_str());
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels, NULL);

  // erase all the coarse mesh levels
  const unsigned erased_levels = N_ERASED_LEVELS;
  mlMsh.EraseCoarseLevels(erased_levels);

  unsigned dim = mlMsh.GetDimension();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("u", LAGRANGE, SECOND, 2);


  mlSol.Initialize("All");
  mlSol.Initialize("u", InitialValueU);

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  // ******* Set boundary conditions *******
  mlSol.GenerateBdc("All");


  // ========= Problem ==========================
  MultiLevelProblem ml_prob(&mlSol);

  ml_prob.SetQuadratureRuleAllGeomElems(fe_quad_rule_2);
  ml_prob.set_all_abstract_fe();


  // ========= System ==========================
  LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ("FracProblem");
  system.AddSolutionToSystemPDE("u");

  // ******* System FEM Assembly *******
  system.SetAssembleFunction(AssembleFracProblem);
  system.SetMaxNumberOfLinearIterations(1);
  //system.SetAssembleFunction(AssembleFEM);
  // ******* set MG-Solver *******
  system.SetMgType(V_CYCLE);

  system.SetAbsoluteLinearConvergenceTolerance(1.e-50);
  //   system.SetNonLinearConvergenceTolerance(1.e-9);
//   system.SetMaxNumberOfNonLinearIterations(20);

  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);

  // ******* Set Preconditioner *******
  system.SetLinearEquationSolverType(FEMuS_DEFAULT);

  system.init();

  //dense =============
  //dense =============
  //dense =============
  const unsigned solType = mlSol.GetSolutionType("u");

  for(int level = 0; level < mlMsh.GetNumberOfLevels(); level++) {

    Mesh*                    msh = mlMsh.GetLevel(level);
    unsigned    nprocs = msh->n_processors();
    unsigned    iproc = msh->processor_id();

    int MM_size = msh->_dofOffset[solType][nprocs];
    int MM_local_size = msh->_dofOffset[solType][iproc + 1] - msh->_dofOffset[solType][iproc];

//   SparseMatrix* CC;
//   CC = SparseMatrix::build().release();
    system._LinSolver[level]->_KK->init(MM_size, MM_size, MM_local_size, MM_local_size, MM_local_size, MM_size - MM_local_size);
    system._LinSolver[level]->_KK->zero();
  }
  //dense =============
  //dense =============
  //dense =============


  // ******* Set Smoother *******
  system.SetSolverFineGrids(GMRES);

  system.SetPreconditionerFineGrids(ILU_PRECOND);

  system.SetTolerances(1.e-20, 1.e-20, 1.e+50, 100);

  system.MGsolve();


  //solve the generalized eigenvalue problem and compute the eigenpairs
  GetHsNorm(numberOfUniformLevels  - erased_levels - 1, ml_prob);


  // ******* Print solution *******
  mlSol.SetWriter(VTK);
  std::vector<std::string> print_vars;
  print_vars.push_back("All");
  mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);

  //ierr = SlepcFinalize();
  //CHKERRQ(ierr);

  return 0;

} //end main




void AssembleFracProblem(MultiLevelProblem& ml_prob)
{
//void GetEigenPair(MultiLevelProblem & ml_prob, Mat &CCSLEPc, Mat &MMSLEPc) {

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("FracProblem");

  const unsigned level = N_UNIFORM_LEVELS - N_ERASED_LEVELS - 1;


  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             MM = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)


  CurrentElem < double > geom_element1(dim, msh);            // must be adept if the domain is moving, otherwise double
  CurrentElem < double > geom_element2(dim, msh);            // must be adept if the domain is moving, otherwise double

  constexpr unsigned int space_dim = 3;
//***************************************************
  //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
//***************************************************
  std::vector < std::vector < double > >  JacI_qp(space_dim);
  std::vector < std::vector < double > >  Jac_qp(dim);
  for(unsigned d = 0; d < dim; d++) {
    Jac_qp[d].resize(space_dim);
  }
  for(unsigned d = 0; d < space_dim; d++) {
    JacI_qp[d].resize(dim);
  }

  double detJac_qp;
  std::vector < std::vector < /*const*/ elem_type_templ_base< double, double > *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
//***************************************************

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned solType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  std::vector < double > solu1;
  std::vector < double > solu2;

  const unsigned soluPdeIndex = mlPdeSys->GetSolPdeIndex("u");    // get the position of "u" in the pdeSys object

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector < vector < double > > x1(dim);    // local coordinates
  vector < vector < double > > x2(dim);    // local coordinates
  for(unsigned k = 0; k < dim; k++) {
    x1[k].reserve(maxSize);
    x2[k].reserve(maxSize);
  }

  vector < double > phi;
  vector < double > phi_x;

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);

  vector< int > l2GMap1; // local to global mapping
  vector< int > l2GMap2; // local to global mapping
  l2GMap1.reserve(maxSize);
  l2GMap2.reserve(maxSize);

//   Local matrices and rhs for laplacian and mass matrix
  vector < double > MMlocal;
  MMlocal.reserve(maxSize * maxSize);
  vector< double >         Res_local;
  Res_local.reserve(maxSize);  // local redidual vector

//   Local matrices and rhs for adaptive quadrature
  vector< double >         Res_local_refined;
  Res_local_refined.reserve(maxSize);  // local redidual vector
  vector < double > CClocal_refined;
  CClocal_refined.reserve(maxSize * maxSize);

//   Non local matrices and vectors for H^s laplacian operator
//   vector< double >         Res_nonlocal;
//   Res_nonlocal.reserve(maxSize);  // local redidual vector
  vector< double >         Res_nonlocalI;
  Res_nonlocalI.reserve(maxSize);  // local redidual vector
  vector< double >         Res_nonlocalJ;
  Res_nonlocalJ.reserve(maxSize);  // local redidual vector
//   vector < double > CClocal;
//   CClocal.reserve(maxSize * maxSize);
  vector < double > CClocalII;
  CClocalII.reserve(maxSize * maxSize);
  vector < double > CClocalIJ;
  CClocalIJ.reserve(maxSize * maxSize);
  vector < double > CClocalJI;
  CClocalJI.reserve(maxSize * maxSize);
  vector < double > CClocalJJ;
  CClocalJJ.reserve(maxSize * maxSize);

  MM->zero(); // Set to zero all the entries of the Global Matrix
  RES->zero();

//   int MM_size = msh->_dofOffset[solType][nprocs];
//   int MM_local_size = msh->_dofOffset[solType][iproc + 1] - msh->_dofOffset[solType][iproc];
//
//   SparseMatrix* CC;
//   CC = SparseMatrix::build().release();
//   CC->init(MM_size, MM_size, MM_local_size, MM_local_size, MM_local_size, MM_size - MM_local_size);
//   CC->zero();



  const double s_frac = S_FRAC;

  const double check_limits = 1.;//1./(1. - s_frac); // - s_frac;

  double C_ns = 2 * (1 - USE_Cns) + USE_Cns * s_frac * pow(2, (2. * s_frac)) * tgamma((dim + 2. * s_frac) / 2.) / (pow(M_PI, dim / 2.) * tgamma(1 -  s_frac)) ;

  for(int kproc = 0; kproc < nprocs; kproc++) {
    for(int jel = msh->_elementOffset[kproc]; jel < msh->_elementOffset[kproc + 1]; jel++) {

      short unsigned ielGeom2;
      unsigned nDof2;
      unsigned nDofx2;
      unsigned nDofu2;

      if(iproc == kproc) {
        ielGeom2 = msh->GetElementType(jel);
        nDof2  = msh->GetElementDofNumber(jel, solType);    // number of solution element dofs
        nDofx2 = msh->GetElementDofNumber(jel, xType);    // number of coordinate element dofs
      }

      MPI_Bcast(&ielGeom2, 1, MPI_UNSIGNED_SHORT, kproc, MPI_COMM_WORLD);
      MPI_Bcast(&nDof2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
      MPI_Bcast(&nDofx2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);

      // resize local arrays
      l2GMap2.resize(nDof2);


      for(int k = 0; k < dim; k++) {
        x2[k].resize(nDofx2);
      }
      solu2.resize(nDof2);

      // local storage of global mapping and solution ********************
      if(iproc == kproc) {
        for(unsigned j = 0; j < nDof2; j++) {
          l2GMap2[j] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, j, jel);  // global to global mapping between solution node and pdeSys dof
        }
      }
      MPI_Bcast(&l2GMap2[0], nDof2, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
      // ******************************************************************

      // local storage of coordinates  #######################################
      if(iproc == kproc) {
        for(unsigned j = 0; j < nDofx2; j++) {
          unsigned xDof  = msh->GetSolutionDof(j, jel, xType);  // global to global mapping between coordinates node and coordinate dof
          for(unsigned k = 0; k < dim; k++) {
            x2[k][j] = (*msh->_topology->_Sol[k])(xDof);  // global extraction and local storage for the element coordinates
          }
        }
        for(unsigned j = 0; j < nDof2; j++) {
          unsigned jDof  = msh->GetSolutionDof(j, jel, solType);  // global to global mapping between coordinates node and coordinate dof
          solu2[j] = (*sol->_Sol[soluIndex])(jDof);  // global extraction and local storage for the element coordinates
        }
      }
      for(unsigned k = 0; k < dim; k++) {
        MPI_Bcast(& x2[k][0], nDofx2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }
      MPI_Bcast(& solu2[0], nDof2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      // ######################################################################

      // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      if(iproc == kproc) {
        geom_element2.set_coords_at_dofs_and_geom_type(jel, xType);
      }
      for(unsigned k = 0; k < dim; k++) {
        MPI_Bcast(& geom_element2.get_coords_at_dofs()[k][0], nDofx2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }
      for(unsigned k = 0; k < space_dim; k++) {
        MPI_Bcast(& geom_element2.get_coords_at_dofs_3d()[k][0], nDofx2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }
      // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

//       const unsigned jgNumber = msh->_finiteElement[ielGeom2][solType]->GetGaussPointNumber();
      const unsigned jgNumber = ml_prob.GetQuadratureRule(ielGeom2).GetGaussPointsNumber();

      vector < vector < double > > xg2(jgNumber);
      vector <double> weight2(jgNumber);
      vector < vector <double> > phi2(jgNumber);  // local test function
      std::vector< double > solY(jgNumber, 0.);

      for(unsigned jg = 0; jg < jgNumber; jg++) {

//         msh->_finiteElement[ielGeom2][solType]->Jacobian(x2, jg, weight2[jg], phi2[jg], phi_x);

        elem_all[ielGeom2][xType]->JacJacInv(/*x2*/geom_element2.get_coords_at_dofs_3d(), jg, Jac_qp, JacI_qp, detJac_qp, space_dim);
        weight2[jg] = detJac_qp * ml_prob.GetQuadratureRule(ielGeom2).GetGaussWeightsPointer()[jg];
        elem_all[ielGeom2][solType]->shape_funcs_current_elem(jg, JacI_qp, phi2[jg], phi_x /*boost::none*/, boost::none /*phi_u_xx*/, space_dim);




        xg2[jg].assign(dim, 0.);
        solY[jg] = 0.;

        for(unsigned j = 0; j < nDof2; j++) {
          solY[jg] += solu2[j] * phi2[jg][j];
          for(unsigned k = 0; k < dim; k++) {
            xg2[jg][k] += x2[k][j] * phi2[jg][j];
          }
        }
      }


      // element loop: each process loops only on the elements that owns
      for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

        short unsigned ielGeom1 = msh->GetElementType(iel);
        unsigned nDof1  = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
        unsigned nDofx1 = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

        // resize local arrays
        l2GMap1.resize(nDof1);
        //std::vector<bool>bdcDirichlet(nDof1);

        for(int k = 0; k < dim; k++) {
          x1[k].resize(nDofx1);
        }
        solu1.resize(nDof1);

        // local storage of global mapping and solution
        for(unsigned i = 0; i < nDof1; i++) {
          l2GMap1[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
          //unsigned solDof = msh->GetSolutionDof(i, iel, solType);    // global to global mapping between solution node and solution dof
          //bdcDirichlet[i] = ( (*sol->_Bdc[soluIndex])(solDof) < 1.5)? false:false;
        }

        // local storage of coordinates
        for(unsigned i = 0; i < nDofx1; i++) {
          unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof
          for(unsigned k = 0; k < dim; k++) {
            x1[k][i] = (*msh->_topology->_Sol[k])(xDof);  // global extraction and local storage for the element coordinates
          }
        }
        for(unsigned i = 0; i < nDof1; i++) {
          unsigned iDof  = msh->GetSolutionDof(i, iel, solType);  // global to global mapping between coordinates node and coordinate dof
          solu1[i] = (*sol->_Sol[soluIndex])(iDof);  // global extraction and local storage for the element coordinates
        }

//         CClocal.assign(nDof1 * nDof2, 0.);   //resize
        CClocalII.assign(nDof1 * nDof2, 0.);   //resize
        CClocalIJ.assign(nDof1 * nDof2, 0.);   //resize
        CClocalJI.assign(nDof1 * nDof2, 0.);   //resize
        CClocalJJ.assign(nDof1 * nDof2, 0.);   //resize
//         Res_nonlocal.assign(nDof1, 0);    //resize
        Res_nonlocalI.assign(nDof1, 0);    //resize
        Res_nonlocalJ.assign(nDof1, 0);    //resize

        if(iel == jel) {
          Res_local.assign(nDof1, 0);    //resize
          MMlocal.assign(nDof1 * nDof1, 0.);
          if(Nsplit != 0) {
//             Vectors and matrices for adaptive quadrature
            Res_local_refined.assign(nDof1, 0);    //resize
            CClocal_refined.assign(nDof1 * nDof1, 0.);
          }
        }

        // *** Gauss point loop ***
        const unsigned igNumber = msh->_finiteElement[ielGeom1][solType]->GetGaussPointNumber();


        double weight1;
        vector < double > phi1;  // local test function

        double weight3;
        vector < double > phi3;  // local test function
        double weight4;
        vector < double > phi4;  // local test function
        double solX = 0.;
        std::vector<double> sol_u_x(space_dim);
        std::fill(sol_u_x.begin(), sol_u_x.end(), 0.);


        std::vector < std::vector < std::vector <double > > > aP(3);
        if(Nsplit > 0) {
          for(unsigned jtype = 0; jtype < solType + 1; jtype++) {
            ProjectNodalToPolynomialCoefficients(aP[jtype], x1, ielGeom1, jtype) ;
          }
        }


        for(unsigned ig = 0; ig < igNumber; ig++) {

          msh->_finiteElement[ielGeom1][solType]->Jacobian(x1, ig, weight1, phi1, phi_x);

          // evaluate the solution, the solution derivatives and the coordinates in the gauss point
          vector < double > xg1(dim, 0.);
          solX = 0.;

          for(unsigned i = 0; i < nDof1; i++) {
            solX += solu1[i] * phi1[i];
            for(unsigned d = 0; d < sol_u_x.size(); d++)   sol_u_x[d] += solu1[i] * phi_x[i * dim + d];
            for(unsigned k = 0; k < dim; k++) {
              xg1[k] += x1[k][i] * phi1[i];
            }
          }

          if(iel == jel) {

            for(unsigned i = 0; i < nDof1; i++) {
              for(unsigned j = 0; j < nDof1; j++) {
                MMlocal[ i * nDof1 + j ] += OP_L2 * phi1[i] * phi1[j] * weight1;
              }
              double mass_res_i = phi1[i] * solX ;
              Res_local[ i ] += OP_L2 * weight1 * mass_res_i ;
              Res_local[ i ] += - RHS_ONE * weight1 * (phi1[i] * (-1.));
            }

//          ---------------------
//          Laplacian assembly
//          Residual
            std::fill(sol_u_x.begin(), sol_u_x.end(), 0.);
            for(unsigned i = 0; i < nDof1; i++) {
              double laplace_res_du_u_i = 0.;
              for(unsigned kdim = 0; kdim < dim; kdim++) {
                laplace_res_du_u_i  +=  phi_x   [i * dim + kdim] * sol_u_x[kdim];
              }
              Res_local[ i ] += - OP_H1 * weight1 * (- laplace_res_du_u_i);

//          Matrix
              for(unsigned j = 0; j < nDof1; j++) {

                double laplace_mat_i_j = 0.;
                for(unsigned kdim = 0; kdim < dim; kdim++) {
                  laplace_mat_i_j    += phi_x   [i * dim + kdim] *
                                        phi_x   [j * dim + kdim];
                }
                MMlocal[ i * nDof1 + j ]  += OP_H1 * weight1 *  laplace_mat_i_j;
              }
            }
//          ---------------------
//          Mixed integral ((Rn-Omega) x Omega) assembly (based on the analytic result of integrals)
            if(dim == 1) {
//               double ex_1 = EX_1;
//               double ex_2 = EX_2;
//               double dist_1 = 0.;
//               double dist_2 = 0.;
//               for(int k = 0; k < dim; k++) {
//                 dist_1 += sqrt((xg1[k] - ex_1) * (xg1[k] - ex_1));
//                 dist_2 += sqrt((xg1[k] - ex_2) * (xg1[k] - ex_2));
//               }
//               double mixed_term = pow(dist_1, -2. * s_frac) + pow(dist_2, - 2. * s_frac);
//
//               for(unsigned i = 0; i < nDof1; i++) {
//                 for(unsigned j = 0; j < nDof1; j++) {
//                   MMlocal[ i * nDof1 + j ] += (C_ns / 2.) * check_limits * (1. / s_frac) * OP_Hhalf * phi1[i] * phi1[j] * weight1 * mixed_term;
//                 }
//                 Res_local[ i ] += (C_ns / 2.) * check_limits * (1. / s_frac) * OP_Hhalf * weight1 * phi1[i] * solX * mixed_term;
//               }
            }
            if(dim == 2) {
              double ex[4] = {EX_1, EX_2, EY_1, EY_2};

//               double CC = 1. / (2 * s_frac * (1 + 2. * s_frac));

              double teta[4], CCC[4];
              teta[0] = atan2((ex[3] - xg1[1]), (ex[1] - xg1[0]));
              teta[1] = atan2((ex[3] - xg1[1]), (ex[0] - xg1[0]));
              teta[2] = atan2((ex[2] - xg1[1]), (ex[0] - xg1[0])) + 2 * M_PI;
              teta[3] = atan2((ex[2] - xg1[1]), (ex[1] - xg1[0])) + 2 * M_PI;

              double mixed_term = 0.;

              for(unsigned qq = 0; qq < 4; qq++) {
                if(qq == 3) teta[3] -= 2. * M_PI ;

                if(qq == 0)
                  mixed_term += 2.* (Antiderivative1(teta[1], s_frac, ex[3] - xg1[1]) -
                                     Antiderivative1(teta[0], s_frac, ex[3] - xg1[1]));
                else if(qq  == 2)
                  mixed_term += 2.* (Antiderivative1(teta[3], s_frac, ex[2] - xg1[1]) -
                                     Antiderivative1(teta[2], s_frac, ex[2] - xg1[1]));
                else if(qq  == 1)
                  mixed_term += 2. * (Antiderivative2(teta[2], s_frac, ex[0] - xg1[0]) -
                                      Antiderivative2(teta[1], s_frac, ex[0] - xg1[0]));
                else
                  mixed_term += 2. * (Antiderivative2(teta[0], s_frac, ex[1] - xg1[0]) -
                                      Antiderivative2(teta[3], s_frac, ex[1] - xg1[0]));


              }
              
              
              
              
              
// //     New approach for numerical integral (START)
// //     -----------------------------------
// //     -----------------------------------              
//               
//     double mixed_term1 = 0;
//     for(int kel = msh->_elementOffset[iproc]; kel < msh->_elementOffset[iproc + 1]; kel++) {
//                   // *** Face Gauss point loop (boundary Integral) ***
//     for ( unsigned jface = 0; jface < msh->GetElementFaceNumber ( kel ); jface++ ) {
//       int faceIndex = el->GetBoundaryIndex(kel, jface);
//       // look for boundary faces
//       if ( faceIndex >= 1 ) {  
//         const unsigned faceGeom = msh->GetElementFaceType ( kel, jface );
//         unsigned faceDofs = msh->GetElementFaceDofNumber (kel, jface, solType);         
//         vector  < vector  <  double> > faceCoordinates ( dim ); // A matrix holding the face coordinates rowwise.
//         for ( int k = 0; k < dim; k++ ) {
//           faceCoordinates[k].resize (faceDofs);
//         }
//         for ( unsigned i = 0; i < faceDofs; i++ ) {
//           unsigned inode = msh->GetLocalFaceVertexIndex ( kel, jface, i ); // face-to-element local node mapping.
//           for ( unsigned k = 0; k < dim; k++ ) {
//             faceCoordinates[k][i] =  x1[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
//           }
//         }
// //         valido per 2D, verifica che funzioni!!
//         double teta2 = atan2((faceCoordinates[1][1] - xg1[1]), ( faceCoordinates[0][1] - xg1[0]));
//         double teta1 = atan2((faceCoordinates[1][0] - xg1[1]), ( faceCoordinates[0][0] - xg1[0])); 
//         
//         double delta_teta = fabs ( teta2 - teta1 );
//         
//         vector <double> mid_sur;
//         mid_sur.resize(dim);
//         for( unsigned k = 0; k < dim; k++ ) {
//           mid_sur[k] = ( faceCoordinates[k][1] + faceCoordinates[k][0] ) * 0.5;
//         }
//         double dist2 = 0;
//         for(int k = 0; k < dim; k++) {
//           dist2 += (xg1[k] - mid_sur[k]) * (xg1[k] - mid_sur[k]);
//         }
//         double dist = sqrt( dist2 );
//         mixed_term1 += (1. / (2. * s_frac )) * pow(dist, -  2. * s_frac) * delta_teta;
//         
//       }
//     } 
//     }
//     
// //     New approach for numerical integral (END)
// //     -----------------------------------
// //     -----------------------------------               
              
              


              for(unsigned i = 0; i < nDof1; i++) {
                for(unsigned j = 0; j < nDof1; j++) {
                  MMlocal[ i * nDof1 + j ] += (C_ns / 2.) * check_limits /** (1. / s_frac)*/ * OP_Hhalf * phi1[i] * phi1[j] * weight1 * mixed_term;
                }
                Res_local[ i ] += (C_ns / 2.) * check_limits /** (1. / s_frac)*/ * OP_Hhalf * weight1 * phi1[i] * solX * mixed_term;
              }
            }

//          ---------------------
//          Adaptive quadrature for iel == jel


            if(Nsplit != 0) {

              std::cout.precision(14);
              std::vector< std::vector<std::vector<double>>> x3;
              
              for(unsigned split = 0; split <= Nsplit; split++) {

//                 unsigned size_part;
//                 if(dim == 1) size_part = 2;
//                 else size_part = (split != Nsplit) ? 12 : 4;
                
                if(dim == 1) GetElementPartition1D(xg1, x1, split, x3);
                else if(dim == 2) {
                  //GetElementPartition2D(xg1, x1, split, x3);
                  GetElementPartitionQuad(xg1, x1, split, Nsplit, x3);
                }
                
                //for(unsigned r = 0; r < size_part; r++) {
                for(unsigned r = 0; r < x3.size(); r++) {
                

                  for(unsigned jg = 0; jg < igNumber; jg++) {


                    msh->_finiteElement[ielGeom1][solType]->Jacobian(x3[r], jg, weight3, phi3, phi_x);
                    
                    vector < double > xg3(dim, 0.);

                    for(unsigned i = 0; i < nDof1; i++) {
                      for(unsigned k = 0; k < dim; k++) {
                        xg3[k] += x3[r][k][i] * phi3[i];
                      }
                    }

                    std::vector<double> xi3(dim, 0.);

                    GetClosestPointInReferenceElement(x1, xg3, ielGeom1, xi3);
                    GetInverseMapping(solType, ielGeom1, aP, xg3, xi3, 1000);

                    msh->_finiteElement[ielGeom1][solType]->GetPhi(phi3, xi3);

                    double solY3 = 0.;
                    for(unsigned i = 0; i < nDof1; i++) {
                      solY3 += solu1[i] * phi3[i];
                    }

                    double dist_xyz3 = 0;
                    for(unsigned k = 0; k < dim; k++) {
                      dist_xyz3 += (xg1[k] - xg3[k]) * (xg1[k] - xg3[k]);
                    }

                    const double denom3 = pow(dist_xyz3, (double)((dim / 2.) + s_frac));
                    
                    for(unsigned i = 0; i < nDof1; i++) {

                      Res_local_refined[ i ]    +=      - (C_ns / 2.) * OP_Hhalf * check_limits *
                                                        ((solX - solY3) * (phi1[i] - phi3[i]) * weight3 / denom3
                                                        ) * weight1 ;

                      for(unsigned j = 0; j < nDof2; j++) {
                        CClocal_refined[ i * nDof2 + j ] += (C_ns / 2.) * OP_Hhalf * check_limits *
                                                            ((phi1[j] - phi3[j]) * (phi1[i] - phi3[i]) * weight3 / denom3
                                                            ) * weight1 ;

                      }
                    }
                    
// //               Mixed integral 1D
// //               ------------------
                    if(ig == 0 && dim == 1 ) {
                      double ex_1 = EX_1;
                      double ex_2 = EX_2;
                      double dist_1 = 0.;
                      double dist_2 = 0.;
                      for(int k = 0; k < dim; k++) {
                        dist_1 += sqrt((xg3[k] - ex_1) * (xg3[k] - ex_1));
                        dist_2 += sqrt((xg3[k] - ex_2) * (xg3[k] - ex_2));
                      }
                      double mixed_term = pow(dist_1, -2. * s_frac) + pow(dist_2, - 2. * s_frac);

                      for(unsigned i = 0; i < nDof1; i++) {
                        for(unsigned j = 0; j < nDof1; j++) {
                          MMlocal[ i * nDof1 + j ] += (C_ns / 2.) * check_limits * (1. / s_frac) * OP_Hhalf * phi3[i] * phi3[j] * weight3 * mixed_term;
                        }
                        Res_local[ i ] += (C_ns / 2.) * check_limits * (1. / s_frac) * OP_Hhalf * weight3 * phi3[i] * solY3 * mixed_term;
                      }
                    }

                  }
                }
              }
            }

          } // end iel == jel loop
          
          
          
          
          
          
          //     New approach for numerical integral (START)
//     -----------------------------------
//     -----------------------------------              
              
    double mixed_term1 = 0;
//     for(int kel = msh->_elementOffset[iproc]; kel < msh->_elementOffset[iproc + 1]; kel++) {
                  // *** Face Gauss point loop (boundary Integral) ***
    for ( unsigned jface = 0; jface < msh->GetElementFaceNumber ( jel ); jface++ ) {
      int faceIndex = el->GetBoundaryIndex(jel, jface);
      // look for boundary faces
      if ( faceIndex >= 1 ) {  
        const unsigned faceGeom = msh->GetElementFaceType ( jel, jface );
        unsigned faceDofs = msh->GetElementFaceDofNumber (jel, jface, solType);         
        vector  < vector  <  double> > faceCoordinates ( dim ); // A matrix holding the face coordinates rowwise.
        for ( int k = 0; k < dim; k++ ) {
          faceCoordinates[k].resize (faceDofs);
        }
        for ( unsigned i = 0; i < faceDofs; i++ ) {
          unsigned inode = msh->GetLocalFaceVertexIndex ( jel, jface, i ); // face-to-element local node mapping.
          for ( unsigned k = 0; k < dim; k++ ) {
            faceCoordinates[k][i] =  x1[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
          }
        }
//         valido per 2D, verifica che funzioni!!
        double teta2 = atan2((faceCoordinates[1][1] - xg1[1]), ( faceCoordinates[0][1] - xg1[0]));
        double teta1 = atan2((faceCoordinates[1][0] - xg1[1]), ( faceCoordinates[0][0] - xg1[0])); 
        
        double delta_teta = fabs ( teta2 - teta1 );
        
        vector <double> mid_sur;
        mid_sur.resize(dim);
        for( unsigned k = 0; k < dim; k++ ) {
          mid_sur[k] = ( faceCoordinates[k][1] + faceCoordinates[k][0] ) * 0.5;
        }
        double dist2 = 0;
        for(int k = 0; k < dim; k++) {
          dist2 += (xg1[k] - mid_sur[k]) * (xg1[k] - mid_sur[k]);
        }
        double dist = sqrt( dist2 );
        mixed_term1 += pow(dist, -  2. * s_frac) * delta_teta;
        
      }
    } 
//     }
    
//     New approach for numerical integral (END)
//     -----------------------------------
//     ----------------------------------- 
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          

          if((Nsplit == 0 || iel != jel) && OP_Hhalf != 0) {
            for(unsigned jg = 0; jg < jgNumber; jg++) {

              double dist_xyz = 0;
              for(unsigned k = 0; k < dim; k++) {
                dist_xyz += (xg1[k] - xg2[jg][k]) * (xg1[k] - xg2[jg][k]);
              }

              const double denom = pow(dist_xyz, (double)((dim / 2.) + s_frac));

              for(unsigned i = 0; i < nDof1; i++) {

//                Res_nonlocal[ i ]         +=      - (C_ns / 2.) * OP_Hhalf *  check_limits * (solX - solY[jg]) * (phi1[i] - phi2[jg][i]) * weight1 * weight2[jg]  / denom;

                Res_nonlocalI[ i ]         +=      - (C_ns / 2.) * OP_Hhalf *  check_limits * (solX - solY[jg]) * (phi1[i]) * weight1 * weight2[jg]  / denom;

                Res_nonlocalJ[ i ]         +=      - (C_ns / 2.) * OP_Hhalf *  check_limits * (solX - solY[jg]) * (- phi2[jg][i]) * weight1 * weight2[jg]  / denom;

                for(unsigned j = 0; j < nDof2; j++) {
//                 CClocal[ i * nDof2 + j ] += (C_ns / 2.) * OP_Hhalf * check_limits * (phi1[j] - phi2[jg][j]) * (phi1[i] - phi2[jg][i]) * weight1 * weight2[jg] / denom;

                  CClocalII[ i * nDof2 + j ] += (C_ns / 2.) * OP_Hhalf * check_limits * phi1[j]  * phi1[i] * weight1 * weight2[jg] / denom;

                  CClocalIJ[ i * nDof2 + j ] += (C_ns / 2.) * OP_Hhalf * check_limits * (- phi2[jg][j]) * phi1[i] * weight1 * weight2[jg] / denom;

                  CClocalJI[ i * nDof2 + j ] += (C_ns / 2.) * OP_Hhalf * check_limits * (phi1[j]) * (- phi2[jg][i]) * weight1 * weight2[jg] / denom;

                  CClocalJJ[ i * nDof2 + j ] += (C_ns / 2.) * OP_Hhalf * check_limits * (- phi2[jg][j]) * (- phi2[jg][i]) * weight1 * weight2[jg] / denom;


                }

              }



            } //endl jg loop

          }

        } //endl ig loop

        if(iel == jel) {
          MM->add_matrix_blocked(MMlocal, l2GMap1, l2GMap1);
          RES->add_vector_blocked(Res_local, l2GMap1);

          if(Nsplit != 0) {
            MM->add_matrix_blocked(CClocal_refined, l2GMap1, l2GMap1);
            RES->add_vector_blocked(Res_local_refined, l2GMap1);
          }
        }
//        MM->add_matrix_blocked(CClocal, l2GMap1, l2GMap2);

        MM->add_matrix_blocked(CClocalII, l2GMap1, l2GMap1);
        MM->add_matrix_blocked(CClocalIJ, l2GMap1, l2GMap2);
        MM->add_matrix_blocked(CClocalJI, l2GMap2, l2GMap1);
        MM->add_matrix_blocked(CClocalJJ, l2GMap2, l2GMap2);

//        RES->add_vector_blocked(Res_nonlocal, l2GMap1);
        RES->add_vector_blocked(Res_nonlocalI, l2GMap1);
        RES->add_vector_blocked(Res_nonlocalJ, l2GMap2);
      } // end iel loop


    } //end jel loop
  } //end kproc loop

  MM->close();
  RES->close();



//   PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
//   PetscObjectSetName((PetscObject)viewer, "FSI matrix");
//   PetscViewerPushFormat(viewer, PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*>(MM))->mat(), viewer);
// //   MatView((static_cast<PetscMatrix*> (MM))->mat(),  PETSC_VIEWER_STDOUT_WORLD );
//   double a;
//   std::cin >> a;


//     PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,900,900,&viewer);
//   PetscObjectSetName((PetscObject)viewer,"FSI matrix");
//   PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*> (MM))->mat(),viewer);
// //   MatView((static_cast<PetscMatrix*> (MM))->mat(),  PETSC_VIEWER_STDOUT_WORLD );
//
//   VecView((static_cast<PetscVector*> (RES))->vec(),  PETSC_VIEWER_STDOUT_WORLD );
//   double a;
//   std::cin>>a;

  // ***************** END ASSEMBLY *******************
}










void GetHsNorm(const unsigned level,  MultiLevelProblem& ml_prob)
{


  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)


  CurrentElem < double > geom_element1(dim, msh);            // must be adept if the domain is moving, otherwise double
  CurrentElem < double > geom_element2(dim, msh);            // must be adept if the domain is moving, otherwise double

  constexpr unsigned int space_dim = 3;
//***************************************************
  //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
//***************************************************
  std::vector < std::vector < double > >  JacI_qp(space_dim);
  std::vector < std::vector < double > >  Jac_qp(dim);
  for(unsigned d = 0; d < dim; d++) {
    Jac_qp[d].resize(space_dim);
  }
  for(unsigned d = 0; d < space_dim; d++) {
    JacI_qp[d].resize(dim);
  }

  double detJac_qp;
  std::vector < std::vector < /*const*/ elem_type_templ_base< double, double > *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
//***************************************************

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned solType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  std::vector < double > solu1;
  std::vector < double > solu2;


  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector < vector < double > > x1(dim);    // local coordinates
  vector < vector < double > > x2(dim);    // local coordinates
  for(unsigned k = 0; k < dim; k++) {
    x1[k].reserve(maxSize);
    x2[k].reserve(maxSize);
  }

  vector < double > phi;
  vector < double > phi_x;

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);



  const double s_frac = S_FRAC;

  double sol_qp = 0.;
  std::vector< double > sol_x_qp(space_dim);
  std::fill(sol_x_qp.begin(), sol_x_qp.end(), 0.);
  double JxWeight = 0.;

  double integral_iproc_L2 = 0.;
  double integral_iproc_H1 = 0.;


  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    geom_element1.set_coords_at_dofs_and_geom_type(iel, xType);

    const short unsigned ielGeom1 = geom_element1.geom_type();

    unsigned nDof_u  = msh->GetElementDofNumber(iel, solType);
    solu1.resize(nDof_u);

    for(unsigned i = 0; i < solu1.size(); i++) {
      unsigned iDof  = msh->GetSolutionDof(i, iel, solType);  // global to global mapping between coordinates node and coordinate dof
      solu1[i] = (*sol->_Sol[soluIndex])(iDof);  // global extraction and local storage for the element coordinates
    }



    const unsigned igNumber = ml_prob.GetQuadratureRule(ielGeom1).GetGaussPointsNumber();


    for(unsigned ig = 0; ig < igNumber; ig++) {

      elem_all[ielGeom1][xType]->JacJacInv(geom_element1.get_coords_at_dofs_3d(), ig, Jac_qp, JacI_qp, detJac_qp, space_dim);

      JxWeight = detJac_qp * ml_prob.GetQuadratureRule(ielGeom1).GetGaussWeightsPointer()[ig];

      elem_all[ielGeom1][solType]->shape_funcs_current_elem(ig, JacI_qp, phi, phi_x /*boost::none*/, boost::none /*phi_xx*/, space_dim);


      sol_qp = 0.;
      std::fill(sol_x_qp.begin(), sol_x_qp.end(), 0.);

      for(unsigned i = 0; i <  solu1.size(); i++) {
        sol_qp += solu1[i] * phi[i];
        for(unsigned d = 0; d < sol_x_qp.size(); d++)   sol_x_qp[d] += solu1[i] * phi_x[i * space_dim + d];

      }

      integral_iproc_L2 += JxWeight * sol_qp * sol_qp;

      for(unsigned d = 0; d < sol_x_qp.size(); d++) integral_iproc_H1 += JxWeight * sol_x_qp[d] * sol_x_qp[d];


    }


  }


  printf("L2 integral on processor %d = %f \n", iproc, integral_iproc_L2);

  double J_L2 = 0.;
  MPI_Allreduce(&integral_iproc_L2, &J_L2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);     //THIS IS THE RIGHT ONE!!

  std::cout << "L2 integral after Allreduce: " << sqrt(J_L2) << std::endl;


  printf("H1 integral on processor %d = %f \n", iproc, integral_iproc_H1);

  double J_H1 = 0.;
  MPI_Allreduce(&integral_iproc_H1, &J_H1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);     //THIS IS THE RIGHT ONE!!

  std::cout << "H1 integral after Allreduce: " << sqrt(J_H1) << std::endl;




  double integral_iproc_Hhalf = 0.;


  for(int kproc = 0; kproc < nprocs; kproc++) {
    for(int jel = msh->_elementOffset[kproc]; jel < msh->_elementOffset[kproc + 1]; jel++) {

      short unsigned ielGeom2;
      unsigned nDof2;
      unsigned nDofx2;
      unsigned nDofu2;

      if(iproc == kproc) {
        ielGeom2 = msh->GetElementType(jel);
        nDof2  = msh->GetElementDofNumber(jel, solType);    // number of solution element dofs
        nDofx2 = msh->GetElementDofNumber(jel, xType);    // number of coordinate element dofs
      }

      MPI_Bcast(&ielGeom2, 1, MPI_UNSIGNED_SHORT, kproc, MPI_COMM_WORLD);
      MPI_Bcast(&nDof2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
      MPI_Bcast(&nDofx2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);



      for(int k = 0; k < dim; k++) {
        x2[k].resize(nDofx2);
      }
      solu2.resize(nDof2);


      // local storage of coordinates  #######################################
      if(iproc == kproc) {
        for(unsigned j = 0; j < nDofx2; j++) {
          unsigned xDof  = msh->GetSolutionDof(j, jel, xType);  // global to global mapping between coordinates node and coordinate dof
          for(unsigned k = 0; k < dim; k++) {
            x2[k][j] = (*msh->_topology->_Sol[k])(xDof);  // global extraction and local storage for the element coordinates
          }
        }
        for(unsigned j = 0; j < nDof2; j++) {
          unsigned jDof  = msh->GetSolutionDof(j, jel, solType);  // global to global mapping between coordinates node and coordinate dof
          solu2[j] = (*sol->_Sol[soluIndex])(jDof);  // global extraction and local storage for the element coordinates
        }
      }
      for(unsigned k = 0; k < dim; k++) {
        MPI_Bcast(& x2[k][0], nDofx2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }
      MPI_Bcast(& solu2[0], nDof2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      // ######################################################################

      // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      if(iproc == kproc) {
        geom_element2.set_coords_at_dofs_and_geom_type(jel, xType);
      }
      for(unsigned k = 0; k < dim; k++) {
        MPI_Bcast(& geom_element2.get_coords_at_dofs()[k][0], nDofx2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }
      for(unsigned k = 0; k < space_dim; k++) {
        MPI_Bcast(& geom_element2.get_coords_at_dofs_3d()[k][0], nDofx2, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }
      // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      const unsigned jgNumber = ml_prob.GetQuadratureRule(ielGeom2).GetGaussPointsNumber();

      vector < vector < double > > xg2(jgNumber);
      vector <double> weight2(jgNumber);
      vector < vector <double> > phi2(jgNumber);  // local test function
      std::vector< double > solY(jgNumber, 0.);

      for(unsigned jg = 0; jg < jgNumber; jg++) {

//          msh->_finiteElement[ielGeom2][solType]->Jacobian(x2, jg, weight2[jg], phi2[jg], phi_x);

        elem_all[ielGeom2][xType]->JacJacInv(/*x2*/geom_element2.get_coords_at_dofs_3d(), jg, Jac_qp, JacI_qp, detJac_qp, space_dim);
        weight2[jg] = detJac_qp * ml_prob.GetQuadratureRule(ielGeom2).GetGaussWeightsPointer()[jg];
        elem_all[ielGeom2][solType]->shape_funcs_current_elem(jg, JacI_qp, phi2[jg], phi_x /*boost::none*/, boost::none /*phi_u_xx*/, space_dim);

        xg2[jg].assign(dim, 0.);
        solY[jg] = 0.;

        for(unsigned j = 0; j < nDof2; j++) {
          solY[jg] += solu2[j] * phi2[jg][j];
          for(unsigned k = 0; k < dim; k++) {
            xg2[jg][k] += x2[k][j] * phi2[jg][j];
          }
        }
      }

      // element loop: each process loops only on the elements that owns
      for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

        short unsigned ielGeom1 = msh->GetElementType(iel);
        unsigned nDof1  = msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
        unsigned nDofx1 = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

        // resize local arrays
//         l2GMap1.resize(nDof1);
        //std::vector<bool>bdcDirichlet(nDof1);

        for(int k = 0; k < dim; k++) {
          x1[k].resize(nDofx1);
        }
        solu1.resize(nDof1);


        // local storage of coordinates
        for(unsigned i = 0; i < nDofx1; i++) {
          unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof
          for(unsigned k = 0; k < dim; k++) {
            x1[k][i] = (*msh->_topology->_Sol[k])(xDof);  // global extraction and local storage for the element coordinates
          }
        }
        for(unsigned i = 0; i < nDof1; i++) {
          unsigned iDof  = msh->GetSolutionDof(i, iel, solType);  // global to global mapping between coordinates node and coordinate dof
          solu1[i] = (*sol->_Sol[soluIndex])(iDof);  // global extraction and local storage for the element coordinates
        }


        // *** Gauss point loop ***
        const unsigned igNumber = msh->_finiteElement[ielGeom1][solType]->GetGaussPointNumber();
//         const unsigned igNumber = ml_prob.GetQuadratureRule(ielGeom1).GetGaussPointsNumber();

        double weight1;
        vector < double > phi1;  // local test function
        double  solX = 0.;

        for(unsigned ig = 0; ig < igNumber; ig++) {

          msh->_finiteElement[ielGeom1][solType]->Jacobian(x1, ig, weight1, phi1, phi_x);

          // evaluate the solution, the solution derivatives and the coordinates in the gauss point
          vector < double > xg1(dim, 0.);
          solX = 0.;

          for(unsigned i = 0; i < nDof1; i++) {
            solX += solu1[i] * phi1[i];
            for(unsigned k = 0; k < dim; k++) {
              xg1[k] += x1[k][i] * phi1[i];
            }
          }


          for(unsigned jg = 0; jg < jgNumber; jg++) {

            double dist_xyz = 0;
            for(unsigned k = 0; k < dim; k++) {
              dist_xyz += (xg1[k] - xg2[jg][k]) * (xg1[k] - xg2[jg][k]);
            }

            const double denom = pow(dist_xyz, (dim / 2.) + s_frac);

            const double sol_diff = (solX - solY[jg]);

//             integral_iproc_Hhalf +=  weight1 *  weight2[jg];
            integral_iproc_Hhalf +=  weight1 * weight2[jg] * (sol_diff * sol_diff) / denom;


          } //endl jg loop
        } //endl ig loop


      } // end iel loop


    } //end jel loop
  } //end kproc loop

  ////////////////////////////////////////

  printf("H-1/2 integral on processor %d = %.40f \n", iproc, integral_iproc_Hhalf);

  double J = 0.;
  MPI_Allreduce(&integral_iproc_Hhalf, &J, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);     //THIS IS THE RIGHT ONE!!

  std::cout << "H-1/2 integral after Allreduce squared: " << J << std::endl;
  std::cout << "H-1/2 integral after Allreduce: " << sqrt(J) << std::endl;

  //return;                                                  //ignore the rest

  // ***************** END ASSEMBLY *******************
}

double Antiderivative1(const double &theta, const double &s, const double &y)
{
  return -(1. / tan(theta) * hypergeometric(0.5, 0.5 - s, 1.5, pow(cos(theta), 2.)) * pow(pow(sin(theta), 2.), 0.5 - s)) /
         (2.* s * pow(y * 1 / sin(theta), 2. * s));
}

double Antiderivative2(const double &theta, const double &s, const double &x)
{
  return (pow(pow(cos(theta), 2.), 0.5 - s) * hypergeometric(0.5, 0.5 - s, 1.5, pow(sin(theta), 2)) * tan(theta)) /
         (2.* s * pow(x * 1. / cos(theta), 2. * s));
}



void GetElementPartition1D(const std::vector <double >  & xg1, const std::vector < std::vector <double > > & x1, const unsigned &split,  std::vector < std::vector < std::vector<double>>> &x)
{
  unsigned dim = 1;
  unsigned left = 0;
  unsigned right = 1;

  if(split == 0) { //init
    x.resize(2);
    x[left].resize(dim);
    x[right].resize(dim);
    for(unsigned k = 0; k < dim; k++) {
      x[left][k].resize(x1[0].size());
      x[right][k].resize(x1[0].size());
//       for(unsigned k = 0; k < dim; k++) {
        x[left][k][0] = x1[k][0];
        x[left][k][1] = 0.5 * (x[left][k][0] + xg1[k]);
        x[left][k][2] = 0.5 * (x[left][k][0] + x[left][k][1]);
        x[right][k][1] = x1[k][1];
        x[right][k][0] = 0.5 * (x[right][k][1] + xg1[k]);
        x[right][k][2] = 0.5 * (x[right][k][0] + x[right][k][1]);
//       }
    }
  }
  else if(split == Nsplit) {
    for(unsigned k = 0; k < dim; k++) {
      x[left][k][0] = x[left][k][1];
      x[left][k][1] = xg1[k];
      x[left][k][2] = 0.5 * (x[left][k][0] + x[left][k][1]);

      x[right][k][1] = x[right][k][0];
      x[right][k][0] = xg1[k];
      x[right][k][2] = 0.5 * (x[right][k][0] + x[right][k][1]);
    }
  }
  else {
    for(unsigned k = 0; k < dim; k++) {
      x[left][k][0] = x[left][k][1];
      x[left][k][1] = 0.5 * (x[left][k][0] + xg1[k]);
      x[left][k][2] = 0.5 * (x[left][k][0] + x[left][k][1]);

      x[right][k][1] = x[right][k][0];
      x[right][k][0] = 0.5 * (x[right][k][1] + xg1[k]);
      x[right][k][2] = 0.5 * (x[right][k][0] + x[right][k][1]);
    }
  }
}

// double x1i[1][3]={{-1.,1.,0.}};
//
// void GetElementPartition1D2(const std::vector <double >  & xig,
//                             const unsigned &split,
//                             const std::vector < std::vector <double > > & x1,
//                             std::vector < std::vector < std::vector<double>>> &xi,
//                             std::vector < std::vector < std::vector<double>>> &x) {
//   unsigned dim = 1;
//   unsigned left = 0;
//   unsigned right = 1;
//
//   if(split == 0) { //init
//     xi.resize(2);
//     xi[left].resize(dim);
//     xi[right].resize(dim);
//     for(unsigned k = 0; k < dim; k++) {
//       xi[left][k].resize(3);
//       xi[right][k].resize(3);
//       for(unsigned k = 0; k < dim; k++) {
//         xi[left][k][0] = x1i[k][0];
//         xi[left][k][1] = 0.5 * (xi[left][k][0] + xig[k]);
//         xi[left][k][2] = 0.5 * (xi[left][k][0] + xi[left][k][1]);
//         xi[right][k][1] = x1i[k][1];
//         xi[right][k][0] = 0.5 * (xi[right][k][1] + xig[k]);
//         xi[right][k][2] = 0.5 * (xi[right][k][0] + xi[right][k][1]);
//       }
//     }
//   }
//   else if(split == Nsplit) {
//     for(unsigned k = 0; k < dim; k++) {
//       xi[left][k][0] = xi[left][k][1];
//       xi[left][k][1] = xig[k];
//       xi[left][k][2] = 0.5 * (xi[left][k][0] + xi[left][k][1]);
//
//       xi[right][k][1] = xi[right][k][0];
//       xi[right][k][0] = xig[k];
//       xi[right][k][2] = 0.5 * (xi[right][k][0] + xi[right][k][1]);
//     }
//   }
//   else {
//     for(unsigned k = 0; k < dim; k++) {
//       xi[left][k][0] = xi[left][k][1];
//       xi[left][k][1] = 0.5 * (xi[left][k][0] + xig[k]);
//       xi[left][k][2] = 0.5 * (xi[left][k][0] + xi[left][k][1]);
//
//       xi[right][k][1] = xi[right][k][0];
//       xi[right][k][0] = 0.5 * (xi[right][k][1] + xig[k]);
//       xi[right][k][2] = 0.5 * (xi[right][k][0] + xi[right][k][1]);
//     }
//   }
//
//   x.resize(2);
//   for(unsigned r = 0; r < 2; r++){
//     x[r].resize(dim);
//     for(unsigned k = 0; k < dim; k++) {
//       x[r][k].assign(3,0.);
//     }
//   }
//
//    for(unsigned r = 0; r < 2; r++){
//      for(unsigned i = 0; i < 3; i++){
//        msh->_finiteElement[ielGeom1][solType]->GetPhi(phi3, xi[r]);
//
//
//     }
//    }
//
//
//
//
//
// }
//














void GetElementPartition2D(const std::vector <double >  & xg1, const std::vector < std::vector <double > > & x1, const unsigned &split,  std::vector < std::vector < std::vector<double>>> &x) {
  unsigned dim = 2;
  unsigned bl = 0; // bottom left
  unsigned br = 1; // bottom right
  unsigned tr = 2; // top left
  unsigned tl = 3; // top right
  unsigned bl1 = 4; // bottom left
  unsigned bl2 = 5; // bottom left
  unsigned br1 = 6; // bottom right
  unsigned br2 = 7; // bottom right
  unsigned tr1 = 8; // top left
  unsigned tr2 = 9; // top left
  unsigned tl1 = 10; // top right
  unsigned tl2 = 11; // top right
  
  double ex_x_1;
  double ex_x_2;
  double ex_y_1;
  double ex_y_2;
  
  unsigned size_part = (split != Nsplit) ? 12 : 4;
  
  if(split == 0) { //init
    x.resize(size_part);
    for(unsigned j = 0; j < size_part; j++) {
      x[j].resize(dim);
      for(unsigned k = 0; k < dim; k++) {
        x[j][k].resize(x1[0].size());
      }
    }
    ex_x_1 = x1[0][0];
    ex_x_2 = x1[0][1];
    ex_y_1 = x1[1][0];
    ex_y_2 = x1[1][3];
  }
  else {
    ex_x_1 = x[bl][0][1];
    ex_x_2 = x[br][0][0];
    ex_y_1 = x[bl][1][2];
    ex_y_2 = x[tl][1][1];
  }
  
  
  //     Prototipo: x[quadrante][dim][numero_nodo]
  x[bl][1][0] = x[bl][1][1] = x[br][1][0] = x[br][1][1] = ex_y_1;
  x[tl][1][2] = x[tl][1][3] = x[tr][1][2] = x[tr][1][3] = ex_y_2;
  x[bl][0][0] = x[bl][0][3] = x[tl][0][0] = x[tl][0][3] = ex_x_1;
  x[br][0][1] = x[br][0][2] = x[tr][0][1] = x[tr][0][2] = ex_x_2;
  
  
  if(split == Nsplit) {
    x[bl][1][2] = x[bl][1][3] = x[br][1][2] = x[br][1][3] = x[tl][1][0] = x[tl][1][1] = x[tr][1][0] = x[tr][1][1] = xg1[1];
    x[bl][0][1] = x[bl][0][2] = x[tl][0][1] = x[tl][0][2] = x[br][0][0] = x[br][0][3] = x[tr][0][0] = x[tr][0][3] = xg1[0];
  }
  else {
    x[bl][1][2] = x[bl][1][3] = x[br][1][2] = x[br][1][3] = 0.5 * (ex_y_1 + xg1[1]);
    x[tl][1][0] = x[tl][1][1] = x[tr][1][0] = x[tr][1][1] = 0.5 * (ex_y_2 + xg1[1]);
    x[bl][0][1] = x[bl][0][2] = x[tl][0][1] = x[tl][0][2] = 0.5 * (ex_x_1 + xg1[0]);
    x[br][0][0] = x[br][0][3] = x[tr][0][0] = x[tr][0][3] = 0.5 * (ex_x_2 + xg1[0]);
  }
  
  if(split != Nsplit) {
    x[bl1][1][0] = x[bl1][1][1] = x[br1][1][0] = x[br1][1][1] = ex_y_1;
    x[tl1][1][2] = x[tl1][1][3] = x[tr1][1][2] = x[tr1][1][3] = ex_y_2;
    x[bl1][1][2] = x[bl1][1][3] = x[br1][1][2] = x[br1][1][3] = 0.5 * (ex_y_1 + xg1[1]);
    x[bl2][1][0] = x[bl2][1][1] = x[br2][1][0] = x[br2][1][1] = 0.5 * (ex_y_1 + xg1[1]);
    x[tl1][1][0] = x[tl1][1][1] = x[tr1][1][0] = x[tr1][1][1] = 0.5 * (ex_y_2 + xg1[1]);
    x[tl2][1][2] = x[tl2][1][3] = x[tr2][1][2] = x[tr2][1][3] = 0.5 * (ex_y_2 + xg1[1]);
    x[bl2][1][2] = x[bl2][1][3] = x[br2][1][2] = x[br2][1][3] = xg1[1];
    x[tl2][1][0] = x[tl2][1][1] = x[tr2][1][0] = x[tr2][1][1] = xg1[1];
    x[bl2][0][0] = x[bl2][0][3] = x[tl2][0][0] = x[tl2][0][3] = ex_x_1;
    x[br2][0][1] = x[br2][0][2] = x[tr2][0][1] = x[tr2][0][2] = ex_x_2;
    x[bl2][0][1] = x[bl2][0][2] = x[tl2][0][1] = x[tl2][0][2] = 0.5 * (ex_x_1 + xg1[0]);
    x[bl1][0][0] = x[bl1][0][3] = x[tl1][0][0] = x[tl1][0][3] = 0.5 * (ex_x_1 + xg1[0]);
    x[br2][0][0] = x[br2][0][3] = x[tr2][0][0] = x[tr2][0][3] = 0.5 * (ex_x_2 + xg1[0]);
    x[br1][0][1] = x[br1][0][2] = x[tr1][0][1] = x[tr1][0][2] = 0.5 * (ex_x_2 + xg1[0]);
    x[bl1][0][1] = x[bl1][0][2] = x[tl1][0][1] = x[tl1][0][2] = xg1[0];
    x[br1][0][0] = x[br1][0][3] = x[tr1][0][0] = x[tr1][0][3] = xg1[0];
  }
  
  
  for(unsigned qq = 0; qq < size_part; qq++) {
    for(unsigned k = 0; k < dim; k++) { //middle point formula
      x[qq][k][4] = 0.5 * (x[qq][k][0] + x[qq][k][1]);
      x[qq][k][5] = 0.5 * (x[qq][k][1] + x[qq][k][2]);
      x[qq][k][6] = 0.5 * (x[qq][k][2] + x[qq][k][3]);
      x[qq][k][7] = 0.5 * (x[qq][k][3] + x[qq][k][0]);
      x[qq][k][8] = 0.5 * (x[qq][k][0] + x[qq][k][2]);
    }
  }
  
}


const unsigned ijndex[2][12][2] = {
  { {0, 0}, {3, 0}, {3, 3}, {0, 3},
  {1, 0}, {0, 1},
  {2, 0}, {3, 1},
  {2, 3}, {3, 2},
  {1, 3}, {0, 2}
  },
  {{0, 0}, {1, 0}, {1, 1}, {0, 1}}
};

void GetElementPartitionQuad(const std::vector <double >  & xg1, const std::vector < std::vector <double > > & xNodes, const unsigned & split, const unsigned & totalNumberofSplits,  std::vector < std::vector < std::vector<double>>> &x) {
  unsigned dim = 2;
  
  unsigned solType;
  unsigned size = xNodes[0].size();
  
  if(size == 4) {
    solType = 0; //lagrange linear
  }
  else if(size == 8) {
    solType = 1; //lagrange serendipity
  }
  else if(size == 9) {
    solType = 2; //lagrange quadratic
  }
  else {
    std::cout << "abort in GetElementPartitionQuad" << std::endl;
    abort();
  }
  
  
  unsigned bl = 0; // bottom left
  unsigned br = 1; // bottom right
  unsigned tr = 2; // top left
  unsigned tl = 3; // top right
  
  std::vector < double > XX;
  std::vector < double > YY;
  
  unsigned size_part = 12;
  unsigned splitType = 0;
  
  if(split < totalNumberofSplits) { //init && update
    
    XX.resize(5);
    YY.resize(5);
    
    if(split == 0) { //init
      
      x.resize(size_part);
      for(unsigned j = 0; j < size_part; j++) {
        x[j].resize(dim);
        for(unsigned k = 0; k < dim; k++) {
          x[j][k].resize(size);
        }
      }
      
      XX[0] = xNodes[0][0];
      XX[4] = xNodes[0][1];
      YY[0] = xNodes[1][0];
      YY[4] = xNodes[1][3];
      
    }
    else { //update
      XX[0] = x[bl][0][1];
      XX[4] = x[br][0][0];
      YY[0] = x[bl][1][2];
      YY[4] = x[tl][1][1];
    }
    XX[2] = xg1[0];
    XX[1] = 0.5 * (XX[0] + XX[2]);
    XX[3] = 0.5 * (XX[2] + XX[4]);
    
    YY[2] = xg1[1];
    YY[1] = 0.5 * (YY[0] + YY[2]);
    YY[3] = 0.5 * (YY[2] + YY[4]);
  }
  else { //close
    
    XX.resize(3);
    YY.resize(3);
    
    XX[0] = x[bl][0][1];
    XX[1] = xg1[0];
    XX[2] = x[br][0][0];
    YY[0] = x[bl][1][2];
    YY[1] = xg1[1];
    YY[2] = x[tl][1][1];
    
    size_part = 4;
    splitType = 1;
    x.resize(size_part);
    for(unsigned j = 0; j < size_part; j++) {
      x[j].resize(dim);
      for(unsigned k = 0; k < dim; k++) {
        x[j][k].resize(size);
      }
    }
  }
  
  for(unsigned qq = 0; qq < size_part; qq++) {
    unsigned i = ijndex[splitType][qq][0];
    x[qq][0][0] = x[qq][0][3] = XX[i];
    x[qq][0][1] = x[qq][0][2] = XX[i + 1];
    
    unsigned j = ijndex[splitType][qq][1];
    x[qq][1][0] = x[qq][1][1] = YY[j];
    x[qq][1][2] = x[qq][1][3] = YY[j + 1];
  }
  if(solType > 0) {
    for(unsigned qq = 0; qq < size_part; qq++) {
      for(unsigned k = 0; k < dim; k++) { //middle point formula
        x[qq][k][4] = 0.5 * (x[qq][k][0] + x[qq][k][1]);
        x[qq][k][5] = 0.5 * (x[qq][k][1] + x[qq][k][2]);
        x[qq][k][6] = 0.5 * (x[qq][k][2] + x[qq][k][3]);
        x[qq][k][7] = 0.5 * (x[qq][k][3] + x[qq][k][0]);
      }
    }
  }
  
  if(solType > 1) {
    for(unsigned qq = 0; qq < size_part; qq++) {
      for(unsigned k = 0; k < dim; k++) { //middle point formula
        x[qq][k][8] = 0.5 * (x[qq][k][0] + x[qq][k][2]);
      }
    }
  }
  
}
