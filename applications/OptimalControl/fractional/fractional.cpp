
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

#include "Assemble_jacobian.hpp"


using namespace femus;

#define N_UNIFORM_LEVELS  4
#define N_ERASED_LEVELS   3
#define S_FRAC 0.25

#define OP_L2       0
#define OP_H1       0
#define OP_Hhalf    1
#define RHS_ONE     1

#define UNBOUNDED   1

#define USE_Cns     1


#define EX_1       -1.
#define EX_2        1.
#define EY_1       -1.
#define EY_2        1.


#include "../fractional_functions.hpp"

#define Nsplit      10

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



void GetHsNorm(const unsigned level, MultiLevelProblem& ml_prob);

void AssembleFracProblem(MultiLevelProblem& ml_prob);


int main(int argc, char** argv)
{


  //quadr rule order
//   const std::string fe_quad_rule_1 = "fifth";
//   const std::string fe_quad_rule_2 = "sixth";
  const std::string fe_quad_rule_1 = "seventh";
  const std::string fe_quad_rule_2 = "eighth";



  // ======= Init ========================
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);

  // ======= Files ========================
  Files files; 
        files.CheckIODirectories();
        files.RedirectCout();

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

  ml_prob.SetFilesHandler(&files);
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
  
  unsigned dimension = pow ( pow(2, numberOfUniformLevels) * 2 + 1, dim );
  system.SetSparsityPatternMinimumSize (dimension, "u");
    
  system.init();

  //dense =============
  //dense =============
  //dense =============
  
  
  
//  const unsigned solType = mlSol.GetSolutionType("u");
//   for(int level = 0; level < mlMsh.GetNumberOfLevels(); level++) {
// 
//     Mesh* msh = mlMsh.GetLevel(level);
//     unsigned nprocs = msh->n_processors();
//     unsigned iproc = msh->processor_id();
// 
//     int KK_size = msh->_dofOffset[solType][nprocs];
//     int KK_local_size = msh->_dofOffset[solType][iproc + 1] - msh->_dofOffset[solType][iproc];
// 
// //   SparseMatrix* CC;
// //   CC = SparseMatrix::build().release();
//     system._LinSolver[level]->_KK->init(KK_size, KK_size, KK_local_size, KK_local_size, KK_local_size, KK_size - KK_local_size);
//     system._LinSolver[level]->_KK->zero();
//   }
  //dense =============
  //dense =============
  //dense =============


  // ******* Set Smoother *******
  system.SetSolverFineGrids(GMRES);

  system.SetPreconditionerFineGrids(ILU_PRECOND);

  system.SetTolerances(1.e-20, 1.e-20, 1.e+50, 100);

  system.MGsolve();
//   system.assemble_call(1);  //to only call the assemble function

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
  SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
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
  unsigned soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned solType   = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

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
  vector < double > KK_local;  KK_local.reserve(maxSize * maxSize);
  vector < double > Res_local; Res_local.reserve(maxSize);

//   Local matrices and rhs for adaptive quadrature
  vector < double > Res_local_refined; Res_local_refined.reserve(maxSize);
  vector < double > CClocal_refined;   CClocal_refined.reserve(maxSize * maxSize);

  vector < double > KK_mixed;   KK_mixed.reserve(maxSize * maxSize);
  vector < double > Res_mixed;  Res_mixed.reserve(maxSize);

//   Non local matrices and vectors for H^s laplacian operator
//   vector< double >         Res_nonlocal;
//   Res_nonlocal.reserve(maxSize);  // local residual vector
  vector< double >         Res_nonlocalI;  Res_nonlocalI.reserve(maxSize);
  vector< double >         Res_nonlocalJ;  Res_nonlocalJ.reserve(maxSize);
//   vector < double > CClocal;
//   CClocal.reserve(maxSize * maxSize);

  vector < double > CC_nonlocal_II;  CC_nonlocal_II.reserve(maxSize * maxSize);
  vector < double > CC_nonlocal_IJ;  CC_nonlocal_IJ.reserve(maxSize * maxSize);
  vector < double > CC_nonlocal_JI;  CC_nonlocal_JI.reserve(maxSize * maxSize);
  vector < double > CC_nonlocal_JJ;  CC_nonlocal_JJ.reserve(maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix
  RES->zero();

//   int KK_size = msh->_dofOffset[solType][nprocs];
//   int KK_local_size = msh->_dofOffset[solType][iproc + 1] - msh->_dofOffset[solType][iproc];
//
//   SparseMatrix* CC;
//   CC = SparseMatrix::build().release();
//   CC->init(KK_size, KK_size, KK_local_size, KK_local_size, KK_local_size, KK_size - KK_local_size);
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
      //unsigned n_face;

      if(iproc == kproc) {
        ielGeom2 = msh->GetElementType(jel);
        nDof2  = msh->GetElementDofNumber(jel, solType);    // number of solution element dofs
        nDofx2 = msh->GetElementDofNumber(jel, xType);    // number of coordinate element dofs

      }

      MPI_Bcast(&ielGeom2, 1, MPI_UNSIGNED_SHORT, kproc, MPI_COMM_WORLD);
      MPI_Bcast(&nDof2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
      MPI_Bcast(&nDofx2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);
      //MPI_Bcast(&n_face, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);

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

      std::vector <int> bd_face(0);
      unsigned nFaces;
      //std::vector <unsigned> faceDofs(n_face, 0);
      //vector < vector <unsigned> > inode(n_face);
      if(iproc == kproc) {
        for(unsigned jface = 0; jface < msh->GetElementFaceNumber(jel); jface++) {
          int faceIndex = el->GetBoundaryIndex(jel, jface);
//           faceDofs[jface] = msh->GetElementFaceDofNumber(jel, jface, solType);
//           inode[jface].resize(faceDofs[jface]);
// //       inode[jface].assign(faceDofs[jface], 0);
//           for(unsigned i = 0; i < faceDofs[jface]; i++) {
//             inode[jface][i] = msh->GetLocalFaceVertexIndex(jel, jface, i);    // face-to-element local node mapping.
//           }
//           MPI_Bcast(& inode[jface][0], faceDofs[jface], MPI_UNSIGNED, kproc, MPI_COMM_WORLD);

          // look for boundary faces
          if(faceIndex >= 1) {
            unsigned i = bd_face.size();
            bd_face.resize(i + 1);
            bd_face[i] = jface;
          }
        }
        nFaces = bd_face.size();
      }

      MPI_Bcast(& nFaces, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD);

      bd_face.resize(nFaces);
      MPI_Bcast(& bd_face[0], nFaces, MPI_INT, kproc, MPI_COMM_WORLD);




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
        CC_nonlocal_II.assign(nDof1 * nDof2, 0.);   //resize
        CC_nonlocal_IJ.assign(nDof1 * nDof2, 0.);   //resize
        CC_nonlocal_JI.assign(nDof1 * nDof2, 0.);   //resize
        CC_nonlocal_JJ.assign(nDof1 * nDof2, 0.);   //resize
//         Res_nonlocal.assign(nDof1, 0);    //resize
        Res_nonlocalI.assign(nDof1, 0);    //resize
        Res_nonlocalJ.assign(nDof1, 0);    //resize

        Res_mixed.assign(nDof1, 0);    //resize
        KK_mixed.assign(nDof1 * nDof1, 0.);

        if(iel == jel) {
          Res_local.assign(nDof1, 0);    //resize
          KK_local.assign(nDof1 * nDof1, 0.);
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
                KK_local[ i * nDof1 + j ] += OP_L2 * phi1[i] * phi1[j] * weight1;
              }
              double mass_res_i = phi1[i] * solX ;
              Res_local[ i ] += OP_L2 * weight1 * mass_res_i ;
              Res_local[ i ] += - RHS_ONE * weight1 * (phi1[i] * (-1.) /** ( sin(2 * acos(0.0) * x1[0][i])) * ( sin(2 * acos(0.0) * x1[1][i]))*/);
            }

//============  Laplacian assembly ==================

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
                KK_local[ i * nDof1 + j ]  += OP_H1 * weight1 *  laplace_mat_i_j;
              }
            }
//============  Laplacian assembly ==================

//============  Mixed integral - Analytical ((Rn-Omega) x Omega) assembly (based on the analytic result of integrals) ==================
//             if(dim == 1 && UNBOUNDED == 1) {
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
//                   KK_local[ i * nDof1 + j ] += (C_ns / 2.) * check_limits * (1. / s_frac) * OP_Hhalf * phi1[i] * phi1[j] * weight1 * mixed_term;
//                 }
//                 Res_local[ i ] += (C_ns / 2.) * check_limits * (1. / s_frac) * OP_Hhalf * weight1 * phi1[i] * solX * mixed_term;
//               }
//             }
//             if(dim == 2 && UNBOUNDED == 1) {
//               double ex[4] = {EX_1 - xg1[0], EX_2 - xg1[0], EY_1 - xg1[1], EY_2 - xg1[1]};
// 
//               double teta[4], CCC[4];
//               teta[0] = atan2(ex[3], ex[1]);
//               teta[1] = atan2(ex[3], ex[0]);
//               teta[2] = atan2(ex[2], ex[0]) + 2 * M_PI;
//               teta[3] = atan2(ex[2], ex[1]) + 2 * M_PI;
// 
//               double mixed_term = 0.;
// 
//               for(unsigned qq = 0; qq < 4; qq++) {
//                 if(qq == 3) teta[3] -= 2. * M_PI ;
// 
//                 if(qq == 0)
//                   mixed_term += 2.* (Antiderivative1(teta[1], s_frac, ex[3]) -
//                                      Antiderivative1(teta[0], s_frac, ex[3]));
//                 else if(qq  == 2)
//                   mixed_term += 2.* (Antiderivative1(teta[3], s_frac, ex[2]) -
//                                      Antiderivative1(teta[2], s_frac, ex[2]));
//                 else if(qq  == 1)
//                   mixed_term += 2. * (Antiderivative2(teta[2], s_frac, ex[0]) -
//                                       Antiderivative2(teta[1], s_frac, ex[0]));
//                 else
//                   mixed_term += 2. * (Antiderivative2(teta[0], s_frac, ex[1]) -
//                                       Antiderivative2(teta[3], s_frac, ex[1]));
// 
// 
//               }

// //               if(iel == 0 && ig == 4)   sum_int += mixed_term;
// //               std::cout<<"sum_int = " << sum_int << "\n";


//               for(unsigned i = 0; i < nDof1; i++) {
//                 for(unsigned j = 0; j < nDof1; j++) {
//                   KK_local[ i * nDof1 + j ] += (C_ns / 2.) * check_limits * OP_Hhalf * phi1[i] * phi1[j] * weight1 * mixed_term;
//                 }
//                 Res_local[ i ] += (C_ns / 2.) * check_limits * OP_Hhalf * weight1 * phi1[i] * solX * mixed_term;
//               }
//             }
//============  Mixed integral - Analytical ((Rn-Omega) x Omega) assembly (based on the analytic result of integrals) ==================

//============ Adaptive quadrature for iel == jel ==================
            if(Nsplit != 0 && OP_Hhalf != 0) {

              std::cout.precision(14);
              std::vector< std::vector<std::vector<double>>> x3;

              for(unsigned split = 0; split <= Nsplit; split++) {

//                 unsigned size_part;
//                 if(dim == 1) size_part = 2;
//                 else size_part = (split != Nsplit) ? 12 : 4;

                if(dim == 1) GetElementPartition1D(xg1, x1, split, Nsplit, x3);
                else if(dim == 2) {
                  //GetElementPartition2D(xg1, x1, split, Nsplit, x3);
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
//============ Adaptive quadrature for iel == jel ==================

              if(UNBOUNDED == 1) {
//============ Mixed integral 1D - Analytical ==================
                if(ig == 0 && dim == 1) {
                  double ex_1 = EX_1;
                  double ex_2 = EX_2;
                  double dist_1 = 0.;
                  double dist_2 = 0.;
                  for(int k = 0; k < dim; k++) {
                    dist_1 += sqrt((xg3[k] - ex_1) * (xg3[k] - ex_1));
                    dist_2 += sqrt((xg3[k] - ex_2) * (xg3[k] - ex_2));
                  }
                  double mixed_term = (pow(dist_1, -2. * s_frac) + pow(dist_2, - 2. * s_frac)) * (1. / s_frac) ;

                  for(unsigned i = 0; i < nDof1; i++) {
                    for(unsigned j = 0; j < nDof1; j++) {
                      KK_local[ i * nDof1 + j ] += (C_ns / 2.) * check_limits * OP_Hhalf * phi3[i] * phi3[j] * weight3 * mixed_term;
                    }
                    Res_local[ i ] += (C_ns / 2.) * check_limits * OP_Hhalf * weight3 * phi3[i] * solY3 * mixed_term;
                  }
                }
//============ Mixed integral 1D - Analytical ==================
//============ Mixed Integral 2D - Numerical ==================
                if(ig == 0 && dim == 2) {
                  double mixed_term1 = 0;
//     for(int kel = msh->_elementOffset[iproc]; kel < msh->_elementOffset[iproc + 1]; kel++) {
            // *** Face Gauss point loop (boundary Integral) ***
                  for(unsigned jj = 0; jj < bd_face.size(); jj++) {

                    int jface = bd_face[jj];
                    
                    // look for boundary faces
                    unsigned faceDofs = el->GetNFACENODES(ielGeom2, jface, solType);


                    vector  < vector  <  double> > faceCoordinates(dim);    // A matrix holding the face coordinates rowwise.
                    for(int k = 0; k < dim; k++) {
                      faceCoordinates[k].resize(faceDofs);
                    }
                    for(unsigned i = 0; i < faceDofs; i++) {
                      unsigned inode = el->GetIG(ielGeom2, jface, i);  // face-to-element local node mapping.
                      for(unsigned k = 0; k < dim; k++) {
                        faceCoordinates[k][i] =  x2[k][inode] - xg3[k]; // We extract the local coordinates on the face from local coordinates on the element.
                      }
                    }
                    const unsigned div = 10;
                    vector  < vector  <  double> > interpCoordinates(dim);
                    for(int k = 0; k < dim; k++) {
                      interpCoordinates[k].resize(div + 1); // set "4" as a parameter
                    }
                    for(unsigned n = 0; n <= div; n++) {
                      for(int k = 0; k < dim; k++) {
                        interpCoordinates[k][n] = faceCoordinates[k][0] + n * (faceCoordinates[k][1] - faceCoordinates[k][0]) /  div ;
                      }
                    }
                    for(unsigned n = 0; n < div; n++) {
                      double teta2 = atan2(interpCoordinates[1][n + 1], interpCoordinates[0][n + 1]);
                      double teta1 = atan2(interpCoordinates[1][n], interpCoordinates[0][n]);

                      if(teta2 < teta1) teta2 += 2. * M_PI;

                      double delta_teta = teta2 - teta1;


                      vector <double> mid_point;
                      mid_point.resize(dim);
                      for(unsigned k = 0; k < dim; k++) {
                        mid_point[k] = (interpCoordinates[k][n + 1] + interpCoordinates[k][n]) * 0.5;
                      }
                      double dist2 = 0;
                      for(int k = 0; k < dim; k++) {
                        dist2 += mid_point[k] * mid_point[k];
                      }
                      double dist = sqrt(dist2);
                      mixed_term1 += 2. * pow(dist, -  2. * s_frac) * (1. / (2. * s_frac)) * delta_teta;
                    }
                  }

                  for(unsigned i = 0; i < nDof1; i++) {
                    for(unsigned j = 0; j < nDof1; j++) {
                      KK_mixed[ i * nDof1 + j ] += (C_ns / 2.) * check_limits * OP_Hhalf * phi3[i] * phi3[j] * weight3 * mixed_term1;
                    }
                    Res_mixed[ i ] += (C_ns / 2.) * check_limits * OP_Hhalf * weight3 * phi3[i] * solX * mixed_term1;
                  }
                }
              }
//============ Mixed Integral 2D - Numerical ==================
            } //end jg
          } //end r
        }  //end split
      }  //end if Nsplit != 0
    } // end iel == jel loop

      if(OP_Hhalf != 0) {
          if(iel != jel || Nsplit == 0) {
            
    //============  Mixed integral 1D  ==================
            if(dim == 1 && UNBOUNDED == 1) {
              double ex_1 = EX_1;
              double ex_2 = EX_2;
              double dist_1 = 0.;
              double dist_2 = 0.;
              for(int k = 0; k < dim; k++) {
                dist_1 += sqrt((xg1[k] - ex_1) * (xg1[k] - ex_1));
                dist_2 += sqrt((xg1[k] - ex_2) * (xg1[k] - ex_2));
              }
              double mixed_term = pow(dist_1, -2. * s_frac) + pow(dist_2, - 2. * s_frac);

              for(unsigned i = 0; i < nDof1; i++) {
                for(unsigned j = 0; j < nDof1; j++) {
                  KK_local[ i * nDof1 + j ] += (C_ns / 2.) * check_limits * (1. / s_frac) * OP_Hhalf * phi1[i] * phi1[j] * weight1 * mixed_term;
                }
                Res_local[ i ] += (C_ns / 2.) * check_limits * (1. / s_frac) * OP_Hhalf * weight1 * phi1[i] * solX * mixed_term;
              }
            }
    //============  Mixed integral 1D  ==================        
    
    //============ Mixed Integral - Numerical ==================      
            if( dim == 2 && UNBOUNDED == 1 ) {
            double mixed_term1 = 0;
//     for(int kel = msh->_elementOffset[iproc]; kel < msh->_elementOffset[iproc + 1]; kel++) {
            // *** Face Gauss point loop (boundary Integral) ***
            for(unsigned jj = 0; jj < bd_face.size(); jj++) {

              int jface = bd_face[jj];
              // look for boundary faces

              unsigned faceDofs = el->GetNFACENODES(ielGeom2, jface, solType);

              vector  < vector  <  double> > faceCoordinates(dim);    // A matrix holding the face coordinates rowwise.
              for(int k = 0; k < dim; k++) {
                faceCoordinates[k].resize(faceDofs);
              }
              for(unsigned i = 0; i < faceDofs; i++) {
                unsigned inode = el->GetIG(ielGeom2, jface, i);  // face-to-element local node mapping.
                for(unsigned k = 0; k < dim; k++) {
                  faceCoordinates[k][i] =  x2[k][inode] - xg1[k]; // We extract the local coordinates on the face from local coordinates on the element.
                }
              }
              const unsigned div = 10;
              vector  < vector  <  double> > interpCoordinates(dim);
              for(int k = 0; k < dim; k++) {
                interpCoordinates[k].resize(div + 1); // set "4" as a parameter
              }
              for(unsigned n = 0; n <= div; n++) {
                for(int k = 0; k < dim; k++) {
                  interpCoordinates[k][n] = faceCoordinates[k][0] + n * (faceCoordinates[k][1] - faceCoordinates[k][0]) /  div ;
                }
              }
              for(unsigned n = 0; n < div; n++) {
                double teta2 = atan2(interpCoordinates[1][n + 1], interpCoordinates[0][n + 1]);
                double teta1 = atan2(interpCoordinates[1][n], interpCoordinates[0][n]);

                if(teta2 < teta1) teta2 += 2. * M_PI;

                double delta_teta = teta2 - teta1;


                vector <double> mid_point;
                mid_point.resize(dim);
                for(unsigned k = 0; k < dim; k++) {
                  mid_point[k] = (interpCoordinates[k][n + 1] + interpCoordinates[k][n]) * 0.5;
                }
                double dist2 = 0;
                for(int k = 0; k < dim; k++) {
                  dist2 += mid_point[k] * mid_point[k];
                }
                double dist = sqrt(dist2);
                mixed_term1 += 2. * pow(dist, -  2. * s_frac) * (1. / (2. * s_frac)) * delta_teta;
              }
            }

            for(unsigned i = 0; i < nDof1; i++) {
              for(unsigned j = 0; j < nDof1; j++) {
                KK_mixed[ i * nDof1 + j ] += (C_ns / 2.) * check_limits * OP_Hhalf * phi1[i] * phi1[j] * weight1 * mixed_term1;
              }
              Res_mixed[ i ] += (C_ns / 2.) * check_limits * OP_Hhalf * weight1 * phi1[i] * solX * mixed_term1;
            }
           }
          
//============ Mixed Integral - Numerical ==================

             
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

                  CC_nonlocal_II[ i * nDof2 + j ] += (C_ns / 2.) * OP_Hhalf * check_limits * phi1[j]  * phi1[i] * weight1 * weight2[jg] / denom;

                  CC_nonlocal_IJ[ i * nDof2 + j ] += (C_ns / 2.) * OP_Hhalf * check_limits * (- phi2[jg][j]) * phi1[i] * weight1 * weight2[jg] / denom;

                  CC_nonlocal_JI[ i * nDof2 + j ] += (C_ns / 2.) * OP_Hhalf * check_limits * (phi1[j]) * (- phi2[jg][i]) * weight1 * weight2[jg] / denom;

                  CC_nonlocal_JJ[ i * nDof2 + j ] += (C_ns / 2.) * OP_Hhalf * check_limits * (- phi2[jg][j]) * (- phi2[jg][i]) * weight1 * weight2[jg] / denom;


                  }
                }
              } //endl jg loop
            } //end if(iel != jel || Nsplit == 0)
          } 
          
        } //endl ig loop

        if(iel == jel) {
          KK->add_matrix_blocked(KK_local, l2GMap1, l2GMap1);
          RES->add_vector_blocked(Res_local, l2GMap1);

          if(Nsplit != 0) {
            KK->add_matrix_blocked(CClocal_refined, l2GMap1, l2GMap1);
            RES->add_vector_blocked(Res_local_refined, l2GMap1);
          }
        }
//        KK->add_matrix_blocked(CClocal, l2GMap1, l2GMap2);
        KK->add_matrix_blocked(KK_mixed, l2GMap1, l2GMap1);
        RES->add_vector_blocked(Res_mixed, l2GMap1);

        KK->add_matrix_blocked(CC_nonlocal_II, l2GMap1, l2GMap1);
        KK->add_matrix_blocked(CC_nonlocal_IJ, l2GMap1, l2GMap2);
        KK->add_matrix_blocked(CC_nonlocal_JI, l2GMap2, l2GMap1);
        KK->add_matrix_blocked(CC_nonlocal_JJ, l2GMap2, l2GMap2);

//        RES->add_vector_blocked(Res_nonlocal, l2GMap1);
        RES->add_vector_blocked(Res_nonlocalI, l2GMap1);
        RES->add_vector_blocked(Res_nonlocalJ, l2GMap2);
        
      } // end iel loop


    } //end jel loop
  } //end kproc loop

  KK->close();
  RES->close();

//--- print matrix on file
  //this shows that the matrix is full, because each row has all columns (see the output file)
const unsigned nonlin_iter = 0/*mlPdeSys->GetNonlinearIt()*/;
    assemble_jacobian< double, double >::print_global_jacobian(/*assemble_matrix*/true, ml_prob, KK, nonlin_iter);
//     assemble_jacobian< double, double >::print_global_residual(ml_prob, RES, nonlin_iter);
//--- print matrix on file


//   PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 900, 900, &viewer);
//   PetscObjectSetName((PetscObject)viewer, "FSI matrix");
//   PetscViewerPushFormat(viewer, PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*>(KK))->mat(), viewer);
// //   MatView((static_cast<PetscMatrix*> (KK))->mat(),  PETSC_VIEWER_STDOUT_WORLD );
//   double a;
//   std::cin >> a;


//     PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,900,900,&viewer);
//   PetscObjectSetName((PetscObject)viewer,"FSI matrix");
//   PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*> (KK))->mat(),viewer);
// //   MatView((static_cast<PetscMatrix*> (KK))->mat(),  PETSC_VIEWER_STDOUT_WORLD );
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






