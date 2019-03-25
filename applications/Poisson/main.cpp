
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelMesh.hpp"
#include "TransientSystem.hpp"
#include "NumericVector.hpp"
#include "Fluid.hpp"
#include "Parameter.hpp"
#include "FemusInit.hpp"
#include "SparseMatrix.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "XDMFWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "SolvertypeEnum.hpp"
#include "FElemTypeEnum.hpp"
#include "ParsedFunction.hpp"
#include "InputParser.hpp"
#include "Files.hpp"


using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int facename, const double time) {
  bool test = 1; //dirichlet
  value = 0.;

  if (3 == facename) { //stress
    test = 0;
    value = 0.2;
  }

  return test;
}

void AssemblePoissonMatrixandRhs(MultiLevelProblem& ml_prob);

void show_usage()
{
  std::cout << "Use --inputfile variable to set the input file" << std::endl;
  std::cout << "e.g.: ./Poisson --inputfile ./input/input.json" << std::endl;
}

ParsedFunction fpsource;


int main(int argc, char** argv) {

  std::string path;

  if (argc < 2)
  {
    std::cout << argv[0] << ": You must specify the input file" << std::endl;
    show_usage();
    return 1;
  }

  for (int count = 1; count < argc; ++count)
  {
    std::string arg = argv[count];

    if ((arg == "-h") || (arg == "--help")) {
      show_usage();
      return 0;
    }
    else if ((arg == "-i") || (arg == "--inputfile"))
    {
      if (count + 1 < argc) {
        path = argv[++count];
      }
      else
      {
        std::cerr << "--input file option requires one argument." << std::endl;
        return 1;
      }
    }

    //         else {
    // 	  std::cerr << argv[count] << " : command line argument not recognized" << std::endl;
    // 	  show_usage();
    // 	  return 1;
    // 	}
  }

  /// Init Petsc-MPI communicator
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);

  //Files files;
  //files.CheckIODirectories();
  //files.RedirectCout();

  // input parser pointer
  std::unique_ptr<InputParser> inputparser = InputParser::build(path);

  /// INIT MESH =================================

  unsigned short nm, nr;
  unsigned int nlevels = inputparser->getValue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.type.multigrid.nlevels", 1);
  nm = nlevels;

  nr = 0;

  int tmp = nm;
  nm += nr;
  nr = tmp;

  //Adimensional quantity (Lref,Uref)
  double Lref = 1.;
  double Uref = 1.;

  //Steadystate NonLinearMultiLevelProblem
  MultiLevelMesh ml_msh;

  if (inputparser->isTrue("multilevel_mesh.first.type", "filename"))
  {
    std::string filename = inputparser->getValue("multilevel_mesh.first.type.filename", "./input/input.neu");
    ml_msh.ReadCoarseMesh(filename.c_str(), "seventh", Lref);
  }
  else if (inputparser->isTrue("multilevel_mesh.first.type", "box"))
  {
    int numelemx = inputparser->getValue("multilevel_mesh.first.type.box.nx", 2);
    int numelemy = inputparser->getValue("multilevel_mesh.first.type.box.ny", 2);
    int numelemz = inputparser->getValue("multilevel_mesh.first.type.box.nz", 0);
    double xa = inputparser->getValue("multilevel_mesh.first.type.box.xa", 0.);
    double xb = inputparser->getValue("multilevel_mesh.first.type.box.xb", 1.);
    double ya = inputparser->getValue("multilevel_mesh.first.type.box.ya", 0.);
    double yb = inputparser->getValue("multilevel_mesh.first.type.box.yb", 1.);
    double za = inputparser->getValue("multilevel_mesh.first.type.box.za", 0.);
    double zb = inputparser->getValue("multilevel_mesh.first.type.box.zb", 0.);
    ElemType elemtype = inputparser->getValue("multilevel_mesh.first.type.box.elem_type", QUAD9);
    ml_msh.GenerateCoarseBoxMesh(numelemx, numelemy, numelemz, xa, xb, ya, yb, za, zb, elemtype, "seventh");
  }
  else
  {
    std::cerr << "Error: no input mesh specified. Please check to have added the keyword mesh in the input json file! " << std::endl;
    return 1;
  }

  ml_msh.RefineMesh(nm, nr, NULL);

  ml_msh.PrintInfo();

  MultiLevelSolution ml_sol(&ml_msh);

  // generate solution vector
  FEOrder fe_order = inputparser->getValue("multilevel_solution.multilevel_mesh.first.variable.first.fe_order", FIRST);
  ml_sol.AddSolution("Sol", LAGRANGE, fe_order);

  //Initialize (update Init(...) function)
  ml_sol.Initialize("Sol");

  std::vector<std::string> facenamearray;
  std::vector<ParsedFunction> parsedfunctionarray;
  std::vector<BDCType> bdctypearray;

  if (inputparser->isTrue("multilevel_mesh.first.type", "box")) {
    //Set Boundary (update Dirichlet(...) function)
    ml_sol.InitializeBdc();

    unsigned int bdcsize = inputparser->getSize("multilevel_solution.multilevel_mesh.first.variable.first.boundary_conditions");

    for (unsigned int index = 0; index < bdcsize; ++index) {
      std::string facename = inputparser->getValueFromArray("multilevel_solution.multilevel_mesh.first.variable.first.boundary_conditions", index, "facename", "top");
      facenamearray.push_back(facename);

      BDCType bdctype = inputparser->getValueFromArray("multilevel_solution.multilevel_mesh.first.variable.first.boundary_conditions", index, "bdc_type", DIRICHLET);
      bdctypearray.push_back(bdctype);

      std::string bdcfuncstr = inputparser->getValueFromArray("multilevel_solution.multilevel_mesh.first.variable.first.boundary_conditions", index, "bdc_func", "0.");
      ParsedFunction pfunc(bdcfuncstr, "x,y,z,t");
      parsedfunctionarray.push_back(pfunc);
    }

    for (int i = 0; i < bdcsize; ++i) {
      ml_sol.SetBoundaryCondition_new("Sol", facenamearray[i], bdctypearray[i], false, &parsedfunctionarray[i]);
    }

    ml_sol.GenerateBdc("All");
  }
  else {
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
    ml_sol.GenerateBdc("All");
  }

  MultiLevelProblem ml_prob(&ml_sol);

  // add fluid material
  Parameter parameter(Lref, Uref);

  //BEGIN Poisson MultiLevel Problem
  std::cout << std::endl;
  std::cout << " PDE problem to solve: Poisson " << std::endl;

  LinearImplicitSystem& system2 = ml_prob.add_system<LinearImplicitSystem>("Poisson");
  system2.AddSolutionToSystemPDE("Sol");

  // reading source function
  std::string function;
  function = inputparser->getValue("multilevel_solution.multilevel_mesh.first.variable.first.func_source", "0.");
  std::string variables = "x";
  variables += ",y";
  variables += ",z";
  variables += ",t";

#ifdef HAVE_FPARSER
  fpsource.SetExpression(function);
  fpsource.SetIndependentVariables(variables);
  fpsource.Parse();
#endif

  // Set MG Options
  system2.SetAssembleFunction(AssemblePoissonMatrixandRhs);

  unsigned int max_number_linear_iteration = inputparser->getValue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.max_number_linear_iteration", 6);
  system2.SetMaxNumberOfLinearIterations(max_number_linear_iteration);

  double abs_conv_tol = inputparser->getValue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.abs_conv_tol", 1.e-08);
  system2.SetAbsoluteLinearConvergenceTolerance(abs_conv_tol);

  MgType mgtype = inputparser->getValue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.type.multigrid.mgtype", V_CYCLE);
  system2.SetMgType(mgtype);

  unsigned int npresmoothing = inputparser->getValue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.type.multigrid.npresmoothing", 1);
  system2.SetNumberPreSmoothingStep(npresmoothing);

  unsigned int npostmoothing = inputparser->getValue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.type.multigrid.npostmoothing", 1);
  system2.SetNumberPostSmoothingStep(npostmoothing);

  if (inputparser->isTrue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.type.multigrid.smoother.type", "gmres")) {
    system2.SetLinearEquationSolverType(FEMuS_DEFAULT);
  }
  else if (inputparser->isTrue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.type.multigrid.smoother.type", "asm")) {
    system2.SetLinearEquationSolverType(FEMuS_ASM);
  }

  system2.init();
  //common smoother option
  system2.SetSolverFineGrids(RICHARDSON);
  system2.SetTolerances(1.e-12, 1.e-20, 1.e+50, 4);
  system2.SetPreconditionerFineGrids(SOR_PRECOND);
  //for Vanka and ASM smoothers
  system2.ClearVariablesToBeSolved();
  system2.AddVariableToBeSolved("All");

  if (inputparser->isTrue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.type.multigrid.smoother.type", "asm")) {
    system2.SetNumberOfSchurVariables(0);
    system2.SetElementBlockNumber(4);
  }

  //for Gmres smoother
  system2.SetDirichletBCsHandling(PENALTY);

  // Solve Temperature system
  //system2.PrintSolverInfo(true);
  system2.MGsolve();
  //END Temperature Multilevel Problem

  /// Print all solutions
  std::vector<std::string> print_vars;
  print_vars.push_back("Sol");

  VTKWriter vtkio(&ml_sol);
  vtkio.Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars);


  GMVWriter gmvio(&ml_sol);
  gmvio.SetDebugOutput(true);
  gmvio.Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars);

  //XDMFWriter xdmfio(&ml_sol);
  //xdmfio.Pwrite(DEFAULT_OUTPUTDIR,"biquadratic",print_vars);

  //Destroy all the new systems
  ml_prob.clear();

  return 0;
}

//-----------------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------------
void AssemblePoissonMatrixandRhs(MultiLevelProblem& ml_prob) {

  //pointers and references

  LinearImplicitSystem& mylin_impl_sys = ml_prob.get_system<LinearImplicitSystem>("Poisson");
  const unsigned level = mylin_impl_sys.GetLevelToAssemble();
  //bool assemble_matrix = mylin_impl_sys.GetAssembleMatrix();

  Solution*      mysolution	       = ml_prob._ml_sol->GetSolutionLevel(level);
  LinearEquationSolver*  mylsyspde     = mylin_impl_sys._LinSolver[level];
  Mesh*          mymsh		       = ml_prob._ml_msh->GetLevel(level);
  elem*          myel		       = mymsh->el;
  SparseMatrix*  myKK		       = mylsyspde->_KK;
  NumericVector* myRES		       = mylsyspde->_RES;
  MultiLevelSolution* ml_sol           = ml_prob._ml_sol;

  //data
  const unsigned	dim	= mymsh->GetDimension();
  unsigned 		nel	= mymsh->GetNumberOfElements();
  unsigned 		igrid	= mymsh->GetLevel();
  unsigned 		iproc	= mymsh->processor_id();

  //solution variable
  unsigned SolIndex;
  unsigned SolPdeIndex;
  SolIndex = ml_sol->GetIndex("Sol");
  SolPdeIndex = mylin_impl_sys.GetSolPdeIndex("Sol");
  //solution order
  unsigned order_ind = ml_sol->GetSolutionType(SolIndex);
  //coordinates
  vector< vector < double> > coordinates(dim);

  // declare
  vector< int > metis_node;
  vector< int > KK_dof;
  vector <double> phi;
  vector <double> gradphi;
  vector <double> nablaphi;
  double weight;
  vector< double > F;
  vector< double > B;
  vector<double> normal(3.0);
  double src_term = 0.;
  vector<double> xyzt(4, 0.);
  ParsedFunction* bdcfunc = NULL;

  // reserve
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));
  metis_node.reserve(max_size);
  KK_dof.reserve(max_size);

  for (int i = 0; i < dim; i++)
    coordinates[i].reserve(max_size);

  phi.reserve(max_size);
  gradphi.reserve(max_size * dim);
  unsigned nabla_dim = (3 * (dim - 1) + !(dim - 1));
  nablaphi.reserve(max_size * nabla_dim);
  F.reserve(max_size);
  B.reserve(max_size * max_size);


  // Set to zeto all the entries of the Global Matrix
  myKK->zero();


  // *** element loop ***
  for (int iel = mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielt = mymsh->GetElementType(iel);
    unsigned nve = mymsh->GetElementDofNumber(iel, order_ind);
    unsigned nve2 = mymsh->GetElementDofNumber(iel, 2);
    // resize
    metis_node.resize(nve);
    KK_dof.resize(nve);

    for (int i = 0; i < dim; i++) {
      coordinates[i].resize(nve);
    }

    // set to zero all the entries of the FE matrices
    F.resize(nve);
    memset(&F[0], 0, nve * sizeof(double));

    B.assign(nve * nve, 0.);

    // get local to global mappings
    for (unsigned i = 0; i < nve2; i++) {
      unsigned inode_coord_metis = mymsh->GetSolutionDof(i, iel, 2);

      for (unsigned ivar = 0; ivar < dim; ivar++) {
        coordinates[ivar][i] = (*mymsh->_topology->_Sol[ivar])(inode_coord_metis);
      }

      if (i < nve) {
        metis_node[i] = mymsh->GetSolutionDof(i, iel, order_ind);
        KK_dof[i] = mylsyspde->GetSystemDof(SolIndex, SolPdeIndex, i, iel);
      }

    }




    // *** Gauss point loop ***

    // Supg stabilization tau evaluation
    double V[3] = {0., 0., 0.};
    double nu = 1.;

    if (dim == 1) {
      V[0] = 1.;
      nu = 0.01;
    }
    else if (dim == 2) {
// 	nu=0.0001;
// 	V[0]=sqrt(2)/2;
// 	V[1]=-sqrt(2)/2;
      nu = 1.;
      V[0] = 0.;
      V[1] = 0.;
    }
    else if (dim == 3) {
      nu = 1.;
      V[0] = V[1] = V[2] = 0.;
    }

    double barNu = 0.;
    double vL2Norm2 = 0.;

    for (int i = 0; i < dim; i++) {
      vL2Norm2 += V[i] * V[i];
      unsigned ip = referenceElementDirection[ielt][i][1];
      unsigned im = referenceElementDirection[ielt][i][0];
      double VxiHxi = 0.;

      for (int j = 0; j < dim; j++) {
        VxiHxi += (coordinates[j][ip] - coordinates[j][im]) * V[j];
      }

      double PeXi = VxiHxi / (2.*nu);
      double barXi = (fabs(PeXi) < 1.0e-10) ? 0. : 1. / tanh(PeXi) - 1. / PeXi;
      barNu += barXi * VxiHxi / 2.;
    }

    double supgTau = (vL2Norm2 > 1.0e-15) ? barNu / vL2Norm2 : 0.;
    // End Stabilization stabilization tau evaluation

    for (unsigned ig = 0; ig < mymsh->_finiteElement[ielt][order_ind]->GetGaussPointNumber(); ig++) {
      // *** get Jacobian and test function and test function derivatives ***
      mymsh->_finiteElement[ielt][order_ind]->Jacobian(coordinates, ig, weight, phi, gradphi, nablaphi);
      //current solution
      double SolT = 0;
      vector < double > gradSolT(dim, 0.);
      vector < double > NablaSolT(dim, 0.);

      xyzt.assign(4, 0.);
      unsigned SolType = ml_sol->GetSolutionType("Sol");

      for (unsigned i = 0; i < nve; i++) {
        double soli = (*mysolution->_Sol[SolIndex])(metis_node[i]);

        for (unsigned ivar = 0; ivar < dim; ivar++) {
          xyzt[ivar] += coordinates[ivar][i] * phi[i];
        }

        SolT += phi[i] * soli;

        for (unsigned ivar2 = 0; ivar2 < dim; ivar2++) {
          gradSolT[ivar2] += gradphi[i * dim + ivar2] * soli;
          NablaSolT[ivar2] += nablaphi[i * nabla_dim + ivar2] * soli;
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nve; i++) {
        //BEGIN RESIDUALS A block ===========================
        double advRhs = 0.;
        double lapRhs = 0.;
        double resRhs = 0.;
        double supgPhi = 0.;

        for (unsigned ivar = 0; ivar < dim; ivar++) {
          lapRhs   +=  nu * gradphi[i * dim + ivar] * gradSolT[ivar];
          advRhs   +=  V[ivar] * gradSolT[ivar] * phi[i];
          resRhs   += -nu * NablaSolT[ivar] + V[ivar] * gradSolT[ivar];
          supgPhi  += (V[ivar] * gradphi[i * dim + ivar] + nu * nablaphi[i * nabla_dim + ivar]) * supgTau;
        }

        src_term = fpsource(&xyzt[0]);

        F[i] += (src_term * phi[i] - lapRhs - advRhs
                 + (src_term - resRhs) * supgPhi) * weight;

        //END RESIDUALS A block ===========================

        // *** phi_j loop ***
        for (unsigned j = 0; j < nve; j++) {
          double lap = 0;
          double adv = 0;

          for (unsigned ivar = 0; ivar < dim; ivar++) {
            lap += nu * (gradphi[i * dim + ivar] * gradphi[j * dim + ivar]
                         - nablaphi[j * nabla_dim + ivar] * supgPhi) * weight;
            adv += V[ivar] * gradphi[j * dim + ivar] * (phi[i] + supgPhi) * weight;
          }

          B[i * nve + j] += lap + adv ;
        } // end phij loop
      } // end phii loop
    } // end gauss point loop

    if (ml_prob._ml_sol->_useParsedBCFunction) {

      //number of faces for each type of element
      unsigned nfaces = mymsh->GetElementFaceNumber(iel);

      // loop on faces
      for (unsigned jface = 0; jface < nfaces; jface++) {
        // look for boundary faces
        if (myel->GetBoundaryIndex(iel, jface) > 0) {
          unsigned int faceIndex =  myel->GetBoundaryIndex(iel, jface);

          if (ml_sol->GetBoundaryCondition("Sol", faceIndex - 1u) == NEUMANN && !ml_sol->Ishomogeneous("Sol", faceIndex - 1u)) {
            bdcfunc = (ParsedFunction*)(ml_sol->GetBdcFunction("Sol", faceIndex - 1u));
            unsigned nve = mymsh->GetElementFaceDofNumber(iel, jface, order_ind);
            const unsigned felt = mymsh->GetElementFaceType(iel, jface);

            for (unsigned i = 0; i < nve; i++) {
              unsigned ilocal = mymsh->GetLocalFaceVertexIndex(iel, jface, i);
              unsigned inode_coord_metis = mymsh->GetSolutionDof(ilocal, iel, 2);

              for (unsigned ivar = 0; ivar < dim; ivar++) {
                coordinates[ivar][i] = (*mymsh->_topology->_Sol[ivar])(inode_coord_metis);
              }
            }

            if (felt != 6) {
              for (unsigned igs = 0; igs < mymsh->_finiteElement[felt][order_ind]->GetGaussPointNumber(); igs++) {
                mymsh->_finiteElement[felt][order_ind]->JacobianSur(coordinates, igs, weight, phi, gradphi, normal);

                xyzt.assign(4, 0.);

                for (unsigned i = 0; i < nve; i++) {
                  for (unsigned ivar = 0; ivar < dim; ivar++) {
                    xyzt[ivar] += coordinates[ivar][i] * phi[i];
                  }
                }

                // *** phi_i loop ***
                for (unsigned i = 0; i < nve; i++) {
                  double surfterm_g = (*bdcfunc)(&xyzt[0]);
                  double bdintegral = phi[i] * surfterm_g * weight;
                  unsigned int ilocalnode = mymsh->GetLocalFaceVertexIndex(iel, jface, i);
                  F[ilocalnode] += bdintegral;
                }
              }
            }
            else { // 1D : the side elems are points and does not still exist the point elem
              // in 1D it is only one point
              xyzt[0] = coordinates[0][0];
              xyzt[1] = 0.;
              xyzt[2] = 0.;
              xyzt[3] = 0.;

              double bdintegral = (*bdcfunc)(&xyzt[0]);
              unsigned int ilocalnode = mymsh->GetLocalFaceVertexIndex(iel, jface, 0);
              F[ilocalnode] += bdintegral;
            }
          }
        }
      }
    }
    else {
      double tau = 0.;
      std::vector< double > xx(dim, 0.);
      vector < double > normal(dim, 0);

      // loop on faces
      for (unsigned jface = 0; jface < mymsh->GetElementFaceNumber(iel); jface++) {
        // look for boundary faces
        if (myel->GetFaceElementIndex(iel, jface) < 0) {
          unsigned int faceIndex =  myel->GetBoundaryIndex(iel, jface);

          if (!ml_sol->GetBdcFunction()(xx, "Sol", tau, faceIndex, 0.) && tau != 0.) {
            unsigned nve = mymsh->GetElementFaceDofNumber(iel, jface, order_ind);
            const unsigned felt = mymsh->GetElementFaceType(iel, jface);

            for (unsigned i = 0; i < nve; i++) {
              unsigned int ilocal = mymsh->GetLocalFaceVertexIndex(iel, jface, i);
              unsigned inode = mymsh->GetSolutionDof(ilocal, iel, 2);

              for (unsigned idim = 0; idim < dim; idim++) {
                coordinates[idim][i] = (*mymsh->_topology->_Sol[idim])(inode);
              }
            }

            for (unsigned igs = 0; igs < mymsh->_finiteElement[felt][order_ind]->GetGaussPointNumber(); igs++) {
              mymsh->_finiteElement[felt][order_ind]->JacobianSur(coordinates, igs, weight, phi, gradphi, normal);

              // *** phi_i loop ***
              for (unsigned i = 0; i < nve; i++) {
                double value = phi[i] * tau * weight;
                unsigned int ilocal = mymsh->GetLocalFaceVertexIndex(iel, jface, i);
                F[ilocal]   += value;
              }
            }
          }
        }
      }
    }

    //--------------------------------------------------------------------------------------------------------
    //Sum the local matrices/vectors into the global Matrix/vector

    myRES->add_vector_blocked(F, KK_dof);

    myKK->add_matrix_blocked(B, KK_dof, KK_dof);
  } //end list of elements loop for each subdomain

  myRES->close();
  myKK->close();

  // ***************** END ASSEMBLY *******************

}


