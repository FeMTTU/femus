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
#include <json/json.h>
#include <json/value.h>
#include "ParsedFunction.hpp"
#include "Files.hpp"
#include <cstdlib>


using std::cout;
using std::endl;

using namespace femus;

void AssemblePoissonMatrixandRhs(MultiLevelProblem& ml_prob);

double GetRelativeError(MultiLevelSolution& ml_sol, const bool& H1);

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;
  return dirichlet;
}

static std::string
readInputTestFile(const char* path)
{
  FILE* file = fopen(path, "rb");

  if (!file)
    return std::string("");

  fseek(file, 0, SEEK_END);
  long size = ftell(file);
  fseek(file, 0, SEEK_SET);
  std::string text;
  char* buffer = new char[size + 1000];
  buffer[size] = 0;

  if (fread(buffer, 1, size, file) == (unsigned long)size)
    text = buffer;

  fclose(file);
  delete[] buffer;
  return text;
}

static void show_usage()
{
  std::cout << "Use --inputfile variable to set the input file" << endl;
  std::cout << "e.g.: ./AMR_ex20 --inputfile ./input/input.json" << std::endl;
}

ParsedFunction fpsource;
ParsedFunction fp_sol;
ParsedFunction fp_dsoldx;
ParsedFunction fp_dsoldy;


int main(int argc, char** argv) {

  std::string path;

  if (argc < 2)
  {
    std::cout << "Error: no input file specified!" << std::endl;
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
  }

  // start reading input from file
  //-----------------------------------------------------------------------------------------------

  std::string input = readInputTestFile(path.c_str());

  if (input.empty())
  {
    printf("Failed to read input or empty input: %s\n", path.c_str());
    return 1;
  }

  Json::Value root;   // will contains the root value after parsing.
  Json::Reader reader;
  bool parsingSuccessful = reader.parse(input, root);

  if (!parsingSuccessful)
  {
    // report to the user the failure and their locations in the document.
    std::cout  << "Failed to parse configuration\n" << reader.getFormatedErrorMessages();
    return 1;
  }

  std::string filename = root["mesh"].get("filename", "").asString();

  int numelemx;
  int numelemy;
  int numelemz;
  double xa, xb, ya, yb, za, zb;
  ElemType elemtype;

  bool isBox = root["mesh"].get("box", false).asBool();

  if (isBox) {
    numelemx = root["mesh"].get("box", "").get("nx", 2).asUInt();
    numelemy = root["mesh"].get("box", "").get("ny", 2).asUInt();
    numelemz = root["mesh"].get("box", "").get("nz", 2).asUInt();
    xa = root["mesh"].get("box", "").get("xa", 0.).asDouble();
    xb = root["mesh"].get("box", "").get("xb", 1.).asDouble();
    ya = root["mesh"].get("box", "").get("ya", 0.).asDouble();
    yb = root["mesh"].get("box", "").get("yb", 1.).asDouble();
    za = root["mesh"].get("box", "").get("za", 0.).asDouble();
    zb = root["mesh"].get("box", "").get("zb", 0.).asDouble();
    std::string elemtypestr = root["mesh"].get("box", "").get("elem_type", "Quad9").asString();

    if (elemtypestr == "Quad9")
    {
      elemtype = QUAD9;
    }
    else if (elemtypestr == "Tri6")
    {
      elemtype = TRI6;
    }
    else if (elemtypestr == "Edge3")
    {
      elemtype = EDGE3;
    }
    else if (elemtypestr == "Hex27")
    {
      elemtype = HEX27;
    }
    else
    {
      elemtype = INVALID_ELEM;
    }
    
  }

  std::string variableName = root["variable"].get("name", "Q").asString();

  FEOrder fe_order;
  std::string fe_order_str = root["variable"].get("fe_order", "first").asString();

  if (!strcmp(fe_order_str.c_str(), "first"))
  {
    fe_order = FIRST;
  }
  else if (!strcmp(fe_order_str.c_str(), "serendipity"))
  {
    fe_order = SERENDIPITY;
  }
  else if (!strcmp(fe_order_str.c_str(), "second"))
  {
    fe_order = SECOND;
  }
  else
  {
    std::cerr << " Error: Lagrange finite element order not supported!" << std::endl;
    exit(1);
  }


  unsigned int nlevels 	= root["mgsolver"].get("nlevels", 1).asInt();
  unsigned int SMRlevels 	= root["mgsolver"].get("SMRlevels", 0).asInt();
  std::string  AMR       	= root["mgsolver"].get("AMR", "no").asString();
  unsigned int maxAMRlevels 	= root["mgsolver"].get("maxAMRlevels", 0).asInt();
  std::string  AMRnorm 	= root["mgsolver"].get("AMRnorm", "l2").asString();
  double       AMRthreshold 	= root["mgsolver"].get("AMRthreshold", "0.001").asDouble();
  unsigned int npresmoothing 	= root["mgsolver"].get("npresmoothing", 1).asUInt();
  unsigned int npostmoothing 	= root["mgsolver"].get("npostsmoothing", 1).asUInt();
  std::string smoother_type  	= root["mgsolver"].get("smoother_type", "gmres").asString();
  std::string mg_type        	= root["mgsolver"].get("mgtype", "V_cycle").asString();
  unsigned int max_number_linear_iteration = root["mgsolver"].get("max_number_linear_iteration", 6).asUInt();
  double abs_conv_tol        = root["mgsolver"].get("abs_conv_tol", 1.e-09).asDouble();

  MgType mgtype;

  if (!strcmp("V_cycle", mg_type.c_str()))
  {
    mgtype = V_CYCLE;
  }
  else if (!strcmp("F_cycle", mg_type.c_str()))
  {
    mgtype = F_CYCLE;
  }
  else if (!strcmp("F_cycle", mg_type.c_str()))
  {
    mgtype = F_CYCLE;
  }
  else {
    cout << "The selected MG cycle does not exist!" << endl;
    exit(1);
  }

  bool Gmres = 0, Asm = 0;
  if (!strcmp("gmres", smoother_type.c_str()))     Gmres = 1;
  else if (!strcmp("asm", smoother_type.c_str()))       Asm = 1;

  if (Gmres + Asm == 0) {
    cout << "The selected MG smoother does not exist!" << endl;
    exit(1);
  }

  // end reading input from file
  //-----------------------------------------------------------------------------------------------


  // reading function
  std::string variables = "x";
  variables += ",y";
  variables += ",z";
  variables += ",t";
  std::string function;

#ifdef HAVE_FPARSER

  function = root["variable"].get("func_source1", "0.").asString();
  fpsource.SetExpression(function);
  fpsource.SetIndependentVariables(variables);
  fpsource.Parse();

  function = root["func_sol1"].get("sol", "0.").asString();
  fp_sol.SetExpression(function);
  fp_sol.SetIndependentVariables(variables);
  fp_sol.Parse();

  function = root["func_sol1"].get("dsoldx", "0.").asString();
  fp_dsoldx.SetExpression(function);
  fp_dsoldx.SetIndependentVariables(variables);
  fp_dsoldx.Parse();

  function = root["func_sol1"].get("dsoldy", "0.").asString();
  fp_dsoldy.SetExpression(function);
  fp_dsoldy.SetIndependentVariables(variables);
  fp_dsoldy.Parse();


#endif

  std::vector<std::string> facenamearray;
  std::vector<ParsedFunction> parsedfunctionarray;
  std::vector<BDCType> bdctypearray;

//   const Json::Value boundary_conditions = root["boundary_conditions"];
// 
//   for (unsigned int index = 0; index < boundary_conditions.size(); ++index) {
// 
//     std::string facename = boundary_conditions[index].get("facename", "to").asString();
//     facenamearray.push_back(facename);
// 
//     std::string bdctypestr = boundary_conditions[index].get("bdc_type", "to").asString();
//     BDCType bdctype = DIRICHLET;
// 
//     if (bdctypestr.compare("dirichlet") == 0) {
//       bdctype = DIRICHLET;
//     }
//     else if (bdctypestr.compare("neumann") == 0) {
//       bdctype = NEUMANN;
//     }
//     else {
//       std::cout << "Boundary condition not implemented! Default one (Dirichlet) is set." << std::endl;
//     }
// 
//     bdctypearray.push_back(bdctype);
// 
//     std::string bdcfuncstr = boundary_conditions[index].get("bdc_func", "0.").asString();
//     ParsedFunction pfunc(bdcfuncstr, "x,y,z,t");
//     parsedfunctionarray.push_back(pfunc);
//   }

  //---------------------------------------------------------------------------


  /// Init Petsc-MPI communicator
  FemusInit mpinit(argc, argv, MPI_COMM_WORLD);

  //Files files;
  //files.CheckIODirectories();
  //files.RedirectCout();

  /// INIT MESH =================================

  unsigned short nm, nr;
  nm = nlevels;
  nr = 0; //SMRlevels;


  int tmp = nm;
  nm += nr;
  nr = tmp;

  //Adimensional quantity (Lref,Uref)
  double Lref = 1.;
  double Uref = 1.;

  //Steadystate NonLinearMultiLevelProblem
  MultiLevelMesh ml_msh;

  if (filename != "")
  {
    ml_msh.ReadCoarseMesh(filename.c_str(), "seventh", Lref);
  }
  else
  {
    ml_msh.GenerateCoarseBoxMesh(numelemx, numelemy, numelemz, xa, xb, ya, yb, za, zb, elemtype, "seventh");
  }

  //ml_msh.RefineMesh(nm,nr, SetRefinementFlag);
  ml_msh.RefineMesh(nm, nr, NULL);

  ml_msh.PrintInfo();

  MultiLevelSolution ml_sol(&ml_msh);

  // generate solution vector
  ml_sol.AddSolution("Sol", LAGRANGE, fe_order);

  //Initialize (update Init(...) function)
  ml_sol.Initialize("Sol");

  //Set Boundary (update Dirichlet(...) function)
 // ml_sol.InitializeBdc();

//     ml_sol.SetBoundaryCondition("Sol","right", NEUMANN, false, false, &bdcfunc);
//     ml_sol.SetBoundaryCondition("Sol","top", NEUMANN);

//   for (int i = 0; i < boundary_conditions.size(); ++i) {
//     ml_sol.SetBoundaryCondition_new("Sol", facenamearray[i], bdctypearray[i], false, &parsedfunctionarray[i]);
//   }

  //ml_sol.GenerateBdc("All");


  ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  ml_sol.GenerateBdc("Sol");


  std::cout<<"AAAAAAAAAAAAAAAAAAAAA"<<std::flush;

  // ml_msh.AddMeshLevel(SetRefinementFlag);
  // ml_sol.AddSolutionLevel();


  MultiLevelProblem ml_prob(&ml_sol);


// ml_prob.parameters.set<func>("func_source") = fpsource;

  // add fluid material
  Parameter parameter(Lref, Uref);

  //BEGIN Poisson MultiLevel Problem
  std::cout << std::endl;
  std::cout << " PDE problem to solve: Poisson " << std::endl;

  LinearImplicitSystem& system2 = ml_prob.add_system<LinearImplicitSystem>("Poisson");
  system2.AddSolutionToSystemPDE("Sol");

  // Set MG Options
  system2.SetAssembleFunction(AssemblePoissonMatrixandRhs);
  system2.SetMaxNumberOfLinearIterations(max_number_linear_iteration);
  system2.SetAbsoluteLinearConvergenceTolerance(abs_conv_tol);
  system2.SetMgType(mgtype);
  system2.SetNumberPreSmoothingStep(npresmoothing);
  system2.SetNumberPostSmoothingStep(npostmoothing);

  //Set Smoother Options
  if (Gmres) 		system2.SetLinearEquationSolverType(FEMuS_DEFAULT);
  else if (Asm) 	system2.SetLinearEquationSolverType(FEMuS_ASM);

  system2.init();

  //system2.SetAMRSetOptions(AMR,AMRlevels,AMRnorm,AMRthreshold,SetRefinementFlag);
  system2.SetAMRSetOptions(AMR, maxAMRlevels, AMRnorm, AMRthreshold);

  //common smoother option
  //system2.SetSolverFineGrids(GMRES);
  system2.SetSolverFineGrids(RICHARDSON);
  system2.SetTolerances(1.e-12, 1.e-20, 1.e+50, 4);
  system2.SetPreconditionerFineGrids(ILU_PRECOND);
  //for Vanka and ASM smoothers
  system2.ClearVariablesToBeSolved();
  system2.AddVariableToBeSolved("All");

  if (Asm )
  {
    system2.SetNumberOfSchurVariables(0);
    system2.SetElementBlockNumber(4);
  }

  //for Gmres smoother
  system2.SetDirichletBCsHandling(PENALTY);

  // Solve Temperature system
  
  system2.MGsolve();
  //END Temperature Multilevel Problem

  /// Print all solutions
  std::vector<std::string> print_vars;
  print_vars.push_back("Sol");

  //VTKWriter vtkio(&ml_sol);
  //vtkio.ParallelWrite(files.GetOutputPath(),"biquadratic",print_vars);

  //GMVWriter gmvio(&ml_sol);
  //gmvio.ParallelWrite(files.GetOutputPath(),"biquadratic",print_vars);

  VTKWriter vtkio(&ml_sol);
  vtkio.Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars);

  //GMVWriter gmvio(&ml_sol);
  //gmvio.ParallelWrite(DEFAULT_OUTPUTDIR,"biquadratic",print_vars);

  //
  //     XDMFWriter xdmfio(ml_sol);
  //     xdmfio.write("biquadratic",print_vars);

  int  iproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  if(AMRnorm == "l2") {
  double l2error = GetRelativeError(ml_sol, false);

  if (iproc == 0) printf("\n||Sol_h-Sol||_L2 / ||Sol||_L2  = %g \n", l2error);
   ml_msh.PrintInfo();
  }
  else if(AMRnorm == "H1") {
  double H1error = GetRelativeError(ml_sol, true);

  if (iproc == 0) printf("\n||Sol_h-Sol||_H1 / ||Sol||_H1  = %g \n", H1error);
   ml_msh.PrintInfo();
      
 }
 else {std::cout << std::endl<< std::endl << " *********AMRnorm use 'H1' or 'l2'************" << std::endl;}

  //Destroy all the new systems
  ml_prob.clear();

  return 0;
}


//------------------------------------------------------------------------------------------------------------
void AssemblePoissonMatrixandRhs(MultiLevelProblem& ml_prob) {

  //pointers and references
  LinearImplicitSystem& mylin_impl_sys = ml_prob.get_system<LinearImplicitSystem>("Poisson");
  const unsigned level = mylin_impl_sys.GetLevelToAssemble();

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
  nablaphi.reserve(max_size * (3 * (dim - 1) + !(dim - 1)));
  F.reserve(max_size);
  B.reserve(max_size * max_size);

  // Set to zeto all the entries of the Global Matrix
  myKK->zero();

  // *** element loop ***
  for (int iel = mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielt = mymsh->GetElementType(iel);
    unsigned nve = mymsh->GetElementDofNumber(iel, order_ind);

    // resize
    metis_node.resize(nve);
    KK_dof.resize(nve);
    phi.resize(nve);
    gradphi.resize(nve * dim);
    nablaphi.resize(nve * (3 * (dim - 1) + !(dim - 1)));

    for (int i = 0; i < dim; i++) {
      coordinates[i].resize(nve);
    }

    // set to zero all the entries of the FE matrices
    F.resize(nve);
    memset(&F[0], 0, nve * sizeof(double));
    B.resize(nve * nve);
    memset(&B[0], 0, nve * nve * sizeof(double));


    // get local to global mappings
    for (unsigned i = 0; i < nve; i++) {
      unsigned inode_coord_metis = mymsh->GetSolutionDof(i, iel, 2);
      metis_node[i] = mymsh->GetSolutionDof(i, iel, order_ind);

      for (unsigned ivar = 0; ivar < dim; ivar++) {
        coordinates[ivar][i] = (*mymsh->_topology->_Sol[ivar])(inode_coord_metis);
      }

      KK_dof[i] = mylsyspde->GetSystemDof(SolIndex, SolPdeIndex, i, iel);
    }


    // *** Gauss poit loop ***
    for (unsigned ig = 0; ig < ml_prob._ml_msh->_finiteElement[ielt][order_ind]->GetGaussPointNumber(); ig++) {
      // *** get Jacobian and test function and test function derivatives ***
      ml_prob._ml_msh->_finiteElement[ielt][order_ind]->Jacobian(coordinates, ig, weight, phi, gradphi, nablaphi);
      //current solution
      double SolT = 0;
      vector < double > gradSolT(dim, 0.);

      for (unsigned ivar = 0; ivar < dim; ivar++) {
        gradSolT[ivar] = 0;
      }

      xyzt.assign(4, 0.);
      unsigned SolType = ml_sol->GetSolutionType("Sol");

      for (unsigned i = 0; i < nve; i++) {
        double soli = (*mysolution->_Sol[SolIndex])(metis_node[i]);

        for (unsigned ivar = 0; ivar < dim; ivar++) {
          xyzt[ivar] += coordinates[ivar][i] * phi[i];
        }

        SolT += phi[i] * soli;

        for (unsigned ivar2 = 0; ivar2 < dim; ivar2++) gradSolT[ivar2] += gradphi[i * dim + ivar2] * soli;
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nve; i++) {
        //BEGIN RESIDUALS A block ===========================
        double Adv_rhs = 0;
        double Lap_rhs = 0;

        for (unsigned ivar = 0; ivar < dim; ivar++) {
          Lap_rhs += gradphi[i * dim + ivar] * gradSolT[ivar];
        }


        //src_term = Source(xyzt);
#ifdef HAVE_FPARSER
        //src_term = fpsource.Eval(xyzt);
        src_term = fpsource(&xyzt[0]);
#endif

        F[i] += (-Lap_rhs + src_term * phi[i]) * weight;

        //END RESIDUALS A block ===========================
        // *** phi_j loop ***
        for (unsigned j = 0; j < nve; j++) {
          double Lap = 0;
          double Adv1 = 0;

          for (unsigned ivar = 0; ivar < dim; ivar++) {
            // Laplacian
            Lap  += gradphi[i * dim + ivar] * gradphi[j * dim + ivar] * weight;
          }

          B[i * nve + j] += Lap;
        } // end phij loop

      } // endif assemble_matrix
    } // end gauss point loop

    //number of faces for each type of element
    //number of faces for each type of element
//     unsigned nfaces = mymsh->GetElementFaceNumber(iel);
// 
//     // loop on faces
//     for (unsigned jface = 0; jface < nfaces; jface++) {
// 
//       // look for boundary faces
//       if (myel->GetBoundaryIndex(iel, jface) > 0) {
// 
//         unsigned int faceIndex =  myel->GetBoundaryIndex(iel, jface) - 1u;
// 
//         if (ml_sol->GetBoundaryCondition("Sol", faceIndex) == NEUMANN && !ml_sol->Ishomogeneous("Sol", faceIndex)) {
// 
//           bdcfunc = (ParsedFunction*)(ml_sol->GetBdcFunction("Sol", faceIndex));
//           unsigned nve = mymsh->GetElementFaceDofNumber(iel, jface, order_ind);
//           const unsigned felt = mymsh->GetElementFaceType(iel, jface);
// 
//           for (unsigned i = 0; i < nve; i++) {
//             unsigned ilocal = mymsh->GetLocalFaceVertexIndex(iel, jface, i);
//             unsigned inode_coord_metis = mymsh->GetSolutionDof(ilocal, iel, 2);
// 
//             for (unsigned ivar = 0; ivar < dim; ivar++) {
//               coordinates[ivar][i] = (*mymsh->_topology->_Sol[ivar])(inode_coord_metis);
//             }
//           }
// 
//           if (felt != 6)
//           {
//             for (unsigned igs = 0; igs < ml_prob._ml_msh->_finiteElement[felt][order_ind]->GetGaussPointNumber(); igs++) {
//               ml_prob._ml_msh->_finiteElement[felt][order_ind]->JacobianSur(coordinates, igs, weight, phi, gradphi, normal);
// 
//               xyzt.assign(4, 0.);
// 
//               for (unsigned i = 0; i < nve; i++) {
//                 for (unsigned ivar = 0; ivar < dim; ivar++) {
//                   xyzt[ivar] += coordinates[ivar][i] * phi[i];
//                 }
//               }
// 
//               // *** phi_i loop ***
//               for (unsigned i = 0; i < nve; i++) {
//                 double surfterm_g = (*bdcfunc)(&xyzt[0]);
//                 double bdintegral = phi[i] * surfterm_g * weight;
//                 unsigned int ilocalnode = mymsh->GetLocalFaceVertexIndex(iel, jface, i);
//                 F[ilocalnode] += bdintegral;
//               }
//             }
//           }
//           else // 1D : the side elems are points and does not still exist the point elem
//           {
//             // in 1D it is only one point
//             xyzt[0] = coordinates[0][0];
//             xyzt[1] = 0.;
//             xyzt[2] = 0.;
//             xyzt[3] = 0.;
// 
//             double bdintegral = (*bdcfunc)(&xyzt[0]);
//             unsigned int ilocalnode = mymsh->GetLocalFaceVertexIndex(iel, jface, 0);
//             F[ilocalnode] += bdintegral;
//           }
//         }
//       }
//     }


    //--------------------------------------------------------------------------------------------------------
    //Sum the local matrices/vectors into the global Matrix/vector

    myRES->add_vector_blocked(F, KK_dof);
    myKK->add_matrix_blocked(B, KK_dof, KK_dof);
  } //end list of elements loop for each subdomain

  myRES->close();
  myKK->close();

  // ***************** END ASSEMBLY *******************

}

// ***********************************************

double GetRelativeError(MultiLevelSolution& ml_sol, const bool& H1) {

  int  iproc, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


  NumericVector* error_vec;
  NumericVector* solution_vec;
  error_vec = NumericVector::build().release();
  solution_vec = NumericVector::build().release();

  if (nprocs == 1) {
    error_vec->init(nprocs, 1, false, SERIAL);
    solution_vec->init(nprocs, 1, false, SERIAL);
  }
  else {
    error_vec->init(nprocs, 1, false, PARALLEL);
    solution_vec->init(nprocs, 1, false, PARALLEL);
  }

  error_vec->zero();
  solution_vec->zero();

  unsigned gridn = ml_sol._mlMesh->GetNumberOfLevels();

  for (int ilevel = gridn - 1; ilevel < gridn; ilevel++) {
    Solution*      solution  = ml_sol.GetSolutionLevel(ilevel);
    Mesh*          msh	     = ml_sol._mlMesh->GetLevel(ilevel);
    unsigned 	   iproc     = msh->processor_id();


    const unsigned	dim	= msh->GetDimension();
    vector< vector < double> > coordinates(dim);

    vector< int > metis_node;

    vector <double> phi;
    vector <double> gradphi;
    vector <double> nablaphi;

    double weight;

    // reserve
    const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));
    metis_node.reserve(max_size);

    for (int i = 0; i < dim; i++)
      coordinates[i].reserve(max_size);

    phi.reserve(max_size);
    gradphi.reserve(max_size * dim);
    nablaphi.reserve(max_size * (3 * (dim - 1) + !(dim - 1)));


    unsigned SolIndex;
    SolIndex = ml_sol.GetIndex("Sol");
    unsigned SolOrder = ml_sol.GetSolutionType(SolIndex);


    for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

      short unsigned ielt = msh->GetElementType(iel);
      unsigned nve = msh->GetElementDofNumber(iel, SolOrder);
      // resize
      metis_node.resize(nve);
      phi.resize(nve);
      gradphi.resize(nve * dim);
      nablaphi.resize(nve * (3 * (dim - 1) + !(dim - 1)));

      for (int i = 0; i < dim; i++) {
        coordinates[i].resize(nve);
      }

      // get local to global mappings
      for (unsigned i = 0; i < nve; i++) {
        unsigned inode_coord_metis = msh->GetSolutionDof(i, iel, 2);
        metis_node[i] = msh->GetSolutionDof(i, iel, SolOrder);

        for (unsigned idim = 0; idim < dim; idim++) {
          coordinates[idim][i] = (*msh->_topology->_Sol[idim])(inode_coord_metis);
        }
      }

      for (unsigned ig = 0; ig < ml_sol._mlMesh->_finiteElement[ielt][SolOrder]->GetGaussPointNumber(); ig++) {
        // *** get Jacobian and test function and test function derivatives ***
        ml_sol._mlMesh->_finiteElement[ielt][SolOrder]->Jacobian(coordinates, ig, weight, phi, gradphi, nablaphi);
        //current solution
        double SolT = 0;
        vector < double > gradSolT(dim, 0.);

        for (unsigned ivar = 0; ivar < dim; ivar++) {
          gradSolT[ivar] = 0;
        }

        double pi = acos(-1.);
        double x[4] = {0., 0., 0., 0.};

        unsigned SolType = ml_sol.GetSolutionType("Sol");

        for (unsigned i = 0; i < nve; i++) {
          double soli = (*solution->_Sol[SolIndex])(metis_node[i]);

          for (unsigned ivar = 0; ivar < dim; ivar++) {
            x[ivar] += coordinates[ivar][i] * phi[i];
          }

          SolT += phi[i] * soli;

          for (unsigned ivar2 = 0; ivar2 < dim; ivar2++) {
            gradSolT[ivar2] += gradphi[i * dim + ivar2] * soli;
          }
        }

        double SolExact, dSolExactdx, dSolExactdy;
#ifdef HAVE_FPARSER
        SolExact    = fp_sol(x);
        dSolExactdx = fp_dsoldx(x);
        dSolExactdy = fp_dsoldy(x);
#endif

        error_vec->add(iproc, ((SolT - SolExact) * (SolT - SolExact) +
                               H1 * ((gradSolT[0] - dSolExactdx) * (gradSolT[0] - dSolExactdx) +
                                     (gradSolT[1] - dSolExactdy) * (gradSolT[1] - dSolExactdy))
                              )*weight);

        solution_vec->add(iproc, (SolExact * SolExact +
                                  H1 * (dSolExactdx * dSolExactdx + dSolExactdy * dSolExactdy))*weight);
      }
    }

  }

  error_vec->close();
  solution_vec->close();

  double l2_error = error_vec->l1_norm();
  double l2_solution = solution_vec->l1_norm();

  delete error_vec;
  delete solution_vec;

  return sqrt(l2_error / l2_solution);

}

