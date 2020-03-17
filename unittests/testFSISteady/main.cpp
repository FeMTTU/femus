#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"
#include "FemusInit.hpp"
#include "SparseMatrix.hpp"
#include "VTKWriter.hpp"
#include "MonolithicFSINonLinearImplicitSystem.hpp"
#include "FElemTypeEnum.hpp"

using std::cout;
using std::endl;
using namespace femus;

void AssembleMatrixResFSI(MultiLevelProblem& ml_prob);

bool SetBoundaryCondition(const std::vector < double >& x, const char name[],
                          double& value, const int FaceName, const double = 0.);

//------------------------------------------------------------------------------------------------------------------

int main(int argc, char** args) {

  /// Init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  unsigned short nm, nr;
  std::cout << "#MULTIGRID levels? (>=1) \n";
  //std::cin>>nm;
  nm = 4;

  std::cout << "#MAX_REFINEMENT levels? (>=0) \n";
  //std::cin>>nr;
  nr = 0;
  int tmp = nm;
  nm += nr;
  nr = tmp;

  char* infile = new char [50];

  sprintf(infile, "./input/fsifirst.neu");

  double Lref = 1.;
  double Uref = 1.;
  double rhof = 1000.;
  double muf = 1.;
  double rhos = 1000;
  double ni = 0.4;
  double E = 1400000;

  MultiLevelMesh ml_msh(nm, nr, infile, "fifth", Lref, NULL);

  MultiLevelSolution ml_sol(&ml_msh);

  //Start System Variables
  ml_sol.AddSolution("DX", LAGRANGE, SECOND, 1);
  ml_sol.AddSolution("DY", LAGRANGE, SECOND, 1);
//   ml_sol.AssociatePropertyToSolution("DX","Displacement"); // Add this line
//   ml_sol.AssociatePropertyToSolution("DY","Displacement"); // Add this line
  ml_sol.AddSolution("U", LAGRANGE, SECOND, 1);
  ml_sol.AddSolution("V", LAGRANGE, SECOND, 1);

  // Since the Pressure is a Lagrange multiplier it is used as an implicit variable
  ml_sol.AddSolution("P", DISCONTINUOUS_POLYNOMIAL, FIRST, 1);
  ml_sol.AssociatePropertyToSolution("P", "Pressure"); // Add this line

  //Initialize (update Init(...) function)
  ml_sol.Initialize("All");

  //Set Boundary (update Dirichlet(...) function)
  ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  ml_sol.GenerateBdc("DX", "Steady");
  ml_sol.GenerateBdc("DY", "Steady");
  ml_sol.GenerateBdc("U", "Steady");
  ml_sol.GenerateBdc("V", "Steady");
  ml_sol.GenerateBdc("P", "Steady");

  MultiLevelProblem ml_prob(&ml_sol);

  Parameter par(Lref, Uref);

  // Generate Solid Object
  Solid solid(par, E, ni, rhos, "Neo-Hookean");
  cout << "Solid properties: " << endl;
  cout << solid << endl;

  // Generate Fluid Object
  Fluid fluid(par, muf, rhof, "Newtonian");
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;

  // Add fluid object
  ml_prob.parameters.set<Fluid>("Fluid") = fluid;

  // Add Solid Object
  ml_prob.parameters.set<Solid>("Solid") = solid;

  //create systems
  // add the system FSI to the MultiLevel problem
  MonolithicFSINonLinearImplicitSystem& system = ml_prob.add_system<MonolithicFSINonLinearImplicitSystem> ("Fluid-Structure-Interaction");
  system.AddSolutionToSystemPDE("DX");
  system.AddSolutionToSystemPDE("DY");
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  system.AddSolutionToSystemPDE("P");

  // init all the systems
  system.init();

  // System Fluid-Structure-Interaction
  system.SetAssembleFunction(AssembleMatrixResFSI);
  system.SetMaxNumberOfLinearIterations(1);
  system.SetAbsoluteLinearConvergenceTolerance(1.e-8);
  system.SetMgType(F_CYCLE);
  system.SetMaxNumberOfNonLinearIterations(4);
  system.SetNonLinearConvergenceTolerance(1.e-5);

  system.SetDirichletBCsHandling(PENALTY);
  system.SetSolverFineGrids(GMRES);
  system.SetPreconditionerFineGrids(ILU_PRECOND);
  system.SetTolerances(1.e-12, 1.e-20, 1.e+50, 20);


  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  VTKWriter vtkio(&ml_sol);
  vtkio.SetMovingMesh(mov_vars);

  // Solving Fluid-Structure-Interaction system
  std::cout << std::endl;
  std::cout << " *********** Fluid-Structure-Interaction ************  " << std::endl;
  system.SetOuterSolver(PREONLY);
  system.MGsolve();

  double l2normvarDX = ml_sol.GetSolutionLevel(3)->GetSolutionName("DX").l2_norm();

  double l2normvarDXStored = 0.00422796021240;

  //std::cout.precision(14);

  std::cout << "Solution DX l2norm: " << l2normvarDX << std::endl;

  if (fabs(l2normvarDX - l2normvarDXStored) > 1.e-07)
  {
    exit(1);
  }

  double l2normvarDY = ml_sol.GetSolutionLevel(3)->GetSolutionName("DY").l2_norm();

  double l2normvarDYStored = 0.06728194901640;

  std::cout << "Solution DY l2norm: " << l2normvarDY << std::endl;

  if (fabs(l2normvarDY - l2normvarDYStored) > 1.e-07)
  {
    exit(1);
  }


  double l2normvarU = ml_sol.GetSolutionLevel(3)->GetSolutionName("U").l2_norm();

  double l2normvarUStored = 43.30221796101648;

  std::cout << "Solution U l2norm: " << l2normvarU << std::endl;

  if (fabs(l2normvarU - l2normvarUStored) > 1.e-06)
  {
    exit(1);
  }

  double l2normvarV = ml_sol.GetSolutionLevel(3)->GetSolutionName("V").l2_norm();

  double l2normvarVStored = 9.83398554915716;

  std::cout << "Solution V l2norm: " << l2normvarV << std::endl;

  if (fabs(l2normvarV - l2normvarVStored) > 1.e-06)
  {
    exit(1);
  }

  double l2normvarP = ml_sol.GetSolutionLevel(3)->GetSolutionName("P").l2_norm();

  double l2normvarPStored = 5.87173860743601;

  std::cout << "Solution P l2norm: " << l2normvarP << std::endl;

  if (fabs(l2normvarP - l2normvarPStored) > 1.e-05)
  {
    exit(1);
  }

//   // print solutions
//   std::vector < std::string > variablesToBePrinted;
//   variablesToBePrinted.push_back("All");
//
//   VTKWriter vtkIO(&ml_sol);
//   vtkIO.write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);


  // Destroy all the new systems
  ml_prob.clear();

  delete [] infile;
  return 0;
}

//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int facename, const double time) {
  bool test = 1; //dirichlet
  value = 0.;

  if (!strcmp(name, "U")) {
    if (1 == facename) { //inflow
      test = 1;
      double um = 0.2;
      value = 1.5 * um * 4.0 / 0.1681 * x[1] * (0.41 - x[1]);
    }
    else if (2 == facename) { //outflow
      test = 0;
//    test=1;
      value = 0.;
    }
    else if (3 == facename) { // no-slip fluid wall
      test = 1;
      value = 0.;
    }
    else if (4 == facename) { // no-slip solid wall
      test = 1;
      value = 0.;
    }
  }
  else if (!strcmp(name, "V")) {
    if (1 == facename) {        //inflow
      test = 1;
      value = 0.;
    }
    else if (2 == facename) {   //outflow
      test = 0;
//    test=1;
      value = 0.;
    }
    else if (3 == facename) {   // no-slip fluid wall
      test = 1;
      value = 0;
    }
    else if (4 == facename) {   // no-slip solid wall
      test = 1;
      value = 0.;
    }
  }
  else if (!strcmp(name, "W")) {
    if (1 == facename) {
      test = 1;
      value = 0.;
    }
    else if (2 == facename) {
      test = 1;
      value = 0.;
    }
    else if (3 == facename) {
      test = 1;
      value = 0.;
    }
    else if (4 == facename) {
      test = 1;
      value = 0.;
    }
  }
  else if (!strcmp(name, "P")) {
    if (1 == facename) {
      test = 0;
      value = 0.;
    }
    else if (2 == facename) {
      test = 0;
      value = 0.;
    }
    else if (3 == facename) {
      test = 0;
      value = 0.;
    }
    else if (4 == facename) {
      test = 0;
      value = 0.;
    }
  }
  else if (!strcmp(name, "DX")) {
    if (1 == facename) {     //inflow
      test = 1;
      value = 0.;
    }
    else if (2 == facename) { //outflow
      test = 1;
      value = 0.;
    }
    else if (3 == facename) { // no-slip fluid wall
      test = 0; //0
      value = 0.;
    }
    else if (4 == facename) { // no-slip solid wall
      test = 1;
      value = 0.;
    }
  }
  else if (!strcmp(name, "DY")) {
    if (1 == facename) {     //inflow
      test = 0; // 0
      value = 0.;
    }
    else if (2 == facename) { //outflow
      test = 0; // 0
      value = 0.;
    }
    else if (3 == facename) { // no-slip fluid wall
      test = 1;
      value = 0.;
    }
    else if (4 == facename) { // no-slip solid wall
      test = 1;
      value = 0.;
    }
  }
  else if (!strcmp(name, "DZ")) {
    if (1 == facename) {     //inflow
      test = 1;
      value = 0.;
    }
    else if (2 == facename) { //outflow
      test = 1;
      value = 0.;
    }
    else if (3 == facename) { // no-slip fluid wall
      test = 1;
      value = 0.;
    }
    else if (4 == facename) { // no-slip solid wall
      test = 1;
      value = 0.;
    }
  }

  return test;
}

void AssembleMatrixResFSI(MultiLevelProblem& ml_prob) {

  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;

  //pointers and references
  MonolithicFSINonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<MonolithicFSINonLinearImplicitSystem>("Fluid-Structure-Interaction");
  const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();

  MultiLevelSolution*	 ml_sol	                      = ml_prob._ml_sol;
  Solution*	 mysolution  	                      = ml_sol->GetSolutionLevel(level);

  LinearEquationSolver*  myLinEqSolver	              = my_nnlin_impl_sys._LinSolver[level];
  Mesh*		mymsh		=  ml_prob._ml_msh->GetLevel(level);
  elem*		myel		=  mymsh->el;
  SparseMatrix*	myKK		=  myLinEqSolver->_KK;
  NumericVector* myRES		=  myLinEqSolver->_RES;
  vector <int>&	myKKIndex	=  myLinEqSolver->KKIndex;

  const unsigned dim = mymsh->GetDimension();
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

  // local objects
  vector<double> SolVAR(2 * dim + 1);
  vector<vector<double> > GradSolVAR(2 * dim);

  for (int i = 0; i < 2 * dim; i++) {
    GradSolVAR[i].resize(dim);
  }

  vector<vector<double> > GradSolhatVAR(dim);

  for (int i = 0; i < dim; i++) {
    GradSolhatVAR[i].resize(dim);
  }

  vector <int> metis_node1;
  vector <int> metis_node2;
  vector <bool> solidmark;

  vector <double > phi;
  vector <double > phi_hat;

  vector <double> gradphi;
  vector <double> gradphi_hat;

  vector <double> nablaphi;
  vector <double> nablaphi_hat;

  metis_node1.reserve(max_size);
  metis_node2.reserve(max_size);
  solidmark.reserve(max_size);
  phi.reserve(max_size);
  phi_hat.reserve(max_size);
  gradphi.reserve(max_size * dim);
  gradphi_hat.reserve(max_size * dim);


  nablaphi.reserve(max_size * (3 * (dim - 1)));
  nablaphi_hat.reserve(max_size * (3 * (dim - 1)));

  const double* phi1;

  double Weight = 0.;
  double Weight_nojac = 0.;
  double Weight_hat = 0.;

  vector <vector < double> > vx(dim);
  vector <vector < double> > vx_hat(dim);

  for (int i = 0; i < dim; i++) {
    vx[i].reserve(max_size);
    vx_hat[i].reserve(max_size);
  }

  vector< vector< double > > Rhs(2 * dim + 1);
  vector< vector< vector< double > > > B(2 * dim + 1);

  for (int i = 0; i < 2 * dim + 1; i++) {
    B[i].resize(2 * dim + 1);
  }

  vector< vector< int > > dofsVAR(2 * dim + 1);

  // ------------------------------------------------------------------------
  // Physical parameters
  double rhof	 	= ml_prob.parameters.get<Fluid>("Fluid").get_density();
  double rhos 		= ml_prob.parameters.get<Solid>("Solid").get_density();
  double mu_lame 	= ml_prob.parameters.get<Solid>("Solid").get_lame_shear_modulus();
  double lambda_lame 	= ml_prob.parameters.get<Solid>("Solid").get_lame_lambda();
  double mus		= mu_lame / rhof;
  double IRe 		= ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();
  double lambda		= lambda_lame / rhof;
  double betafsi	= rhos / rhof;
  double betans		= 1.;
  int    solid_model	= ml_prob.parameters.get<Solid>("Solid").get_physical_model();

  //physical quantity
  double Jnp1_hat;
  double Jn_hat;
  double I_bleft;
  double I_e;
  double Cauchy[3][3];
  double Cauchy_old[3][3];
  double tg_stiff_matrix[3][3];
  //initialization C tensor: Saint-Venaint Kirchoff model : solid_model==1;
  const double Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};
  double C_mat[3][3][3][3];

  for (int I = 0; I < 3; ++I) {
    for (int J = 0; J < 3; ++J) {
      for (int K = 0; K < 3; ++K) {
        for (int L = 0; L < 3; ++L) {
          C_mat[I][J][K][L] = 2.*mus * Id2th[I][K] * Id2th[J][L];
        }
      }
    }
  }

  // ale map
  double _lambda_map = 0.;
  double _mu_ale[3] = {1., 1., 1.};

  // gravity
  double _gravity[3] = {0., 0., 0.};

  // newton algorithm
  bool nwtn_alg = false;

  // -----------------------------------------------------------------
  // space discretization parameters
  unsigned order_ind2 = ml_sol->GetSolutionType(ml_sol->GetIndex("U"));
  unsigned order_ind1 = ml_sol->GetSolutionType(ml_sol->GetIndex("P"));

  // mesh and procs
  unsigned nel    = mymsh->GetNumberOfElements();
  unsigned igrid  = mymsh->GetLevel();
  unsigned iproc  = mymsh->processor_id();

  //----------------------------------------------------------------------------------
  //variable-name handling
  const char varname[7][3] = {"DX", "DY", "DZ", "U", "V", "W", "P"};
  vector <unsigned> indexVAR(2 * dim + 1);
  vector <unsigned> indVAR(2 * dim + 1);
  vector <unsigned> SolType(2 * dim + 1);

  for (unsigned ivar = 0; ivar < dim; ivar++) {
    indVAR[ivar] = ml_sol->GetIndex(&varname[ivar][0]);
    indVAR[ivar + dim] = ml_sol->GetIndex(&varname[ivar + 3][0]);
    SolType[ivar] = ml_sol->GetSolutionType(&varname[ivar][0]);
    SolType[ivar + dim] = ml_sol->GetSolutionType(&varname[ivar + 3][0]);
    indexVAR[ivar] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);
    indexVAR[ivar + dim] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar + 3][0]);
  }

  indexVAR[2 * dim] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[6][0]);
  indVAR[2 * dim] = ml_sol->GetIndex(&varname[6][0]);
  SolType[2 * dim] = ml_sol->GetSolutionType(&varname[6][0]);
  //----------------------------------------------------------------------------------

  start_time = clock();

  myKK->zero();

  /// *** element loop ***
  for (int iel = mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc + 1]; iel++) {

    unsigned kel        = iel;
    short unsigned kelt = mymsh->GetElementType(kel);
    unsigned nve        = mymsh->GetElementDofNumber(kel, order_ind2);
    unsigned nve1       = mymsh->GetElementDofNumber(kel, order_ind1);
    int flag_mat        = mymsh->GetElementMaterial(kel);

    //*******************************************************************************************************

    //initialization of everything is in common fluid and solid

    //Rhs
    for (int i = 0; i < 2 * dim; i++) {
      dofsVAR[i].resize(nve);
      Rhs[indexVAR[i]].resize(nve);
      memset(&Rhs[indexVAR[i]][0], 0, nve * sizeof(double));
    }

    dofsVAR[2 * dim].resize(nve1);
    Rhs[indexVAR[2 * dim]].resize(nve1);
    memset(&Rhs[indexVAR[2 * dim]][0], 0, nve1 * sizeof(double));

    //Kinematic relation (solid) and ALE Map (fluid)
    for (int i = 0; i < dim; i++) {
      B[indexVAR[i]][indexVAR[i]].resize(nve * nve);
      memset(&B[indexVAR[i]][indexVAR[i]][0], 0, nve * nve * sizeof(double));
    }

    //Stiffness Matrix (solid)
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++) {
        B[indexVAR[dim + i]][indexVAR[j]].resize(nve * nve);
        memset(&B[indexVAR[dim + i]][indexVAR[j]][0], 0, nve * nve * sizeof(double));

// 	B[indexVAR[dim+i]][indexVAR[dim+j]].resize(nve*nve);
// 	memset(&B[indexVAR[dim+i]][indexVAR[dim+j]][0],0,nve*nve*sizeof(double));
      }
    }

    //Mass Matrix (solid and fluid) and Diffusion Matrix (fluid only)
    for (int i = 0; i < dim; i++) {
      B[indexVAR[dim + i]][indexVAR[dim + i]].resize(nve * nve);
      memset(&B[indexVAR[dim + i]][indexVAR[dim + i]][0], 0, nve * nve * sizeof(double));

      if (nwtn_alg == true) {
        for (int idim2 = 1; idim2 < dim; idim2++) {
          B[indexVAR[dim + i]][indexVAR[dim + (i + idim2) % dim]].resize(nve * nve);
          memset(&B[indexVAR[dim + i]][indexVAR[dim + (i + idim2) % dim]][0], 0, nve * nve * sizeof(double));
        }
      }
    }

    //Pressure gradient (fluid and solid)
    for (int i = 0; i < dim; i++) {
      B[indexVAR[dim + i]][indexVAR[2 * dim]].resize(nve * nve1);
      memset(&B[indexVAR[dim + i]][indexVAR[2 * dim]][0], 0, nve * nve1 * sizeof(double));
    }

    // Pressure Mass Matrix
    B[indexVAR[2 * dim]][indexVAR[2 * dim]].resize(nve1 * nve1);
    memset(&B[indexVAR[2 * dim]][indexVAR[2 * dim]][0], 0, nve1 * nve1 * sizeof(double));

    if (flag_mat == 2) { //initialization for fluid only
      // Fluid Continuity Matrix: divergence of the velocity
      for (int i = 0; i < dim; i++) {
        B[indexVAR[2 * dim]][indexVAR[dim + i]].resize(nve1 * nve);
        memset(&B[indexVAR[2 * dim]][indexVAR[dim + i]][0], 0, nve1 * nve * sizeof(double));
      }
    }
    else { // initialization for solid only
      // Kinematic relation
      for (int i = 0; i < dim; i++) {
        B[indexVAR[i]][indexVAR[dim + i]].resize(nve * nve);
        memset(&B[indexVAR[i]][indexVAR[dim + i]][0], 0, nve * nve * sizeof(double));
      }

      // Solid Continuity Matrix: divergence of the displacemnet
      for (int i = 0; i < dim; i++) {
        B[indexVAR[2 * dim]][indexVAR[i]].resize(nve1 * nve);
        memset(&B[indexVAR[2 * dim]][indexVAR[i]][0], 0, nve1 * nve * sizeof(double));
      }
    }

    // ----------------------------------------------------------------------------------------
    // coordinates, displacement, velocity dofs

    metis_node2.resize(nve);
    metis_node1.resize(nve1);
    solidmark.resize(nve);
    phi.resize(nve);
    phi_hat.resize(nve);
    gradphi.resize(nve * dim);
    gradphi_hat.resize(nve * dim);
    nablaphi.resize(nve * (3 * (dim - 1)));
    nablaphi_hat.resize(nve * (3 * (dim - 1)));

    for (int i = 0; i < dim; i++) {
      vx[i].resize(nve);
      vx_hat[i].resize(nve);
    }

    for (unsigned i = 0; i < nve; i++) {
      // gambit nodes
     
      // dof metis
      unsigned inode_Metis = mymsh->GetSolutionDof(i, kel, 2);
      metis_node2[i] = inode_Metis;

      //unsigned inode_Metis=mymsh->GetSolutionDof(inode,2);
      // flag to know if the node "inode" lays on the fluid-solid interface
      solidmark[i] = mymsh->GetSolidMark(inode_Metis); // to check

      for (int j = 0; j < dim; j++) {
        //Updated coordinates (Moving frame)
        vx[j][i] = (*mymsh->_topology->_Sol[j])(inode_Metis) + (*mysolution->_Sol[indVAR[j]])(inode_Metis);
        //Fixed coordinates (Reference frame)
        vx_hat[j][i] = (*mymsh->_topology->_Sol[j])(inode_Metis);
        // displacement dofs
        dofsVAR[j][i] = myLinEqSolver->GetSystemDof(indVAR[j], indexVAR[j], i, kel);
        // velocity dofs
        dofsVAR[j + dim][i] = myLinEqSolver->GetSystemDof(indVAR[j + dim], indexVAR[j + dim], i, kel);
      }
    }

    // pressure dofs
    for (unsigned i = 0; i < nve1; i++) {
      metis_node1[i] = mymsh->GetSolutionDof(i, kel, SolType[2 * dim]);
      dofsVAR[2 * dim][i] = myLinEqSolver->GetSystemDof(indVAR[2 * dim], indexVAR[2 * dim], i, kel);
    }

    // ----------------------------------------------------------------------------------------


    /// *** Gauss point loop ***
    for (unsigned ig = 0; ig < ml_prob._ml_msh->_finiteElement[kelt][order_ind2]->GetGaussPointNumber(); ig++) {

      // *** get Jacobian and test function and test function derivatives in the moving frame***
      ml_prob._ml_msh->_finiteElement[kelt][order_ind2]->Jacobian(vx, ig, Weight, phi, gradphi, nablaphi);
      ml_prob._ml_msh->_finiteElement[kelt][order_ind2]->Jacobian(vx_hat, ig, Weight_hat, phi_hat, gradphi_hat, nablaphi_hat);
      phi1 = ml_prob._ml_msh->_finiteElement[kelt][order_ind1]->GetPhi(ig);

      if (flag_mat == 2) Weight_nojac = ml_prob._ml_msh->_finiteElement[kelt][order_ind2]->GetGaussWeight(ig);

      // ---------------------------------------------------------------------------
      // displacement and velocity
      for (int i = 0; i < 2 * dim; i++) {
        SolVAR[i] = 0.;

        for (int j = 0; j < dim; j++) {
          GradSolVAR[i][j] = 0.;

          if (i < dim) {
            GradSolhatVAR[i][j] = 0.;
          }
        }

        for (unsigned inode = 0; inode < nve; inode++) {
          unsigned sol_dof = metis_node2[inode];

          double soli = (*mysolution->_Sol[indVAR[i]])(sol_dof);
          SolVAR[i] += phi[inode] * soli;

          for (int j = 0; j < dim; j++) {
            GradSolVAR[i][j] += gradphi[inode * dim + j] * soli;

            if (i < dim) {
              GradSolhatVAR[i][j]   += gradphi_hat[inode * dim + j] * soli;
            }
          }
        }
      }

      // pressure
      SolVAR[2 * dim] = 0.;

      for (unsigned inode = 0; inode < nve1; inode++) {
        double soli = (*mysolution->_Sol[indVAR[2 * dim]])(metis_node1[inode]);
        SolVAR[2 * dim] += phi1[inode] * soli;
      }

      // ---------------------------------------------------------------------------

      //BEGIN FLUID ASSEMBLY ============

      if (flag_mat == 2) {

        //divergence of the velocity
        double div_vel = 0.;
        double div_w = 0.;

        for (int i = 0; i < dim; i++) {
          div_vel += GradSolVAR[dim + i][i];
// 	    div_w+=(GradSolVAR[i][i]-GradSolOldVAR[i][i])*(1./dt);
        }

        { // Laplace operator + adection operator + Mass operator

          const double* gradfi = &gradphi[0];
          const double* fi = &phi[0];

          // *** phi_i loop ***
          for (unsigned i = 0; i < nve; i++, gradfi += dim, fi++) {

            //BEGIN RESIDUALS A block ===========================
            double LapmapVAR[3] = {0., 0., 0.};

            for (int idim = 0; idim < dim; idim++) {
              for (int idim2 = 0; idim2 < dim; idim2++) {
                LapmapVAR[idim] += (_mu_ale[idim2] * gradphi_hat[i * dim + idim2] * GradSolhatVAR[idim][idim2]);
              }
            }

            // Residual ALE equations
            for (int idim = 0; idim < dim; idim++) {
              Rhs[indexVAR[idim]][i] += (!solidmark[i]) * (-LapmapVAR[idim] * Weight_nojac);
            }

            //-------------------------------------------------------------------------------

            double LapvelVAR[3] = {0., 0., 0.};
            double AdvaleVAR[3] = {0., 0., 0.};

            for (int idim = 0.; idim < dim; idim++) {
              for (int idim2 = 0.; idim2 < dim; idim2++) {
                LapvelVAR[idim] += gradphi[i * dim + idim2] * GradSolVAR[dim + idim][idim2];
                AdvaleVAR[idim] += SolVAR[dim + idim2] * GradSolVAR[dim + idim][idim2] * phi[i];;
              }
            }

            // Residual Momentum equations
            for (int idim = 0; idim < dim; idim++) {
              Rhs[indexVAR[dim + idim]][i] += (
                                                -AdvaleVAR[idim] * Weight
// 					      -0.5*div_vel*SolVAR[dim+idim]*phi[i]*(Weight)
// 					      +div_w*SolVAR[dim+idim]*phi[i]*(Weight)
                                                - IRe * LapvelVAR[idim] * Weight
                                                + SolVAR[2 * dim] * gradphi[i * dim + idim] * Weight
                                              );
            }

            //END RESIDUALS A block ===========================

            const double* gradfj = &gradphi[0];
            const double* fj = &phi[0];

            //  *** phi_j loop ***
            for (unsigned j = 0; j < nve; j++, gradfj += dim, fj++) {

              //Laplacian
              double Lap = 0.;

              for (int idim = 0; idim < dim; idim++) {
                Lap += (*(gradfi + idim)) * (*(gradfj + idim));
              }

              double LapXweight = Lap * Weight;

              //advection term I
              double Adv1 = 0.;
              double Adv2 = ((*fi)) * ((*fj)) * Weight;

              for (int idim = 0; idim < dim; idim++) {
                Adv1 += SolVAR[dim + idim] * (*(gradfj + idim)) * (*(fi)) * Weight;
              }

              double div_stab = 0.*0.5 * div_vel * ((*fi)) * ((*fj)) * (Weight);

              double div_ale = 0.*div_w * ((*fi)) * ((*fj)) * (Weight);

              for (int idim = 0; idim < dim; idim++) {
                B[indexVAR[dim + idim]][indexVAR[dim + idim]][i * nve + j] += IRe * LapXweight + Adv1 + div_stab - div_ale;

                // Advection term II
                if (nwtn_alg == true) {
                  B[indexVAR[dim + idim]][indexVAR[dim + idim]][i * nve + j]  += Adv2 * GradSolVAR[dim + idim][idim];

                  for (unsigned idim2 = 1; idim2 < dim; idim2++) {
                    B[indexVAR[dim + idim]][indexVAR[dim + (idim + idim2) % dim]][i * nve + j] += Adv2 * GradSolVAR[dim + idim][(idim + idim2) % dim];
                  }
                }
              }


              double Lap_ale = 0.;

              for (int idim = 0; idim < dim; idim++) {
                Lap_ale += _mu_ale[idim] * (*(gradfi + idim)) * (*(gradfj + idim));
              }

              // Laplacian ALE map
              for (int idim = 0; idim < dim; idim++) {
                B[indexVAR[0 + idim]][indexVAR[idim]][i * nve + j] += (!solidmark[i]) * Lap_ale * Weight_nojac;
              }

            } // end phi_j loop
          } // end phi loop
        } // end A + Bt

        { //Gradient of Pressure operator

          const double* gradfi = &gradphi[0];
          const double* fi = phi1;

          /// *** phi_i loop ***
          for (unsigned i = 0; i < nve; i++, gradfi += dim, fi++) {

            const double* fj = phi1;

            /// *** phi_j loop ***
            for (unsigned j = 0; j < nve1; j++, fj++) {
              for (int idim = 0; idim < dim; idim++) {
                B[indexVAR[dim + idim]][indexVAR[2 * dim]][i * nve1 + j] -= ((*(gradfi + idim)) * (*fj)) * Weight;
              }
            } // end phi_j loop
          } // end phi_i loop
        } // End Bt

        { // Divergence of the Velocity operator
          const double* fi = phi1;

          // *** phi_i loop ***
          for (unsigned i = 0; i < nve1; i++, fi++) {

            //BEGIN RESIDUALS B block ===========================
            Rhs[indexVAR[2 * dim]][i] += -(-((*fi)) * div_vel) * Weight;
            //END RESIDUALS  B block ===========================

            const double* gradfj = &gradphi[0];

            // *** phi_j loop ***
            for (unsigned j = 0; j < nve; j++, gradfj += dim) {
              for (int idim = 0; idim < dim; idim++) {
                B[indexVAR[2 * dim]][indexVAR[idim + dim]][i * nve + j] -= ((*fi) * (*(gradfj + idim))) * Weight;
              }
            }
          }
        }
      }
      //END FLUID ASSEMBLY ============
      //*******************************************************************************************************
      //BEGIN SOLID ASSEMBLY ============

      else {

        //------------------------------------------------------------------------------------------------------------
        if (solid_model == 0) {
          double e[3][3];
          double e_old[3][3];

          //computation of the stress tensor
          for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
              e[i][j] = 0.5 * (GradSolhatVAR[i][j] + GradSolhatVAR[j][i]);
            }
          }

          I_e = 0;

          for (int i = 0; i < dim; i++) {
            I_e += e[i][i];
          }

          for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
              //incompressible
              Cauchy[i][j] = 2 * mus * e[i][j];
              Cauchy_old[i][j] = 2 * mus * e_old[i][j];
            }
          }
        }

        else if (solid_model == 1) {
          double F[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
          double b_left[3][3];

          for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
              F[i][j] += GradSolhatVAR[i][j];
            }
          }

          Jnp1_hat =  F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
                      - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];

          // computation of the the three deformation tensor b
          for (int I = 0; I < 3; ++I) {
            for (int J = 0; J < 3; ++J) {
              b_left[I][J] = 0.;

              for (int K = 0; K < 3; ++K) {
                //left Cauchy-Green deformation tensor or Finger tensor (b = F*F^T)
                b_left[I][J] += F[I][K] * F[J][K];
              }

              Cauchy[I][J] = (mus / Jnp1_hat) * (b_left[I][J] - Id2th[I][J]);
            }
          }

          I_bleft = b_left[0][0] + b_left[1][1] + b_left[2][2];

          //compressible case
          //             for (int ii=0; ii<3; ++ii) {
          //               for (int jj=0; jj<3; ++jj) {
          //                 //Cauchy stress tensor
          //                 Cauchy[ii][jj] = (_mus/Jnp1_hat)*(b_left[ii][jj] - Id2th[ii][jj]) + (_lambda/Jnp1_hat)*log(Jnp1_hat)*Id2th[ii][jj];
          //                 for (int k=0; k<3; ++k) {
          //                   for (int l=0; l<3; ++l) {
          //                   }
          //                 }
          //               }
          //             }

          //for the incompressible(nearly incompressible) case
          for (int ii = 0; ii < 3; ++ii) {
            for (int jj = 0; jj < 3; ++jj) {
              for (int kk = 0; kk < 3; ++kk) {
                for (int ll = 0; ll < 3; ++ll) {
                  C_mat[ii][jj][kk][ll] = 2.*mus * pow(Jnp1_hat, -1.6666666666666) * (
                                            0.333333333333 * I_bleft * Id2th[ii][kk] * Id2th[jj][ll]        //1/3*I_c*i
                                            // 	                        +0.111111111111*I_C*Id2th[i][j]*Id2th[k][l]             //1/9*I_b*IxI
                                            // 				-0.333333333333*b_left[i][j]*Id2th[k][l]                //-1/3*b*I
                                            // 				-0.333333333333*Id2th[i][j]*b_left[k][l]                //-1/3*b*I
                                          )
                                          - SolVAR[2 * dim] * (Id2th[ii][jj] * Id2th[kk][ll] - 2.*Id2th[ii][kk] * Id2th[jj][ll]); // -p(IxI-2i)
                }
              }
            }
          }

          //Old deformation gradient
          double F_old[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

          for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
              F_old[i][j] += GradSolhatVAR[i][j];
            }
          }

          Jn_hat =  F_old[0][0] * F_old[1][1] * F_old[2][2] + F_old[0][1] * F_old[1][2] * F_old[2][0] + F_old[0][2] * F_old[1][0] * F_old[2][1]
                    - F_old[2][0] * F_old[1][1] * F_old[0][2] - F_old[2][1] * F_old[1][2] * F_old[0][0] - F_old[2][2] * F_old[1][0] * F_old[0][1] ;

          // computation of the the three deformation tensor b
          for (int I = 0; I < 3; ++I) {
            for (int J = 0; J < 3; ++J) {
              b_left[I][J] = 0.;

              for (int k = 0; k < 3; ++k) {
                //left Cauchy-Green deformation tensor or F_oldinger tensor (b = F_old*F_old^T)
                b_left[I][J] += F_old[I][k] * F_old[J][k];
              }

              Cauchy_old[I][J] = (mus / Jn_hat) * (b_left[I][J] - Id2th[I][J]);
            }
          }
        }

        //----------------------------------------------------------------------------------------------------------------------------

        /////////////

        ///Mass + Stiffness operator
        {
          const double* gradfi = &gradphi[0];
          const double* fi = &phi[0];

          /// *** phi_i loop ***
          for (unsigned i = 0; i < nve; i++, gradfi += dim, fi++) {


            //BEGIN RESIDUALS A + Bt block ===========================

            // Residual ALE equations
            for (int idim = 0; idim < dim; idim++) {
              Rhs[indexVAR[idim]][i] += (-phi[i] * (-SolVAR[dim + idim])) * Weight_hat;
            }

            double CauchyDIR[3] = {0., 0., 0.};

            for (int idim = 0.; idim < dim; idim++) {
              for (int idim2 = 0.; idim2 < dim; idim2++) {
                CauchyDIR[idim] += gradphi[i * dim + idim2] * Cauchy[idim][idim2];
              }
            }

            // Residual Momentum equations
            for (int idim = 0; idim < dim; idim++) {
              Rhs[indexVAR[dim + idim]][i] += (
                                                phi[i] * _gravity[idim] * Weight_hat
                                                - CauchyDIR[idim] * Weight
                                                + SolVAR[2 * dim] * gradphi[i * dim + idim] * Weight
                                              );

            }

            //---------------------------------------------------------------------------------------------------------------------------------

            //END RESIDUALS A + Bt block ===========================

            const double* gradfj = &gradphi[0];
            const double* fj = &phi[0];

            // *** phi_j loop ***
            for (unsigned j = 0; j < nve; j++, gradfj += dim, fj++) {

              for (int icount = 0; icount < dim; ++icount) {
                for (int jcount = 0; jcount < dim; ++jcount) {
                  tg_stiff_matrix[icount][jcount] = 0.;

                  for (int kcount = 0; kcount < dim; ++kcount) {
                    for (int lcount = 0; lcount < dim; ++lcount) {
                      tg_stiff_matrix[icount][jcount] += (*(gradfi + kcount)) * 0.25 * (
                                                           C_mat[icount][kcount][jcount][lcount] + C_mat[icount][kcount][lcount][jcount]
                                                           + C_mat[kcount][icount][jcount][lcount] + C_mat[kcount][icount][lcount][jcount]
                                                         ) * (*(gradfj + lcount));
                    }
                  }
                }
              }

              //geometric tangent stiffness matrix
              double geom_tg_stiff_matrx = 0.;

              for (int kcount = 0; kcount < dim; ++kcount) {
                for (int lcount = 0; lcount < dim; ++lcount) {
                  geom_tg_stiff_matrx += (*(gradfi + kcount)) * Cauchy[kcount][lcount] * (*(gradfj + lcount));
                }
              }

              /// Stiffness operator -- Elasticity equation (Linear or not)
              for (int idim = 0; idim < dim; idim++) {
                B[indexVAR[dim + idim]][indexVAR[0 + idim]][i * nve + j] += geom_tg_stiff_matrx * Weight;

                for (int idim2 = 0; idim2 < dim; idim2++) {
                  B[indexVAR[dim + idim]][indexVAR[0 + idim2]][i * nve + j] += tg_stiff_matrix[0 + idim][0 + idim2] * Weight;
                }
              }

              /// Kinematic equation v = du/dt --> In the steady state we write \deltau^n+1 - \deltav^n+1 = v - 0
              //
              for (int idim = 0; idim < dim; idim++) {
                // -(v_n+1,eta)
                B[indexVAR[0 + idim]][indexVAR[dim + idim]][i * nve + j] -= (*(fi)) * (*(fj)) * Weight_hat;
                //  (u_n+1,eta)
                B[indexVAR[0 + idim]][indexVAR[0 + idim]][i * nve + j] += (*(fi)) * (*(fj)) * Weight_hat;
              }
            }
          }
        }
        ////////////
        { ///Gradient of Pressure
          const double* gradfi = &gradphi[0];

          // *** phi_i loop ***
          for (unsigned i = 0; i < nve; i++, gradfi += dim) {
            const double* fj = phi1;

            // *** phi_j loop ***
            for (unsigned j = 0; j < nve1; j++, fj++) {
              for (int idim = 0; idim < dim; idim++) {
                B[indexVAR[dim + idim]][indexVAR[2 * dim]][i * nve1 + j] -= ((*(gradfi + idim)) * (*fj)) * Weight;
              }
            }
          }
        }
        ////////////
        { ///Divergence of the Displacement
          const double* fi = phi1;

          // *** phi_i loop ***
          for (unsigned i = 0; i < nve1; i++, fi++) {

            //BEGIN RESIDUALS B block ===========================

            if (solid_model == 0) {
              Rhs[indexVAR[2 * dim]][i] += -(-((*fi)) * (I_e + (1. / lambda) * SolVAR[2 * dim])) * Weight_hat;
            }
            else if (solid_model == 1) {
              Rhs[indexVAR[2 * dim]][i] += -(-((*fi)) * (log(Jnp1_hat) / Jnp1_hat + (1. / lambda) * SolVAR[2 * dim])) * Weight_hat;
            }

            //END RESIDUALS B block ===========================

            const double* gradfj = &gradphi[0];

            // *** phi_j loop ***
            for (unsigned j = 0; j < nve; j++, gradfj += dim) {
              for (int idim = 0; idim < dim; idim++) {
                B[indexVAR[2 * dim]][indexVAR[idim]][i * nve + j] -= ((*fi) * (*(gradfj + idim))) * Weight;
              }
            }
          }
        }
        //  /////////////
        { ///Pressure Mass term
          const double* fi = phi1;

          // *** phi_i loop ***
          for (unsigned i = 0; i < nve1; i++, fi++) {
            const double* fj = phi1;

            // *** phi_j loop ***
            for (unsigned j = 0; j < nve1; j++, fj++) {
              B[indexVAR[2 * dim]][indexVAR[2 * dim]][i * nve1 + j] -= (1. / lambda) * ((*fi) * (*fj)) * Weight_hat;
            }
          }
        }  //end pressure mass term
        //---------------------------------------------------------------------------------------------------------------------------------
      }

      //END SOLID ASSEMBLY ============
    }


    //BEGIN local to global assembly
    // ALE mapping
    for (int i = 0; i < dim; i++) {
      myRES->add_vector_blocked(Rhs[indexVAR[i]], dofsVAR[i]);
      myKK ->add_matrix_blocked(B[indexVAR[i]][indexVAR[i]], dofsVAR[i], dofsVAR[i]);

      if (flag_mat != 2) { //Solid only
        myKK->add_matrix_blocked(B[indexVAR[i]][indexVAR[i + dim]], dofsVAR[i], dofsVAR[i + dim]);
      }
    }

    // Momentum equation
    for (int i = 0; i < dim; i++) {
      myRES->add_vector_blocked(Rhs[indexVAR[dim + i]], dofsVAR[dim + i]);
      myKK->add_matrix_blocked(B[indexVAR[dim + i]][indexVAR[dim + i]], dofsVAR[dim + i], dofsVAR[dim + i]);

      if (nwtn_alg == true) {
        for (unsigned idim2 = 1; idim2 < dim; idim2++) {
          myKK->add_matrix_blocked(B[indexVAR[dim + i]][indexVAR[dim + (i + idim2) % dim]], dofsVAR[dim + i], dofsVAR[dim + (i + idim2) % dim]);
        }
      }

      myKK->add_matrix_blocked(B[indexVAR[dim + i]][indexVAR[2 * dim]], dofsVAR[dim + i], dofsVAR[2 * dim]);

      for (int j = 0; j < dim; j++) {
        myKK->add_matrix_blocked(B[indexVAR[dim + i]][indexVAR[j]], dofsVAR[dim + i], dofsVAR[j]);
      }
    }

    //P-continuity equation
    myRES->add_vector_blocked(Rhs[indexVAR[2 * dim]], dofsVAR[2 * dim]);

    for (int i = 0; i < dim; i++) {
      if (flag_mat == 2) { //Fluid only
        myKK->add_matrix_blocked(B[indexVAR[2 * dim]][indexVAR[dim + i]], dofsVAR[2 * dim], dofsVAR[dim + i]);
      }
      else { //Solid only
        myKK->add_matrix_blocked(B[indexVAR[2 * dim]][indexVAR[i]], dofsVAR[2 * dim], dofsVAR[i]);
      }
    }

    myKK->add_matrix_blocked(B[indexVAR[2 * dim]][indexVAR[2 * dim]], dofsVAR[2 * dim], dofsVAR[2 * dim]);
    //END local to global assembly

  } //end list of elements loop

  // close residual vector and matrix

  myKK->close();
  myRES->close();

  // *************************************
  end_time = clock();
  AssemblyTime += (end_time - start_time);
  // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

}
