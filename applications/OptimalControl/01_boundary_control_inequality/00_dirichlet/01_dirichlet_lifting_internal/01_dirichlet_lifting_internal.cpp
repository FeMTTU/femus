#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
// #include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"
#include "NumericVector.hpp"
#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"



#include "../../../param.hpp"


using namespace femus;



 //Unknown definition  ==================
 const std::vector< Unknown >  provide_list_of_unknowns(const unsigned int dimension) {
     
     
  std::vector< FEFamily > feFamily;
  std::vector< FEOrder >   feOrder;

                        feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE);
 
                        feOrder.push_back(/*FIRST*/SECOND);  //same
                        feOrder.push_back(/*FIRST*/SECOND);  //same
                        feOrder.push_back(/*FIRST*/SECOND);
                        feOrder.push_back(/*FIRST*/SECOND);  //same
 

  assert( feFamily.size() == feOrder.size() );
 
 std::vector< Unknown >  unknowns(feFamily.size());

   unknowns[0]._name      = "state";
   unknowns[1]._name      = "control";
   unknowns[2]._name      = "adjoint";
   unknowns[3]._name      = "mu";

   unknowns[0]._is_sparse = true;
   unknowns[1]._is_sparse = true;
   unknowns[2]._is_sparse = true;
   unknowns[3]._is_sparse = true;
   
     for (unsigned int u = 0; u < unknowns.size(); u++) {
         
              unknowns[u]._fe_family  = feFamily[u];
              unknowns[u]._fe_order   = feOrder[u];
              unknowns[u]._time_order = 0;
              unknowns[u]._is_pde_unknown = true;
              
     }
 
 
   return unknowns;
     
}





double Solution_set_initial_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {

    double value = 0.;

    if(!strcmp(name,"state")) {
        value = 0.;
    }
    else if(!strcmp(name,"control")) {
        value = 0.;
    }
    else if(!strcmp(name,"adjoint")) {
        value = 0.;
    }
    else if(!strcmp(name,"mu")) {
        value = 0.;
    }
    
    
    else if(!strcmp(name,"TargReg")) {
        value = ctrl::ElementTargetFlag(x);
    }
    else if(!strcmp(name,"ContReg")) {
        value = ctrl::ControlDomainFlag_internal_restriction(x);
    }
    else if(!strcmp(name,"act_flag")) {
        value = 0.;
    }


    return value;
}



bool Solution_set_boundary_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int faceName, const double time) {

  bool dirichlet = true; //dirichlet
  value = 0.;
  
  if(!strcmp(name,"control")) {
      value = 0.;
    if (faceName == FACE_FOR_CONTROL) {
        if (x[ ctrl::axis_direction_Gamma_control(faceName) ] > GAMMA_CONTROL_LOWER - 1.e-5 && x[ ctrl::axis_direction_Gamma_control(faceName) ] < GAMMA_CONTROL_UPPER + 1.e-5)    
            dirichlet = false;
    }
  }
  
  
  //MU
  if(!strcmp(name,"mu")) {
//       value = 0.;
//   if (faceName == FACE_FOR_CONTROL)
    dirichlet = false;
  }
  
  return dirichlet;
}



void AssembleLiftRestrProblem(MultiLevelProblem& ml_prob);



int main(int argc, char** args) {

  // ======= Init ========================
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
  // ======= Problem ========================
  MultiLevelProblem ml_prob;
  
  // ======= Files - BEGIN  ========================
  Files files; 
        files.CheckIODirectories(use_output_time_folder);
        files.RedirectCout(redirect_cout_to_file);

  // ======= Problem, Files ========================
  ml_prob.SetFilesHandler(&files);
  // ======= Files - END  ========================

  // ======= Problem, Quad Rule - BEGIN  ========================
    std::vector< std::string > fe_quad_rule_vec;
  fe_quad_rule_vec.push_back("seventh");
  
  ml_prob.SetQuadratureRuleAllGeomElemsMultiple(fe_quad_rule_vec);
  ml_prob.set_all_abstract_fe_multiple();
  // ======= Problem, Quad Rule - END  ========================

  
  // ======= Mesh, Coarse reading - BEGIN ==================
  MultiLevelMesh ml_mesh;
   
//   std::string input_file = "parametric_square_1x1.med";
//   std::string input_file = "square_parametric.med";
//   std::string input_file = "Mesh_3_groups_with_bdry_nodes.med";
  std::string input_file = "Mesh_3_groups_with_bdry_nodes_coarser.med";
  std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << input_file;
  const std::string infile = mystream.str();
  const double Lref = 1.;

  const bool read_groups = true;
  const bool read_boundary_groups = true;
    
  ml_mesh.ReadCoarseMeshFileReadingBeforePartitioning(infile.c_str(), Lref, read_groups, read_boundary_groups);
    
  ml_mesh.GetLevelZero(0)->build_dofmap_all_fe_families_and_elem_and_node_structures();
 

  ml_mesh.BuildFETypesBasedOnExistingCoarseMeshGeomElements();
  
  ml_mesh.PrepareNewLevelsForRefinement();


  ml_mesh.InitializeQuadratureWithFEEvalsOnExistingCoarseMeshGeomElements(fe_quad_rule_vec[0].c_str()); ///@todo keep it only for compatibility with old ElemType, because of its destructor 
  // ======= Mesh, Coarse reading - END ==================

  
  // ======= Mesh: Refinement - BEGIN ==================
  unsigned numberOfUniformLevels = N_UNIFORM_LEVELS;
  unsigned numberOfSelectiveLevels = 0;
  const unsigned erased_levels = N_ERASED_LEVELS;
  
  ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  // ======= Mesh: Refinement - END ==================

  // ======= Mesh: COARSE ERASING - BEGIN  ========================
  ml_mesh.EraseCoarseLevels(erased_levels/*numberOfUniformLevels - 1*/);
  ml_mesh.PrintInfo();
  // ======= Mesh: COARSE ERASING - END  ========================

  // ======= Solution - BEGIN ==================
  MultiLevelSolution ml_sol(&ml_mesh);

  ml_sol.SetWriter(VTK);
  ml_sol.GetWriter()->SetDebugOutput(true);
  
 // ======= Problem, Mesh and Solution  ==================
 ml_prob.SetMultiLevelMeshAndSolution(& ml_sol);
  // ======= Solution - END ==================


 // ======= Solutions that are Unknowns - BEGIN ==================
  std::vector< Unknown > unknowns = provide_list_of_unknowns( ml_mesh.GetDimension() );
 
  for (unsigned int u = 0; u < unknowns.size(); u++) { ml_sol.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order, unknowns[u]._time_order, unknowns[u]._is_pde_unknown); }
   
 for (unsigned int u = 0; u < unknowns.size(); u++)  { ml_sol.Initialize(unknowns[u]._name.c_str(), Solution_set_initial_conditions, & ml_prob); }
  
  ml_sol.AttachSetBoundaryConditionFunction(Solution_set_boundary_conditions);
   for (unsigned int u = 0; u < unknowns.size(); u++)  {  ml_sol.GenerateBdc(unknowns[u]._name.c_str(), (unknowns[u]._time_order == 0) ? "Steady" : "Time_dependent", & ml_prob);  }
 // ======= Solutions that are Unknowns - END ==================

 
  // ======= Solutions that are not Unknowns - BEGIN  ==================
  ml_sol.AddSolution("TargReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  ml_sol.Initialize("TargReg",     Solution_set_initial_conditions, & ml_prob);

  ml_sol.AddSolution("ContReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  ml_sol.Initialize("ContReg",     Solution_set_initial_conditions, & ml_prob);

  // ******** active flag - BEGIN 
  //MU
  const unsigned int  n_components_ctrl = 1;
  
  const unsigned int act_set_fake_time_dep_flag = 2;  //this is needed to be able to use _SolOld
  const bool         act_flag_is_an_unknown_of_a_pde = false;

  unsigned int index_control = 0;
    for (unsigned int u = 0; u < unknowns.size(); u++) {
        if ( !(unknowns[u]._name.compare("control")) ) index_control = u;
    }
   std::vector<std::string> act_set_flag_name(n_components_ctrl);
   act_set_flag_name[0] = "act_flag";

   ml_sol.AddSolution(act_set_flag_name[0].c_str(), unknowns[index_control]._fe_family, unknowns[index_control]._fe_order, act_set_fake_time_dep_flag, act_flag_is_an_unknown_of_a_pde);               
   ml_sol.Initialize(act_set_flag_name[0].c_str(), Solution_set_initial_conditions, & ml_prob);
  // ******** active flag - END 

  // ======= Solutions that are not Unknowns - END  ==================

  
  //== Solution: CHECK SOLUTION FE TYPES between Unknowns and Not Unknowns - BEGIN  ==
  if ( ml_sol.GetSolutionType("control") != ml_sol.GetSolutionType(act_set_flag_name[0].c_str())) abort();
  //== Solution: CHECK SOLUTION FE TYPES between Unknowns and Not Unknowns - END ==
  

  // ======= Problem, System - BEGIN ========================
  NonLinearImplicitSystemWithPrimalDualActiveSetMethod & system_opt = ml_prob.add_system < NonLinearImplicitSystemWithPrimalDualActiveSetMethod > ("LiftRestr"); //MU
  
  system_opt.SetActiveSetFlagName(act_set_flag_name); //MU

  for (unsigned int u = 0; u < unknowns.size(); u++)  { 
  system_opt.AddSolutionToSystemPDE(unknowns[u]._name.c_str());
  }
  
  
  system_opt.SetAssembleFunction(AssembleLiftRestrProblem);
  
// *****************
  const unsigned n_components_state = n_components_ctrl;
  std::vector<std::string> state_vars(n_components_state);  state_vars[0] = "state";
  std::vector<std::string> ctrl_vars(n_components_ctrl);     ctrl_vars[0] = "control";
  
  system_opt.set_state_vars(state_vars);
  system_opt.set_ctrl_vars(ctrl_vars);
  
  system_opt.SetDebugNonlinear(true);
  system_opt.SetDebugFunction(ctrl::compute_cost_functional_regularization_lifting_internal);
// *****************
//   system_opt.SetMaxNumberOfNonLinearIterations(2);

  // initialize and solve the system
  system_opt.init();
  
  system_opt.MGsolve();
//   system.assemble_call_before_boundary_conditions(1);
  // ======= Problem, System  - END ========================
  
  // ======= Print - BEGIN  ========================
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("all");
  ml_sol.GetWriter()->Write(files.GetOutputPath(), "biquadratic", variablesToBePrinted);
  // ======= Print - END  ========================

  return 0;
  
}


void AssembleLiftRestrProblem(MultiLevelProblem& ml_prob) {

  NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystemWithPrimalDualActiveSetMethod> ("LiftRestr");   // pointer to the linear implicit system named "LiftRestr"
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*             JAC = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

  unsigned    iproc = msh->processor_id();

  constexpr bool print_algebra_global = false;
  constexpr bool print_algebra_local = false;
  
  //=============== Geometry ========================================
  unsigned solType_coords = FE_DOMAIN;
 
  CurrentElem < double > geom_element_iel(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
   //***************************************************  
  
  // quadratures ********************************
  double weight;
  

 //********************* state *********************** 
 //***************************************************  
  vector <double> phi_u; 
  vector <double> phi_u_x;

  phi_u.reserve(max_size);
  phi_u_x.reserve(max_size * dim);
  
 
  unsigned solIndex_u =  solIndex_u = ml_sol->GetIndex("state");
  unsigned solType_u = ml_sol->GetSolutionType(solIndex_u);

  unsigned solPdeIndex_u  =  solPdeIndex_u = mlPdeSys->GetSolPdeIndex("state");

  vector < double >  sol_u;
  sol_u.reserve(max_size);
  vector< int > l2GMap_u;
  l2GMap_u.reserve(max_size);
 //***************************************************  
 //***************************************************  

  
 //******************** control ********************** 
 //***************************************************   
  vector <double> phi_ctrl;
  vector <double> phi_ctrl_x;

  phi_ctrl.reserve(max_size);
  phi_ctrl_x.reserve(max_size * dim);
  
  unsigned solIndex_ctrl;
  solIndex_ctrl = ml_sol->GetIndex("control");
  unsigned solType_ctrl = ml_sol->GetSolutionType(solIndex_ctrl);

  unsigned solPdeIndex_ctrl;
  solPdeIndex_ctrl = mlPdeSys->GetSolPdeIndex("control");

  vector < double >  sol_ctrl;
  sol_ctrl.reserve(max_size);
  vector< int > l2GMap_ctrl;
  l2GMap_ctrl.reserve(max_size);
 //***************************************************  
 //***************************************************  
  
  
 //********************* adjoint ********************* 
 //***************************************************  
  vector <double> phi_adj;
  vector <double> phi_adj_x;

  phi_adj.reserve(max_size);
  phi_adj_x.reserve(max_size * dim);
 
  
  unsigned solIndex_adj;
  solIndex_adj = ml_sol->GetIndex("adjoint");
  unsigned solType_adj = ml_sol->GetSolutionType(solIndex_adj);

  unsigned solPdeIndex_adj;
  solPdeIndex_adj = mlPdeSys->GetSolPdeIndex("adjoint");

  vector < double >  sol_adj;
    sol_adj.reserve(max_size);
  vector< int > l2GMap_adj;
    l2GMap_adj.reserve(max_size);
 //***************************************************  
 //***************************************************  

 //****************** mu ******************************  
 //***************************************************  
  vector <double> phi_mu;
  vector <double> phi_mu_x;

  phi_mu.reserve(max_size);
  phi_mu_x.reserve(max_size * dim);

  
  unsigned solIndex_mu;
  solIndex_mu = ml_sol->GetIndex("mu");
   
  unsigned solPdeIndex_mu;
  solPdeIndex_mu = mlPdeSys->GetSolPdeIndex("mu");
  
  unsigned solType_mu = ml_sol->GetSolutionType(solIndex_mu);
  vector < double >  sol_mu;   sol_mu.reserve(max_size);
  vector < int > l2GMap_mu;   l2GMap_mu.reserve(max_size);

  
    const unsigned int n_components_ctrl = 1;
    const unsigned int first_loc_comp_ctrl = 0;

  //************** variables for ineq constraints: act flag ****************************   
    std::vector <unsigned int> solIndex_act_flag_sol(n_components_ctrl);
  
  ctrl_inequality::store_act_flag_in_old(mlPdeSys, ml_sol, sol, solIndex_act_flag_sol);
    
  
  
  //********* variables for ineq constraints *****************
     std::vector <std::vector < double/*int*/ > > sol_actflag(n_components_ctrl);    //flag for active set
     std::vector <std::vector < double > > ctrl_lower(n_components_ctrl);  
     std::vector <std::vector < double > > ctrl_upper(n_components_ctrl);  
  
     for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) {
          sol_actflag[kdim].reserve(max_size);
           ctrl_lower[kdim].reserve(max_size);
           ctrl_upper[kdim].reserve(max_size);
     }
     
  const int ineq_flag = INEQ_FLAG;
  const double c_compl = C_COMPL;
  //***************************************************  

 //***************************************************  
 //********* WHOLE SET OF VARIABLES ****************** 

  const int n_unknowns =  mlPdeSys->GetSolPdeIndex().size();
 
  vector< int > l2GMap_AllVars;
  l2GMap_AllVars.reserve(n_unknowns*max_size);
  
  vector< double > Res;
  Res.reserve(n_unknowns*max_size);

  vector < double > Jac;
  Jac.reserve( n_unknowns*max_size * n_unknowns*max_size);
  

//***************************************************
    enum Pos_in_matrix {pos_mat_state = 0, pos_mat_ctrl, pos_mat_adj, pos_mat_mu}; //these are known at compile-time 
                    ///@todo these are the positions in the MlSol object or in the Matrix? I'd say the matrix, but we have to check where we use it...

                    
    assert(pos_mat_state   == mlPdeSys->GetSolPdeIndex("state"));
    assert(pos_mat_ctrl    == mlPdeSys->GetSolPdeIndex("control"));
    assert(pos_mat_adj     == mlPdeSys->GetSolPdeIndex("adjoint"));
    assert(pos_mat_mu      == mlPdeSys->GetSolPdeIndex("mu"));

 std::vector < unsigned >   pos_ctrl_in_mat(n_components_ctrl);  pos_ctrl_in_mat[0] = pos_mat_ctrl;
 std::vector < unsigned >   pos_mu_in_mat(n_components_ctrl);    pos_mu_in_mat[0] = pos_mat_mu;
    
    //***************************************************

  vector < std::string > Solname(n_unknowns);
  Solname[0] = "state";
  Solname[1] = "control";
  Solname[2] = "adjoint";
  Solname[3] = "mu";
  

 //***************************************************
   enum Pos_in_Sol {pos_sol_state = 0, pos_sol_ctrl, pos_sol_adj, pos_sol_mu, pos_sol_targreg, pos_sol_contreg, pos_sol_actflag}; //these are known at compile-time 

        assert(pos_sol_actflag == solIndex_act_flag_sol[0]);
//***************************************************
    
  


  vector < unsigned > SolPdeIndex(n_unknowns);
  vector < unsigned > SolIndex(n_unknowns);  
  vector < unsigned > SolFEType(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolPdeIndex[ivar]	= mlPdeSys->GetSolPdeIndex(Solname[ivar].c_str());
    SolIndex[ivar]	= ml_sol->GetIndex        (Solname[ivar].c_str());
    SolFEType[ivar]	= ml_sol->GetSolutionType(SolIndex[ivar]);
  }

    vector < unsigned > Sol_n_el_dofs_Mat_vol(n_unknowns);
 //***************************************************  

//***************************************************
    vector < vector < double > >  sol_eldofs_Mat(n_unknowns);  //should have Mat order
    for(int k = 0; k < n_unknowns; k++) {        sol_eldofs_Mat[k].reserve(max_size);    }


    //----------- quantities (at dof objects) ------------------------------
    std::vector< int >       L2G_dofmap_Mat_AllVars;
    L2G_dofmap_Mat_AllVars.reserve( n_unknowns * max_size );
    vector < vector < int > >     L2G_dofmap_Mat(n_unknowns);     //should have Mat order
    for(int i = 0; i < n_unknowns; i++) {
        L2G_dofmap_Mat[i].reserve(max_size);
    }

    
    
 //********************* DATA ************************ 
  double u_des = ctrl::DesiredTarget();
  double alpha = ALPHA_CTRL_VOL;
  double beta  = BETA_CTRL_VOL;
  double penalty_outside_control_domain = PENALTY_OUTSIDE_CONTROL_DOMAIN;         // penalty for zero control outside
 //***************************************************  

 
 
      // ***************** ADD PART - BEGIN *******************
  RES->zero();
  if (assembleMatrix)  JAC->zero();

    

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
    const short unsigned ielGeom = geom_element_iel.geom_type();
    
   geom_element_iel.set_elem_center_3d(iel, solType_coords);

  //************* set target domain flag **************
   int target_flag = 0;
   target_flag = ctrl::ElementTargetFlag(geom_element_iel.get_elem_center_3d());
 //*************************************************** 
   
    
 //**************** state **************************** 
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u);
    sol_u    .resize(nDof_u);
    l2GMap_u.resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
     unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u);
      sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);
      l2GMap_u[i] = pdeSys->GetSystemDof(solIndex_u, solPdeIndex_u, i, iel);
    }
 //***************************************************  
 
 //***************** control ************************* 
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);
    sol_ctrl    .resize(nDof_ctrl);
    l2GMap_ctrl.resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl);
      sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);
      l2GMap_ctrl[i] = pdeSys->GetSystemDof(solIndex_ctrl, solPdeIndex_ctrl, i, iel);
    } 
 //*************************************************** 
 

 //************** adjoint **************************** 
    unsigned nDof_adj  = msh->GetElementDofNumber(iel, solType_adj);
        sol_adj    .resize(nDof_adj);
        l2GMap_adj.resize(nDof_adj);
    for (unsigned i = 0; i < sol_adj.size(); i++) {
      unsigned solDof_adj = msh->GetSolutionDof(i, iel, solType_adj);
      sol_adj[i] = (*sol->_Sol[solIndex_adj])(solDof_adj);
      l2GMap_adj[i] = pdeSys->GetSystemDof(solIndex_adj, solPdeIndex_adj, i, iel);
    } 
 //***************************************************  
 
 //************** mu **************************** 
    unsigned nDof_mu  = msh->GetElementDofNumber(iel, solType_mu);
    sol_mu   .resize(nDof_mu);
    l2GMap_mu.resize(nDof_mu);
    for (unsigned i = 0; i < sol_mu.size(); i++) {
      unsigned solDof_mu = msh->GetSolutionDof(i, iel, solType_mu);
      sol_mu[i] = (*sol->_Sol[solIndex_mu])(solDof_mu);
      l2GMap_mu[i] = pdeSys->GetSystemDof(solIndex_mu, solPdeIndex_mu, i, iel);
    }
 //***************************************************  
    
 
    
 //******************** ALL VARS ********************* 
    unsigned nDof_AllVars = nDof_u + nDof_ctrl + nDof_adj + nDof_mu; 
    int nDof_max    =  nDof_u;   // TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
    
    if(nDof_adj > nDof_max) 
    {
      nDof_max = nDof_adj;
      }
    
    if(nDof_ctrl > nDof_max)
    {
      nDof_max = nDof_ctrl;
    }
    
    
    Res.resize(nDof_AllVars);
    std::fill(Res.begin(), Res.end(), 0.);

    Jac.resize(nDof_AllVars * nDof_AllVars);
    std::fill(Jac.begin(), Jac.end(), 0.);
    
    l2GMap_AllVars.resize(0);
    l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_u.begin(),l2GMap_u.end());
    l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_ctrl.begin(),l2GMap_ctrl.end());
    l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_adj.begin(),l2GMap_adj.end());
    l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_mu.begin(),l2GMap_mu.end());
    
  for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
	Sol_n_el_dofs_Mat_vol[k] = ndofs_unk;
  }
    
 //*************************************************** 
    
 //***** set control flag ****************************
  int control_el_flag = 0;
  control_el_flag = ctrl::ControlDomainFlag_internal_restriction(geom_element_iel.get_elem_center_3d());
  std::vector<int> control_node_flag(nDof_ctrl, 0);
  if (control_el_flag == 1) std::fill(control_node_flag.begin(), control_node_flag.end(), 1);
 //*************************************************** 
  
 //========= gauss value quantities ==================   
	double sol_u_gss = 0.;
	double sol_adj_gss = 0.;
	double sol_ctrl_gss = 0.;
	double sol_mu_gss = 0.;
	std::vector<double> sol_u_x_gss(dim);       std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	std::vector<double> sol_adj_x_gss(dim);     std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
	std::vector<double> sol_ctrl_x_gss(dim);    std::fill(sol_ctrl_x_gss.begin(), sol_ctrl_x_gss.end(), 0.);
	std::vector<double> sol_mu_x_gss(dim);      std::fill(sol_mu_x_gss.begin(), sol_mu_x_gss.end(), 0.);
 //===================================================   

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
	msh->_finiteElement[ielGeom][solType_u]   ->Jacobian(geom_element_iel.get_coords_at_dofs_3d(), ig, weight, phi_u, phi_u_x, boost::none);
    msh->_finiteElement[ielGeom][solType_ctrl]->Jacobian(geom_element_iel.get_coords_at_dofs_3d(), ig, weight, phi_ctrl, phi_ctrl_x, boost::none);
    msh->_finiteElement[ielGeom][solType_adj] ->Jacobian(geom_element_iel.get_coords_at_dofs_3d(), ig, weight, phi_adj, phi_adj_x, boost::none);
	msh->_finiteElement[ielGeom][solType_mu]  ->Jacobian(geom_element_iel.get_coords_at_dofs_3d(), ig, weight, phi_mu, phi_mu_x, boost::none);

	
    sol_u_gss = 0.;
	sol_adj_gss = 0.;
	sol_ctrl_gss = 0.;
	sol_mu_gss = 0.;
    
	std::fill(sol_u_x_gss.begin(),sol_u_x_gss.end(), 0.);
	std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
	std::fill(sol_ctrl_x_gss.begin(), sol_ctrl_x_gss.end(), 0.);
	std::fill(sol_mu_x_gss.begin(), sol_mu_x_gss.end(), 0.);
	
	for (unsigned i = 0; i < nDof_u; i++) {
	                                                sol_u_gss      += sol_u[i] * phi_u[i];
                   for (unsigned d = 0; d < dim; d++)   sol_u_x_gss[d] += sol_u[i] * phi_u_x[i * dim + d];
          }
	
	for (unsigned i = 0; i < nDof_adj; i++) {
	                                                sol_adj_gss      += sol_adj[i] * phi_adj[i];
                   for (unsigned d = 0; d < dim; d++)   sol_adj_x_gss[d] += sol_adj[i] * phi_adj_x[i * dim + d];
        }
	
	for (unsigned i = 0; i < nDof_ctrl; i++) {
	                                                sol_ctrl_gss      += sol_ctrl[i] * phi_ctrl[i];
                   for (unsigned d = 0; d < dim; d++)   sol_ctrl_x_gss[d] += sol_ctrl[i] * phi_ctrl_x[i * dim + d];
        }
        
    for (unsigned i = 0; i < nDof_mu; i++) {
	                                                sol_mu_gss      += sol_mu[i] * phi_mu[i];
		   for (unsigned d = 0; d < dim; d++)   sol_mu_x_gss[d] += sol_mu[i] * phi_mu_x[i * dim + d];
	}
        
//==========FILLING WITH THE EQUATIONS ===========
	// *** phi_i loop ***
        for (unsigned i = 0; i < nDof_max; i++) {
	  
	      double laplace_rhs_du_adj_i = 0.; //
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_u )         laplace_rhs_du_adj_i             +=  (phi_u_x   [i * dim + kdim] * sol_adj_x_gss[kdim]);
	      }
	      
              double laplace_rhs_dctrl_ctrl_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_ctrl )         laplace_rhs_dctrl_ctrl_i      +=  (phi_ctrl_x   [i * dim + kdim] * sol_ctrl_x_gss[kdim]);
	      }
	      
	      double laplace_rhs_dctrl_adj_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_ctrl )         laplace_rhs_dctrl_adj_i       +=  (phi_ctrl_x   [i * dim + kdim] * sol_adj_x_gss[kdim]);
	      }
	      
	      double laplace_rhs_dadj_u_i = 0.;  //
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_adj )         laplace_rhs_dadj_u_i           +=  (phi_adj_x   [i * dim + kdim] * sol_u_x_gss[kdim]);
	      }
	      
	      double laplace_rhs_dadj_ctrl_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_adj )         laplace_rhs_dadj_ctrl_i        +=  (phi_adj_x   [i * dim + kdim] * sol_ctrl_x_gss[kdim]);
	      }
//======================Residuals - BEGIN =======================
          // FIRST ROW
	  if (i < nDof_u)                      Res[0      + i] += - weight * (target_flag * phi_u[i] * ( sol_u_gss + sol_ctrl_gss - u_des) - laplace_rhs_du_adj_i );
          // SECOND ROW
	  if (i < nDof_ctrl)  {
	     if ( control_el_flag == 1)    {    Res[nDof_u + i] +=  /*(control_node_flag[i]) **/ - weight * (target_flag * phi_ctrl[i] * ( sol_u_gss + sol_ctrl_gss - u_des) 
													                                                  + alpha * phi_ctrl[i] * sol_ctrl_gss
		                                                                                              + beta * laplace_rhs_dctrl_ctrl_i
		                                                                                              - laplace_rhs_dctrl_adj_i /// @todo this is here because variations of restricted functions are also restricted
													      /*+ ineq_flag * sol_mu_gss*/ ); }
	      else if ( control_el_flag == 0) { Res[nDof_u + i] +=   (- penalty_outside_control_domain)  *  (1 - control_node_flag[i]) * (sol_ctrl[i] - 0.); }
	  }
          // THIRD ROW
          if (i < nDof_adj)        Res[nDof_u + nDof_ctrl + i] += /*-weight * phi_adj[i] * sol_adj_gss - 6.;*/- weight *  ( - laplace_rhs_dadj_u_i - laplace_rhs_dadj_ctrl_i ) ;

          // FOURTH ROW
       if (i < nDof_mu)  {
            if ( control_el_flag == 0) {  Res[nDof_u + nDof_ctrl + nDof_adj + i] +=   (- penalty_outside_control_domain) *  (1 - control_node_flag[i]) * (sol_mu[i] - 0.); }   //MU
       }
       
           
//====================== Residuals - END =======================
	      
//====================== Jacobians - BEGIN =======================
          if (assembleMatrix) {
	    
            // *** phi_j loop ***
            for (unsigned j = 0; j < nDof_max; j++) {
              //double laplace_mat_du_u = 0.;
              double laplace_mat_du_adj = 0.;
              double laplace_mat_dadj_u = 0.;
              double laplace_mat_dctrl_adj = 0.;
              double laplace_mat_dadj_ctrl = 0.;
              double laplace_mat_dctrl_ctrl = 0.;

              for (unsigned kdim = 0; kdim < dim; kdim++) {
              //if ( i < nDof_u && j < nDof_u )           laplace_mat_du_u           += (phi_u_x   [i * dim + kdim] * phi_u_x   [j * dim + kdim]);
              if ( i < nDof_u && j < nDof_adj )         laplace_mat_du_adj         += (phi_u_x   [i * dim + kdim] * phi_adj_x [j * dim + kdim]);
              if ( i < nDof_adj && j < nDof_u )         laplace_mat_dadj_u         += (phi_adj_x [i * dim + kdim] * phi_u_x   [j * dim + kdim]);  //equal to the previous
              if ( i < nDof_ctrl && j < nDof_adj )      laplace_mat_dctrl_adj      += (phi_ctrl_x[i * dim + kdim] * phi_adj_x [j * dim + kdim]);
              if ( i < nDof_adj && j < nDof_ctrl )      laplace_mat_dadj_ctrl      += (phi_adj_x [i * dim + kdim] * phi_ctrl_x[j * dim + kdim]);  //equal to the previous
              if ( i < nDof_ctrl && j < nDof_ctrl )     laplace_mat_dctrl_ctrl     += (phi_ctrl_x  [i * dim + kdim] * phi_ctrl_x  [j * dim + kdim]);
	      }

              //============ delta_state row ============================
              //DIAG BLOCK delta_state - state
	      if ( i < nDof_u && j < nDof_u )       
		Jac[ (0 + i) * nDof_AllVars   +
		     (0 + j)                            ]  += weight * target_flag * phi_u[j] *  phi_u[i];
              
	      // BLOCK  delta_state - control
              if ( i < nDof_u && j < nDof_ctrl )   
		Jac[ (0 + i) * nDof_AllVars   +
                     (nDof_u + j)                       ]  += weight * target_flag  * phi_ctrl[j] *  phi_u[i];
	      
              // BLOCK  delta_state - adjoint
              if ( i < nDof_u && j < nDof_adj )  
		Jac[  (0 + i) * nDof_AllVars  +
                      (nDof_u + nDof_ctrl + j)          ]  += weight * (-1.) * laplace_mat_du_adj;
              
	      //=========== delta_control row ===========================     
	      if ( control_el_flag == 1)  {

	      //BLOCK delta_control - state
              if ( i < nDof_ctrl   && j < nDof_u   ) 
		Jac[ (nDof_u + i) * nDof_AllVars  +
		     (0 + j)                            ]  += ( control_node_flag[i]) * weight * target_flag * phi_u[j] * phi_ctrl[i];
		
	      //BLOCK delta_control - control
              if ( i < nDof_ctrl   && j < nDof_ctrl   )
		Jac[ (nDof_u + i) * nDof_AllVars +
		     (nDof_u + j)                       ]  += ( control_node_flag[i]) * weight * ( beta * /*control_el_flag  **/ laplace_mat_dctrl_ctrl 
		                                                                                + alpha * /*control_el_flag **/ phi_ctrl[i] * phi_ctrl[j] 
		                                                                                            + target_flag * phi_ctrl[i] * phi_ctrl[j] );
              
	      //BLOCK delta_control - adjoint
              if ( i < nDof_ctrl   && j < nDof_adj  ) 
		Jac[ (nDof_u + i) * nDof_AllVars  + 
		     (nDof_u + nDof_ctrl + j)           ]  += ( control_node_flag[i]) * weight * (-1.) * laplace_mat_dctrl_adj;
	      	      
	        }
	      
	      else if ( control_el_flag == 0)  {  
		
              //BLOCK delta_control - control
              if ( i < nDof_ctrl   && j < nDof_ctrl   &&  i==j ) {
		 Jac[ (nDof_u + i) * nDof_AllVars +
		      (nDof_u + j)                      ]  +=  penalty_outside_control_domain * (1 - control_node_flag[i]);
		}
	      
	      }
	      
		     
	      //=========== delta_adjoint row ===========================
              // BLOCK delta_adjoint - state	      
              if ( i < nDof_adj && j < nDof_u )   
		Jac[ (nDof_u + nDof_ctrl + i) * nDof_AllVars +
		     (0 + j)                            ]  += weight * (-1) * laplace_mat_dadj_u;   
	      
              // BLOCK delta_adjoint - control   
              if ( i < nDof_adj && j < nDof_ctrl )  
		Jac[ (nDof_u + nDof_ctrl + i)  * nDof_AllVars +
		     (nDof_u  + j)                      ]  += weight * (-1) * laplace_mat_dadj_ctrl; 
		     
		     
// 	      // BLOCK delta_adjoint - adjoint   
//               if ( i < nDof_adj && j < nDof_adj )  
// 		Jac[ (nDof_u + nDof_ctrl + i)  * nDof_AllVars +
// 		     (nDof_u + nDof_ctrl + j)                      ]  += weight * phi_adj[j] *  phi_adj[i]; 
    
	      
	      //============= delta_mu row ===============================
                if ( i < nDof_mu && j < nDof_mu && i==j )  {   
	        if ( control_el_flag == 0)  {  
		  Jac[   (nDof_u + nDof_ctrl + nDof_adj + i) * nDof_AllVars +  (nDof_u + nDof_ctrl + nDof_adj + j) ]  += penalty_outside_control_domain * (1 - control_node_flag[i]);    //MU
                }
	      }
	      
	      
            } // end phi_j loop
          } // endif assemble_matrix
//====================== Jacobians - END =======================

        } // end phi_i loop
        
      } // end gauss point loop

      
      
    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
   // std::cout << " ************* Element ************** " << iel << " **************************************** " << std::endl;     

// // //     if (control_el_flag == 0) {  //elements that should have zero control
// //          for (unsigned i_unk = 0; i_unk < n_unknowns; i_unk++) {
// //     std::cout << " ======= Row === " << i_unk << " =================================================== " << std::endl;     
// //         unsigned int row_block_offset = 0;
// // 	         for (unsigned k = 0; k < i_unk; k++) row_block_offset += Sol_n_el_dofs_Mat_vol[k];
// //          for (unsigned j_unk = 0; j_unk < n_unknowns; j_unk++) {
// //     std::cout << " ======= Column === " << j_unk << " ================== " << std::endl;     
// //         unsigned int column_block_offset = 0;
// // 	         for (unsigned k = 0; k < j_unk; k++) column_block_offset += Sol_n_el_dofs_Mat_vol[k];
// // 	  
// //          for (unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[i_unk]; i++) {
// // // 	      std::cout << Res[nDof_u + nDof_ctrl + nDof_adj + i ] << " " << std::endl;
// // 	   for (unsigned j = 0; j < Sol_n_el_dofs_Mat_vol[j_unk]; j++) {
// // 	      std::cout <<  " " << std::setfill(' ') << std::setw(10) << Jac[ (row_block_offset + i) * nDof_AllVars + ( column_block_offset + j) ] << " ";
// // 	    }
// // 	      std::cout << std::endl;
// // 	 }
// // 
// // 	 } //j_unk
// // 	} //i_unk
	 
	 
// // // 	}


 //========== end of integral-based part

    //copy the value of the adept::adoube aRes in double Res and store
    RES->add_vector_blocked(Res, l2GMap_AllVars);
      if (assembleMatrix)  JAC->add_matrix_blocked(Jac, l2GMap_AllVars, l2GMap_AllVars);

      
     if (print_algebra_local) {
         assemble_jacobian<double,double>::print_element_residual(iel, Res, Sol_n_el_dofs_Mat_vol, 10, 5);
         assemble_jacobian<double,double>::print_element_jacobian(iel, Jac, Sol_n_el_dofs_Mat_vol, 10, 5);
     }

       
  } //end element loop for each process
      
  //MU
  
  
// // // 
// // //   unsigned int global_ctrl_size = pdeSys->KKoffset[ctrl_index+1][iproc] - pdeSys->KKoffset[ctrl_index][iproc];
// // //   
// // //   std::vector<double>  one_times_mu(global_ctrl_size, 0.);
// // //   std::vector<int>    positions(global_ctrl_size);
// // // //  double position_mu_i;
// // //   for (unsigned i = 0; i < positions.size(); i++) {
// // //     positions[i] = pdeSys->KKoffset[ctrl_index][iproc] + i;
// // // //     position_mu_i = pdeSys->KKoffset[mu_index][iproc] + i;
// // // //     std::cout << position_mu_i << std::endl;
// // //     one_times_mu[i] = ineq_flag * 1. * (*sol->_Sol[solIndex_mu])(i/*position_mu_i*/) ;
// // //   }
// // //     RES->add_vector_blocked(one_times_mu, positions);
  std::vector<unsigned int> ctrl_index(1);  ctrl_index[0] = mlPdeSys->GetSolPdeIndex("control");
  std::vector<unsigned int>   mu_index(1);    mu_index[0] = mlPdeSys->GetSolPdeIndex("mu");
    
ctrl_inequality::add_one_times_mu_res_ctrl(iproc,
                               ineq_flag,
                               ctrl_index,
                               mu_index,
                               SolIndex,
                               sol,
                               mlPdeSys,
                               pdeSys,
                               RES);
    
    
    
     
RES->close();
if (assembleMatrix) JAC->close();  /// This is needed for the parallel, when splitting the add part from the insert part!!!
      // ***************** ADD PART - END  *******************
    
    
 


//   ***************** INSERT PART - BEGIN (must go AFTER the sum, clearly) *******************
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

// -------
   geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
      
   geom_element_iel.set_elem_center_3d(iel, solType_coords);
// -------
   
// -------
    el_dofs_unknowns_vol(sol, msh, pdeSys, iel,
                        SolFEType,
                        SolIndex,
                        SolPdeIndex,
                        Sol_n_el_dofs_Mat_vol,
                        sol_eldofs_Mat,
                        L2G_dofmap_Mat);
      
    //***** set control flag ****************************
  int control_el_flag = 0;
  control_el_flag = ctrl::ControlDomainFlag_internal_restriction(geom_element_iel.get_elem_center_3d());
 
    
    if (control_el_flag == 1) {

        
  ctrl_inequality::update_active_set_flag_for_current_nonlinear_iteration
  (msh,
   sol,
   iel,
   geom_element_iel.get_coords_at_dofs/*_3d*/(),
   sol_eldofs_Mat,
   Sol_n_el_dofs_Mat_vol,
   c_compl,
   pos_mu_in_mat,
   pos_ctrl_in_mat,
   solIndex_act_flag_sol,
   ctrl_lower,
   ctrl_upper,
   sol_actflag);
  
      


    ctrl_inequality::node_insertion(iel,
                   msh,
                   L2G_dofmap_Mat,
                   pos_mu_in_mat,
                   pos_ctrl_in_mat,
                   sol_eldofs_Mat,
                   Sol_n_el_dofs_Mat_vol,
                   sol_actflag,
                   ctrl_lower,
                   ctrl_upper,
                   ineq_flag,
                   c_compl,
                   JAC,
                   RES,
                   assembleMatrix
                   );
   

       }

 
    
    
 //============= delta_ctrl-delta_mu row ===============================
  if (assembleMatrix) { JAC->matrix_set_off_diagonal_values_blocked(L2G_dofmap_Mat[pos_mat_ctrl], L2G_dofmap_Mat[pos_mat_mu], ineq_flag * 1.); }
  
  
  }
//   ***************** INSERT PART - END (must go AFTER the sum, clearly) *******************
  
  RES->close();

  if (assembleMatrix) JAC->close();
  
  // ***************** END ASSEMBLY *******************
  
  
  print_global_residual_jacobian(print_algebra_global,
                                 ml_prob,
                                 mlPdeSys,
                                 pdeSys,
                                 RES,
                                 JAC,
                                 iproc,
                                 assembleMatrix);


  return;
}



