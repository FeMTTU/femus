#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"
#include "NumericVector.hpp"

#include "./elliptic_nonlin_param.hpp"

using namespace femus;


 
  // ./elliptic_nonlin -mat_view draw -draw_save myfile.ppm


double  nonlin_term_function(const double& v) {
    
   return 1.;
//    return 1./( (1. - v) );
//    return 0.01*1./( (1. - v)*(1. - v) );
//     return exp(v);
 }


double  nonlin_term_derivative(const double& v) {
    
    return 0.;
//    return  +2. * 1./( (1. - v)*(1. - v) ); 
//    return 0.01* (+2.) * 1./( (1. - v)*(1. - v)*(1. - v) ); 
//     return exp(v);
 }


 inline unsigned int res_row_index(const std::vector<unsigned int>& _Sol_n_el_dofs, const int my_row_pos, const int i) {

    assert(i < _Sol_n_el_dofs[my_row_pos]); 
    
    unsigned int pos_previous = 0;
    for (unsigned k = 0; k < my_row_pos; k++) pos_previous += _Sol_n_el_dofs[k];

    return pos_previous + i;
  }
  

 inline unsigned int jac_row_col_index(const std::vector<unsigned int>& _Sol_n_el_dofs, const int nDof_AllVars, const int my_row_pos, const int my_col_pos, const int i, const int j) {

     assert(i < _Sol_n_el_dofs[my_row_pos]); 
     assert(j < _Sol_n_el_dofs[my_col_pos]); 
     
    unsigned int pos_previous_row = 0;
    unsigned int pos_previous_col = 0;
    for (unsigned k = 0; k < my_row_pos; k++) pos_previous_row += _Sol_n_el_dofs[k];
    for (unsigned k = 0; k < my_col_pos; k++) pos_previous_col += _Sol_n_el_dofs[k];

    return (pos_previous_row + i) * nDof_AllVars + (pos_previous_col + j);
  }

  
 inline double laplacian_row(const vector < unsigned > & SolFEType, const vector < vector < double > > & phi_x_fe_qp, const vector < vector < double > > & sol_grad_qp, const int my_row_pos, const int my_col_pos, const int i, const unsigned dim_in) {
  
              double laplace_rhs_i_j = 0.;

              for (unsigned kdim = 0; kdim < dim_in; kdim++) laplace_rhs_i_j  += phi_x_fe_qp[SolFEType[my_row_pos]][i * dim_in + kdim] *
                                                                                           sol_grad_qp[my_col_pos][kdim];

     return laplace_rhs_i_j;
}
  
  
 inline double laplacian_row_col(const vector < unsigned > & SolFEType, const vector < vector < double > > & phi_x_fe_qp, const int my_row_pos, const int my_col_pos, const int i, const int j, const unsigned dim_in) {
  
              double laplace_mat_i_j = 0.;

              for (unsigned kdim = 0; kdim < dim_in; kdim++) laplace_mat_i_j  += phi_x_fe_qp[SolFEType[my_row_pos]][i * dim_in + kdim] *
                                                                                 phi_x_fe_qp[SolFEType[my_col_pos]][j * dim_in + kdim];

     return laplace_mat_i_j;
}


double SetInitialCondition (const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {
         
           double value = 0.;

             if(!strcmp(name,"control")) {
                 value = 0.;
             }
           
             if(!strcmp(name,"mu")) {
                 value = 0.;
             }
             if(!strcmp(name,"state")) {
                 value = 0.;
             }
           
             if(!strcmp(name,"adjoint")) {
                 value = 0.;
             }
             if(!strcmp(name,"TargReg")) {
                 value = ElementTargetFlag(x);
             }
           
             if(!strcmp(name,"ContReg")) {
                 value = ControlDomainFlag_internal_restriction(x);
             }
             if(!strcmp(name,"act_flag")) {
                 value = 0.;
             }
           
      return value;   
}
  

  
bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int faceName, const double time) {

  bool dirichlet = true; //dirichlet
  value = 0.;
  
  if(!strcmp(name,"control")) {
      value = 0.;
//   if (faceName == 3)
//     dirichlet = false;
  
  }
  
  if(!strcmp(name,"mu")) {
//       value = 0.;
//   if (faceName == 3)
    dirichlet = false;
  }
  
  return dirichlet;
}


void ComputeIntegral(const MultiLevelProblem& ml_prob);

void AssembleProblem(MultiLevelProblem& ml_prob);



int main(int argc, char** args) {

    // ======= Init ========================
  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
    // ======= Files ========================
  Files files; 
        files.CheckIODirectories();
        files.RedirectCout();

    // ======= Quad Rule ========================
  std::string fe_quad_rule("seventh");
 /* "seventh" is the order of accuracy that is used in the gauss integration scheme
    In the future it is not going to be an argument of the mesh function   */
  
    // ======= Mesh ========================
  MultiLevelMesh mlMsh;
  mlMsh.GenerateCoarseBoxMesh(NSUB_X,NSUB_Y,0,0.,1.,0.,1.,0.,0.,QUAD9,fe_quad_rule.c_str());
  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  mlMsh.PrintInfo();

    // ======= Solution ========================
  MultiLevelSolution mlSol(&mlMsh);  // define the multilevel solution and attach the mlMsh object to it

  // add variables to mlSol
  mlSol.AddSolution("state",   LAGRANGE, FIRST);
  mlSol.AddSolution("control", LAGRANGE, FIRST);
  mlSol.AddSolution("adjoint", LAGRANGE, FIRST);
  mlSol.AddSolution("mu",      LAGRANGE, FIRST);  
  mlSol.AddSolution("TargReg", DISCONTINOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  mlSol.AddSolution("ContReg", DISCONTINOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  const unsigned int fake_time_dep_flag = 2;  //this is needed to be able to use _SolOld
  const std::string act_set_flag_name = "act_flag";
  mlSol.AddSolution(act_set_flag_name.c_str(), LAGRANGE, FIRST,fake_time_dep_flag);               

    // ======= Problem ========================
  MultiLevelProblem ml_prob(&mlSol);  // define the multilevel problem attach the mlSol object to it

  ml_prob.SetFilesHandler(&files);
  
    // ======= Initial values ========================
  mlSol.Initialize("All");    // initialize all variables to zero

//   mlSol.Initialize("All", SetInitialCondition, &ml_prob); //unfortunately if I do this it sets all to zero //I would like to do an attach function similar to the BC
  mlSol.Initialize("state",   SetInitialCondition, &ml_prob);
  mlSol.Initialize("control", SetInitialCondition, &ml_prob);
  mlSol.Initialize("adjoint", SetInitialCondition, &ml_prob);
  mlSol.Initialize("mu",      SetInitialCondition, &ml_prob);
  mlSol.Initialize("TargReg", SetInitialCondition, &ml_prob);
  mlSol.Initialize("ContReg", SetInitialCondition, &ml_prob);
  mlSol.Initialize(act_set_flag_name.c_str(),  SetInitialCondition, &ml_prob);

    // ======= Boundary Conditions ========================
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);  // attach the boundary condition function and generate boundary data

//   mlSol.GenerateBdc("All");  //this would do it also for the non-equation-related variables
  mlSol.GenerateBdc("state");
  mlSol.GenerateBdc("control");
  mlSol.GenerateBdc("adjoint");
  mlSol.GenerateBdc("mu");  //we need this for all Pde variables to make the matrix iterations work... but this should be related to the matrix and not to the sol... The same for the initial condition

    // ======= System ========================
  NonLinearImplicitSystemWithPrimalDualActiveSetMethod& system = ml_prob.add_system < NonLinearImplicitSystemWithPrimalDualActiveSetMethod > ("OptSys");    // add system in ml_prob

  system.SetActiveSetFlagName(act_set_flag_name);

  //here we decide the order in the matrix!
  const std::vector < std::string > sol_matrix_pos = {"state","control","adjoint","mu"};
  for (unsigned k = 0; k < sol_matrix_pos.size(); k++)  system.AddSolutionToSystemPDE(sol_matrix_pos[k].c_str());  
  
  // attach the assembling function to system
  system.SetAssembleFunction(AssembleProblem);
  
  mlSol.SetWriter(VTK);   //need to move this here for the DebugNonlinear function
  mlSol.GetWriter()->SetDebugOutput(true);
  
  system.SetDebugNonlinear(true);
  system.SetDebugFunction(ComputeIntegral);
  //   system.SetMaxNumberOfNonLinearIterations(2);

    // ======= Solve ========================
  system.init();    // initialize and solve the system
  system.MGsolve();
  
//   ComputeIntegral(ml_prob);
 
    // ======= Final Print ========================
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("all");
  mlSol.GetWriter()->Write(files.GetOutputPath()/*DEFAULT_OUTPUTDIR*/, "biquadratic", variablesToBePrinted);    // print solutions


  return 0;
}


void AssembleProblem(MultiLevelProblem& ml_prob) {

  // ************** J dx = f - J x_old ******************  
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystemWithPrimalDualActiveSetMethod> ("OptSys");
  const unsigned          level      = mlPdeSys->GetLevelToAssemble();
  const bool          assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                          msh = ml_prob._ml_msh->GetLevel(level);

  MultiLevelSolution*          mlSol = ml_prob._ml_sol;
  Solution*                      sol = ml_prob._ml_sol->GetSolutionLevel(level);

  LinearEquationSolver*       pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*                   KK = pdeSys->_KK;
  NumericVector*                 RES = pdeSys->_RES;

  const unsigned     dim = msh->GetDimension();                                 // get the domain dimension of the problem
  const unsigned    dim2 = (3 * (dim - 1) + !(dim - 1));                        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  const unsigned   iproc = msh->processor_id(); 

 //***************************************************  
  const unsigned int n_unknowns = mlPdeSys->GetSolPdeIndex().size();
  
  enum Sol_pos{pos_state=0,pos_ctrl,pos_adj,pos_mu};  //these are known at compile-time 
  //right now this is the same as Sol_pos = SolPdeIndex, but now I want to have flexible rows and columns
  // the ROW index corresponds to the EQUATION we want to solve
  // the COLUMN index corresponds to the column variables 
  
  // Now, first of all we have to make sure that the Sparsity pattern is correctly built...
  // This means that we should write the Sparsity pattern with COLUMNS independent of ROWS!
  // But that implies that the BOUNDARY CONDITIONS will be enforced in OFF-DIAGONAL BLOCKS!
  // In fact, having the row order be different from the column order implies that the diagonal blocks would be RECTANGULAR.  
  // Although this is in principle possible, I would avoid doing that for now.
  
//   const unsigned int pos_state = mlPdeSys->GetSolPdeIndex("state");   //these are known at run-time, so they are slower
//   const unsigned int pos_ctrl  = mlPdeSys->GetSolPdeIndex("control"); //these are known at run-time, so they are slower
//   const unsigned int pos_adj   = mlPdeSys->GetSolPdeIndex("adjoint"); //these are known at run-time, so they are slower
//   const unsigned int pos_mu    = mlPdeSys->GetSolPdeIndex("mu");      //these are known at run-time, so they are slower
  
  assert(pos_state == mlPdeSys->GetSolPdeIndex("state"));
  assert(pos_ctrl  == mlPdeSys->GetSolPdeIndex("control"));
  assert(pos_adj   == mlPdeSys->GetSolPdeIndex("adjoint"));
  assert(pos_mu    == mlPdeSys->GetSolPdeIndex("mu"));
 //***************************************************  

 //***************************************************  
  vector < vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
  for (unsigned i = 0; i < dim; i++) {
    coordX[i].reserve(maxSize);
  }
 //***************************************************   

 //***************************************************  
  double weight_qp; // gauss point weight
  
  vector < vector < double > > phi_fe_qp(NFE_FAMS);
  vector < vector < double > > phi_x_fe_qp(NFE_FAMS);
  vector < vector < double > > phi_xx_fe_qp(NFE_FAMS);
 
  for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_fe_qp[fe].reserve(maxSize);
      phi_x_fe_qp[fe].reserve(maxSize*dim);
     phi_xx_fe_qp[fe].reserve(maxSize*(3*(dim-1)));
   }

 //***************************************************  
 //********* WHOLE SET OF VARIABLES ****************** 
  const int solType_max = 2;  //biquadratic

  vector < std::string > Solname(n_unknowns);
  Solname[pos_state] = "state";
  Solname[pos_ctrl]  = "control";
  Solname[pos_adj]   = "adjoint";
  Solname[pos_mu]    = "mu";
  
 int m_b_f[n_unknowns][n_unknowns];
     m_b_f[pos_state][pos_state] = 1; //nonzero
     m_b_f[pos_state][pos_ctrl]  = 0; //THIS IS ZERO IN NON-LIFTING APPROACHES (there are also pieces that go on top on blocks that are already meant to be different from zero)
     m_b_f[pos_state][pos_adj]   = 1; //nonzero
     m_b_f[pos_state][pos_mu]    = 0;  //this is zero without state constraints
     
     m_b_f[pos_ctrl][pos_state]  = 0;//THIS IS ZERO IN NON-LIFTING APPROACHES
     m_b_f[pos_ctrl][pos_ctrl]   = 1;//nonzero
     m_b_f[pos_ctrl][pos_adj]    = 1;//nonzero
     m_b_f[pos_ctrl][pos_mu]     = 1;//nonzero
     
     m_b_f[pos_adj][pos_state]  = 1; //nonzero
     m_b_f[pos_adj][pos_ctrl]   = 1; //nonzero
     m_b_f[pos_adj][pos_adj]    = 0; //this must always be zero 
     m_b_f[pos_adj][pos_mu]     = 0; //this must always be zero 
     
     m_b_f[pos_mu][pos_state]  = 0; //this is zero without state constraints
     m_b_f[pos_mu][pos_ctrl]   = 1; //nonzero
     m_b_f[pos_mu][pos_adj]    = 0; //this must always be zero 
     m_b_f[pos_mu][pos_mu]     = 1; //nonzero
 
  //***************************************************  
  vector < vector < int > > L2G_dofmap(n_unknowns);     for(int i = 0; i < n_unknowns; i++) { L2G_dofmap[i].reserve(maxSize); }
            vector< int >   L2G_dofmap_AllVars; L2G_dofmap_AllVars.reserve( n_unknowns*maxSize );
            vector< double >         Res;                      Res.reserve( n_unknowns*maxSize );
            vector< double >         Jac;                      Jac.reserve( n_unknowns*maxSize * n_unknowns*maxSize);
  //***************************************************  
  
  //***************************************************  
  vector < unsigned > SolPdeIndex(n_unknowns);  //index as in the row/column of the matrix (diagonal blocks are square)
  vector < unsigned > SolIndex(n_unknowns);     //index as in the MultilevelSolution vector
  vector < unsigned > SolFEType(n_unknowns);    //FEtype of each MultilevelSolution       

  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolPdeIndex[ivar] = mlPdeSys->GetSolPdeIndex(  Solname[ivar].c_str() );
       SolIndex[ivar] = mlSol->GetIndex         (  Solname[ivar].c_str() );
      SolFEType[ivar] = mlSol->GetSolutionType  ( SolIndex[ivar]);
  }
  //***************************************************  

  vector < unsigned > Sol_n_el_dofs(n_unknowns);
   
  //----------- quantities at dof objects ------------------------------
  vector < vector < double > >     sol_eldofs(n_unknowns);
  for(int k=0; k<n_unknowns; k++)  sol_eldofs[k].reserve(maxSize);
  
  
  //------------ quantities at quadrature points ---------------------
            vector<double>        sol_qp(n_unknowns);
    vector< vector<double> > sol_grad_qp(n_unknowns);
    
      std::fill(sol_qp.begin(), sol_qp.end(), 0.);
    for (unsigned  k = 0; k < n_unknowns; k++) {
        sol_grad_qp[k].resize(dim);
        std::fill(sol_grad_qp[k].begin(), sol_grad_qp[k].end(), 0.);
    }

//************** act flag **************************** 
   std::string act_flag_name = "act_flag";
   unsigned int solIndex_act_flag = mlSol->GetIndex(act_flag_name.c_str());
   unsigned int solFEType_act_flag = mlSol->GetSolutionType(solIndex_act_flag); 
      if(sol->GetSolutionTimeOrder(solIndex_act_flag) == 2) {
        *(sol->_SolOld[solIndex_act_flag]) = *(sol->_Sol[solIndex_act_flag]);
      }
    
    
  //********* variables for ineq constraints *****************
  const int ineq_flag = INEQ_FLAG;
  const double c_compl = C_COMPL;
  vector < double/*int*/ >  sol_actflag;   sol_actflag.reserve(maxSize); //flag for active set
  vector < double >  ctrl_lower;   ctrl_lower.reserve(maxSize);
  vector < double >  ctrl_upper;   ctrl_upper.reserve(maxSize);
 //***************************************************  

  //********************* DATA ************************ 
  double u_des = DesiredTarget();
  double alpha = ALPHA_CTRL_VOL;
  double beta  = BETA_CTRL_VOL;
  double penalty_strong = 10e+14;
 //***************************************************  

  RES->zero();
  if (assembleMatrix)  KK->zero();
    
  // element loop: each process loops only on the elements that it owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);    // element geometry type

 //******************** GEOMETRY ********************* 
    unsigned nDofx = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs
    for (int i = 0; i < dim; i++)  coordX[i].resize(nDofx);
    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, coordXType);  // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        coordX[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

   // elem average point 
    vector < double > elem_center(dim);   
    for (unsigned j = 0; j < dim; j++) {  elem_center[j] = 0.;  }
    for (unsigned j = 0; j < dim; j++) {  
      for (unsigned i = 0; i < nDofx; i++) {
         elem_center[j] += coordX[j][i];
       }
    }
    
   for (unsigned j = 0; j < dim; j++) { elem_center[j] = elem_center[j]/nDofx; }
 //***************************************************  
  
 //****** set target domain flag ********************* 
   int target_flag = 0;
   target_flag = ElementTargetFlag(elem_center);
 //*************************************************** 
   
   //all vars###################################################################  
  for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned  ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
	   Sol_n_el_dofs[k] = ndofs_unk;
          sol_eldofs[k].resize(ndofs_unk);
          L2G_dofmap[k].resize(ndofs_unk); 
    for (unsigned i = 0; i < ndofs_unk; i++) {
           unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);                        // global to global mapping between solution node and solution dof // via local to global solution node
           sol_eldofs[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);                            // global extraction and local storage for the solution
           L2G_dofmap[k][i] = pdeSys->GetSystemDof(SolIndex[k], SolPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }
  //all vars###################################################################
    
    
 //************** update active set flag for current nonlinear iteration **************************** 
 // 0: inactive; 1: active_a; 2: active_b
   assert(Sol_n_el_dofs[pos_mu] == Sol_n_el_dofs[pos_ctrl]);
   sol_actflag.resize(Sol_n_el_dofs[pos_mu]);
   ctrl_lower.resize(Sol_n_el_dofs[pos_mu]);
   ctrl_upper.resize(Sol_n_el_dofs[pos_mu]);
     std::fill(sol_actflag.begin(), sol_actflag.end(), 0);
     std::fill(ctrl_lower.begin(), ctrl_lower.end(), 0.);
     std::fill(ctrl_upper.begin(), ctrl_upper.end(), 0.);
   
    for (unsigned i = 0; i < sol_actflag.size(); i++) {
        std::vector<double> node_coords_i(dim,0.);
        for (unsigned d = 0; d < dim; d++) node_coords_i[d] = coordX[d][i];
        ctrl_lower[i] = InequalityConstraint(node_coords_i,false);
        ctrl_upper[i] = InequalityConstraint(node_coords_i,true);
        
    if      ( (sol_eldofs[pos_mu][i] + c_compl * (sol_eldofs[pos_ctrl][i] - ctrl_lower[i] )) < 0 )  sol_actflag[i] = 1;
    else if ( (sol_eldofs[pos_mu][i] + c_compl * (sol_eldofs[pos_ctrl][i] - ctrl_upper[i] )) > 0 )  sol_actflag[i] = 2;
    }
 
 //************** act flag **************************** 
    unsigned nDof_act_flag  = msh->GetElementDofNumber(iel, solFEType_act_flag);    // number of solution element dofs
    
    for (unsigned i = 0; i < nDof_act_flag; i++) {
      unsigned solDof_mu = msh->GetSolutionDof(i, iel, solFEType_act_flag); 
      (sol->_Sol[solIndex_act_flag])->set(solDof_mu,sol_actflag[i]);     
    }    

    //******************** ALL VARS ********************* 
    unsigned nDof_AllVars = 0;
    for (unsigned  k = 0; k < n_unknowns; k++) { nDof_AllVars += Sol_n_el_dofs[k]; }
 // TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
    int nDof_max    =  0;   
      for (unsigned  k = 0; k < n_unknowns; k++)     {
          if(Sol_n_el_dofs[k] > nDof_max)    nDof_max = Sol_n_el_dofs[k];
      }
  
    Res.resize(nDof_AllVars);                   std::fill(Res.begin(), Res.end(), 0.);
    Jac.resize(nDof_AllVars * nDof_AllVars);    std::fill(Jac.begin(), Jac.end(), 0.);
    
    L2G_dofmap_AllVars.resize(0);
      for (unsigned  k = 0; k < n_unknowns; k++)     L2G_dofmap_AllVars.insert(L2G_dofmap_AllVars.end(),L2G_dofmap[k].begin(),L2G_dofmap[k].end());
 //*************************************************** 


 //***** set control flag ****************************
  int control_el_flag = 0;
      control_el_flag = ControlDomainFlag_internal_restriction(elem_center);
  std::vector<int> control_node_flag(Sol_n_el_dofs[pos_ctrl],0);
  if (control_el_flag == 1) std::fill(control_node_flag.begin(), control_node_flag.end(), 1);
 //*************************************************** 
  



      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType_max]->GetGaussPointNumber(); ig++) {
	
      // *** get gauss point weight, test function and test function partial derivatives ***
      for(int fe=0; fe < NFE_FAMS; fe++) {
         ml_prob._ml_msh->_finiteElement[ielGeom][fe]->Jacobian(coordX,ig,weight_qp,phi_fe_qp[fe],phi_x_fe_qp[fe],phi_xx_fe_qp[fe]);
      }
   //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
    ml_prob._ml_msh->_finiteElement[ielGeom][BIQUADR_FE]->Jacobian(coordX,ig,weight_qp,phi_fe_qp[BIQUADR_FE],phi_x_fe_qp[BIQUADR_FE],phi_xx_fe_qp[BIQUADR_FE]);
    
 //========= fill gauss value quantities ==================   
   std::fill(sol_qp.begin(), sol_qp.end(), 0.);
   for (unsigned  k = 0; k < n_unknowns; k++) { std::fill(sol_grad_qp[k].begin(), sol_grad_qp[k].end(), 0.); }
    
    for (unsigned  k = 0; k < n_unknowns; k++) {
	for (unsigned i = 0; i < Sol_n_el_dofs[k]; i++) {
	                                                         sol_qp[k]    += sol_eldofs[k][i] *   phi_fe_qp[SolFEType[k]][i];
                   for (unsigned d = 0; d < dim; d++)   sol_grad_qp[k][d] += sol_eldofs[k][i] * phi_x_fe_qp[SolFEType[k]][i * dim + d];
       }        
    }
 //========= fill gauss value quantities ==================   

 
//==========FILLING THE EQUATIONS ===========
	// *** phi_i loop ***
        for (unsigned i = 0; i < nDof_max; i++) {
	  
//======================Residuals=======================
          // ===============
        Res[ res_row_index(Sol_n_el_dofs,pos_state,i) ] += - weight_qp * (target_flag * phi_fe_qp[SolFEType[pos_state]][i] * (
                                                                                              m_b_f[pos_state][pos_state] * sol_qp[pos_state] 
                                                                                            - u_des ) 
                                                                                          + m_b_f[pos_state][pos_adj]   * ( laplacian_row(SolFEType, phi_x_fe_qp, sol_grad_qp, pos_state, pos_adj, i, dim)  
                                                                                                                            - sol_qp[pos_adj] *  nonlin_term_derivative(phi_fe_qp[SolFEType[pos_state]][i]) *  sol_qp[pos_ctrl]) );
        
          

          // ==============
	     if ( control_el_flag == 1)        Res[ res_row_index(Sol_n_el_dofs,pos_ctrl,i) ] +=  /*(control_node_flag[i]) **/ - weight_qp * (
                                                                                  + m_b_f[pos_ctrl][pos_ctrl] * alpha * phi_fe_qp[SolFEType[pos_ctrl]][i] * sol_qp[pos_ctrl]
		                                                                          + m_b_f[pos_ctrl][pos_ctrl] *  beta * laplacian_row(SolFEType, phi_x_fe_qp, sol_grad_qp, pos_ctrl, pos_ctrl, i, dim) 
		                                                                          - m_b_f[pos_ctrl][pos_adj]  *         phi_fe_qp[SolFEType[pos_ctrl]][i] * sol_qp[pos_adj] * nonlin_term_function(sol_qp[pos_state])
                                                                                                                         );
	      else if ( control_el_flag == 0)  Res[ res_row_index(Sol_n_el_dofs,pos_ctrl,i) ] +=  /*(1 - control_node_flag[i]) **/ m_b_f[pos_ctrl][pos_ctrl] * (- penalty_strong) * (sol_eldofs[pos_ctrl][i]);

          // =============
        Res[ res_row_index(Sol_n_el_dofs,pos_adj,i) ] += - weight_qp *  ( + m_b_f[pos_adj][pos_state] * laplacian_row(SolFEType, phi_x_fe_qp, sol_grad_qp, pos_adj, pos_state, i, dim) 
                                                                          - m_b_f[pos_adj][pos_ctrl]  * phi_fe_qp[SolFEType[pos_adj]][i] * sol_qp[pos_ctrl] * nonlin_term_function(sol_qp[pos_state]) 
                                                                        );

//======================End Residuals=======================
	      
          if (assembleMatrix) {
	    
            // *** phi_j loop ***
            for (unsigned j = 0; j < nDof_max; j++) {
                
              //============ delta_state row ============================
		Jac[ jac_row_col_index(Sol_n_el_dofs, nDof_AllVars, pos_state, pos_state, i, j) ]  += 
                                                      m_b_f[pos_state][pos_state] * weight_qp * target_flag * phi_fe_qp[SolFEType[pos_state]][j] *  phi_fe_qp[SolFEType[pos_state]][i];
              
		Jac[ jac_row_col_index(Sol_n_el_dofs, nDof_AllVars, pos_state, pos_adj, i, j)  ]  +=
		                                              m_b_f[pos_state][pos_adj] * weight_qp * (
                                                                          laplacian_row_col(SolFEType, phi_x_fe_qp, pos_state, pos_adj, i, j, dim)
                                                                                                       - phi_fe_qp[SolFEType[pos_adj]][j] *
                                                                                                         nonlin_term_derivative(phi_fe_qp[SolFEType[pos_state]][i]) * 
                                                                                                         sol_qp[pos_ctrl]
                                                                                              );
              
	      //=========== delta_control row ===========================     
	      if ( control_el_flag == 1)  {
		Jac[ jac_row_col_index(Sol_n_el_dofs, nDof_AllVars, pos_ctrl, pos_ctrl, i, j) ]  += 
		                                              m_b_f[pos_ctrl][pos_ctrl] * ( control_node_flag[i]) * weight_qp * (
                                                                                  beta * control_el_flag  * laplacian_row_col(SolFEType, phi_x_fe_qp, pos_ctrl, pos_ctrl, i, j, dim) 
                                                                               + alpha * control_el_flag  * phi_fe_qp[SolFEType[pos_ctrl]][i] * phi_fe_qp[SolFEType[pos_ctrl]][j] 
		                                                                                    );
              
		Jac[ jac_row_col_index(Sol_n_el_dofs, nDof_AllVars, pos_ctrl, pos_adj, i, j) ]  += m_b_f[pos_ctrl][pos_adj] * control_node_flag[i] * weight_qp * (
                                                                                                                                      - phi_fe_qp[SolFEType[pos_ctrl]][i] * 
                                                                                                                                        phi_fe_qp[SolFEType[pos_adj]][j] * 
                                                                                                                                        nonlin_term_function(sol_qp[pos_state])
                                                                                                    );
	        }
	      
	      else if ( control_el_flag == 0 && i == j)  {  
        Jac[ jac_row_col_index(Sol_n_el_dofs, nDof_AllVars, pos_ctrl, pos_ctrl, i, j) ]  += m_b_f[pos_ctrl][pos_ctrl] * (1 - control_node_flag[i]) * penalty_strong;
	      }
	      
	      //=========== delta_adjoint row ===========================
		Jac[ jac_row_col_index(Sol_n_el_dofs, nDof_AllVars, pos_adj, pos_state, i, j) ]  += m_b_f[pos_adj][pos_state] * weight_qp * (
                                                                                                                                      laplacian_row_col(SolFEType, phi_x_fe_qp, pos_adj, pos_state, i, j, dim)
                                                                                                                                      - phi_fe_qp[SolFEType[pos_adj]][i] *
                                                                                                                                        nonlin_term_derivative(phi_fe_qp[SolFEType[pos_state]][j]) * 
                                                                                                                                        sol_qp[pos_ctrl]
                                                                                                                                    );

        Jac[ jac_row_col_index(Sol_n_el_dofs, nDof_AllVars, pos_adj, pos_ctrl, i, j)  ]  += m_b_f[pos_adj][pos_ctrl] * weight_qp * ( 
                                                                                                                                      - phi_fe_qp[SolFEType[pos_adj]][i] * 
                                                                                                                                        phi_fe_qp[SolFEType[pos_ctrl]][j] * 
                                                                                                                                        nonlin_term_function(sol_qp[pos_state])
                                                                                                                                   ); 
		     
	      
            } // end phi_j loop
          } // endif assemble_matrix

        } // end phi_i loop
        
      } // end gauss point loop

      
      
    //std::vector<double> Res_ctrl (Sol_n_el_dofs[pos_ctrl]); std::fill(Res_ctrl.begin(),Res_ctrl.end(), 0.);
    for (unsigned i = 0; i < sol_eldofs[pos_ctrl].size(); i++) {
       unsigned n_els_that_node = 1;
     if ( control_el_flag == 1) {
// 	Res[Sol_n_el_dofs[pos_state] + i] += - ( + n_els_that_node * ineq_flag * sol_eldofs[pos_mu][i] /*- ( 0.4 + sin(M_PI * x[0][i]) * sin(M_PI * x[1][i]) )*/ );
// 	Res_ctrl[i] =  Res[Sol_n_el_dofs[pos_state] + i];
      }
    }
    

 //========== end of integral-based part

                          RES->add_vector_blocked(Res, L2G_dofmap_AllVars);
      if (assembleMatrix)  KK->add_matrix_blocked(Jac, L2G_dofmap_AllVars, L2G_dofmap_AllVars);
    
    
 //========== dof-based part, without summation
 
 //============= delta_mu row ===============================
      std::vector<double> Res_mu (Sol_n_el_dofs[pos_mu]); std::fill(Res_mu.begin(),Res_mu.end(), 0.);
      
    for (unsigned i = 0; i < sol_actflag.size(); i++) {
      if (sol_actflag[i] == 0) {  //inactive
         Res_mu [i] = (- ineq_flag) * ( 1. * sol_eldofs[pos_mu][i] ); 
// 	 Res_mu [i] = Res[res_row_index(Sol_n_el_dofs,pos_mu,i)]; 
      }
      else if (sol_actflag[i] == 1) {  //active_a 
	     Res_mu [i] = (- ineq_flag) * c_compl * ( sol_eldofs[pos_ctrl][i] - ctrl_lower[i]);
      }
      else if (sol_actflag[i] == 2) {  //active_b 
	     Res_mu [i] = (- ineq_flag) * c_compl * ( sol_eldofs[pos_ctrl][i] - ctrl_upper[i]);
      }
    }
    
    RES->insert(Res_mu, L2G_dofmap[pos_mu]);
    
 //============= delta_ctrl - mu ===============================
 KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap[pos_ctrl], L2G_dofmap[pos_mu], m_b_f[pos_ctrl][pos_mu] * ineq_flag * 1.);
  
 //============= delta_mu - ctrl row ===============================
 for (unsigned i = 0; i < sol_actflag.size(); i++) if (sol_actflag[i] != 0 ) sol_actflag[i] = m_b_f[pos_mu][pos_ctrl] * ineq_flag * c_compl;    
  
 KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap[pos_mu], L2G_dofmap[pos_ctrl], sol_actflag);

 //============= delta_mu - mu row ===============================
  for (unsigned i = 0; i < sol_actflag.size(); i++) sol_actflag[i] =  m_b_f[pos_mu][pos_mu] * ( ineq_flag * (1 - sol_actflag[i]/c_compl)  + (1-ineq_flag) * 1. ); //can do better to avoid division, maybe use modulo operator 

  KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap[pos_mu], L2G_dofmap[pos_mu], sol_actflag );
  
  } //end element loop for each process
  
  RES->close();

  if (assembleMatrix) KK->close();
//  std::ostringstream mat_out; mat_out << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "jacobian" << mlPdeSys->GetNonlinearIt()  << ".txt";
//   KK->print_matlab(mat_out.str(),"ascii"); 
//    KK->print();
  
  // ***************** END ASSEMBLY *******************
  
  unsigned int global_ctrl_size = pdeSys->KKoffset[pos_ctrl+1][iproc] - pdeSys->KKoffset[pos_ctrl][iproc];
  
  std::vector<double>  one_times_mu(global_ctrl_size, 0.);
  std::vector<int>        positions(global_ctrl_size);
//  double position_mu_i;
  for (unsigned i = 0; i < positions.size(); i++) {
    positions[i] = pdeSys->KKoffset[pos_ctrl][iproc] + i;
//     position_mu_i = pdeSys->KKoffset[pos_mu][iproc] + i;
//     std::cout << position_mu_i << std::endl;
    one_times_mu[i] =  m_b_f[pos_ctrl][pos_mu] * ineq_flag * 1. * (*sol->_Sol[SolIndex[pos_mu]])(i/*position_mu_i*/) ;
  }
    RES->add_vector_blocked(one_times_mu, positions);
//     RES->print();
    
  return;
}



void ComputeIntegral(const MultiLevelProblem& ml_prob)    {
  
  const NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystemWithPrimalDualActiveSetMethod> ("OptSys");
  const unsigned          level      = mlPdeSys->GetLevelToAssemble();

  Mesh*                          msh = ml_prob._ml_msh->GetLevel(level);

  MultiLevelSolution*          mlSol = ml_prob._ml_sol;
  Solution*                      sol = ml_prob._ml_sol->GetSolutionLevel(level);

  const unsigned     dim = msh->GetDimension();                                 // get the domain dimension of the problem
  const unsigned    dim2 = (3 * (dim - 1) + !(dim - 1));                        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  const unsigned   iproc = msh->processor_id(); 
  
 //*************************************************** 
  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
  for (unsigned i = 0; i < dim; i++) {
    x[i].reserve(maxSize);
  }
 //*************************************************** 

 //*************************************************** 
  double weight_qp; // gauss point weight
  
 //***************************************************  
  double alpha = ALPHA_CTRL_VOL;
  double beta  = BETA_CTRL_VOL;

 //******************** state ************************ 
 //*************************************************** 
  vector <double> phi_u;
  vector <double> phi_u_x; 
  vector <double> phi_u_xx;

  phi_u.reserve(maxSize);
  phi_u_x.reserve(maxSize * dim);
  phi_u_xx.reserve(maxSize * dim2);
  
 
  unsigned solIndex_u = mlSol->GetIndex("state");    // get the position of "state" in the ml_sol object
  unsigned solType_u  = mlSol->GetSolutionType(solIndex_u);    // get the finite element type for "state"

  vector < double >  sol_u; // local solution
  sol_u.reserve(maxSize);
  
  double u_gss = 0.;
 //*************************************************** 
 //*************************************************** 

 //******************** control ********************** 
 //*************************************************** 
  vector <double> phi_ctrl;
  vector <double> phi_ctrl_x;
  vector <double> phi_ctrl_xx;

  phi_ctrl.reserve(maxSize);
  phi_ctrl_x.reserve(maxSize * dim);
  phi_ctrl_xx.reserve(maxSize * dim2);
  
  unsigned solIndex_ctrl = mlSol->GetIndex("control");
  unsigned solType_ctrl  = mlSol->GetSolutionType(solIndex_ctrl);

  vector < double >  sol_ctrl;
  sol_ctrl.reserve(maxSize);
  
  double ctrl_gss = 0.;
  double ctrl_x_gss = 0.;
 //*************************************************** 
 //***************************************************  

  
 //********************* desired ********************* 
 //*************************************************** 
  vector <double> phi_udes;
  vector <double> phi_udes_x;
  vector <double> phi_udes_xx;

    phi_udes.reserve(maxSize);
    phi_udes_x.reserve(maxSize * dim);
    phi_udes_xx.reserve(maxSize * dim2);
 
  vector < double >  sol_udes; // local solution
  sol_udes.reserve(maxSize);
  
   double udes_gss = 0.;
 //*************************************************** 
 //*************************************************** 

 //*************************************************** 
 //********* WHOLE SET OF VARIABLES ****************** 
  const int solType_max = 2;  //biquadratic
 //*************************************************** 
  
 //********************* DATA ************************ 
  double u_des = DesiredTarget();
 //*************************************************** 
  
  double integral_target = 0.;
  double integral_alpha  = 0.;
  double integral_beta   = 0.;

    
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);    // element geometry type
 
 //***************** GEOMETRY ************************ 
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
    for (int i = 0; i < dim; i++)  x[i].resize(nDofx);
    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

   // elem average point 
    vector < double > elem_center(dim);   
   for (unsigned j = 0; j < dim; j++) {  elem_center[j] = 0.;  }
   for (unsigned j = 0; j < dim; j++) {  
      for (unsigned i = 0; i < nDofx; i++) {
         elem_center[j] += x[j][i];
       }
    }
    
   for (unsigned j = 0; j < dim; j++) { elem_center[j] = elem_center[j]/nDofx; }
 //*************************************************** 
  
 //************** set target domain flag *************
   int target_flag = 0;
   target_flag = ElementTargetFlag(elem_center);
 //*************************************************** 

   
 //**************** state **************************** 
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u);
        sol_u    .resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
     unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u);
            sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);
    }
 //*************************************************** 


 //************** control **************************** 
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);
    sol_ctrl    .resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl);
      sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);
    } 
 //***************************************************  
 
 
 //**************** u_des **************************** 
    unsigned nDof_udes  = msh->GetElementDofNumber(iel, solType_u);
    sol_udes    .resize(nDof_udes);
    for (unsigned i = 0; i < sol_udes.size(); i++) {
      sol_udes[i] = u_des;  //dof value
    } 
 //*************************************************** 

 
 //******************* ALL VARS ********************** 
    int nDof_max    =  nDof_u;   // TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
    
    if(nDof_udes > nDof_max) 
    {
      nDof_max = nDof_udes;
      }
    
    if(nDof_ctrl > nDof_max)
    {
      nDof_max = nDof_ctrl;
    }
    
 //*************************************************** 
   
      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType_max]->GetGaussPointNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
	msh->_finiteElement[ielGeom][solType_u]               ->Jacobian(x, ig, weight_qp, phi_u, phi_u_x, phi_u_xx);
    msh->_finiteElement[ielGeom][solType_u/*solTypeudes*/]->Jacobian(x, ig, weight_qp, phi_udes, phi_udes_x, phi_udes_xx);
    msh->_finiteElement[ielGeom][solType_ctrl]            ->Jacobian(x, ig, weight_qp, phi_ctrl, phi_ctrl_x, phi_ctrl_xx);

	      u_gss = 0.; for (unsigned i = 0; i < nDof_u; i++)       u_gss  += sol_u[i]    * phi_u[i];		
	   ctrl_gss = 0.; for (unsigned i = 0; i < nDof_ctrl; i++) ctrl_gss  += sol_ctrl[i] * phi_ctrl[i];  
	  udes_gss  = 0.; for (unsigned i = 0; i < nDof_udes; i++)  udes_gss += sol_udes[i] * phi_udes[i]; 
    ctrl_x_gss  = 0.; for (unsigned i = 0; i < nDof_ctrl; i++)  {
                      for (unsigned idim = 0; idim < dim; idim ++) ctrl_x_gss  += sol_ctrl[i] * phi_ctrl_x[i + idim * nDof_ctrl];
                        }

               integral_target += target_flag * weight_qp * (u_gss - udes_gss) * (u_gss - udes_gss);
               integral_alpha  +=       alpha * weight_qp * ctrl_gss * ctrl_gss;
               integral_beta   +=        beta * weight_qp * ctrl_x_gss * ctrl_x_gss;
	  
      } // end gauss point loop
  } //end element loop

  std::cout << "The value of the integral_target is " << std::setw(11) << std::setprecision(10) << integral_target << std::endl;
  std::cout << "The value of the integral_alpha  is " << std::setw(11) << std::setprecision(10) << integral_alpha << std::endl;
  std::cout << "The value of the integral_beta   is " << std::setw(11) << std::setprecision(10) << integral_beta << std::endl;
  std::cout << "The value of the total integral  is " << std::setw(11) << std::setprecision(10) << integral_target + 0.5 * alpha * integral_alpha + 0.5 * beta * integral_beta << std::endl;
  
return;
  
}
  
  

