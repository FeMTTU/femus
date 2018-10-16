#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "NumericVector.hpp"

#include "./elliptic_nonlin_param.hpp"

using namespace femus;

double InitialValueContReg(const std::vector < double >& x) {
  return ControlDomainFlag_internal_restriction(x);
}

double InitialValueTargReg(const std::vector < double >& x) {
  return ElementTargetFlag(x);
}

double InitialValueState(const std::vector < double >& x) {
  return 0.;
}

double InitialValueAdjoint(const std::vector < double >& x) {
  return 0.;
}

double InitialValueMu(const std::vector < double >& x) {
  return 0.;
}

double InitialValueControl(const std::vector < double >& x) {
  return 0.;
}

bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int faceName, const double time) {

  bool dirichlet = true; //dirichlet
  value = 0.;
  
  if(!strcmp(name,"control")) {
      value = 0.;
  if (faceName == 3)
    dirichlet = false;
  
  }
  
  if(!strcmp(name,"mu")) {
//       value = 0.;
//   if (faceName == 3)
    dirichlet = false;
  }
  
  return dirichlet;
}


double ComputeIntegral(MultiLevelProblem& ml_prob);

void AssembleProblem(MultiLevelProblem& ml_prob);



int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
    // ======= Files ========================
  Files files; 
        files.CheckIODirectories();
	files.RedirectCout();

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;

  mlMsh.GenerateCoarseBoxMesh(NSUB_X,NSUB_Y,0,0.,1.,0.,1.,0.,0.,QUAD9,"seventh");
 /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the furure it is not going to be an argument of this function   */
  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("state", LAGRANGE, FIRST);
  mlSol.AddSolution("control", LAGRANGE, FIRST);
  mlSol.AddSolution("adjoint", LAGRANGE, FIRST);
  mlSol.AddSolution("mu", LAGRANGE, FIRST);  
  mlSol.AddSolution("TargReg",  DISCONTINOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  mlSol.AddSolution("ContReg",  DISCONTINOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field

  
  mlSol.Initialize("All");    // initialize all varaibles to zero

  mlSol.Initialize("state", InitialValueState);
  mlSol.Initialize("control", InitialValueControl);
  mlSol.Initialize("adjoint", InitialValueAdjoint);
  mlSol.Initialize("mu", InitialValueMu);
  mlSol.Initialize("TargReg", InitialValueTargReg);
  mlSol.Initialize("ContReg", InitialValueContReg);

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("state");
  mlSol.GenerateBdc("control");
  mlSol.GenerateBdc("adjoint");
  mlSol.GenerateBdc("mu");  //we need add this to make the matrix iterations work...

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);
  
  mlProb.SetFilesHandler(&files);

 // add system  in mlProb as a Linear Implicit System
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("OptSys");

  system.AddSolutionToSystemPDE("state");  
  system.AddSolutionToSystemPDE("control");  
  system.AddSolutionToSystemPDE("adjoint");  
  system.AddSolutionToSystemPDE("mu");  
  
  // attach the assembling function to system
  system.SetAssembleFunction(AssembleProblem);
  
  mlSol.SetWriter(VTK);
  mlSol.GetWriter()->SetDebugOutput(true);
  
  system.SetDebugNonlinear(true);
//   system.SetMaxNumberOfNonLinearIterations(2);

  // initialize and solve the system
  system.init();
  system.MGsolve();
  
  ComputeIntegral(mlProb);
 
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("all");
  mlSol.GetWriter()->Write(files.GetOutputPath()/*DEFAULT_OUTPUTDIR*/, "biquadratic", variablesToBePrinted);

  return 0;
}


void AssembleProblem(MultiLevelProblem& ml_prob) {

  // ************** J dx = f - J x_old ******************  
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  NonLinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystem> ("OptSys");
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
  vector < vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
  for (unsigned i = 0; i < dim; i++) {
    coordX[i].reserve(maxSize);
  }
 //***************************************************   

 //***************************************************  
  double weight; // gauss point weight
  
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

  const int n_unknowns = 4;
 
  enum Sol_pos{pos_state=0,pos_ctrl,pos_adj,pos_mu};
  
  vector < std::string > Solname(n_unknowns);
  Solname[pos_state] = "state";
  Solname[pos_ctrl]  = "control";
  Solname[pos_adj]   = "adjoint";
  Solname[pos_mu]    = "mu";
  
 int m_b_f[n_unknowns][n_unknowns];   //m_b_f
     m_b_f[0][0] = 1;
     m_b_f[0][1] = 1;
     m_b_f[0][2] = 1;
     m_b_f[0][3] = 1;
     m_b_f[1][0] = 1;
     m_b_f[1][1] = 1;
     m_b_f[1][2] = 1;
     m_b_f[1][3] = 1;
     m_b_f[2][0] = 1;
     m_b_f[2][1] = 1;
     m_b_f[2][2] = 1;
     m_b_f[2][3] = 1;
     m_b_f[3][0] = 1;
     m_b_f[3][1] = 1;
     m_b_f[3][2] = 1;
     m_b_f[3][3] = 1;
 
  //***************************************************  
  vector < vector < int > > L2G_dofmap(n_unknowns);     for(int i = 0; i < n_unknowns; i++) { L2G_dofmap[i].reserve(maxSize); }
            vector< int >   L2G_dofmap_AllVars; L2G_dofmap_AllVars.reserve( n_unknowns*maxSize );
            vector< double >         Res;                      Res.reserve( n_unknowns*maxSize );
            vector< double >         Jac;                      Jac.reserve( n_unknowns*maxSize * n_unknowns*maxSize);
  //***************************************************  
  
  //***************************************************  
  vector < unsigned > SolPdeIndex(n_unknowns);
  vector < unsigned > SolIndex(n_unknowns);  
  vector < unsigned > SolFEType(n_unknowns);  

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
  
  //********* variables for ineq constraints *****************
  const int ineq_flag = INEQ_FLAG;
  const double ctrl_lower =  CTRL_BOX_LOWER;
  const double ctrl_upper =  CTRL_BOX_UPPER;
  assert(ctrl_lower < ctrl_upper);
  const double c_compl = C_COMPL;
  vector < double/*int*/ >  sol_actflag;   sol_actflag.reserve(maxSize); //flag for active set
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
   
   //STATE###################################################################  
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
  //CTRL###################################################################
    
    
 //************** update active set flag for current nonlinear iteration **************************** 
 // 0: inactive; 1: active_a; 2: active_b
   assert(Sol_n_el_dofs[pos_mu] == Sol_n_el_dofs[pos_ctrl]);
   sol_actflag.resize(Sol_n_el_dofs[pos_mu]);
     std::fill(sol_actflag.begin(), sol_actflag.end(), 0);
   
    for (unsigned i = 0; i < sol_actflag.size(); i++) {  
    if      ( (sol_eldofs[pos_mu][i] + c_compl * (sol_eldofs[pos_ctrl][i] - ctrl_lower )) < 0 )  sol_actflag[i] = 1;
    else if ( (sol_eldofs[pos_mu][i] + c_compl * (sol_eldofs[pos_ctrl][i] - ctrl_upper )) > 0 )  sol_actflag[i] = 2;
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
         ml_prob._ml_msh->_finiteElement[ielGeom][fe]->Jacobian(coordX,ig,weight,phi_fe_qp[fe],phi_x_fe_qp[fe],phi_xx_fe_qp[fe]);
      }
   //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
    ml_prob._ml_msh->_finiteElement[ielGeom][BIQUADR_FE]->Jacobian(coordX,ig,weight,phi_fe_qp[BIQUADR_FE],phi_x_fe_qp[BIQUADR_FE],phi_xx_fe_qp[BIQUADR_FE]);
    
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
	      double laplace_rhs_dstate_adj_i = 0.; //
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < Sol_n_el_dofs[pos_state] )         laplace_rhs_dstate_adj_i             +=  (phi_x_fe_qp[SolFEType[pos_state]][i * dim + kdim] * sol_grad_qp[pos_adj][kdim]);
	      }
	      
              double laplace_rhs_dctrl_ctrl_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < Sol_n_el_dofs[pos_ctrl] )         laplace_rhs_dctrl_ctrl_i      +=  (phi_x_fe_qp[SolFEType[pos_ctrl]]   [i * dim + kdim] * sol_grad_qp[pos_ctrl][kdim]);
	      }
	      
	      double laplace_rhs_dctrl_adj_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < Sol_n_el_dofs[pos_ctrl] )         laplace_rhs_dctrl_adj_i       +=  (phi_x_fe_qp[SolFEType[pos_ctrl]]   [i * dim + kdim] * sol_grad_qp[pos_adj][kdim]);
	      }
	      
	      double laplace_rhs_dadj_state_i = 0.;  //
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < Sol_n_el_dofs[pos_adj] )         laplace_rhs_dadj_state_i           +=  (phi_x_fe_qp[SolFEType[pos_adj]][i * dim + kdim] * sol_grad_qp[pos_state][kdim]);
	      }
	      
	      double laplace_rhs_dadj_ctrl_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < Sol_n_el_dofs[pos_adj] )         laplace_rhs_dadj_ctrl_i        +=  (phi_x_fe_qp[SolFEType[pos_adj]][i * dim + kdim] * sol_grad_qp[pos_ctrl][kdim]);
	      }
	      
          // FIRST ROW
	  if (i < Sol_n_el_dofs[pos_state])                      Res[0      + i] += - weight * (target_flag * phi_fe_qp[SolFEType[pos_state]][i] * (
                                                                                              m_b_f[pos_state][pos_state] * sol_qp[pos_state] 
                                                                                            + m_b_f[pos_state][pos_ctrl]  * sol_qp[pos_ctrl] 
                                                                                            - u_des ) 
                                                                                          - m_b_f[pos_state][pos_adj]   * laplace_rhs_dstate_adj_i );
          // SECOND ROW
	  if (i < Sol_n_el_dofs[pos_ctrl])  {
	     if ( control_el_flag == 1)        Res[Sol_n_el_dofs[pos_state] + i] +=  /*(control_node_flag[i]) **/ - weight * (
                                                                                 target_flag * phi_fe_qp[SolFEType[pos_ctrl]][i] * (
                                                                                      m_b_f[pos_ctrl][pos_state] * sol_qp[pos_state] 
                                                                                    + m_b_f[pos_ctrl][pos_ctrl]  * sol_qp[pos_ctrl] 
                                                                                    - u_des ) 
                                                                                  + m_b_f[pos_ctrl][pos_ctrl] * alpha * phi_fe_qp[SolFEType[pos_ctrl]][i] * sol_qp[pos_ctrl]
		                                                                          + m_b_f[pos_ctrl][pos_ctrl] *  beta * laplace_rhs_dctrl_ctrl_i 
		                                                                          - m_b_f[pos_ctrl][pos_adj]  * laplace_rhs_dctrl_adj_i 
                                                                                                                         );
	      else if ( control_el_flag == 0)  Res[Sol_n_el_dofs[pos_state] + i] +=  /*(1 - control_node_flag[i]) **/ m_b_f[pos_ctrl][pos_ctrl] * (- penalty_strong) * (sol_eldofs[pos_ctrl][i]);
	  }
          // THIRD ROW
          if (i < Sol_n_el_dofs[pos_adj])  Res[Sol_n_el_dofs[pos_state] + Sol_n_el_dofs[pos_ctrl] + i] += /*-weight * phi_adj[i] * sol_qp[pos_adj] - 6.;*/
                                                                                                          - weight *  ( - m_b_f[pos_adj][pos_state] * laplace_rhs_dadj_state_i 
                                                                                                                        - m_b_f[pos_adj][pos_ctrl]  * laplace_rhs_dadj_ctrl_i );

//======================Residuals=======================
	      
          if (assembleMatrix) {
	    
            // *** phi_j loop ***
            for (unsigned j = 0; j < nDof_max; j++) {
                
              double laplace_mat_dstate_adj = 0.;
              double laplace_mat_dadj_state = 0.;
              double laplace_mat_dctrl_adj = 0.;
              double laplace_mat_dadj_ctrl = 0.;
              double laplace_mat_dctrl_ctrl = 0.;

              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < Sol_n_el_dofs[pos_state] && 
                   j < Sol_n_el_dofs[pos_adj]  )     laplace_mat_dstate_adj     += (phi_x_fe_qp[SolFEType[pos_state]]   [i * dim + kdim] * 
                                                                                    phi_x_fe_qp[SolFEType[pos_adj]][j * dim + kdim]);
              if ( i < Sol_n_el_dofs[pos_adj]  &&
                   j < Sol_n_el_dofs[pos_state] )    laplace_mat_dadj_state     += (phi_x_fe_qp[SolFEType[pos_adj]][i * dim + kdim] *
                                                                                    phi_x_fe_qp[SolFEType[pos_state]] [j * dim + kdim]);  //equal to the previous
              if ( i < Sol_n_el_dofs[pos_ctrl] &&
                   j < Sol_n_el_dofs[pos_adj] )      laplace_mat_dctrl_adj      += (phi_x_fe_qp[SolFEType[pos_ctrl]][i * dim + kdim] *
                                                                                    phi_x_fe_qp[SolFEType[pos_adj]][j * dim + kdim]);
              if ( i < Sol_n_el_dofs[pos_adj]  && 
                   j < Sol_n_el_dofs[pos_ctrl] )     laplace_mat_dadj_ctrl      += (phi_x_fe_qp[SolFEType[pos_adj]][i * dim + kdim] *
                                                                                    phi_x_fe_qp[SolFEType[pos_ctrl]][j * dim + kdim]);  //equal to the previous
              if ( i < Sol_n_el_dofs[pos_ctrl] &&
                   j < Sol_n_el_dofs[pos_ctrl] )     laplace_mat_dctrl_ctrl     += (phi_x_fe_qp[SolFEType[pos_ctrl]]  [i * dim + kdim] *
                                                                                    phi_x_fe_qp[SolFEType[pos_ctrl]]  [j * dim + kdim]);
	      }

              //============ delta_state row ============================
              //DIAG BLOCK delta_state - state
	      if ( i < Sol_n_el_dofs[pos_state] && 
               j < Sol_n_el_dofs[pos_state] )       
		Jac[ (0 + i) * nDof_AllVars   +
		     (0 + j)                            ]  += m_b_f[pos_state][pos_state] * weight * target_flag * phi_fe_qp[SolFEType[pos_state]][j] *  phi_fe_qp[SolFEType[pos_state]][i];
              
	      // BLOCK  delta_state - control
              if ( i < Sol_n_el_dofs[pos_state] && 
                   j < Sol_n_el_dofs[pos_ctrl] )   
		Jac[ (0 + i) * nDof_AllVars   +
                     (Sol_n_el_dofs[pos_state] + j)                       ]  += m_b_f[pos_state][pos_ctrl] * weight * target_flag * phi_fe_qp[SolFEType[pos_state]][i] * phi_fe_qp[SolFEType[pos_ctrl]][j];
	      
              // BLOCK  delta_state - adjoint
              if ( i < Sol_n_el_dofs[pos_state] && 
                   j < Sol_n_el_dofs[pos_adj] )  
		Jac[  (0 + i) * nDof_AllVars  +
                      (Sol_n_el_dofs[pos_state] + Sol_n_el_dofs[pos_ctrl] + j)   ]  += m_b_f[pos_state][pos_adj] * weight * (-1) * laplace_mat_dstate_adj;
              
	      //=========== delta_control row ===========================     
	      if ( control_el_flag == 1)  {

	      //BLOCK delta_control - state
              if ( i < Sol_n_el_dofs[pos_ctrl]   && 
                   j < Sol_n_el_dofs[pos_state]   ) 
		Jac[ (Sol_n_el_dofs[pos_state] + i) * nDof_AllVars  +
		     (0 + j)                            ]  += m_b_f[pos_ctrl][pos_state] * ( control_node_flag[i]) * weight * target_flag * phi_fe_qp[SolFEType[pos_ctrl]][i] * phi_fe_qp[SolFEType[pos_state]][j];
		
	      //BLOCK delta_control - control
              if ( i < Sol_n_el_dofs[pos_ctrl]   &&
                   j < Sol_n_el_dofs[pos_ctrl]   )
		Jac[ (Sol_n_el_dofs[pos_state] + i) * nDof_AllVars +
		     (Sol_n_el_dofs[pos_state] + j)                       ]  += m_b_f[pos_ctrl][pos_ctrl] * ( control_node_flag[i]) * weight * (
                                                                                           beta * control_el_flag  * laplace_mat_dctrl_ctrl 
		                                                                                + alpha * control_el_flag * phi_fe_qp[SolFEType[pos_ctrl]][i] * phi_fe_qp[SolFEType[pos_ctrl]][j] 
		                                                                                            + target_flag * phi_fe_qp[SolFEType[pos_ctrl]][i] * phi_fe_qp[SolFEType[pos_ctrl]][j] );
              
	      //BLOCK delta_control - adjoint
              if ( i < Sol_n_el_dofs[pos_ctrl]  &&
                   j < Sol_n_el_dofs[pos_adj]  ) 
		Jac[ (Sol_n_el_dofs[pos_state] + i) * nDof_AllVars  + 
		     (Sol_n_el_dofs[pos_state] + Sol_n_el_dofs[pos_ctrl] + j)           ]  += m_b_f[pos_ctrl][pos_adj] * ( control_node_flag[i]) * weight * (-1) * laplace_mat_dctrl_adj;
	      	      
	        }
	      
	      else if ( control_el_flag == 0)  {  
		
              //BLOCK delta_control - control
              if ( i < Sol_n_el_dofs[pos_ctrl]  &&
                   j < Sol_n_el_dofs[pos_ctrl]  &&  
                   i == j ) {
		 Jac[ (Sol_n_el_dofs[pos_state] + i) * nDof_AllVars +
		      (Sol_n_el_dofs[pos_state] + j)                      ]  += m_b_f[pos_ctrl][pos_ctrl] * (1 - control_node_flag[i]) * penalty_strong;
		    }
	      
	      }
	      
		     
	      //=========== delta_adjoint row ===========================
              // BLOCK delta_adjoint - state	      
              if ( i < Sol_n_el_dofs[pos_adj] && 
                   j < Sol_n_el_dofs[pos_state] )   
		Jac[ (Sol_n_el_dofs[pos_state] + Sol_n_el_dofs[pos_ctrl] + i) * nDof_AllVars +
		     (0 + j)                            ]  += m_b_f[pos_adj][pos_state] * weight * (-1) * laplace_mat_dadj_state;   
	      
              // BLOCK delta_adjoint - control   
              if ( i < Sol_n_el_dofs[pos_adj] &&
                   j < Sol_n_el_dofs[pos_ctrl] )  
		Jac[ (Sol_n_el_dofs[pos_state] + Sol_n_el_dofs[pos_ctrl] + i)  * nDof_AllVars +
		     (Sol_n_el_dofs[pos_state]  + j)                      ]  += m_b_f[pos_adj][pos_ctrl] * weight * (-1) * laplace_mat_dadj_ctrl; 
		     
		     
// 	      // BLOCK delta_adjoint - adjoint   
//               if ( i < Sol_n_el_dofs[pos_adj] && j < Sol_n_el_dofs[pos_adj] )  
// 		Jac[ (Sol_n_el_dofs[pos_state] + Sol_n_el_dofs[pos_ctrl] + i)  * nDof_AllVars +
// 		     (Sol_n_el_dofs[pos_state] + Sol_n_el_dofs[pos_ctrl] + j)                      ]  += weight * phi_adj[j] *  phi_adj[i]; 
    
	      
	      //============= delta_mu row ===============================
//	      if (sol_actflag[i] == 0) //inactive
//	      { // BLOCK delta_mu - mu	      
// 	        if ( i < Sol_n_el_dofs[pos_mu] && j < Sol_n_el_dofs[pos_mu] && i==j )   
// 		  Jac[ (Sol_n_el_dofs[pos_state] + Sol_n_el_dofs[pos_ctrl] + Sol_n_el_dofs[pos_adj] + i) * nDof_AllVars +
// 		       (Sol_n_el_dofs[pos_state] + Sol_n_el_dofs[pos_ctrl] + Sol_n_el_dofs[pos_adj] + j)]  = 1. ;  
// 	     // }
// 	      else //active
// 	      { // BLOCK delta_mu - ctrl	      
//                 if ( i < Sol_n_el_dofs[pos_mu] && j < Sol_n_el_dofs[pos_ctrl] && i==j )   
// 		  Jac[ (Sol_n_el_dofs[pos_state] + Sol_n_el_dofs[pos_ctrl] + Sol_n_el_dofs[pos_adj] + i) * nDof_AllVars +
// 		       (Sol_n_el_dofs[pos_state] + j)                       ]  = c_compl * 1. ; 
	     // }
	      
            } // end phi_j loop
          } // endif assemble_matrix

        } // end phi_i loop
        
      } // end gauss point loop

      
      
    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
   // std::cout << " ************* Element ************** " << iel << " **************************************** " << std::endl;     

// // //     if (control_el_flag == 0) {  //elements that should have zero control
// //          for (unsigned i_unk = 0; i_unk < n_unknowns; i_unk++) {
// //     std::cout << " ======= Row === " << i_unk << " =================================================== " << std::endl;     
// //         unsigned int row_block_offset = 0;
// // 	         for (unsigned k = 0; k < i_unk; k++) row_block_offset += Sol_n_el_dofs[k];
// //          for (unsigned j_unk = 0; j_unk < n_unknowns; j_unk++) {
// //     std::cout << " ======= Column === " << j_unk << " ================== " << std::endl;     
// //         unsigned int column_block_offset = 0;
// // 	         for (unsigned k = 0; k < j_unk; k++) column_block_offset += Sol_n_el_dofs[k];
// // 	  
// //          for (unsigned i = 0; i < Sol_n_el_dofs[i_unk]; i++) {
// // // 	      std::cout << Res[Sol_n_el_dofs[pos_state] + Sol_n_el_dofs[pos_ctrl] + Sol_n_el_dofs[pos_adj] + i ] << " " << std::endl;
// // 	   for (unsigned j = 0; j < Sol_n_el_dofs[j_unk]; j++) {
// // 	      std::cout <<  " " << std::setfill(' ') << std::setw(10) << Jac[ (row_block_offset + i) * nDof_AllVars + ( column_block_offset + j) ] << " ";
// // 	    }
// // 	      std::cout << std::endl;
// // 	 }
// // 
// // 	 } //j_unk
// // 	} //i_unk
	 
	 
// // // 	}
    //std::vector<double> Res_ctrl (Sol_n_el_dofs[pos_ctrl]); std::fill(Res_ctrl.begin(),Res_ctrl.end(), 0.);
    for (unsigned i = 0; i < sol_eldofs[pos_ctrl].size(); i++){
       unsigned n_els_that_node = 1;
     if ( control_el_flag == 1) {
// 	Res[Sol_n_el_dofs[pos_state] + i] += - ( + n_els_that_node * ineq_flag * sol_eldofs[pos_mu][i] /*- ( 0.4 + sin(M_PI * x[0][i]) * sin(M_PI * x[1][i]) )*/ );
// 	Res_ctrl[i] =  Res[Sol_n_el_dofs[pos_state] + i];
      }
    }//------------------------------->>>>>>
    
//     std::vector<double> Res_u (Sol_n_el_dofs[pos_state]); std::fill(Res_u.begin(),Res_u.end(), 0.);
//     for (unsigned i = 0; i < sol_eldofs[pos_state].size(); i++){
// 	Res[0 + i] = - ( sol_eldofs[pos_state][i] - 8. );
// 	Res_u[i] = Res[0 + i];
//     }
    
//     std::vector<double> Res_adj (Sol_n_el_dofs[pos_adj]); std::fill(Res_adj.begin(),Res_adj.end(), 0.);
//     for (unsigned i = 0; i < sol_eldofs[pos_adj].size(); i++){
// 	Res[Sol_n_el_dofs[pos_state] + Sol_n_el_dofs[pos_ctrl] + i] = - (sol_eldofs[pos_adj][i] - 7.);
// 	Res_adj[i] = Res[Sol_n_el_dofs[pos_state] + Sol_n_el_dofs[pos_ctrl] + i];
//     }
    

 //========== end of integral-based part

    //copy the value of the adept::adoube aRes in double Res and store
                          RES->add_vector_blocked(Res, L2G_dofmap_AllVars);
      if (assembleMatrix)  KK->add_matrix_blocked(Jac, L2G_dofmap_AllVars, L2G_dofmap_AllVars);
    
    
 //========== dof-based part, without summation
 
 //============= delta_mu row ===============================
      std::vector<double> Res_mu (Sol_n_el_dofs[pos_mu]); std::fill(Res_mu.begin(),Res_mu.end(), 0.);
      
    for (unsigned i = 0; i < sol_actflag.size(); i++) {
      if (sol_actflag[i] == 0) {  //inactive
         Res_mu [i] = (- ineq_flag) * ( 1. * sol_eldofs[pos_mu][i] ); 
// 	 Res_mu [i] = Res[Sol_n_el_dofs[pos_state] + Sol_n_el_dofs[pos_ctrl] + Sol_n_el_dofs[pos_adj] + i]; 
      }
      else if (sol_actflag[i] == 1) {  //active_a 
	     Res_mu [i] = (- ineq_flag) * c_compl * ( sol_eldofs[pos_ctrl][i] - ctrl_lower);
      }
      else if (sol_actflag[i] == 2) {  //active_b 
	     Res_mu [i] = (- ineq_flag) * c_compl * ( sol_eldofs[pos_ctrl][i] - ctrl_upper);
      }
    }

    
    RES->insert(Res_mu, L2G_dofmap[pos_mu]);
//     RES->insert(Res_ctrl, l2GMap_ctrl);
//     RES->insert(Res_u, l2GMap_u);
//     RES->insert(Res_adj, l2GMap_adj);
    
//  //============= delta_state-delta_state row ===============================
//  KK->matrix_set_off_diagonal_values_blocked(l2GMap_u, l2GMap_u, 1.);

//  //============= delta_ctrl-delta_ctrl row ===============================
//  KK->matrix_set_off_diagonal_values_blocked(l2GMap_ctrl, l2GMap_ctrl, 1.);
 
//  //============= delta_adj-delta_adj row ===============================
//  KK->matrix_set_off_diagonal_values_blocked(l2GMap_adj, l2GMap_adj, 1.);
  
 //============= delta_ctrl-delta_mu row ===============================
 KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap[pos_ctrl], L2G_dofmap[pos_mu], m_b_f[pos_ctrl][pos_mu] * ineq_flag * 1.);//------------------------------->>>>>>
  
 //============= delta_mu-delta_ctrl row ===============================
 for (unsigned i = 0; i < sol_actflag.size(); i++) if (sol_actflag[i] != 0 ) sol_actflag[i] = m_b_f[pos_mu][pos_ctrl] * ineq_flag * c_compl;    
  
 KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap[pos_mu], L2G_dofmap[pos_ctrl], sol_actflag);

 //============= delta_mu-delta_mu row ===============================
  for (unsigned i = 0; i < sol_actflag.size(); i++) sol_actflag[i] =  m_b_f[pos_mu][pos_mu] * ( ineq_flag * (1 - sol_actflag[i]/c_compl)  + (1-ineq_flag) * 1. ); //can do better to avoid division, maybe use modulo operator 

  KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap[pos_mu], L2G_dofmap[pos_mu], sol_actflag );
  
  } //end element loop for each process
  
  RES->close();

  if (assembleMatrix) KK->close();
//  std::ostringstream mat_out; mat_out << "matrix" << mlPdeSys->_nonliniteration  << ".txt";
//   KK->print_matlab(mat_out.str(),"ascii"); //  KK->print();
  
  // ***************** END ASSEMBLY *******************
  unsigned int ctrl_index = mlPdeSys->GetSolPdeIndex("control");
  unsigned int   mu_index = mlPdeSys->GetSolPdeIndex("mu");

  unsigned int global_ctrl_size = pdeSys->KKoffset[ctrl_index+1][iproc] - pdeSys->KKoffset[ctrl_index][iproc];
  
  std::vector<double>  one_times_mu(global_ctrl_size, 0.);
  std::vector<int>    positions(global_ctrl_size);
//  double position_mu_i;
  for (unsigned i = 0; i < positions.size(); i++) {
    positions[i] = pdeSys->KKoffset[ctrl_index][iproc] + i;
//     position_mu_i = pdeSys->KKoffset[mu_index][iproc] + i;
//     std::cout << position_mu_i << std::endl;
    one_times_mu[i] =  m_b_f[pos_ctrl][pos_mu] * ineq_flag * 1. * (*sol->_Sol[SolIndex[pos_mu]])(i/*position_mu_i*/) ;
  }
    RES->add_vector_blocked(one_times_mu, positions);
//     RES->print();
    
  return;
}



double ComputeIntegral(MultiLevelProblem& ml_prob)    {
  
  NonLinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystem> ("OptSys");
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
  double weight; // gauss point weight
  
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
  
 
  unsigned solIndex_u;
  solIndex_u = mlSol->GetIndex("state");    // get the position of "state" in the ml_sol object
  unsigned solType_u = mlSol->GetSolutionType(solIndex_u);    // get the finite element type for "state"

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
  
  unsigned solIndex_ctrl;
  solIndex_ctrl = mlSol->GetIndex("control");
  unsigned solType_ctrl = mlSol->GetSolutionType(solIndex_ctrl);

  vector < double >  sol_ctrl;
  sol_ctrl.reserve(maxSize);
//   vector< int > l2GMap_ctrl;
//   l2GMap_ctrl.reserve(maxSize);
  
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
 
//  unsigned solIndex_udes;
//   solIndex_udes = mlSol->GetIndex("Tdes");    // get the position of "state" in the ml_sol object
//   unsigned solType_udes = mlSol->GetSolutionType(solIndex_udes);    // get the finite element type for "state"

  vector < double >  sol_udes; // local solution
  sol_udes.reserve(maxSize);
//   vector< int > l2GMap_udes;
//   l2GMap_udes.reserve(maxSize);
   double udes_gss = 0.;
 //*************************************************** 
 //*************************************************** 


 //*************************************************** 
 //********* WHOLE SET OF VARIABLES ****************** 
  const int solType_max = 2;  //biquadratic

  const int n_unknowns = 4;
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
	msh->_finiteElement[ielGeom][solType_u]               ->Jacobian(x, ig, weight, phi_u, phi_u_x, phi_u_xx);
    msh->_finiteElement[ielGeom][solType_u/*solTypeudes*/]->Jacobian(x, ig, weight, phi_udes, phi_udes_x, phi_udes_xx);
    msh->_finiteElement[ielGeom][solType_ctrl]            ->Jacobian(x, ig, weight, phi_ctrl, phi_ctrl_x, phi_ctrl_xx);

	      u_gss = 0.; for (unsigned i = 0; i < nDof_u; i++)       u_gss  += sol_u[i]    * phi_u[i];		
	   ctrl_gss = 0.; for (unsigned i = 0; i < nDof_ctrl; i++) ctrl_gss  += sol_ctrl[i] * phi_ctrl[i];  
	  udes_gss  = 0.; for (unsigned i = 0; i < nDof_udes; i++)  udes_gss += sol_udes[i] * phi_udes[i]; 
    ctrl_x_gss  = 0.; for (unsigned i = 0; i < nDof_ctrl; i++)  {
                      for (unsigned idim = 0; idim < dim; idim ++) ctrl_x_gss  += sol_ctrl[i] * phi_ctrl_x[i + idim * nDof_ctrl];
                        }

               integral_target += target_flag * weight * (u_gss +  ctrl_gss - udes_gss) * (u_gss +  ctrl_gss - udes_gss);
               integral_alpha  += target_flag * alpha * weight * ctrl_gss * ctrl_gss;
               integral_beta   += target_flag * beta * weight * ctrl_x_gss * ctrl_x_gss;
	  
      } // end gauss point loop
  } //end element loop

  std::cout << "The value of the integral_target is " << std::setw(11) << std::setprecision(10) << integral_target << std::endl;
  std::cout << "The value of the integral_alpha  is " << std::setw(11) << std::setprecision(10) << integral_alpha << std::endl;
  std::cout << "The value of the integral_beta   is " << std::setw(11) << std::setprecision(10) << integral_beta << std::endl;
  std::cout << "The value of the total integral  is " << std::setw(11) << std::setprecision(10) << integral_target + integral_alpha + integral_beta << std::endl;
  
return integral_target + integral_alpha + integral_beta;
  
}
  
  

