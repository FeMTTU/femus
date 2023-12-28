 
#ifndef OPT_SYSTEMS_NEUMANN_HPP
#define OPT_SYSTEMS_NEUMANN_HPP



#include  "03_opt_system_inequalities.hpp"

#include "SparseMatrix.hpp"


namespace femus  {

namespace elliptic  {

class pure_boundary  {

  public:
    
static void assemble_elliptic_neumann_control(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data

  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use

  NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystemWithPrimalDualActiveSetMethod> ("LiftRestr");   // pointer to the linear implicit system named "LiftRestr"
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->GetMeshElements();  // pointer to the elem object in msh (level)

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

 //*************************************************** 
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
  std::vector < std::vector < double > > x(dim);
  std::vector < std::vector < double> >  x_bdry(dim);
  for (unsigned i = 0; i < dim; i++) {
         x[i].reserve(maxSize);
	 x_bdry[i].reserve(maxSize);
  }

 //***************************************************  

 //***************************************************  
  double weight = 0.; // gauss point weight
  double weight_bdry = 0.; // gauss point weight on the boundary


 //******************** state ************************ 
 //*************************************************** 
  std::vector <double> phi_u;  // local test function
  std::vector <double> phi_u_x; // local test function first order partial derivatives
  std::vector <double> phi_u_xx; // local test function second order partial derivatives

  phi_u.reserve(maxSize);
  phi_u_x.reserve(maxSize * dim);
  phi_u_xx.reserve(maxSize * dim2);
  
  unsigned solIndex_u    = ml_sol->GetIndex("state");    // get the position of "state" in the ml_sol object
  unsigned solType_u     = ml_sol->GetSolutionType(solIndex_u);    // get the finite element type for "state"
  unsigned solPdeIndex_u = mlPdeSys->GetSolPdeIndex("state");    // get the position of "state" in the pdeSys object

  std::vector < double >  sol_u;     sol_u.reserve(maxSize);
  std::vector < int > l2GMap_u;    l2GMap_u.reserve(maxSize);
 //***************************************************  
 //***************************************************  

  
 //********************* adjoint *********************
 //*************************************************** 
  std::vector <double> phi_adj;  // local test function
  std::vector <double> phi_adj_x; // local test function first order partial derivatives
  std::vector <double> phi_adj_xx; // local test function second order partial derivatives

  phi_adj.reserve(maxSize);
  phi_adj_x.reserve(maxSize * dim);
  phi_adj_xx.reserve(maxSize * dim2);
 
  unsigned solIndex_adj    = ml_sol->GetIndex("adjoint");    // get the position of "state" in the ml_sol object
  unsigned solType_adj     = ml_sol->GetSolutionType(solIndex_adj);    // get the finite element type for "state"
  unsigned solPdeIndex_adj = mlPdeSys->GetSolPdeIndex("adjoint");    // get the position of "state" in the pdeSys object

  std::vector < double >  sol_adj;   sol_adj.reserve(maxSize);
  std::vector < int > l2GMap_adj; l2GMap_adj.reserve(maxSize);

  //boundary adjoint shape functions  
  std::vector <double> phi_adj_bdry;  
  std::vector <double> phi_adj_x_bdry; 

  phi_adj_bdry.reserve(maxSize);
  phi_adj_x_bdry.reserve(maxSize * dim);

 //*************************************************** 
 //***************************************************  

  
 //********************* bdry cont *******************
 //*************************************************** 
  std::vector <double> phi_ctrl_bdry;  
  std::vector <double> phi_ctrl_x_bdry; 

  phi_ctrl_bdry.reserve(maxSize);
  phi_ctrl_x_bdry.reserve(maxSize * dim);
  
  unsigned solIndex_ctrl = ml_sol->GetIndex("control");
  unsigned solType_ctrl = ml_sol->GetSolutionType(solIndex_ctrl);

  unsigned solPdeIndex_ctrl = mlPdeSys->GetSolPdeIndex("control");

  std::string ctrl_name = "control";
  
 std::vector < double >  sol_ctrl;   sol_ctrl.reserve(maxSize);
 std::vector < int > l2GMap_ctrl;   l2GMap_ctrl.reserve(maxSize);
 //*************************************************** 
 //*************************************************** 
  
 
     const unsigned int n_components_ctrl = 1;

     
 //****************** mu ******************************  
 //***************************************************  
  unsigned solIndex_mu;
  solIndex_mu = ml_sol->GetIndex("mu");    // get the position of "mu" in the ml_sol object
   
  unsigned solPdeIndex_mu;
  solPdeIndex_mu = mlPdeSys->GetSolPdeIndex("mu");
  
  unsigned solType_mu = ml_sol->GetSolutionType(solIndex_mu);    // get the finite element type for "mu"
  std::vector < double >  sol_mu;   sol_mu.reserve(maxSize);
  std::vector < int > l2GMap_mu;   l2GMap_mu.reserve(maxSize);

  
  //************** act flag **************************** 
  std::string act_flag_name = "act_flag";
  unsigned int solIndex_act_flag = ml_sol->GetIndex(act_flag_name.c_str());
  unsigned int solFEType_act_flag = ml_sol->GetSolutionType(solIndex_act_flag); 
     if(sol->GetSolutionTimeOrder(solIndex_act_flag) == TIME_DEPENDENT) {
       *(sol->_SolOld[solIndex_act_flag]) = *(sol->_Sol[solIndex_act_flag]);
     }
  
  //********* variables for ineq constraints *****************
  const int ineq_flag = INEQ_FLAG;
  const double c_compl = C_COMPL;
  std::vector < double/*int*/ >  sol_actflag;   sol_actflag.reserve(maxSize); //flag for active set
  std::vector < double >  ctrl_lower;   ctrl_lower.reserve(maxSize);
  std::vector < double >  ctrl_upper;   ctrl_upper.reserve(maxSize);
  //***************************************************  

 
 //***************************************************  
  //********* WHOLE SET OF VARIABLES *****************
  const int solType_max = 2;  //biquadratic

  const int n_unknowns = 4;
 
  std::vector < int > l2GMap_AllVars; // local to global mapping
  l2GMap_AllVars.reserve(n_unknowns*maxSize);
  
  std::vector < double > Res; // local redidual vector
  Res.reserve(n_unknowns*maxSize);

  std::vector < double > Jac;
  Jac.reserve( n_unknowns*maxSize * n_unknowns*maxSize);
 //*************************************************** 

  
 //********************* DATA ************************ 
  double u_des = femus::ctrl::square_or_cube :: cost_functional_without_regularization::DesiredTargetVec()[0];
  double alpha = ALPHA_CTRL_BDRY;
  double beta  = BETA_CTRL_BDRY;
  double penalty_outside_control_boundary = 1.e50;       // penalty for zero control outside Gamma_c
  double penalty_strong_bdry = 1.e20;  // penalty for boundary equation on Gamma_c
  double penalty_ctrl = 1.e10;         //penalty for u=q
 //*************************************************** 
  
  
  if (assembleMatrix)  KK->zero();

    
 // element loop: each process loops only on the elements that owns
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    short unsigned kelGeom = msh->GetElementType(iel);    // element geometry type

 //******************** GEOMETRY **********************
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
    for (int i = 0; i < dim; i++)  x[i].resize(nDofx);
    
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof // via local to global solution node

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->GetTopology()->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

   // elem average point 
    std::vector < double > elem_center(dim);   
    for (unsigned j = 0; j < dim; j++) {  elem_center[j] = 0.;  }
    for (unsigned j = 0; j < dim; j++) {  
      for (unsigned i = 0; i < nDofx; i++) {
         elem_center[j] += x[j][i];
       }
    }
    
   for (unsigned j = 0; j < dim; j++) { elem_center[j] = elem_center[j]/nDofx; }
 //*************************************************** 
  
 //************* set target domain flag **************
   int target_flag = 0;
   target_flag = femus::ctrl::square_or_cube :: cost_functional_without_regularization::ElementTargetFlag(elem_center);
 //*************************************************** 
   
 //************** set control flag *******************
        const bool does_iel_contain_Gamma_c = ctrl::square_or_cube :: pure_boundary< ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >  ::volume_elem_contains_a_Gamma_control_face( sol, msh, iel );
  int control_el_flag = does_iel_contain_Gamma_c? 1 : 0;
        
  std::vector<int> control_node_flag(nDofx,0);
//   if (control_el_flag == 0) std::fill(control_node_flag.begin(), control_node_flag.end(), 0);
 //*************************************************** 
    
 //********************* state ***********************
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u);    // number of solution element dofs
    sol_u   .resize(nDof_u);
    l2GMap_u.resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
      unsigned solDof_u = msh->GetSolutionDof(i,iel, solType_u);    // global to global mapping between solution node and solution dof
      sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);      // global extraction and local storage for the solution
      l2GMap_u[i] = pdeSys->GetSystemDof(solIndex_u, solPdeIndex_u, i, iel);    // global to global mapping between solution node and pdeSys dof
    }
 //********************* state ***********************

 //******************* bdry cont *********************
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);    // number of solution element dofs
    sol_ctrl    .resize(nDof_ctrl);
    l2GMap_ctrl.resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl);    // global to global mapping between solution node and solution dof
      sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);      // global extraction and local storage for the solution
      l2GMap_ctrl[i] = pdeSys->GetSystemDof(solIndex_ctrl, solPdeIndex_ctrl, i, iel);    // global to global mapping between solution node and pdeSys dof
    } 
 //***************************************************  

 //**************** adjoint **************************
    unsigned nDof_adj  = msh->GetElementDofNumber(iel, solType_adj);    
    sol_adj    .resize(nDof_adj);
    l2GMap_adj.resize(nDof_adj);
    for (unsigned i = 0; i < sol_adj.size(); i++) {
      unsigned solDof_adj = msh->GetSolutionDof(i, iel, solType_adj);
      sol_adj[i] = (*sol->_Sol[solIndex_adj])(solDof_adj);
      l2GMap_adj[i] = pdeSys->GetSystemDof(solIndex_adj, solPdeIndex_adj, i, iel);
    } 
 //*************** adjoint **************************** 

 //*************** mu ********************************* 
    unsigned nDof_mu  = msh->GetElementDofNumber(iel, solType_mu);    // number of solution element dofs
    sol_mu   .resize(nDof_mu);
    l2GMap_mu.resize(nDof_mu);
    for (unsigned i = 0; i < sol_mu.size(); i++) {
      unsigned solDof_mu = msh->GetSolutionDof(i, iel, solType_mu);   // global to global mapping between solution node and solution dof
      sol_mu[i] = (*sol->_Sol[solIndex_mu])(solDof_mu);      // global extraction and local storage for the solution 
      l2GMap_mu[i] = pdeSys->GetSystemDof(solIndex_mu, solPdeIndex_mu, i, iel);   // global to global mapping between solution node and pdeSys dof
    }
 
 //********************* ALL VARS ********************* 
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
 //***************************************************

    
 //===================================================   

	// Perform face loop over elements that contain some control face
	if (control_el_flag == 1) {
	  
	  double tau=0.;
	  std::vector <double> normal(dim,0);
	       
	  // loop on faces of the current element

	  for(unsigned jface=0; jface < msh->GetElementFaceNumber(iel); jface++) {
            std::vector < double > xyz_bdc(3,0.);  //not being used, because the boundaries are identified by the face numbers
	    // look for boundary faces
	    if(el->GetFaceElementIndex(iel,jface) < 0) {
	      unsigned int face = -( msh->GetMeshElements()->GetFaceElementIndex(iel,jface)+1);
	      
		
// 	      if( !ml_sol->_SetBoundaryConditionFunction(xx,"U",tau,face,0.) && tau!=0.){
	      if(  face == 3) { //control face

 //===================================================   
		//we use the dirichlet flag to say: if dirichlet = true, we set 1 on the diagonal. if dirichlet = false, we put the boundary equation
	      bool  dir_bool = ml_sol->GetBdcFunction()(xyz_bdc,ctrl_name.c_str(),tau,face,0.);

 //===================================================   

		
		unsigned nve_bdry = msh->GetElementFaceDofNumber(iel,jface,solType_ctrl);
            for (unsigned idim = 0; idim < dim; idim++) {  x_bdry[idim].resize(nve_bdry); }
		const unsigned felt_bdry = msh->GetElementFaceType(iel, jface);    
		for(unsigned i=0; i < nve_bdry; i++) {
		  unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i);
                  unsigned iDof = msh->GetSolutionDof(i_vol, iel, xType);
		  for(unsigned idim=0; idim<dim; idim++) {
		      x_bdry[idim][i]=(*msh->GetTopology()->_Sol[idim])(iDof);
		  }
		}

//========================================================================================================== 

         //************** update active set flag for current nonlinear iteration **************************** 
         // 0: inactive; 1: active_a; 2: active_b
            assert(nDof_mu == nDof_ctrl);
            sol_actflag.resize(nve_bdry/*nDof_mu*/);
            ctrl_lower.resize(nve_bdry/*nDof_mu*/);
            ctrl_upper.resize(nve_bdry/*nDof_mu*/);
           std::fill(sol_actflag.begin(), sol_actflag.end(), 0);
           std::fill(ctrl_lower.begin(), ctrl_lower.end(), 0.);
           std::fill(ctrl_upper.begin(), ctrl_upper.end(), 0.);
   
		      for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
//             for (unsigned i = 0; i < sol_actflag.size(); i++) {
        std::vector<double> node_coords_i(dim,0.);
        for (unsigned d = 0; d < dim; d++) node_coords_i[d] = x_bdry[d][i_bdry];
        ctrl_lower[i_bdry] = ctrl::square_or_cube::mixed_state_or_ctrl_inequality::InequalityConstraint(n_components_ctrl, node_coords_i, false)[0];
        ctrl_upper[i_bdry] = ctrl::square_or_cube::mixed_state_or_ctrl_inequality::InequalityConstraint(n_components_ctrl, node_coords_i, true)[0];

        if      ( (sol_mu[i_vol] + c_compl * (sol_ctrl[i_vol] - ctrl_lower[i_bdry] )) < 0 )  sol_actflag[i_bdry] = 1;
        else if ( (sol_mu[i_vol] + c_compl * (sol_ctrl[i_vol] - ctrl_upper[i_bdry] )) > 0 )  sol_actflag[i_bdry] = 2;
            }
            
        //************** act flag **************************** 
      for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
	    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
      unsigned solDof_actflag = msh->GetSolutionDof(i_vol, iel, solFEType_act_flag); 
      (sol->_Sol[solIndex_act_flag])->set(solDof_actflag,sol_actflag[i_bdry]);     
    }
    

 // ===================================================
 //node-based insertion on the boundary ===============
 // ===================================================
    
 //============= delta_mu row ===============================
      std::vector<double> Res_mu (nDof_mu); std::fill(Res_mu.begin(),Res_mu.end(), 0.);
      
      for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
	    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
//     for (unsigned i = 0; i < sol_actflag.size(); i++) {
      if (sol_actflag[i_bdry] == 0){  //inactive
         Res_mu [i_vol] = - ineq_flag * ( 1. * sol_mu[i_vol] - 0. ); 
// 	 Res_mu [i] = Res[nDof_u + nDof_ctrl + nDof_adj + i]; 
      }
      else if (sol_actflag[i_bdry] == 1){  //active_a 
	 Res_mu [i_vol] = - ineq_flag * ( c_compl *  sol_ctrl[i_vol] - c_compl * ctrl_lower[i_bdry]);
      }
      else if (sol_actflag[i_bdry] == 2){  //active_b 
	Res_mu [i_vol]  =  - ineq_flag * ( c_compl *  sol_ctrl[i_vol] - c_compl * ctrl_upper[i_bdry]);
      }
    }

    
    RES->insert(Res_mu, l2GMap_mu);    
 //============= delta_mu row - end ===============================
    
 //============= delta_mu-delta_ctrl row ===============================
 //auxiliary volume vector for act flag
 unsigned nDof_actflag_vol  = msh->GetElementDofNumber(iel, solFEType_act_flag);
 std::vector<double> sol_actflag_vol(nDof_actflag_vol); 


 for (unsigned i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++) if (sol_actflag[i_bdry] != 0 ) sol_actflag[i_bdry] = ineq_flag * c_compl;    
 
 std::fill(sol_actflag_vol.begin(), sol_actflag_vol.end(), 0.);
    for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
       unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
       sol_actflag_vol[i_vol] = sol_actflag[i_bdry];
    }
 
 KK->matrix_set_off_diagonal_values_blocked(l2GMap_mu, l2GMap_ctrl, sol_actflag_vol);
 //============= delta_mu-delta_ctrl row - end ===============================

 //============= delta_mu-delta_mu row ===============================
  for (unsigned i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++) sol_actflag[i_bdry] =  ineq_flag * (1 - sol_actflag[i_bdry]/c_compl)  + (1-ineq_flag) * 1.;  //can do better to avoid division, maybe use modulo operator 

 std::fill(sol_actflag_vol.begin(), sol_actflag_vol.end(), 0.);
    for (int i_bdry = 0; i_bdry < sol_actflag.size(); i_bdry++)  {
       unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
       sol_actflag_vol[i_vol] = sol_actflag[i_bdry];
    }
  
  KK->matrix_set_off_diagonal_values_blocked(l2GMap_mu, l2GMap_mu, sol_actflag_vol );
 //============= delta_mu-delta_mu row - end ===============================
    
 // =========================================================
 //node-based insertion on the boundary - end ===============
 // =========================================================
		
//========= initialize gauss quantities on the boundary =============================================
                double sol_ctrl_bdry_gss = 0.;
                double sol_adj_bdry_gss = 0.;
                std::vector<double> sol_ctrl_x_bdry_gss(dim);   std::fill(sol_ctrl_x_bdry_gss.begin(), sol_ctrl_x_bdry_gss.end(), 0.);

//========= initialize gauss quantities on the boundary ============================================
		
		for(unsigned ig_bdry=0; ig_bdry < msh->_finiteElement[felt_bdry][solType_ctrl]->GetGaussPointNumber(); ig_bdry++) {
		  
		  msh->_finiteElement[felt_bdry][solType_ctrl]->JacobianSur(x_bdry,ig_bdry,weight_bdry,phi_ctrl_bdry,phi_ctrl_x_bdry,normal);
		  msh->_finiteElement[felt_bdry][solType_adj]->JacobianSur(x_bdry,ig_bdry,weight_bdry,phi_adj_bdry,phi_adj_x_bdry,normal);

//========== temporary soln for surface gradient on a face parallel to the X axis ===================
		  double dx_dxi = 0.;
		 const elem_type_1D * myeltype = static_cast<const elem_type_1D*>(msh->_finiteElement[felt_bdry][solType_ctrl]);
		 const double * myptr = myeltype->GetDPhiDXi(ig_bdry);
		      for (int inode = 0; inode < nve_bdry/*_nc*/; inode++) dx_dxi += myptr[inode] * x_bdry[0][inode];
  
		      for (int inode = 0; inode < nve_bdry/*_nc*/; inode++) {
                            for (int d = 0; d < dim; d++) {
                              if (d==0 ) phi_ctrl_x_bdry[inode + d*nve_bdry/*_nc*/] = myptr[inode]* (1./ dx_dxi);
                              else  phi_ctrl_x_bdry[inode + d*nve_bdry/*_nc*/] = 0.;
                         }
                     }
//========== temporary soln for surface gradient on a face parallel to the X axis ====================
		  
//========== compute gauss quantities on the boundary ================================================
		  sol_ctrl_bdry_gss = 0.;
		  sol_adj_bdry_gss = 0.;
                  std::fill(sol_ctrl_x_bdry_gss.begin(), sol_ctrl_x_bdry_gss.end(), 0.);
		      for (int i_bdry = 0; i_bdry < nve_bdry/*_nc*/; i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
			
			sol_adj_bdry_gss  +=  sol_adj[i_vol] * phi_adj_bdry[i_bdry];
			sol_ctrl_bdry_gss +=  sol_ctrl[i_vol] * phi_ctrl_bdry[i_bdry];
                            for (int d = 0; d < dim; d++) {
			      sol_ctrl_x_bdry_gss[d] += sol_ctrl[i_vol] * phi_ctrl_x_bdry[i_bdry + d*nve_bdry];
			    }
		      }  

//========== compute gauss quantities on the boundary =================================================

		  // *** phi_i loop ***
		  for(unsigned i_bdry=0; i_bdry < nve_bdry; i_bdry++) {
		    
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);

                 double lap_rhs_dctrl_ctrl_bdry_gss_i = 0.;
                 for (unsigned d = 0; d < dim; d++) {
                       if ( i_vol < nDof_ctrl )  lap_rhs_dctrl_ctrl_bdry_gss_i +=  phi_ctrl_x_bdry[i_bdry + d*nve_bdry] * sol_ctrl_x_bdry_gss[d];
                 }
                 
//=============== construct control node flag field on the go  =========================================  
	      /* (control_node_flag[i])       picks nodes on \Gamma_c
	         (1 - control_node_flag[i])   picks nodes on \Omega \setminus \Gamma_c
	       */
	      if (dir_bool == false) { 
		std::cout << " found boundary control nodes ==== " << std::endl;
			for(unsigned k=0; k<control_node_flag.size(); k++) {
				  control_node_flag[i_vol] = 1;
			}
              }
//=============== construct control node flag field on the go  =========================================

//============ Bdry Residuals ==================	

        if (i_vol < nDof_u)     Res[ (0 + i_vol) ]                    += 0.;

        if (i_vol < nDof_ctrl)  Res[ (nDof_u + i_vol) ]               +=  - control_node_flag[i_vol] *  weight_bdry *
                                                              (     alpha * phi_ctrl_bdry[i_bdry] * sol_ctrl_bdry_gss
							         +  beta * lap_rhs_dctrl_ctrl_bdry_gss_i 
							         +  phi_ctrl_bdry[i_bdry]*sol_adj_bdry_gss
							        );  //boundary optimality condition

        if (i_vol < nDof_adj)   Res[ (nDof_u + nDof_ctrl +  i_vol) ]  += - control_node_flag[i_vol] *  weight_bdry * phi_adj_bdry[i_bdry]*sol_ctrl_bdry_gss; 
//============ Bdry Residuals ==================    
		    
		    for(unsigned j_bdry=0; j_bdry < nve_bdry; j_bdry ++) {
		    unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, jface, j_bdry);

		    
//============ boundary control eqns ===========
		    
// SECOND BLOCK ROW
//========block delta_control / control ========
   	           if ( i_vol < nDof_ctrl && j_vol < nDof_ctrl ) {
                     Jac[  
		        (nDof_u + i_vol) * nDof_AllVars +
	                (nDof_u + j_vol)             ]  += control_node_flag[i_vol] *  weight_bdry * (alpha * phi_ctrl_bdry[i_bdry] * phi_ctrl_bdry[j_bdry]);

                     double  lap_mat_dctrl_ctrl_bdry_gss = 0.;
		     for (unsigned d = 0; d < dim; d++) {  lap_mat_dctrl_ctrl_bdry_gss += phi_ctrl_x_bdry[i_bdry + d*nve_bdry] * phi_ctrl_x_bdry[j_bdry + d*nve_bdry];   }
		     
	             Jac[
		        (nDof_u + i_vol) * nDof_AllVars +
	                (nDof_u + j_vol)             ]  += control_node_flag[i_vol] * weight_bdry * beta *  lap_mat_dctrl_ctrl_bdry_gss;
                }
//========== block delta_control/adjoint ========
		   if ( i_vol < nDof_ctrl    && j_vol < nDof_adj )   
		     Jac[ 
			(nDof_u + i_vol) * nDof_AllVars +
		        (nDof_u + nDof_ctrl + j_vol) ]  += control_node_flag[i_vol] * (weight_bdry * phi_adj_bdry[j_bdry] * phi_ctrl_bdry[i_bdry]);      
		      

// THIRD BLOCK ROW
//========= block delta_adjoint/control =========
		   if ( i_vol < nDof_adj    && j_vol < nDof_ctrl )   
		     Jac[ 
			(nDof_u + nDof_ctrl + i_vol) * nDof_AllVars  +
		        (nDof_u + j_vol)             ]  += control_node_flag[i_vol] * (weight_bdry * phi_ctrl_bdry[j_bdry] * phi_adj_bdry[i_bdry]);      
			
//============ boundary control eqn =============			 
		   				
		    }  //end j loop
		  } //end i loop
		}
	      }
	    }
	  }    
	  
	} //end if control element flag
	
	else { //here we set the diagonal to 1 and the rhs to 0
	  
	  
	}
    
//========= gauss value quantities ==============   
	double sol_u_gss = 0.;
	double sol_adj_gss = 0.;
	std::vector<double> sol_u_x_gss(dim);     std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	std::vector<double> sol_adj_x_gss(dim);   std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
//===============================================   
 
 
      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[kelGeom][solType_max]->GetGaussPointNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
        //  ==== state 
	msh->_finiteElement[kelGeom][solType_u]  ->Jacobian(x, ig, weight, phi_u, phi_u_x, phi_u_xx);
        //  ==== adj 
        msh->_finiteElement[kelGeom][solType_adj]->Jacobian(x, ig, weight, phi_adj, phi_adj_x, phi_adj_xx);
          
	sol_u_gss = 0.;
	sol_adj_gss = 0.;
	std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
	
	for (unsigned i = 0; i < nDof_u; i++) {
	                                                sol_u_gss      += sol_u[i] * phi_u[i];
                   for (unsigned d = 0; d < dim; d++)   sol_u_x_gss[d] += sol_u[i] * phi_u_x[i * dim + d];
          }
	
	for (unsigned i = 0; i < nDof_adj; i++) {
	                                                sol_adj_gss      += sol_adj[i] * phi_adj[i];
                   for (unsigned d = 0; d < dim; d++)   sol_adj_x_gss[d] += sol_adj[i] * phi_adj_x[i * dim + d];
        }

//==========FILLING WITH THE EQUATIONS ===========
	// *** phi_i loop ***
        for (unsigned i = 0; i < nDof_max; i++) {
	  
              double laplace_rhs_du_adj_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_u )         laplace_rhs_du_adj_i      +=  (phi_u_x   [i * dim + kdim] * sol_adj_x_gss[kdim]);
	      }
	      
              double laplace_rhs_dadj_u_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_adj )         laplace_rhs_dadj_u_i      +=  (phi_adj_x   [i * dim + kdim] * sol_u_x_gss[kdim]);
	      }
	      
//============ Volume residuals ==================	    
        if (i < nDof_u)     Res[0                  + i] += - weight * ( target_flag * phi_u[i] * ( sol_u_gss - u_des)
	                                                                  - laplace_rhs_du_adj_i); 
        if (i < nDof_ctrl)  Res[nDof_u             + i] += - penalty_outside_control_boundary * ( (1 - control_node_flag[i]) * (  sol_ctrl[i] - 0.)  );
	      
        if (i < nDof_adj)   Res[nDof_u + nDof_ctrl + i] += - weight * (-1) * (laplace_rhs_dadj_u_i);
            
        if (i < nDof_mu)    Res[nDof_u + nDof_ctrl + nDof_adj + i] += - penalty_outside_control_boundary * ( (1 - control_node_flag[i]) * (  sol_mu[i] - 0.)  );
//============  Volume Residuals ==================	    
	      
	      
          if (assembleMatrix) {
	    
            // *** phi_j loop ***
            for (unsigned j = 0; j < nDof_max; j++) {
              double laplace_mat_dadj_u = 0.;
              double laplace_mat_du_adj = 0.;

              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_adj && j < nDof_u )     laplace_mat_dadj_u        +=  (phi_adj_x [i * dim + kdim] * phi_u_x   [j * dim + kdim]);
              if ( i < nDof_u   && j < nDof_adj )   laplace_mat_du_adj        +=  (phi_u_x   [i * dim + kdim] * phi_adj_x [j * dim + kdim]);
		
	      }

              //============ delta_state row ============================
              // BLOCK delta_state / state	      
              if ( i < nDof_u && j < nDof_u )   
		Jac[ i * nDof_AllVars +
		   (0 + j)                                ]  += weight * target_flag *  phi_u[i] * phi_u[j];   
		   
	      //BLOCK delta_state / adjoint
              if ( i < nDof_u && j < nDof_adj )   
		Jac[ i * nDof_AllVars +
		   (nDof_u + nDof_ctrl + j)               ]  += weight * (-1) * laplace_mat_du_adj;
	      
	      
              //=========== delta_control row ===========================
              //enforce control zero outside the control boundary
	      if ( i < nDof_ctrl && j < nDof_ctrl && i==j)
		Jac[    
		   (nDof_u + i) * nDof_AllVars  +
		   (nDof_u + j)                           ]  += penalty_outside_control_boundary * ( (1 - control_node_flag[i]));

		
	      //=========== delta_adjoint row ===========================
	      // BLOCK delta_adjoint / state
	      if ( i < nDof_adj && j < nDof_u ) {
		Jac[    
		   (nDof_u + nDof_ctrl + i) * nDof_AllVars  +
		   (0 + j)                                ]  += weight * (-1) * laplace_mat_dadj_u;
           
           //============= delta_mu row ===============================
	        if ( i < nDof_mu && j < nDof_mu && i==j )   
		  Jac[ (nDof_u + nDof_ctrl + nDof_adj + i) * nDof_AllVars +
		       (nDof_u + nDof_ctrl + nDof_adj + j)]  += penalty_outside_control_boundary * ( (1 - control_node_flag[i]));  

	      }
	    
	    } // end phi_j loop
	    
          } // endif assemble_matrix

        } // end phi_i loop
        
      } // end gauss point loop

/*	if (control_el_flag == 1) {
	  
         for (unsigned i = 0; i < nDof_max; i++) {
            for (unsigned j = 0; j < nDof_max; j++) {
	      std::cout << " " << std::setfill(' ') << std::setw(10) << Jac[ (nDof_u + i) * nDof_AllVars + (nDof_u + j) ];
// 	      std::cout << " " << std::setfill(' ') << std::setw(10) << Jac[ (nDof_u + i) * nDof_AllVars + (nDof_u + nDof_ctrl + j) ];
	    }
	      std::cout << std::endl;
	 }

    std::cout << " ========== " << iel << " ================== " << std::endl;      
         for (unsigned i = 0; i < nDof_max; i++) {
            for (unsigned j = 0; j < nDof_max; j++) {
// 	      std::cout << " " << std::setfill(' ') << std::setw(10) << Jac[ (nDof_u + i) * nDof_AllVars + (nDof_u + j) ];
	      std::cout << " " << std::setfill(' ') << std::setw(10) << Jac[ (nDof_u + i) * nDof_AllVars + (nDof_u + nDof_ctrl + j) ];
	    }
	      std::cout << std::endl;
	 }
	 
	 
	}
    std::cout << " ========== " << iel << " ================== " << std::endl;   */   
    
    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    RES->add_vector_blocked(Res, l2GMap_AllVars);

    if (assembleMatrix) {
      //store K in the global matrix KK
      KK->add_matrix_blocked(Jac, l2GMap_AllVars, l2GMap_AllVars);
    }
    
    
    //========== dof-based part, without summation
 


    
 //============= delta_ctrl-delta_mu row ===============================
 KK->matrix_set_off_diagonal_values_blocked(l2GMap_ctrl, l2GMap_mu, ineq_flag * 1.);//------------------------------->>>>>>
  
  
  } //end element loop for each process
  

  
  unsigned int ctrl_index = mlPdeSys->GetSolPdeIndex("control");
  unsigned int mu_index = mlPdeSys->GetSolPdeIndex("mu");

  unsigned int global_ctrl_size = pdeSys->KKoffset[ctrl_index+1][iproc] - pdeSys->KKoffset[ctrl_index][iproc];
  
  std::vector<double>  one_times_mu(global_ctrl_size, 0.);
  std::vector<int>    positions(global_ctrl_size);

  for (unsigned i = 0; i < positions.size(); i++) {
    positions[i] = pdeSys->KKoffset[ctrl_index][iproc] + i;
    one_times_mu[i] = ineq_flag * 1. * (*sol->_Sol[solIndex_mu])(i/*position_mu_i*/) ;
  }
    RES->add_vector_blocked(one_times_mu, positions);
    
  // ***************** END ASSEMBLY *******************

    if (assembleMatrix) KK->close();
    std::ostringstream mat_out; mat_out << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "matrix_" << mlPdeSys->GetNonlinearIt()  << ".txt";
    KK->print_matlab(mat_out.str(),"ascii"); //  KK->print();

    RES->close();
    std::ostringstream res_out; res_out << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "res_" << mlPdeSys->GetNonlinearIt()  << ".txt";
    std::filebuf res_fb;
    res_fb.open (res_out.str().c_str(),std::ios::out);
    std::ostream  res_file_stream(&res_fb);
    RES->print(res_file_stream);

  return;
}

};


class lifting_internal  {

#define SERVICE 1.


  public:
    
static void assemble_elliptic_neumann_control(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data

  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("LiftRestr");   // pointer to the linear implicit system named "LiftRestr"
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->GetMeshElements();  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

 //***************************************************  
  std::vector < std::vector < double > > x(dim);    // local coordinates
  std::vector < std::vector < double> >  x_bdry(dim);
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
  for (unsigned i = 0; i < dim; i++) {
    x[i].reserve(maxSize);
    x_bdry[i].reserve(maxSize);
  }
 //***************************************************   

 //***************************************************  
  double weight; // gauss point weight
  double weight_bdry; // gauss point weight on the boundary

 //********************* state *********************** 
 //***************************************************  
  std::vector <double> phi_u;  // local test function
  std::vector <double> phi_u_x; // local test function first order partial derivatives
  std::vector <double> phi_u_xx; // local test function second order partial derivatives

  phi_u.reserve(maxSize);
  phi_u_x.reserve(maxSize * dim);
  phi_u_xx.reserve(maxSize * dim2);
  
 
  unsigned solIndex_u;
  solIndex_u = mlSol->GetIndex("state");    // get the position of "state" in the ml_sol object
  unsigned solType_u = mlSol->GetSolutionType(solIndex_u);    // get the finite element type for "state"

  unsigned solPdeIndexThom;
  solPdeIndexThom = mlPdeSys->GetSolPdeIndex("state");    // get the position of "state" in the pdeSys object

  std::vector < double >  sol_u; // local solution
  sol_u.reserve(maxSize);
  std::vector < int > l2GMap_u;
  l2GMap_u.reserve(maxSize);
 //***************************************************  
 //***************************************************  

  
 //******************** control ********************** 
 //***************************************************   
  std::vector <double> phi_ctrl;  // local test function
  std::vector <double> phi_ctrl_x; // local test function first order partial derivatives
  std::vector <double> phi_ctrl_xx; // local test function second order partial derivatives

  phi_ctrl.reserve(maxSize);
  phi_ctrl_x.reserve(maxSize * dim);
  phi_ctrl_xx.reserve(maxSize * dim2);
  
  unsigned solIndex_ctrl;
  solIndex_ctrl = mlSol->GetIndex("control");
  unsigned solType_ctrl = mlSol->GetSolutionType(solIndex_ctrl);

  unsigned solPdeIndex_ctrl;
  solPdeIndex_ctrl = mlPdeSys->GetSolPdeIndex("control");

  std::vector < double >  sol_ctrl; // local solution
  sol_ctrl.reserve(maxSize);
  std::vector < int > l2GMap_ctrl;
  l2GMap_ctrl.reserve(maxSize);
 
  std::string ctrl_name = "control";
  
 //*************************************************** 

  std::vector <double> sol_ctrl_x_vol_at_bdry_gss;
  sol_ctrl_x_vol_at_bdry_gss.reserve(dim);
 //*************************************************** 
   
  std::vector <double> phi_ctrl_vol_at_bdry;  // local test function
  phi_ctrl_vol_at_bdry.reserve(maxSize);
  
  std::vector <double> phi_ctrl_x_vol_at_bdry; // local test function first order partial derivatives
  phi_ctrl_x_vol_at_bdry.reserve(maxSize * dim);

 
 //***************************************************  
 //***************************************************  
  
  
 //********************* adjoint ********************* 
 //***************************************************  
  std::vector <double> phi_adj;  // local test function
  std::vector <double> phi_adj_x; // local test function first order partial derivatives
  std::vector <double> phi_adj_xx; // local test function second order partial derivatives

  phi_adj.reserve(maxSize);
  phi_adj_x.reserve(maxSize * dim);
  phi_adj_xx.reserve(maxSize * dim2);
  
  unsigned solIndex_adj;
  solIndex_adj = mlSol->GetIndex("adjoint");    // get the position of "state" in the ml_sol object
  unsigned solType_adj = mlSol->GetSolutionType(solIndex_adj);    // get the finite element type for "state"

  unsigned solPdeIndex_adj;
  solPdeIndex_adj = mlPdeSys->GetSolPdeIndex("adjoint");    // get the position of "state" in the pdeSys object

   //boundary adjoint shape functions  
  std::vector <double> phi_adj_bdry;  
  std::vector <double> phi_adj_x_bdry; 

  phi_adj_bdry.reserve(maxSize);
  phi_adj_x_bdry.reserve(maxSize * dim);
  
  //***************************************************  
  std::vector < double >  sol_adj; // local solution
    sol_adj.reserve(maxSize);
  std::vector < int > l2GMap_adj;
    l2GMap_adj.reserve(maxSize);
    
    std::vector < double > sol_adj_bdry;
     sol_adj.reserve(maxSize);
 //***************************************************  
 //***************************************************  

    
    

  
 //***************************************************
 //********* WHOLE SET OF VARIABLES ****************** 
  const int solType_max = 2;  //biquadratic

  const int n_vars = 3;
 
  std::vector < int > l2GMap_AllVars; // local to global mapping
  l2GMap_AllVars.reserve(n_vars*maxSize);
  
  std::vector < double > Res; // local redidual vector
  Res.reserve(n_vars*maxSize);

  std::vector < double > Jac;
  Jac.reserve( n_vars*maxSize * n_vars*maxSize);
 //***************************************************  

  
 //********************* DATA ************************ 
  double u_des = femus::ctrl::square_or_cube :: cost_functional_without_regularization::DesiredTargetVec()[0];
  double alpha = ALPHA_CTRL_VOL;
  double beta  = BETA_CTRL_VOL;
  double penalty_strong = 10e+14;
  //double penalty_ctrl = 1.e10;         //penalty for u=q
 //***************************************************  
  
  
  if (assembleMatrix)  KK->zero();

    
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    short unsigned kelGeom = msh->GetElementType(iel);    // element geometry type

 //******************** GEOMETRY ********************* 
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
    for (int i = 0; i < dim; i++)  x[i].resize(nDofx);
    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);  // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->GetTopology()->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

   // elem average point 
    std::vector < double > elem_center(dim);   
    for (unsigned j = 0; j < dim; j++) {  elem_center[j] = 0.;  }
    for (unsigned j = 0; j < dim; j++) {  
      for (unsigned i = 0; i < nDofx; i++) {
         elem_center[j] += x[j][i];
       }
    }
    
   for (unsigned j = 0; j < dim; j++) { elem_center[j] = elem_center[j]/nDofx; }
 //***************************************************  
  
 //****** set target domain flag ********************* 
   int target_flag = 0;
   target_flag = femus::ctrl::square_or_cube :: cost_functional_without_regularization::ElementTargetFlag(elem_center);
 //*************************************************** 
   
    
 //**************** state **************************** 
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u);    // number of solution element dofs
    sol_u    .resize(nDof_u);
    l2GMap_u.resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
     unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u);  // global to global mapping between solution node and solution dof
      sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);            // global extraction and local storage for the solution
      l2GMap_u[i] = pdeSys->GetSystemDof(solIndex_u, solPdeIndexThom, i, iel);  // global to global mapping between solution node and pdeSys dof
    }
 //***************************************************  
 
 //***************** control ************************* 
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);    // number of solution element dofs
    sol_ctrl    .resize(nDof_ctrl);
    l2GMap_ctrl.resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl);    // global to global mapping between solution node and solution dof
      sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);      // global extraction and local storage for the solution
      l2GMap_ctrl[i] = pdeSys->GetSystemDof(solIndex_ctrl, solPdeIndex_ctrl, i, iel);    // global to global mapping between solution node and pdeSys dof
    } 
 //*************************************************** 
 

 //************** adjoint **************************** 
    unsigned nDof_adj  = msh->GetElementDofNumber(iel, solType_adj);    // number of solution element dofs
        sol_adj    .resize(nDof_adj);
        l2GMap_adj.resize(nDof_adj);
    for (unsigned i = 0; i < sol_adj.size(); i++) {
      unsigned solDof_adj = msh->GetSolutionDof(i, iel, solType_adj);   // global to global mapping between solution node and solution dof
      sol_adj[i] = (*sol->_Sol[solIndex_adj])(solDof_adj);      // global extraction and local storage for the solution
      l2GMap_adj[i] = pdeSys->GetSystemDof(solIndex_adj, solPdeIndex_adj, i, iel);   // global to global mapping between solution node and pdeSys dof
    } 
 //***************************************************  

 //******************** ALL VARS ********************* 
    unsigned nDof_AllVars = nDof_u + nDof_ctrl + nDof_adj; 
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
 //*************************************************** 
   
    
 //***** set control flag ****************************
  int control_el_flag = 0;
  control_el_flag = femus::ctrl:: square_or_cube:: lifting_internal< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >::ControlDomainFlag_internal_restriction(elem_center);
  std::vector<int> control_node_flag(nDof_ctrl,0);
  if (control_el_flag == 1) std::fill(control_node_flag.begin(), control_node_flag.end(), 1);
 //*************************************************** 


  
//////////////////Begin Boundary Integral////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if (control_el_flag == 1) {
	  
	  double tau=0.;
	  std::vector <double> normal(dim,0);
// 	  double normal_fixed[3] = {0.,1.,0};
	       
	  // loop on faces of the current element

	  for(unsigned jface=0; jface < msh->GetElementFaceNumber(iel); jface++) {
            std::vector < double > xyz_bdc(3,0.);  //not being used, because the boundaries are identified by the face numbers
	    // look for boundary faces
	    if(el->GetFaceElementIndex(iel,jface) < 0) {
	      unsigned int face = -( msh->GetMeshElements()->GetFaceElementIndex(iel,jface)+1);
	      
		
// 	      if( !ml_sol->_SetBoundaryConditionFunction(xx,"U",tau,face,0.) && tau!=0.){
	      if(  face == 3) { //control face

 //=================================================== 
		//we use the dirichlet flag to say: if dirichlet = true, we set 1 on the diagonal. if dirichlet = false, we put the boundary equation
	      bool  dir_bool = mlSol->GetBdcFunction()(xyz_bdc,ctrl_name.c_str(),tau,face,0.);

 //=================================================== 

		
		unsigned nve_bdry = msh->GetElementFaceDofNumber(iel,jface,solType_ctrl);
		const unsigned felt_bdry = msh->GetElementFaceType(iel, jface);

 //*****************x dofs on the boundary ********************************** 
		for(unsigned i=0; i < nve_bdry; i++) {
		  unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i);
                  unsigned iDof = msh->GetSolutionDof(i_vol, iel, xType);
		  for(unsigned idim=0; idim<dim; idim++) {
		      x_bdry[idim][i] = (*msh->GetTopology()->_Sol[idim])(iDof);
		  }
		}
 //*****************x dofs on the boundary, end********************************** 
		
 //***************** adj dofs on the boundary ********************************** 
		sol_adj_bdry    .resize(nve_bdry/*nDof_adj_bdry*/);
                 for (unsigned i = 0; i < sol_adj_bdry.size(); i++) {
		  unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i);
                     unsigned solDof_adj_vol = msh->GetSolutionDof(i_vol, iel, solType_adj);   // global to global mapping between solution node and solution dof
                          sol_adj_bdry[i] = (*sol->_Sol[solIndex_adj])(solDof_adj_vol);      // global extraction and local storage for the solution
		 }
 //***************** adj dofs on the boundary, end ********************************** 
      
		for(unsigned ig_bdry=0; ig_bdry < msh->_finiteElement[felt_bdry][solType_ctrl]->GetGaussPointNumber(); ig_bdry++) {
		  
		  msh->_finiteElement[felt_bdry][solType_adj]->JacobianSur(x_bdry,ig_bdry,weight_bdry,phi_adj_bdry,phi_adj_x_bdry,normal);
		  msh->_finiteElement[kelGeom][solType_ctrl]->fill_volume_shape_funcs_at_boundary_quadrature_points_on_current_elem(x,x_bdry,jface,ig_bdry,phi_ctrl_vol_at_bdry,phi_ctrl_x_vol_at_bdry);
		  
		      
//=============== grad dot n for residual ========================================= 
//     compute gauss quantities on the boundary through VOLUME interpolation
//            std::fill(sol_ctrl_x_vol_at_bdry_gss.begin(), sol_ctrl_x_vol_at_bdry_gss.end(), 0.);
		  for (int d = 0; d < dim; d++) {sol_ctrl_x_vol_at_bdry_gss[d]=0.;}
		      for (int iv = 0; iv < nDof_ctrl; iv++)  {
			
                            for (int d = 0; d < dim; d++) {
//    std::cout << " ivol " << iv << std::endl;
//   std::cout << " ctrl dofs " << sol_ctrl[iv] << std::endl;
			      sol_ctrl_x_vol_at_bdry_gss[d] += sol_ctrl[iv] * phi_ctrl_x_vol_at_bdry[iv * dim + d];//notice that the convention of the orders x y z is different from vol to bdry
			    }
		      }  
		      
    double grad_ctrl_dot_n_res = 0.;
        for(unsigned d=0; d<dim; d++) {
	  grad_ctrl_dot_n_res += sol_ctrl_x_vol_at_bdry_gss[d]*normal[d];  
	}
		
//=============== grad dot n  for residual =========================================       

//========== compute gauss quantities on the boundary ================================================
	double sol_adj_bdry_gss=0.;
            for (int ib = 0; ib < nve_bdry/*nDof_adj_bdry*/; ib++){ sol_adj_bdry_gss  +=     sol_adj_bdry[ib]*phi_adj_bdry[ib]    ;}

		  // *** phi_i loop ***
		  for(unsigned i_bdry=0; i_bdry < nve_bdry; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);

 //=============== grad phi dot n  for residual =========================================       
      double grad_phi_ctrl_dot_n_res = 0.;
        for(unsigned d=0; d<dim; d++) {
	  grad_phi_ctrl_dot_n_res += phi_ctrl_x_vol_at_bdry[i_vol * dim + d]*normal[d];  
	}
 //=============== grad phi dot n  for residual =========================================       

//=============== construct control node flag field on the go  =========================================    
	      /* (control_node_flag[i])       picks nodes on \Gamma_c
	         (1 - control_node_flag[i])   picks nodes on \Omega \setminus \Gamma_c
	       */
	      if (dir_bool == false) { 
// 		std::cout << " found boundary control nodes ==== " << std::endl;
			for(unsigned k=0; k<control_node_flag.size(); k++) {
				  control_node_flag[i_vol] = 1;
			}
              }
//=============== construct control node flag field on the go  =========================================    

//   std::cout << " graddotn_res " << grad_ctrl_dot_n_res << std::endl;
//    std::cout << " normal " << sol_adj[0] << std::endl;
//    std::cout << " sol_ctrl_x_vol_at_bdry_gss_0 " << sol_ctrl_x_vol_at_bdry_gss[0] << std::endl;
//     std::cout << " sol_ctrl_x_vol_at_bdry_gss_1 " << sol_ctrl_x_vol_at_bdry_gss[1] << std::endl;
// 
//    std::cout << " phi_ctrl_x_vol_at_bdry " << phi_ctrl_x_vol_at_bdry[2] << std::endl;
//    std::cout << " sol_ctrl " << sol_ctrl[0] << std::endl;
  
   std::cout << " sol_adj " << sol_adj_bdry_gss << std::endl;
   std::cout << " grad_phi_ctrl_dot_n_res " << grad_phi_ctrl_dot_n_res << std::endl;
		 
//============ Bdry Residuals ==================	
                if (i_vol < nDof_u)     Res[ (0 + i_vol) ]                    +=  0.; 
		
                if (i_vol < nDof_ctrl)  Res[ (nDof_u + i_vol) ]               +=   - control_node_flag[i_vol] *  weight_bdry *
                                                                                (    grad_phi_ctrl_dot_n_res * sol_adj_bdry_gss
// wrong                                                                                grad_ctrl_dot_n_res * phi_adj_bdry[i_bdry]/*phi_adj_vol_at_bdry[i_vol??]*/
							                         )*SERVICE;    
		
                if (i_vol < nDof_adj)   Res[ (nDof_u + nDof_ctrl +  i_vol) ]  +=  - control_node_flag[i_vol] *  weight_bdry *
                                                                                ( grad_ctrl_dot_n_res * phi_adj_bdry[i_bdry]
							                         ); 
//============ Bdry Residuals ==================    
		    
		    for(unsigned j_bdry=0; j_bdry < nve_bdry; j_bdry ++) {
		         unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, jface, j_bdry);

//============ Bdry Jacobians ==================	
//============ Bdry Jacobians ==================	


// FIRST BLOCK ROW
//============ u = q ===========================	    
// block delta_state/state =====================
// 		if (i_vol < nDof_u && j_vol < nDof_u && i_vol == j_vol)  {
// 		  Jac[    
// 		(0 + i_vol) * nDof_AllVars  +
// 		(0 + j_vol)                                ]  +=  penalty_ctrl * ( control_node_flag[i_vol]);
// 		  
// 		}

// block delta_state/control ===================
// 	      if ( i_vol < nDof_u && j_vol < nDof_ctrl && i_vol == j_vol) {
// 		Jac[    
// 		(0     + i_vol) * nDof_AllVars  +
// 		(nDof_u + j_vol)                           ]  += penalty_ctrl * ( control_node_flag[i_vol]) * (-1.);
// 	
// 	      }
//============ u = q ===========================

		    

// SECOND BLOCK ROW

/*//==========block delta_control/adjoint ========
		   if ( i_vol < nDof_ctrl    && j_vol < nDof_adj)   
		     Jac[ 
			(nDof_u + i_vol) * nDof_AllVars  +
		        (nDof_u + nDof_ctrl + j_vol)             ]  += control_node_flag[i_vol] * (-1) *
		        (
			  weight_bdry * grad_ctrl_dot_n_mat * phi_ctrl_bdry[i_bdry]
// 			  weight_bdry * phi_adj_bdry[j_bdry] * phi_ctrl_bdry[i_bdry]  // for neumann boundary condition
			  
			);    */  
		    

		   
//============ End Bdry Jacobians ==================	
//============ End Bdry Jacobians ==================	
				
	      }  //end j loop
	      
//===================loop over j in the VOLUME (while i is in the boundary)	      
	for(unsigned j=0; j < nDof_max; j ++) {
		      
  //=============== grad dot n  =========================================    
    double grad_ctrl_dot_n_mat = 0.;
        for(unsigned d=0; d<dim; d++) {
	                                       grad_ctrl_dot_n_mat += phi_ctrl_x_vol_at_bdry[j * dim + d]*normal[d];  //notice that the convention of the orders x y z is different from vol to bdry
	}
//=============== grad dot n  =========================================    
// //=============== grad phi dot n  =======================================edit==========
//     double grad_phi_ctrl_dot_n_mat = 0.;
//         for(unsigned d=0; d<dim; d++) {
// 	                                       grad_phi_ctrl_dot_n_mat += phi_ctrl_x_vol_at_bdry[j * dim + d]*normal[d];  //notice that the convention of the orders x y z is different from vol to bdry
// 	}
//=============== grad phi dot n  =================================================
//  std::cout << " gradcontroldotn " << grad_ctrl_dot_n_mat << std::endl;
  

//==========block delta_control\adjoint==================================edit=======
              if ( i_vol < nDof_ctrl   && j < nDof_adj  ) 
		Jac[ 
		   (nDof_u + i_vol) * nDof_AllVars  + 
		   (nDof_u + nDof_ctrl + j)         ]  +=  control_node_flag[i_vol] * 
		         weight_bdry * grad_ctrl_dot_n_mat * phi_adj_bdry[i_bdry]*SERVICE;
			 
  
//==========block delta_adjoint\control ========
		   if ( i_vol < nDof_adj    && j < nDof_ctrl)   
		     Jac[ 
			(nDof_u+ nDof_ctrl + i_vol) * nDof_AllVars  +
		        (nDof_u  + j)             ]  += control_node_flag[i_vol] *
		        (
			  weight_bdry * grad_ctrl_dot_n_mat * phi_adj_bdry[i_bdry]
			);    		      		      
		      
		    }   //end loop i_bdry // j_vol
	      



//========= debugging ==========================    
//   std::cout << "====== phi values on boundary for gauss point " << ig_bdry << std::endl;
//   
//      for(unsigned i=0; i < nve_bdry; i ++) {
//      std::cout << phi_ctrl_bdry[i] << " ";
//      }
//    std::cout << std::endl;
 
//   std::cout << "====== phi derivatives on boundary for gauss point " << ig_bdry << std::endl;
//   
//   for (unsigned d = 0; d < dim; d++) {
//      for(unsigned i_bdry=0; i_bdry < nve_bdry; i_bdry ++) {
//      std::cout << phi_ctrl_x_bdry[i_bdry + d*nve_bdry] << " ";
//      }
//   }
//========== debugging ==========================    

		  }  //end i loop
		}  //end ig_bdry loop
	      }    //end if control face
	      
	    }  //end if boundary faces
	  }    //end loop over faces
	  
	} //end if control element flag
  
  
  
			
///////////////////////End Boundary Integral///////////////////////////////////////////////////////////////////////////////////////////////////////////////

			
			
			
			
			
			
			
			
 //========= gauss value quantities ==================   
	double sol_u_gss = 0.;
	double sol_adj_gss = 0.;
	double sol_ctrl_gss = 0.;
	std::vector<double> sol_u_x_gss(dim);        std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	std::vector<double> sol_adj_x_gss(dim);      std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
	std::vector<double> sol_ctrl_x_gss(dim);     std::fill(sol_ctrl_x_gss.begin(), sol_ctrl_x_gss.end(), 0.);
 //===================================================   

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[kelGeom][solType_max]->GetGaussPointNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
	msh->_finiteElement[kelGeom][solType_u]   ->Jacobian(x, ig, weight, phi_u, phi_u_x, phi_u_xx);
        msh->_finiteElement[kelGeom][solType_ctrl]->Jacobian(x, ig, weight, phi_ctrl, phi_ctrl_x, phi_ctrl_xx);
        msh->_finiteElement[kelGeom][solType_adj] ->Jacobian(x, ig, weight, phi_adj, phi_adj_x, phi_adj_xx);
	
	std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
	std::fill(sol_ctrl_x_gss.begin(), sol_ctrl_x_gss.end(), 0.);
	
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
        
//==========FILLING WITH THE EQUATIONS ===========
	// *** phi_i loop ***
        for (unsigned i = 0; i < nDof_max; i++) {
	  
	      double laplace_rhs_du_adj_i = 0.;
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
	      
	      double laplace_rhs_dadj_u_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_adj )         laplace_rhs_dadj_u_i           +=  (phi_adj_x   [i * dim + kdim] * sol_u_x_gss[kdim]);
	      }
	      
	      double laplace_rhs_dadj_ctrl_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_adj )         laplace_rhs_dadj_ctrl_i        +=  (phi_adj_x   [i * dim + kdim] * sol_ctrl_x_gss[kdim]);
	      }
//======================Residuals=======================
          // FIRST ROW
	  if (i < nDof_u)                      Res[0      + i] += - weight * (target_flag * phi_u[i] * ( sol_u_gss + sol_ctrl_gss - u_des) - laplace_rhs_du_adj_i - 0.);
          // SECOND ROW
	  if (i < nDof_ctrl)  {
	      if ( control_el_flag == 1)       Res[nDof_u + i] +=  /*(control_node_flag[i]) **/ - weight *  (target_flag * phi_ctrl[i] * ( sol_u_gss*SERVICE + sol_ctrl_gss - u_des) 
													      + alpha * phi_ctrl[i] * sol_ctrl_gss
		                                                                                              - laplace_rhs_dctrl_adj_i * SERVICE
		                                                                                              + beta * laplace_rhs_dctrl_ctrl_i);
	      else if ( control_el_flag == 0)  Res[nDof_u + i] +=  /*(1 - control_node_flag[i]) **/ (- penalty_strong) * (sol_ctrl[i] - 0.);
	  }
          // THIRD ROW
          if (i < nDof_adj) Res[nDof_u + nDof_ctrl + i] += - weight * ( - laplace_rhs_dadj_u_i - laplace_rhs_dadj_ctrl_i - 0.) ;
//======================Residuals=======================
	      
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
		(0 + j)                         ]  += weight * target_flag * phi_u[j] *  phi_u[i];
              
	      // BLOCK  delta_state - control
              if ( i < nDof_u && j < nDof_ctrl )   
		Jac[ (0 + i) * nDof_AllVars   +
                (nDof_u + j)                    ]  += weight * target_flag  * phi_ctrl[j] *  phi_u[i];
	      
              // BLOCK  delta_state - adjoint
              if ( i < nDof_u && j < nDof_adj )  
		Jac[  (0 + i) * nDof_AllVars  +
                (nDof_u + nDof_ctrl + j)        ]  += weight * (-1) * laplace_mat_du_adj;
              
	      //=========== delta_control row ===========================     
	      if ( control_el_flag == 1)  {
	      
              //BLOCK delta_control - control
              if ( i < nDof_ctrl   && j < nDof_ctrl   )
		Jac[ (nDof_u + i) * nDof_AllVars +
		(nDof_u + j)                     ]  += ( control_node_flag[i]) * weight * ( beta * control_el_flag  * laplace_mat_dctrl_ctrl + alpha * control_el_flag * phi_ctrl[i] * phi_ctrl[j] + target_flag  * phi_ctrl[i] * phi_ctrl[j]);
              
	      //BLOCK delta_control - state
              if ( i < nDof_ctrl   && j < nDof_u   ) 
		Jac[ (nDof_u + i) * nDof_AllVars  +
		(0 + j)                          ]  += ( control_node_flag[i]) * weight * target_flag * phi_u[j] * phi_ctrl[i]*SERVICE;
	      
	      //BLOCK delta_control - adjoint
              if ( i < nDof_ctrl   && j < nDof_adj  ) 
		Jac[ (nDof_u + i) * nDof_AllVars  + 
		(nDof_u + nDof_ctrl + j)         ]  +=  ( control_node_flag[i]) * weight * (-1) * laplace_mat_dctrl_adj*SERVICE;
	      }
	      
	      else if ( control_el_flag == 0)  {  
		
              //BLOCK delta_control - control
               if ( i < nDof_ctrl   && j < nDof_ctrl &&  i==j ) {
		 Jac[ (nDof_u + i) * nDof_AllVars +
		 (nDof_u + j)                     ] += (1-control_node_flag[i]) * penalty_strong;
		}
	      
	   }
	      
	      //=========== delta_adjoint row ===========================
              // BLOCK delta_adjoint - state	      
              if ( i < nDof_adj && j < nDof_u )   
		Jac[ (nDof_u + nDof_ctrl + i) * nDof_AllVars +
		(0 + j)                           ] += weight * (-1) * laplace_mat_dadj_u;   
	      
              // BLOCK delta_adjoint - control   
              if ( i < nDof_adj && j < nDof_ctrl )  
		Jac[ (nDof_u + nDof_ctrl + i)  * nDof_AllVars +
		(nDof_u  + j)                     ] += weight * (-1) * laplace_mat_dadj_ctrl; 

	      
            } // end phi_j loop
          } // endif assemble_matrix

        } // end phi_i loop
        
      } // end gauss point loop

      
      
    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
	if (control_el_flag == 0) {  //elements that should have zero control
	  
         for (unsigned i = 0; i < nDof_max; i++) {
            for (unsigned j = 0; j < nDof_max; j++) {
// 	      std::cout << Jac[ i * nDof_AllVars +j ] << " " << std::endl;
// 	      std::cout << Jac[ (nDof_u + i) * nDof_AllVars +j ] << " " << std::endl;
// 	      std::cout << Jac[ (nDof_u + nDof_ctrl + i) * nDof_AllVars +j ] << " " << std::endl;
// 	      std::cout << Jac[ i * nDof_AllVars + (nDof_u + j) ] << " " << std::endl;
// 	      std::cout << " " << std::setfill(' ') << std::setw(10) << Jac[ (nDof_u + i) * nDof_AllVars + (nDof_u + j) ];
	      std::cout << " " << std::setfill(' ') << std::setw(10) << Jac[ (nDof_u + nDof_adj + i) * nDof_AllVars + (nDof_u + nDof_adj + j) ];
// 	      std::cout << Jac[ (nDof_u + nDof_ctrl + i) * nDof_AllVars + (nDof_u + j) ] << " " << std::endl;
// 	      std::cout << Jac[ i * nDof_AllVars + (nDof_u + nDof_ctrl + j) ] << " " << std::endl;
//       std::cout << Jac[ (nDof_u + i) * nDof_AllVars + (nDof_u + nDof_ctrl + j) ] << " " << std::endl;
// 	      std::cout << Jac[ (nDof_u + nDof_ctrl + i) * nDof_AllVars + (nDof_u + nDof_ctrl + j) ] << " " << std::endl;
	    }
	      std::cout << std::endl;
	 }
	 
	}
    std::cout << " ========== " << iel << " ================== " << std::endl;     
    
    //copy the value of the adept::adoube aRes in double Res and store
    RES->add_vector_blocked(Res, l2GMap_AllVars);

    if (assembleMatrix) {
      //store K in the global matrix KK
      KK->add_matrix_blocked(Jac, l2GMap_AllVars, l2GMap_AllVars);
    }
  } //end element loop for each process

  RES->close();

  if (assembleMatrix) KK->close();

  // ***************** END ASSEMBLY *******************

  return;
}


};


} //end namespace elliptic

} //end namespace femus



#endif
