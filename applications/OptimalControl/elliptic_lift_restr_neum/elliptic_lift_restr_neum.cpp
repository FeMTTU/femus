#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "VTKWriter.hpp"
#include "LinearImplicitSystem.hpp"
#include "NumericVector.hpp"

#define FACE_FOR_CONTROL 3  //we do control on the right (=2) face
#define AXIS_DIRECTION_CONTROL_SIDE  0  //change this accordingly to the other variable above
#include "../elliptic_param.hpp"

#define SERVICE 1.

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

double InitialValueControl(const std::vector < double >& x) {
  return 0.; 
}

bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int faceName, const double time) {

  bool dirichlet = true; //dirichlet
  value = 0;
  
  if(!strcmp(name,"state")) {
  if (faceName == 3)
    dirichlet = false;
  }
    
   if(!strcmp(name,"adjoint")) {
  if (faceName == 3) { value = 0.;   dirichlet = false; }
  }
  
  return dirichlet;
}


double ComputeIntegral(MultiLevelProblem& ml_prob);

void AssembleLiftRestrProblem(MultiLevelProblem& ml_prob);


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

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
  mlSol.AddSolution("TargReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  mlSol.AddSolution("ContReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field

  
  mlSol.Initialize("All");    // initialize all varaibles to zero

  mlSol.Initialize("state", InitialValueState);
  mlSol.Initialize("control", InitialValueControl);
  mlSol.Initialize("adjoint", InitialValueAdjoint);
  mlSol.Initialize("TargReg", InitialValueTargReg);
  mlSol.Initialize("ContReg", InitialValueContReg);

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("state");
  mlSol.GenerateBdc("control");
  mlSol.GenerateBdc("adjoint");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);
  
 // add system  in mlProb as a Linear Implicit System
  LinearImplicitSystem& system = mlProb.add_system < LinearImplicitSystem > ("LiftRestr");
 
  system.AddSolutionToSystemPDE("state");  
  system.AddSolutionToSystemPDE("control");  
  system.AddSolutionToSystemPDE("adjoint");  
  
  // attach the assembling function to system
  system.SetAssembleFunction(AssembleLiftRestrProblem);

  // initilaize and solve the system
  system.init();
  system.MGsolve();
  
  ComputeIntegral(mlProb);
 
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("state");
  variablesToBePrinted.push_back("control");
  variablesToBePrinted.push_back("adjoint");
  variablesToBePrinted.push_back("TargReg");
  variablesToBePrinted.push_back("ContReg");

    // ******* Print solution *******
  mlSol.SetWriter(VTK);
  mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  return 0;
}


void AssembleLiftRestrProblem(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data

  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("LiftRestr");   // pointer to the linear implicit system named "LiftRestr"
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

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
  vector < vector < double > > x(dim);    // local coordinates
  vector < vector < double> >  x_bdry(dim);
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
  vector <double> phi_u;  // local test function
  vector <double> phi_u_x; // local test function first order partial derivatives
  vector <double> phi_u_xx; // local test function second order partial derivatives

  phi_u.reserve(maxSize);
  phi_u_x.reserve(maxSize * dim);
  phi_u_xx.reserve(maxSize * dim2);
  
 
  unsigned solIndex_u;
  solIndex_u = mlSol->GetIndex("state");    // get the position of "state" in the ml_sol object
  unsigned solType_u = mlSol->GetSolutionType(solIndex_u);    // get the finite element type for "state"

  unsigned solPdeIndexThom;
  solPdeIndexThom = mlPdeSys->GetSolPdeIndex("state");    // get the position of "state" in the pdeSys object

  vector < double >  sol_u; // local solution
  sol_u.reserve(maxSize);
  vector< int > l2GMap_u;
  l2GMap_u.reserve(maxSize);
 //***************************************************  
 //***************************************************  

  
 //******************** control ********************** 
 //***************************************************   
  vector <double> phi_ctrl;  // local test function
  vector <double> phi_ctrl_x; // local test function first order partial derivatives
  vector <double> phi_ctrl_xx; // local test function second order partial derivatives

  phi_ctrl.reserve(maxSize);
  phi_ctrl_x.reserve(maxSize * dim);
  phi_ctrl_xx.reserve(maxSize * dim2);
  
  unsigned solIndex_ctrl;
  solIndex_ctrl = mlSol->GetIndex("control");
  unsigned solType_ctrl = mlSol->GetSolutionType(solIndex_ctrl);

  unsigned solPdeIndex_ctrl;
  solPdeIndex_ctrl = mlPdeSys->GetSolPdeIndex("control");

  vector < double >  sol_ctrl; // local solution
  sol_ctrl.reserve(maxSize);
  vector< int > l2GMap_ctrl;
  l2GMap_ctrl.reserve(maxSize);
 
  std::string ctrl_name = "control";
  
 //*************************************************** 

  vector <double> sol_ctrl_x_vol_at_bdry_gss;
  sol_ctrl_x_vol_at_bdry_gss.reserve(dim);
 //*************************************************** 
   
  vector <double> phi_ctrl_vol_at_bdry;  // local test function
  phi_ctrl_vol_at_bdry.reserve(maxSize);
  
  vector <double> phi_ctrl_x_vol_at_bdry; // local test function first order partial derivatives
  phi_ctrl_x_vol_at_bdry.reserve(maxSize * dim);

 
 //***************************************************  
 //***************************************************  
  
  
 //********************* adjoint ********************* 
 //***************************************************  
  vector <double> phi_adj;  // local test function
  vector <double> phi_adj_x; // local test function first order partial derivatives
  vector <double> phi_adj_xx; // local test function second order partial derivatives

  phi_adj.reserve(maxSize);
  phi_adj_x.reserve(maxSize * dim);
  phi_adj_xx.reserve(maxSize * dim2);
  
  unsigned solIndex_adj;
  solIndex_adj = mlSol->GetIndex("adjoint");    // get the position of "state" in the ml_sol object
  unsigned solType_adj = mlSol->GetSolutionType(solIndex_adj);    // get the finite element type for "state"

  unsigned solPdeIndex_adj;
  solPdeIndex_adj = mlPdeSys->GetSolPdeIndex("adjoint");    // get the position of "state" in the pdeSys object

   //boundary adjoint shape functions  
  vector <double> phi_adj_bdry;  
  vector <double> phi_adj_x_bdry; 

  phi_adj_bdry.reserve(maxSize);
  phi_adj_x_bdry.reserve(maxSize * dim);
  
  //***************************************************  
vector < double >  sol_adj; // local solution
    sol_adj.reserve(maxSize);
  vector< int > l2GMap_adj;
    l2GMap_adj.reserve(maxSize);
    
    vector < double > sol_adj_bdry;
     sol_adj.reserve(maxSize);
 //***************************************************  
 //***************************************************  

    
    

  
 //***************************************************
 //********* WHOLE SET OF VARIABLES ****************** 
  const int solType_max = 2;  //biquadratic

  const int n_vars = 3;
 
  vector< int > l2GMap_AllVars; // local to global mapping
  l2GMap_AllVars.reserve(n_vars*maxSize);
  
  vector< double > Res; // local redidual vector
  Res.reserve(n_vars*maxSize);

  vector < double > Jac;
  Jac.reserve( n_vars*maxSize * n_vars*maxSize);
 //***************************************************  

  
 //********************* DATA ************************ 
  double u_des = DesiredTarget();
  double alpha = ALPHA_CTRL_VOL;
  double beta  = BETA_CTRL_VOL;
  double penalty_strong = 10e+14;
  //double penalty_ctrl = 1.e10;         //penalty for u=q
 //***************************************************  
  
  
  if (assembleMatrix)  KK->zero();

    
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned kelGeom = msh->GetElementType(iel);    // element geometry type

 //******************** GEOMETRY ********************* 
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
    for (int i = 0; i < dim; i++)  x[i].resize(nDofx);
    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);  // global to global mapping between coordinates node and coordinate dof

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
  
 //****** set target domain flag ********************* 
   int target_flag = 0;
   target_flag = ElementTargetFlag(elem_center);
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
  control_el_flag = ControlDomainFlag_internal_restriction(elem_center);
  std::vector<int> control_node_flag(nDof_ctrl,0);
  if (control_el_flag == 1) std::fill(control_node_flag.begin(), control_node_flag.end(), 1);
 //*************************************************** 


  
//////////////////Begin Boundary Integral////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if (control_el_flag == 1) {
	  
	  double tau=0.;
	  vector<double> normal(dim,0);
// 	  double normal_fixed[3] = {0.,1.,0};
	       
	  // loop on faces of the current element

	  for(unsigned jface=0; jface < msh->GetElementFaceNumber(iel); jface++) {
            std::vector < double > xyz_bdc(3,0.);  //not being used, because the boundaries are identified by the face numbers
	    // look for boundary faces
	    if(el->GetFaceElementIndex(iel,jface) < 0) {
	      unsigned int face = -( msh->el->GetFaceElementIndex(iel,jface)+1);
	      
		
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
		      x_bdry[idim][i] = (*msh->_topology->_Sol[idim])(iDof);
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
		  msh->_finiteElement[kelGeom][solType_ctrl]->VolumeShapeAtBoundary(x,x_bdry,jface,ig_bdry,phi_ctrl_vol_at_bdry,phi_ctrl_x_vol_at_bdry);
		  
		      
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



double ComputeIntegral(MultiLevelProblem& ml_prob)    {
  
  
    LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("LiftRestr");   // pointer to the linear implicit system named "LiftRestr"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

 //*************************************************** 
  vector < vector < double > > x(dim);    // local coordinates
  vector < vector < double> >  x_bdry(dim);
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
  for (unsigned i = 0; i < dim; i++) {
    x[i].reserve(maxSize);
    x_bdry[i].reserve(maxSize);
  }
 //*************************************************** 

 //*************************************************** 
  double weight; // gauss point weight
  double weight_bdry = 0.; // gauss point weight on the boundary
  
 //***************************************************  
  double alpha = ALPHA_CTRL_VOL;
  double beta  = BETA_CTRL_VOL;

 //******************** state ************************ 
 //*************************************************** 
  vector <double> phi_u;  // local test function
  vector <double> phi_u_x; // local test function first order partial derivatives
  vector <double> phi_u_xx; // local test function second order partial derivatives

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
  vector <double> phi_ctrl;  // local test function
  vector <double> phi_ctrl_x; // local test function first order partial derivatives
  vector <double> phi_ctrl_xx; // local test function second order partial derivatives

  phi_ctrl.reserve(maxSize);
  phi_ctrl_x.reserve(maxSize * dim);
  phi_ctrl_xx.reserve(maxSize * dim2);
  
 //**************** control bdry**********************  
  vector <double> phi_ctrl_bdry;  
  vector <double> phi_ctrl_x_bdry; 

  phi_ctrl_bdry.reserve(maxSize);
  phi_ctrl_x_bdry.reserve(maxSize * dim);
  
 //***************************************************
  unsigned solIndex_ctrl;
  solIndex_ctrl = mlSol->GetIndex("control");
  unsigned solType_ctrl = mlSol->GetSolutionType(solIndex_ctrl);

  vector < double >  sol_ctrl; // local solution
  sol_ctrl.reserve(maxSize);
 vector< int > l2GMap_ctrl;
  l2GMap_ctrl.reserve(maxSize);
  
  double ctrl_gss = 0.;
  double ctrl_x_gss = 0.;
  

 //*************************************************** 
 //***************************************************  

  
 //********************* desired ********************* 
 //*************************************************** 
  vector <double> phi_udes;  // local test function
  vector <double> phi_udes_x; // local test function first order partial derivatives
  vector <double> phi_udes_xx; // local test function second order partial derivatives

    phi_udes.reserve(maxSize);
    phi_udes_x.reserve(maxSize * dim);
    phi_udes_xx.reserve(maxSize * dim2);
 
  
//  unsigned solIndexTdes;
//   solIndexTdes = mlSol->GetIndex("Tdes");    // get the position of "state" in the ml_sol object
//   unsigned solTypeTdes = mlSol->GetSolutionType(solIndexTdes);    // get the finite element type for "state"

  vector < double >  sol_udes; // local solution
  sol_udes.reserve(maxSize);
 vector< int > l2GMap_Tdes;
    l2GMap_Tdes.reserve(maxSize);
  double udes_gss = 0.;
 //*************************************************** 
 //*************************************************** 


 //*************************************************** 
 //********* WHOLE SET OF VARIABLES ****************** 
  const int solType_max = 2;  //biquadratic

  const int n_vars = 3;
 
  vector< int > l2GMap_AllVars; // local to global mapping
  l2GMap_AllVars.reserve(n_vars*maxSize);
  
  vector< double > Res; // local redidual vector
  Res.reserve(n_vars*maxSize);

  vector < double > Jac;
  Jac.reserve( n_vars*maxSize * n_vars*maxSize);
 //*************************************************** 

  
 //********************* DATA ************************ 
  double u_des = DesiredTarget();
 //*************************************************** 
  
  double integral_target = 0.;
  double integral_alpha  = 0.;
  double integral_beta   = 0.;

    
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned kelGeom = msh->GetElementType(iel);    // element geometry type
 
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
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u);    // number of solution element dofs
        sol_u    .resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
     unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u);  // global to global mapping between solution node and solution dof
            sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);            // global extraction and local storage for the solution
    }
 //*************************************************** 


 //************** control **************************** 
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);    // number of solution element dofs
    sol_ctrl    .resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl);    // global to global mapping between solution node and solution dof
      sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);      // global extraction and local storage for the solution
    } 
 //***************************************************  
 
 
 //**************** u_des **************************** 
    unsigned nDof_udes  = msh->GetElementDofNumber(iel, solType_u);    // number of solution element dofs
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
      for (unsigned ig = 0; ig < msh->_finiteElement[kelGeom][solType_max]->GetGaussPointNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
        //  ==== State 
	msh->_finiteElement[kelGeom][solType_u]   ->Jacobian(x, ig, weight, phi_u, phi_u_x, phi_u_xx);
        //  ==== Adjoint 
        msh->_finiteElement[kelGeom][solType_u/*solTypeTdes*/]->Jacobian(x, ig, weight, phi_udes, phi_udes_x, phi_udes_xx);
        //  ==== Control 
        msh->_finiteElement[kelGeom][solType_ctrl]  ->Jacobian(x, ig, weight, phi_ctrl, phi_ctrl_x, phi_ctrl_xx);

	u_gss = 0.;  for (unsigned i = 0; i < nDof_u; i++) u_gss += sol_u[i] * phi_u[i];		
	ctrl_gss = 0.; for (unsigned i = 0; i < nDof_ctrl; i++) ctrl_gss += sol_ctrl[i] * phi_ctrl[i];  
	udes_gss  = 0.; for (unsigned i = 0; i < nDof_udes; i++)  udes_gss  += sol_udes[i]  * phi_udes[i]; 
        ctrl_x_gss  = 0.; for (unsigned i = 0; i < nDof_ctrl; i++)  
        {
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
  
  

